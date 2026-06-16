#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matcher sweep driver for the 6-hour product -- SERVER paths (/opt/atlas).

Same as sweep_matcher_6hr.py but with the binary/output (and sample geom)
pointed at the /opt/atlas compute server, where the run is heavy enough to
crash a workstation. Requires the noise-reduction build on the server.

Runs the displacement tool over CONSECUTIVE 6-hour scan pairs (no hour
filtering) for several matcher configurations, then samples the center point
and tabulates MATCH_COUNT / detrended-sigma / RMS_RESIDUAL per config, split
day vs night. A/B for the mutual-NN matcher knobs (zgate / zgate-plane /
matchmin / matchmax).

Structure mirrors standard_seperated_displacement_NCC.py / _V2.py.

@author: laserglaciers
"""

import subprocess
import os
import re
import glob
import logging
from datetime import datetime
from multiprocessing import Pool

import numpy as np
import pandas as pd
from tqdm import tqdm
from obstore.store import S3Store
import obstore

# ---------------------------------------------------------------------------
# Config -- SERVER paths
# ---------------------------------------------------------------------------
ATLAS_BIN = '/opt/atlas/helheim-atlas-lidar/displacement/build/atlas'
OUT_BASE  = '/opt/atlas/helheim-atlas-lidar/output/sweep_6hr'  # per-config subdirs
LOG_FILE  = f"./log/sweep_matcher_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

bucket = "atlas-lidar-helheim"
region = "us-west-2"

# Summer window = the daytime-penalty regime we want to fix (freeze-up hasn't
# happened yet). Stage 1: a few days to check direction. Stage 2: ~2-3 weeks.
# Stage-2 window: ~3 weeks, so noise_resid_std and the cross-tod floor are
# trustworthy. (Stage-1 direction check used 3 days: DATE_END = 2019-07-18.)
DATE_START = datetime(2019, 7, 15)
DATE_END   = datetime(2019, 8, 5)

xshift = 6.4   # rule of thumb for 6 hour scan difference
yshift = -.51  # rule of thumb for 6 hour scan difference

NPROC = 16

# Point + buffer for sampling. Must exist on the server for analyze() to run;
# otherwise sync OUT_BASE back and run analyze() on the workstation copy.
POINT_PATH = '/opt/atlas/helheim-atlas-lidar/geoms/test_center_point_utm24N.gpkg'
BUFFER = 200

# Tag appended to the analysis CSVs so successive sweeps don't overwrite each
# other (e.g. 'stage1', 'stage2'). Per-config tif dirs under OUT_BASE are reused
# across runs (skip-if-exists), so widening the window only adds new pairs.
SWEEP_TAG = 'stage2'

# Path 2: shape-GENERATION sweep (find more matches by detecting more shapes).
# minshape is already at its floor (default 1 = no filter; production doesn't
# pass --minshape), so the free lever is topfrac: raising it admits more points
# into the slice -> more occupied cells -> more/larger blobs in sparse daytime
# tiles. gridlen is left out (it rescales the cell-unit match threshold, a
# confound). NCC omitted to isolate the matcher and run faster. If a single
# topfrac band beats 0.5, multi-slice (pooling several bands) should beat it
# further -- that's the code change to do next.


# SWEEP_CONFIGS = {
#     'baseline':          [],                                                  # mutual-NN, no gate (committed defaults)                                                                    
#     'zgate3':            ['--zgate=3'],                                        # absolute-Z gate                                                                                           
#     'zgate3_mm15':       ['--zgate=3', '--matchmin=1.5'],                      # abs-Z gate + tighten interior                                                                             
#     'zgate3_plane':      ['--zgate=3', '--zgate-plane'],                       # plane-residual gate                                                                                       
#     'zgate3_plane_mm15': ['--zgate=3', '--zgate-plane', '--matchmin=1.5'],     # plane gate + tighten interior                                                                             
#     'zgate3_plane_mm10': ['--zgate=3', '--zgate-plane', '--matchmin=1.0'],     # plane gate + tighten more
# }

# Stage-1 generation sweep (commented; kept for the record). Result: slices2
# raised matches ~29% (daytime 11.8->15.1, up to baseline's nighttime level)
# with RMS flat; slices3 added more matches but noise crept up; topfrac was not
# the lever (topfrac09 even lowered daytime match). See memory.
# SWEEP_CONFIGS = {
#     'gen_base':     [],                              # topfrac 0.5 (current default)
#     'topfrac07':    ['--topfrac=0.7'],               # admit more points (denser slice)
#     'topfrac09':    ['--topfrac=0.9'],               # almost the whole cloud
#     'slices2':      ['--slices=2'],                  # split top 0.5 into 2 height bands
#     'slices3':      ['--slices=3'],                  # split top 0.5 into 3 height bands
#     'slices3_tf09': ['--slices=3', '--topfrac=0.9'], # 3 bands over a deeper slice
# }

# Stage-2: run the slice candidates head-to-head over a 3-week window where
# noise_resid_std is trustworthy. slices3_tf09 had the most matches AND the
# lowest formal sigma=RMS/sqrt(N) at stage 1, but its RMS rose above baseline --
# a hint its extra matches may be correlated (dedup at 1 cell didn't fully
# remove them), which would make that formal sigma optimistic. The empirical
# noise_resid_std here adjudicates.
SWEEP_CONFIGS = {
    'gen_base':     [],                              # baseline: top slice (topfrac 0.5)
    'slices2':      ['--slices=2'],                  # +29% matches, RMS flat/better
    'slices3':      ['--slices=3'],                  # best daytime match, slight noise creep
    'slices3_tf09': ['--slices=3', '--topfrac=0.9'], # most matches + lowest formal sigma; RMS up
}

# Displacement tif band layout (no --ncc): 6=MEDIAN_X, 7=MEDIAN_Y,
# 9=MATCH_COUNT, 10=RMS_RESIDUAL.
BAND_MEDX, BAND_MEDY, BAND_MATCH, BAND_RMS = 6, 7, 9, 10

os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)

# ---------------------------------------------------------------------------
# Scan listing / consecutive 6 h pairing (no hour filter)
# ---------------------------------------------------------------------------
store = S3Store(bucket, region=region, skip_signature=True)
listing = obstore.list(store, prefix="copc/south/").collect()
keys = sorted(obj["path"] for obj in listing if obj["path"].endswith(".copc.laz"))


def parse_date(key):
    m = re.search(r"(\d{6})_(\d{6})\.copc\.laz", key)
    if m:
        try:
            return datetime.strptime(m.group(1) + m.group(2), "%y%m%d%H%M%S")
        except ValueError:
            return None


def aws_http(key, bucket=bucket, region=region):
    return f'https://{bucket}.s3.{region}.amazonaws.com/{key}'


filtered = sorted(
    [(k, parse_date(k)) for k in keys
     if parse_date(k) and DATE_START <= parse_date(k) <= DATE_END],
    key=lambda x: x[1],
)
print(f"{len(filtered)} scans: {filtered[0][1]} -> {filtered[-1][1]}")

# Consecutive pairs (like _NCC.py). shift_factor scales xshift/yshift by the
# actual gap so a ~6 h pair gets factor 1, a 12 h gap gets 2, etc.
diff_tups = []
for i in range(len(filtered) - 1):
    key_i, time_i = filtered[i]
    key_j, time_j = filtered[i + 1]
    try:
        timeDiff = (time_j - time_i) / np.timedelta64(1, 'h')
        shift_factor = round(timeDiff) / 6
        diff_tups.append((aws_http(key_i), aws_http(key_j), shift_factor))
    except Exception as e:
        print(f"An error occurred: {e}")
print(f"{len(diff_tups)} consecutive pairs")


# ---------------------------------------------------------------------------
# Build the run list: every (pair x config), skipping existing outputs
# ---------------------------------------------------------------------------
def out_path(out_dir, copc_1, copc_2):
    stem1 = os.path.basename(copc_1).removesuffix('.copc.laz')
    stem2 = os.path.basename(copc_2).removesuffix('.copc.laz')
    return os.path.join(out_dir, f"{stem1}_{stem2}.tif")


args_list = []
for name, extra in SWEEP_CONFIGS.items():
    out_dir = os.path.join(OUT_BASE, name)
    os.makedirs(out_dir, exist_ok=True)
    for copc_1, copc_2, shift_factor in diff_tups:
        if os.path.exists(out_path(out_dir, copc_1, copc_2)):
            continue
        args = [ATLAS_BIN, copc_1, copc_2,
                '--xshift=' + str(round(xshift * shift_factor, 4)),
                '--yshift=' + str(round(yshift * shift_factor, 4)),
                *extra,
                '--tiff', out_dir]
        args_list.append(args)
print(f"{len(args_list)} atlas runs queued across {len(SWEEP_CONFIGS)} configs")


def worker(cmd):
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(LOG_FILE)],
    )
    cmd_str = " ".join(cmd)
    logging.info("START: %s", cmd_str)
    proc = subprocess.Popen(cmd, stderr=subprocess.PIPE)
    _, stderr = proc.communicate()
    if proc.returncode != 0:
        logging.error("FAILED (rc=%d): %s\n%s", proc.returncode, cmd_str,
                      stderr.decode().strip())
    else:
        logging.info("DONE (rc=%d): %s", proc.returncode, cmd_str)


# ---------------------------------------------------------------------------
# Analysis: sample the point, tabulate per config (day vs night)
# ---------------------------------------------------------------------------
def sample_config(name):
    """Per-pair samples for one config: horiz displacement, match, rms."""
    import geopandas as gpd
    import rasterio
    from rasterstats import zonal_stats

    pt = gpd.read_file(POINT_PATH).buffer(BUFFER)
    out_dir = os.path.join(OUT_BASE, name)
    rows = []
    for f in sorted(glob.glob(os.path.join(out_dir, '*.tif'))):
        m = re.search(r"(\d{6})_(\d{6})_(\d{6})_(\d{6})\.tif", os.path.basename(f))
        if not m:
            continue
        t2 = pd.to_datetime('20' + m.group(3) + m.group(4), format='%Y%m%d%H%M%S')
        with rasterio.open(f) as src:
            nod = src.nodata
        def samp(band):
            return zonal_stats(pt, f, band=band, stats=['mean'], nodata=nod)[0]['mean']
        mx, my = samp(BAND_MEDX), samp(BAND_MEDY)
        horiz = np.hypot(mx, my) if (mx is not None and my is not None) else np.nan
        rows.append({'config': name, 't2': t2, 't2_hour': t2.hour,
                     'horiz': horiz, 'match': samp(BAND_MATCH), 'rms': samp(BAND_RMS)})
    return rows


def analyze():
    rows = []
    for name in SWEEP_CONFIGS:
        rows += sample_config(name)
    df = pd.DataFrame(rows).dropna(subset=['horiz'])
    if df.empty:
        print("No samples found to analyze.")
        return

    df = df.sort_values(['config', 't2'])
    # detrended residual of horizontal displacement (noise proxy), per config
    df['resid'] = df.groupby('config')['horiz'].transform(
        lambda s: s - s.rolling(21, center=True, min_periods=7).median())
    df['period'] = np.where(df['t2_hour'].isin([13, 19]), 'day',
                    np.where(df['t2_hour'].isin([1, 7]), 'night', 'other'))

    df.to_csv(os.path.join(OUT_BASE, f'sweep_samples_{SWEEP_TAG}.csv'), index=False)

    # summary: per config x {all, day, night}
    summ = []
    for name, g in df.groupby('config'):
        for period, gg in [('all', g), ('day', g[g.period == 'day']),
                           ('night', g[g.period == 'night'])]:
            if len(gg) < 3:
                continue
            summ.append({'config': name, 'period': period, 'n': len(gg),
                         'match_med': round(gg['match'].median(), 1),
                         'rms_med': round(gg['rms'].median(), 3),
                         'noise_resid_std': round(gg['resid'].std(), 3)})
    sdf = pd.DataFrame(summ)
    sdf.to_csv(os.path.join(OUT_BASE, f'sweep_summary_{SWEEP_TAG}.csv'), index=False)
    print("\n=== sweep summary (match up + noise down, without rms blowup = win) ===")
    print(sdf.to_string(index=False))
    print(f"\nwrote {OUT_BASE}/sweep_summary_{SWEEP_TAG}.csv and sweep_samples_{SWEEP_TAG}.csv")


if __name__ == "__main__":
    if args_list:
        with Pool(processes=NPROC) as pool:
            list(tqdm(pool.imap_unordered(worker, args_list), total=len(args_list)))
    analyze()
