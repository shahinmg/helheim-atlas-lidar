#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matcher sweep driver for the 6-hour product.

Runs the displacement tool over CONSECUTIVE 6-hour scan pairs (no hour
filtering) for several matcher configurations, then samples the center point
and tabulates MATCH_COUNT / detrended-sigma / RMS_RESIDUAL per config, split
day vs night. Use it to A/B the mutual-NN matcher knobs (zgate / zgate-plane /
matchmin / matchmax) added on the noise-reduction branch.

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
# Config
# ---------------------------------------------------------------------------
ATLAS_BIN = '/media/m484s199/qaanaaq/helheim-atlas-lidar/displacement/build/atlas'
OUT_BASE  = '/media/m484s199/qaanaaq/helheim_vels/displacements/sweep_6hr'  # per-config subdirs
LOG_FILE  = f"./log/sweep_matcher_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

bucket = "atlas-lidar-helheim"
region = "us-west-2"

# Summer window = the daytime-penalty regime we want to fix (freeze-up hasn't
# happened yet). Stage 1: a few days to check direction. Stage 2: ~2-3 weeks.
# Stage-1 window: 3 days (01/07/13/19 each) -> ~12 scans, ~11 consecutive pairs,
# all four time-of-day windows covered. Widen to ~2-3 weeks for Stage 2.
DATE_START = datetime(2019, 7, 15)
DATE_END   = datetime(2019, 7, 18)

xshift = 6.4   # rule of thumb for 6 hour scan difference
yshift = -.51  # rule of thumb for 6 hour scan difference

NPROC = 16

# Point + buffer for sampling (same geom/buffer as the velocity CSVs).
POINT_PATH = '/media/m484s199/qaanaaq/helheim_vels/geoms/test_center_point_utm24N.gpkg'
BUFFER = 200

# Each config = extra atlas args beyond the positional/xshift/yshift. NCC is
# omitted on purpose: the matcher change lives in bands 1-13, so dropping --ncc
# isolates it and runs faster. Edit configs 4/5 once stage 1 picks abs-Z vs plane.
SWEEP_CONFIGS = {
    'baseline':          [],                                                  # mutual-NN, no gate (committed defaults)
    'zgate3':            ['--zgate=3'],                                        # absolute-Z gate
    'zgate3_mm15':       ['--zgate=3', '--matchmin=1.5'],                      # abs-Z gate + tighten interior
    'zgate3_plane':      ['--zgate=3', '--zgate-plane'],                       # plane-residual gate
    'zgate3_plane_mm15': ['--zgate=3', '--zgate-plane', '--matchmin=1.5'],     # plane gate + tighten interior
    'zgate3_plane_mm10': ['--zgate=3', '--zgate-plane', '--matchmin=1.0'],     # plane gate + tighten more
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

    df.to_csv(os.path.join(OUT_BASE, 'sweep_samples.csv'), index=False)

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
    sdf.to_csv(os.path.join(OUT_BASE, 'sweep_summary.csv'), index=False)
    print("\n=== sweep summary (match up + noise down, without rms blowup = win) ===")
    print(sdf.to_string(index=False))
    print(f"\nwrote {OUT_BASE}/sweep_summary.csv and sweep_samples.csv")


if __name__ == "__main__":
    if args_list:
        with Pool(processes=NPROC) as pool:
            list(tqdm(pool.imap_unordered(worker, args_list), total=len(args_list)))
    analyze()
