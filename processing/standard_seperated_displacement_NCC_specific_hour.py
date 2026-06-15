#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 19:33:14 2020

@author: laserglaciers
"""

import subprocess
import os
import logging
import numpy as np
from multiprocessing import Pool
#from subProcessUtils import date_sort
import pandas as pd
from obstore.store import S3Store
import obstore
from datetime import datetime
import re
from tqdm import tqdm

LOG_FILE = f"./log/atlas_runs_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
"""
Command for displacement calcs

./atlas /home/laserglaciers/may_laz/laz/200501_010203.laz /home/laserglaciers/may_laz/laz/200501_072310.laz
 --xshift=6.4 --yshift=-.51

 For COPC files in S3
  ./atlas \
      https://atlas-lidar-helheim.s3.us-west-2.amazonaws.com/copc/south/230108_131655.copc.laz \
      https://atlas-lidar-helheim.s3.us-west-2.amazonaws.com/copc/south/230109_131631.copc.laz \
      --xshift=25.6 --yshift=-2.04 \
      --ncc \
      --tiff /media/m484s199/qaanaaq/helheim-atlas-lidar/velocities

"""


# get atlas-lidar-helheim bucket
# skip_signature is like --no-sign-request
bucket = "atlas-lidar-helheim"
region = "us-west-2"
store = S3Store(bucket, region=region, skip_signature=True)
listing = obstore.list(store, prefix="copc/south/").collect()
keys = sorted(obj["path"] for obj in listing if obj["path"].endswith(".copc.laz"))

DATE_START   = datetime(2019, 1, 1)
DATE_END     = datetime(2019, 12, 31)
out_dir = '/opt/atlas/helheim-atlas-lidar/output/ncc_output' # output for the displacement tif

def parse_date(key):
    m = re.search(r"(\d{6})_(\d{6})\.copc\.laz", key)
    if m:
        try:
            return datetime.strptime(m.group(1) + m.group(2), "%y%m%d%H%M%S")
        except ValueError:
            return None

filtered = sorted(
    [(k, parse_date(k)) for k in keys
     if parse_date(k) and DATE_START <= parse_date(k) <= DATE_END],
    key=lambda x: x[1],
)

def hour_match(scans, target_hours=24, tolerance_hours=3):
    """
    Pair each scan with the *later* scan whose separation is closest to
    ``target_hours`` (default 24 h). ``scans`` is the sorted list of
    ``(key, datetime)`` tuples (i.e. ``filtered``).

    A match is kept only when its separation lands within ``tolerance_hours``
    of the target, so gaps in the record don't produce a bogus pair spanning
    several days. Returns a list of ``(key_1_url, key_2_url, shift_factor)``
    tuples in the same shape as the old ``diff_tups``, where ``shift_factor``
    scales the rule-of-thumb 6-hour xshift/yshift (24 h -> factor 4.0).
    """
    hour = np.timedelta64(1, 'h')
    matches = []
    n = len(scans)
    for i in range(n):
        key_i, time_i = scans[i]
        best_j = None
        best_err = None
        best_diff = None
        for j in range(i + 1, n):
            key_j, time_j = scans[j]
            diff_h = (time_j - time_i) / hour
            err = abs(diff_h - target_hours)
            if best_err is None or err < best_err:
                best_err, best_j, best_diff = err, j, diff_h
            # scans are time-sorted: once past the target, error only grows
            if diff_h >= target_hours:
                break

        if best_j is None or best_err > tolerance_hours:
            continue

        key_j, time_j = scans[best_j]
        shift_factor = round(best_diff) / 6
        print(f'time 1: {time_i} time 2: {time_j} '
              f'diff: {best_diff:.2f} h shift factor: {shift_factor}')
        matches.append((aws_http(key_i), aws_http(key_j), shift_factor))

    return matches


keys_filtered  = [k for k, _ in filtered]
dates_filtered = [d for _, d in filtered]
print(f"{len(keys_filtered)} scans: {filtered[0][1]} → {filtered[-1][1]}")


def aws_http(key, bucket=bucket, region=region):
    """
    example: https://atlas-lidar-helheim.s3.us-west-2.amazonaws.com/copc/south/230108_131655.copc.laz
    """

    return f'https://{bucket}.s3.{region}.amazonaws.com/{key}'

## Run atlas command like commented command below import statements and rename the svg file
## to the same base name as laz and move
args_list = []
xshift = 6.4 # rule of thumb for 6 hour scan difference
yshift = -.51 # rule of thumb for 6 hour scan difference

# Pair scans ~24 h apart (closest match within tolerance) instead of
# consecutive scans. Change target_hours to retarget the separation.
diff_tups = hour_match(filtered, target_hours=24)
print(f"{len(diff_tups)} matched pairs (~24 h apart)")


def output_exists(copc_1, copc_2):
    stem1 = os.path.basename(copc_1).removesuffix('.copc.laz')
    stem2 = os.path.basename(copc_2).removesuffix('.copc.laz')
    out_file = os.path.join(out_dir, f"{stem1}_{stem2}.tif")
    return os.path.exists(os.path.join(out_dir, f"{stem1}_{stem2}.tif")), out_file

args_list = []
for i,laz in enumerate(diff_tups):

    copc_1, copc_2, shift_factor = laz
    if output_exists(copc_1, copc_2)[0]:
        print(f"Skipping (already done): {output_exists(copc_1, copc_2)[1]}")
        continue
#    print(f'shift factor: {shift_factor}')
    # print('Calculating displacements between '+' '.join([laz[0],sortDates[i+1][0]]))
    args = ['/opt/atlas/helheim-atlas-lidar/displacement/build/atlas', copc_1, copc_2,
            '--xshift='+str(round(xshift*shift_factor,4)),
            '--yshift='+str(round(yshift*shift_factor,4)),
            '--ncc',
            '--tiff', out_dir]
    args_list.append(args)


def worker(cmd):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[
            logging.FileHandler(LOG_FILE),
        ],
    )
    cmd_str = " ".join(cmd)
    logging.info("START: %s", cmd_str)
    proc = subprocess.Popen(cmd, stderr=subprocess.PIPE)
    _, stderr = proc.communicate()
    if proc.returncode != 0:
        logging.error("FAILED (rc=%d): %s\n%s", proc.returncode, cmd_str, stderr.decode().strip())
    else:
        logging.info("DONE (rc=%d): %s", proc.returncode, cmd_str)



#with Pool(processes=16) as pool:
#    results = pool.map(worker, args_list)


if __name__ == "__main__":
    with Pool(processes=16) as pool:
          results = list(tqdm(pool.imap_unordered(worker, args_list), total=len(args_list)))
# procs = [ subprocess.Popen(i) for i in args_list[:30]]
# for p in procs:
# #    p.wait()
