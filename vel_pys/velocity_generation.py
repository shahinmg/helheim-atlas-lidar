"""
Generate velocity from atlas displacement data
"""

import os
import logging
import traceback
from datetime import datetime
from multiprocessing import Pool

import rasterio
import numpy as np
import pandas as pd
from tqdm import tqdm

LOG_FILE = f"./log/velocity_generation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)

vel_path = "/media/m484s199/qaanaaq/helheim_vels/ncc_output/"
outpath = "/media/m484s199/qaanaaq/helheim_vels/velocity_unsmoothed/"



# file list
displacement_tifs = sorted( [file for file in os.listdir(vel_path) if file.endswith('.tif')] )

# bands to compute velocity
displacement_bands = [1,2,3, 6,7,8, 14, 15]


def setup_logging():
    """
    Configure root logging to write only to the log file
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(process)d] %(levelname)s %(message)s",
        handlers=[
            logging.FileHandler(LOG_FILE),
        ],
        force=True,
    )


def get_time_delta(file_name):
    """
    compute the time difference in hours between scan 1 and scan 2 time

    Parameters
    ----------
    file_name: str
        ATLAS displacement file time
        YYMMDD_HHMMSS_YYMMDD_HHMMSS (time_1_time_2)


    Returns:
        time_delta: float
            The difference between time_2 - time_1 in hours

    """
    time_1 = pd.to_datetime(f'20{file_name[:13]}'.replace('_',''))
    time_2 = pd.to_datetime(f'20{file_name[14:-4]}'.replace('_',''))

    time_delta = time_2 - time_1
    time_delta = time_delta/np.timedelta64(1,'h')

    return time_delta

def calc_velocity(file, in_path=vel_path, out_path=outpath):
    in_file = f'{in_path}{file}'
    out_file = f'{out_path}{file[:-4]}_vel.tif'

    if os.path.exists(out_file):
        logging.info("SKIP: output already exists: %s", out_file)
        return (file, "skipped", None)

    logging.info("START: %s", in_file)

    try:
        with rasterio.open(in_file) as src:

            band_dict = {name: index for index, name in zip(src.indexes, src.descriptions)}

            time_delta = get_time_delta(file)
            logging.info("%s: time_delta = %.4f hours", file, time_delta)

            # read mean bands
            x_mean = src.read(1, masked=True)
            y_mean = src.read(2, masked=True)
            z_mean = src.read(3, masked=True)

            # read median bands
            x_median = src.read(6, masked=True)
            y_median = src.read(7, masked=True)
            z_median = src.read(8, masked=True)

            # read ncc bands
            ncc_x = src.read(14, masked=True)
            ncc_y = src.read(15, masked=True)

            # compute velocity in meters per day
            # compute mean vel
            X_vel_mean = ( x_mean / time_delta ) * 24.0
            Y_vel_mean = ( y_mean / time_delta ) * 24.0
            Z_vel_mean = ( z_mean / time_delta ) * 24.0

            absH_mean = np.sqrt( np.square(X_vel_mean) + np.square(Y_vel_mean) ) #  Magnitude of horizontal velocity


            # compute median velocity
            X_vel_median = ( x_median / time_delta ) * 24.0
            Y_vel_median = ( y_median / time_delta ) * 24.0
            Z_vel_median = ( z_median / time_delta ) * 24.0

            absH_median = np.sqrt( np.square(X_vel_median) + np.square(Y_vel_median) ) #  Magnitude of horizontal velocity

            # compute ncc velocity
            X_vel_ncc = ( ncc_x / time_delta ) * 24.0
            Y_vel_ncc = ( ncc_y / time_delta ) * 24.0
            absH_ncc = np.sqrt( np.square(X_vel_ncc) + np.square(Y_vel_ncc) ) #  Magnitude of horizontal velocity

            # make dictionary of bands I want to write
            bands_out = {
                        'X_vel_mean': X_vel_mean,
                        'Y_vel_mean': Y_vel_mean,
                        'Z_vel_mean': Z_vel_mean,
                        'BEFORE': src.read(band_dict['BEFORE'], masked=True),
                        'AFTER': src.read(band_dict['AFTER'], masked=True),
                        'X_vel_median': X_vel_median,
                        'Y_vel_median': Y_vel_median,
                        'Z_vel_median': Z_vel_median,
                        'MATCH_COUNT': src.read(band_dict['MATCH_COUNT'], masked=True),
                        'RMS_RESIDUAL': src.read(band_dict['RMS_RESIDUAL'], masked=True),
                        'MAD_X': src.read(band_dict['MAD_X'], masked=True),
                        'MAD_Y': src.read(band_dict['MAD_Y'], masked=True),
                        'MAD_Z': src.read(band_dict['MAD_Z'], masked=True),
                        'X_vel_ncc': X_vel_ncc,
                        'Y_vel_ncc': Y_vel_ncc,
                        'NCC_PEAK': src.read(band_dict['NCC_PEAK'], masked=True),
                        'Absolute_horizontal_vel_mean': absH_mean,
                        'Absolute_horizontal_vel_median': absH_median,
                        'Absolute_horizontal_vel_ncc': absH_ncc
                        }



            profile = src.profile.copy()
            profile['count'] = len(bands_out)

        # fill masked (nodata) cells with the raster's declared nodata value so
        # downstream tools treat them as nodata rather than as real velocities
        nodata = profile.get('nodata')
        if nodata is None:
            nodata = -9999.0
            profile['nodata'] = nodata

        with rasterio.open(out_file, mode='w', **profile) as dst:
            for band_num, (name, band_data) in enumerate(bands_out.items(), start=1):
                if np.ma.isMaskedArray(band_data):
                    band_data = band_data.filled(nodata)
                dst.write(band_data.astype(profile['dtype']), band_num)
                dst.set_band_description(band_num, name)

        logging.info("DONE: wrote %d bands to %s", len(bands_out), out_file)
        return (file, "ok", None)

    except Exception as exc:
        logging.error(
            "FAILED: %s -> %s\n%s",
            in_file,
            out_file,
            traceback.format_exc(),
        )
        return (file, "failed", str(exc))


if __name__ == "__main__":
    setup_logging()
    logging.info(
        "Processing %d displacement tifs from %s -> %s",
        len(displacement_tifs), vel_path, outpath,
    )
    os.makedirs(outpath, exist_ok=True)

    with Pool(processes=os.cpu_count(), initializer=setup_logging) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(calc_velocity, displacement_tifs),
                total=len(displacement_tifs),
            )
        )

    succeeded = [f for f, status, _ in results if status == "ok"]
    skipped = [f for f, status, _ in results if status == "skipped"]
    failed = [(f, err) for f, status, err in results if status == "failed"]

    logging.info(
        "COMPLETE: %d succeeded, %d skipped, %d failed",
        len(succeeded), len(skipped), len(failed),
    )
    for f, err in failed:
        logging.error("  failed: %s (%s)", f, err)
