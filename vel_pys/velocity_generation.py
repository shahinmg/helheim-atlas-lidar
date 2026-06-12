"""
Generate velocity from atlas displacement data
"""

import os
import rasterio
import numpy as np
import pandas as pd
from multiprocessing import Pool

ncc_path = "/media/m484s199/qaanaaq/helheim_vels/ncc_output/"
out_path = "/media/m484s199/qaanaaq/helheim_vels/velocity_unsmoothed/"

# file list
displacement_tifs = sorted( [file for file in os.listdir(ncc_path) if file.endswith('.tif')] )

# bands to compute velocity
displacement_bands = [1,2,3, 6,7,8, 14, 15]

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


for file in displacement_tifs:
    with rasterio.open(f'{ncc_path}{file}') as src:
        
        band_dict = {name: index for index, name in zip(src.indexes, src.descriptions)}

        time_delta = get_time_delta(file)
        
        
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

        absH_mean = np.sqrt( np.square(X_vel_mean) + np.square(Y_vel_mean) ) #  Magnitude of displacement 


        # compute median velocity
        X_vel_median = ( x_median / time_delta ) * 24.0
        Y_vel_median = ( y_median / time_delta ) * 24.0
        Z_vel_median = ( z_median / time_delta ) * 24.0

        absH_median = np.sqrt( np.square(X_vel_median) + np.square(Y_vel_median) ) #  Magnitude of displacement 

        # compute ncc velocity
        X_vel_ncc = ( ncc_x / time_delta ) * 24.0
        Y_vel_ncc = ( ncc_y / time_delta ) * 24.0
        absH_ncc = np.sqrt( np.square(X_vel_ncc) + np.square(Y_vel_ncc) ) #  Magnitude of displacement 

        # make dictionary of bands I want to write 
        bands_out = {
                    'X_vel_mean': X_vel_mean,
                    'Y_vel_mean': Y_vel_mean,
                    'Z_vel_mean': Z_vel_mean,
                    'BEFORE': src.read(band_dict['BEFORE']),
                    'AFTER': src.read(band_dict['AFTER']),
                    'X_vel_median': X_vel_median,
                    'Y_vel_median': Y_vel_median,
                    'Z_vel_median': Z_vel_median,
                    'MATCH_COUNT': src.read(band_dict['MATCH_COUNT']),
                    'RMS_RESIDUAL': src.read(band_dict['RMS_RESIDUAL']),
                    'MAD_X': src.read(band_dict['MAD_X']),
                    'MAD_Y': src.read(band_dict['MAD_Y']),
                    'MAD_Z': src.read(band_dict['MAD_Z']),
                    'X_vel_ncc': X_vel_ncc,
                    'Y_vel_ncc': Y_vel_ncc,
                    'NCC_PEAK': src.read(band_dict['NCC_PEAK']),
                    'Absolute_horizontal_vel_mean': absH_mean,
                    'Absolute_horizontal_vel_median': absH_median,
                    'Absolute_horizontal_vel_ncc': absH_ncc
                    }



        profile = src.profile.copy()
        profile['count'] = len(bands_out) 

        out_file = f'{out_path}{file[:-4]}_vel.tif'
        print('Writing file:', out_file)

        for name, index in band_dict.items():
            with rasterio.open(out_file, mode='w', **profile) as dst:

                for band_num, (name, band_data) in enumerate(bands_out.items(), start=1):

                    dst.write(band_data, band_num)
                    dst.set_band_description(band_num, name)

                    



            
        
