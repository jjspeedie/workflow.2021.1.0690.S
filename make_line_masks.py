# python3 -m venv modularcasa
# source modularcasa/bin/activate
# (modularcasa) pip install --upgrade pip wheel
# (modularcasa) pip install casatasks==6.4.3.27
# (modularcasa) pip install casatools==6.4.3.27
# (modularcasa) pip install casadata==2022.9.5
# (modularcasa) pip install numba
# (modularcasa) pip install astropy
# (modularcasa) pip install shutils
# (modularcasa) pip install pandas
'''
ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

Script to experiment with making masks for cleaning of the lines.

To run this script, do:
source modularcasa/bin/activate
(modularcasa) python make_line_masks.py
'''
import os
import sys
import numpy as np
import pandas as pd
import casatools
tb = casatools.table()
import casatasks
from casatasks import impbcor
from casatasks import exportfits
import dictionary_data as ddata # contains data_dict
import dictionary_disk as ddisk # contains disk_dict
import dictionary_mask as dmask # contains mask_dict
import dictionary_lines as dlines # contains line_dict

from JvM_correction_casa6 import do_JvM_correction_and_get_epsilon
from keplerian_mask import make_mask

# """Guidelines for optimizing auto-multithresh parameters:
# Save the masks for each major cycle. To save the intermediate masks,
# type the following on the casa command line:"""
# os.environ['SAVE_ALL_AUTOMASKS']="true"


molecules       = ['SO']#['13CO', '12CO', 'C18O', 'SO']
vres_version    = 'v1' # 14-Jan-2023

for line in molecules:
    for robust in [0.5]:
        for cont in ['']:#, '_wcont']:
            os.system('mkdir '+ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont)
            vis             = ddata.data_dict[line+cont]
            robust          = robust
            imagename       = ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont+'/ABAur_'+line




            imagename +='.dirty'
            linefreq    = dlines.line_dict[line]['freq']
            print('linefreq = ', linefreq)

            if line=='SO': # Make a keplerian_mask
                print("Making Keplerian mask...")
                rms = make_mask(image           = imagename+'.image',
                                inc             = ddisk.disk_dict['incl'],
                                PA              = ddisk.disk_dict['PA_gofish'],
                                mstar           = ddisk.disk_dict['M_star'], # ideally this would be dynamical, measured with gofish/eddy
                                dist            = ddisk.disk_dict['distance'],
                                vlsr            = ddisk.disk_dict['v_sys']*1000., # needs m/s
                                restfreqs       = linefreq,
                                export_FITS     = True,
                                estimate_rms    = True, # prints and returns rms (Jy/beam) outside mask
                                **dmask.mask_dict[line])









# sys.exit()
