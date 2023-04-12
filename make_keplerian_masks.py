# Run this script in a virtual environment, made with the following:
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

Quick script to make Keplerian masks for the data.

To run this script, do:
source modularcasa/bin/activate
(modularcasa) python make_keplerian_masks.py
'''
import os
import sys
import numpy as np
import pandas as pd
import casatools
tb = casatools.table()
ia = casatools.image()
import casatasks
from casatasks import impbcor
from casatasks import exportfits
import dictionary_data as ddata # contains data_dict
import dictionary_disk as ddisk # contains disk_dict
import dictionary_mask as dmask # contains mask_dict
import dictionary_lines as dlines # contains line_dict

from JvM_correction_casa6 import do_JvM_correction_and_get_epsilon
from keplerian_mask import make_keplerian_mask

def get_kep_mask_wrapper(vis,
                        imagename,
                        # maskname,
                        line,
                        imsize=None,
                        cellsize=None,
                        robust=0.5,
                        vres_version='v11',
                        uvtaper=[]):
    """
    """

    imagename +='.clean'
    linefreq   = dlines.line_dict[line]['freq']

    make_keplerian_mask(image           = imagename+'.image',
                        inc             = ddisk.disk_dict['incl'],
                        PA              = ddisk.disk_dict['PA_gofish'],
                        mstar           = ddisk.disk_dict['M_star'], # ideally this would be dynamical, measured with gofish/eddy
                        dist            = ddisk.disk_dict['distance'],
                        vlsr            = ddisk.disk_dict['v_sys']*1000., # needs m/s
                        restfreqs       = linefreq,
                        export_FITS     = True,
                        estimate_rms    = False,
                        **dmask.mask_dict[line+'_keplerian'])



"""
######################################################
################## MASK THE LINES ####################
######################################################
"""

molecules       = ['12CO']#'13CO', 'SO', '12CO']
vres_version    = '' # REMOVED _ FROM ROBUST PREFIX AND ADDED 'DECEMBER' - PUT IT BACK

for line in molecules:
    for robust in [1.5]:
        for cont in ['']:#, '_wcont']:
            os.system('mkdir '+ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'robust'+str(robust)+'_december'+cont)
            vis             = ddata.data_dict[line+cont]
            robust          = robust
            imagename       = ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'robust'+str(robust)+'_december'+cont+'/ABAur_'+line

            print("################################################")
            print("###### About to start imaging measurement set: ", vis)
            print("###### Creating files whose names will start with: ", imagename)
            print("###### Image cubes will have spectral resolution defined by version: ", vres_version)
            print("###### And Briggs robust weighting: ", robust)
            print("################################################")

            get_kep_mask_wrapper(vis           = vis,
                                imagename      = imagename,
                                line           = line,
                                robust         = robust,
                                vres_version   = vres_version,
                                # maskname       = maskname
                                )
