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
# (modularcasa) pip install astropy
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

from keplerian_mask import make_keplerian_mask
from astropy.io import fits
from scipy.ndimage import gaussian_filter

def get_kep_mask_wrapper(imagename,
                        # maskname,
                        line,
                        imsize=None,
                        cellsize=None,
                        robust=0.5,
                        vres_version='v11',
                        tag='.keplerian_mask',
                        anti_tag='.anti_keplerian_mask',
                        smooth_tag='.smooth_keplerian_mask',
                        sigma_channels=5,
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
                        tag             = tag,
                        **dmask.mask_dict[line+'_keplerian'])

    # Create a second mask that has ones where the Keplerian mask has zeros, and vice versa
    hdul = fits.open(imagename+tag+'.fits')
    hdul[0].data -= 1.
    hdul[0].data = np.abs(hdul[0].data)
    hdul.writeto(imagename+anti_tag+'.fits')
    hdul.close() # Does not over-write tag.fits data

    # Create a third mask that is smoothed along the spectral axis (to get rid of streaks in moment maps)
    hdul = fits.open(imagename+tag+'.fits') # Opens original tag.fits data
    # Perform the spectral smoothing; clunky but it works
    for y_ind in np.arange(0, hdul[0].data.shape[1], 1):
        hdul[0].data[:, y_ind, :] = gaussian_filter(hdul[0].data[:, y_ind, :], sigma=sigma_channels)
    for x_ind in np.arange(0, hdul[0].data.shape[2], 1):
        hdul[0].data[:, :, x_ind] = gaussian_filter(hdul[0].data[:, :, x_ind], sigma=sigma_channels)
    hdul.writeto(imagename+smooth_tag+'.fits')
    hdul.close()



"""
######################################################
################## MASK THE LINES ####################
######################################################
"""

molecules       = ['12CO']#'13CO', 'SO', '12CO']
vres_version    = 'v11'
mask_version    = '_m8'
sigma_channels  = 5

for line in molecules:
    for robust in [0.5]:
        for cont in ['']:#, '_wcont']:
            os.system('mkdir '+ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont)
            imagename       = ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont+'/ABAur_'+line

            print("################################################")
            print("###### Creating files whose names will start with: ", imagename)
            print("###### Image cubes will have spectral resolution defined by version: ", vres_version)
            print("###### And Briggs robust weighting: ", robust)
            print("################################################")

            get_kep_mask_wrapper(imagename     = imagename,
                                line           = line,
                                robust         = robust,
                                vres_version   = vres_version,
                                tag            = '.keplerian_mask'+mask_version,
                                anti_tag       = '.anti_keplerian_mask'+mask_version,
                                smooth_tag     = '.smooth_keplerian_mask'+mask_version,
                                sigma_channels = sigma_channels
                                )
