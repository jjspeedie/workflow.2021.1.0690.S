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

Script to image the lines with tclean.

To run this script, do:
source modularcasa/bin/activate
(modularcasa) python image_lines.py
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
from keplerian_mask import make_mask as make_keplerian_mask
from keplerian_mask import make_mask_for_diffuse_emission
from keplerian_mask import make_mask_from_model
# from calc_uvtaper import calc_taper
from shutil import copytree
# from shutil import rmtree

"""Guidelines for optimizing auto-multithresh parameters:
Save the masks for each major cycle. To save the intermediate masks,
type the following on the casa command line:"""
os.environ['SAVE_ALL_AUTOMASKS']="true"


def tclean_wrapper_line(vis,
                        imagename,
                        line,
                        imsize=None,
                        cellsize=None,
                        robust=0.5,
                        vres_version='v2',
                        uvtaper=[]):
    """
    Master tclean wrapper function to image a line.

    Breakdown of the procedure:
        - Make a dirty image cube
        - Estimate the rms noise in the dirty image cube in line free channels
            (makes a mask to do this if one does not already exist)
        - Cleans down to a threshold of 4x the estimated rms noise, using auto-multithresh
            (saves the intermediate masks, saves a summary log file of tclean's iterations)
        - Performs JvM correction and primary beam correction
        - Saves a csv of metrics used for imaging, saves a csv of metrics of the resulting images
    """

    tclean_rm_extensions          = ['.image', '.mask', '.model', '.pb', '.psf', '.residual', '.sumwt']
    tclean_sv_dirty_extensions    = ['.image', '.psf']
    tclean_rm_dirty_extensions    = ['.residual', '.model', '.pb', '.sumwt']
    tclean_sv_clean_extensions    = ['.residual', '.psf', '.model'] #  .image and .JvM_convolved_model will be added to this
    tclean_rm_clean_extensions    = ['.sumwt']
    JvM_extensions                = ['.image','.JvMcorr.image','.JvMcorr_lowres.image']
    JvM_pbcor_extensions          = [ext+'.pbcor' for ext in JvM_extensions]

    """ Make a dirty image for determining the cleaning threshold """
    imagename +='.dirty'

    for ext in tclean_rm_extensions:
        print('Deleting previous file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)

    print("Making a dirty image of "+(line)+"...")
    uvtaper = ''

    imsize      = 2048,                  # to image to FWHM of primary beam; FOV 40 arcsec (diameter)
    cell        = '0.02arcsec',          # samples bmin ~9 times
    scales      = [5, 15, 30, 50,100] # x0.02arcsec = 0.1, 0.3, 0.6, 1., 2. arcsec
    # Dec 2022: scales  = [0, 5, 15, 25] # x0.04arcsec = 0, 0.2, 0.6, 1. arcsec
            # MAPS scales: [0, 5, 15, 25] x0.02arcsec = 0., 0.1, 0.3, 0.5 arcsec
                    # MAPS Huang on GM Aur: [0, 0.4, 1, 2] arcsec

    start       = dlines.line_dict[line][vres_version+'_start']
    width       = dlines.line_dict[line][vres_version+'_width']
    nchan       = dlines.line_dict[line][vres_version+'_nchan']
    linefreq    = dlines.line_dict[line]['freq']
    print('linefreq = ', linefreq)

    casatasks.tclean(vis                    = vis,              # msfile to image
                     imagename              = imagename,        # file names preceding .image, .residual, etc.
                     specmode               = 'cube',
                     restfreq               = linefreq,
                     start                  = start,
                     width                  = width,
                     nchan                  = nchan,
                     outframe               = 'lsrk',
                     veltype                = 'radio',
                     deconvolver            = 'multiscale',
                     scales                 = scales,
                     weighting              = 'briggsbwtaper',  # Ryan's suggestion over 'briggs' (used by MAPS)
                     robust                 = robust,
                     imsize                 = imsize,
                     cell                   = cell,
                     spw                    = '',
                     niter                  = 0,                # to make a dirty image
                     interactive            = False,
                     perchanweightdensity   = True,             # Ryan's suggestion over False (used by MAPS)
                     restoringbeam          = 'common',
                     cycleniter             = 300,
                     cyclefactor            = 1,
                     uvtaper                = uvtaper,
                     savemodel              = 'none'
                     )

    # print("Exporting dirty image products of the "+line+" line to FITS...")
    # for ext in tclean_sv_dirty_extensions:
    #     exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)
    for ext in tclean_rm_dirty_extensions:
        print('(Space saving) Deleting this dirty file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)


    print("Dirty image complete, estimating rms noise in the dirty image...")
    rms = casatasks.imstat(imagename=imagename+'.image', chans='0~9')['rms'][0]
    print("Estimated rms noise in the full FOV of the first 10 channels: %.2f mJy/beam"%(rms*1e3))
    threshold   = "%.8f" %(10.*rms*1e3)+'mJy'

    """ Clean down to the cleaning threshold """
    imagename = imagename.replace('.dirty', '.initial')

    for ext in tclean_rm_extensions:
        print('Deleting previous file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)


    print("Doing tclean to generate the model for the initial mask (using pb mask)...")
    casatasks.tclean(vis                    = vis,              # msfile to image
                     imagename              = imagename,        # file names preceding .image, .residual, etc.
                     specmode               = 'cube',
                     restfreq               = linefreq,
                     start                  = start,
                     width                  = width,
                     nchan                  = nchan,
                     outframe               = 'lsrk',
                     veltype                = 'radio',
                     deconvolver            = 'multiscale',
                     scales                 = scales,
                     weighting              = 'briggsbwtaper',  # Ryan's suggestion over 'briggs' (used by MAPS)
                     robust                 = robust,
                     imsize                 = imsize,
                     cell                   = cell,
                     spw                    = '',
                     threshold              = threshold,
                     perchanweightdensity   = True,             # Ryan's suggestion over False (used by MAPS)
                     restoringbeam          = 'common',
                     cyclefactor            = 1,
                     uvtaper                = uvtaper,
                     savemodel              = 'none',

                     niter                  = 100000,            # to reach the threshold
                     cycleniter             = 300,
                     interactive            = False,
                     usemask                = 'pb',
                     pbmask                 = 0.2)


    # print("Exporting clean image products of the "+line+" line to FITS...")
    # for ext in tclean_sv_clean_extensions:
    #     exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)
    for ext in tclean_rm_clean_extensions:
        print('(Space saving) Deleting this clean file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)


    print("Saving imaging metrics and image metrics...")
    imaging_metrics = {'.dirty rms (mJy/beam)' : rms*1e3,
                       'tclean threshold' : threshold,
                       # 'JvM epsilon' : JvM_epsilon,
                       }
    imaging_info = pd.DataFrame(data=imaging_metrics, index=[0])
    imaging_info.to_csv(imagename.replace('.initial', '.imaging_info')+'.csv')
    print("Saved imaging info csv.")

    # rows = []
    # for ext in JvM_extensions + JvM_pbcor_extensions:
    #     image_metrics = {}
    #     image_metrics['peak intensity (mJy/beam)'] = estimate_peak_intensity(imagename=imagename+ext, mask=mask)
    #     image_metrics['disk flux (mJy/beam)'] = estimate_disk_flux(imagename=imagename+ext, mask=mask)
    #     image_metrics['rms noise (uJy/beam)'] = estimate_rms(imagename=imagename+ext, region=region)
    #     image_metrics['SNR'] = estimate_SNR(imagename=imagename+ext, mask=mask, region=region)
    #     rows.append(pd.Series(image_metrics, name=ext))
    #
    # image_info = pd.concat(rows, axis=1)
    # image_info.transpose().to_csv(imagename.replace('.clean', '.image_info')+'.csv')
    # print("Saved image info csv.")

    print("***** Now for the important part --> make_mask_from_model...")

    level = 10. * casatasks.imstat(imagename+'.model')['mean'][0]
    print('10*mean = level = ', level)

    make_mask_from_model(imagename+'.model', level)

    print("Done processing the "+line+" line! (for robust %.1f)"%(robust))




"""
######################################################
################## IMAGE THE LINES ###################
######################################################
"""

molecules       = ['C18O', '13CO', '12CO', 'SO']#'C18O', 'SO']#'SO', '13CO', '12CO', 'C18O']#, 'SO']
vres_version    = 'v2' # 17-Jan-2023

for line in molecules:
    for robust in [1.5]:
        for cont in ['']:#, '_wcont']:
            os.system('mkdir '+ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_initial_mask'+cont)
            vis             = ddata.data_dict[line+cont]
            robust          = robust
            imagename       = ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_initial_mask'+cont+'/ABAur_'+line

            print("################################################")
            print("###### About to start imaging measurement set: ", vis)
            print("###### Creating files whose names will start with: ", imagename)
            print("###### Image cubes will have spectral resolution defined by version: ", vres_version)
            print("###### And Briggs robust weighting: ", robust)
            print("################################################")

            tclean_wrapper_line(vis            = vis,
                                imagename      = imagename,
                                line           = line,
                                robust         = robust,
                                vres_version   = vres_version
                                )





            # sys.exit()




# sys.exit()
