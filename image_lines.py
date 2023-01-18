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
from keplerian_mask import make_mask_of_linefree_channels
# from calc_uvtaper import calc_taper
# from shutil import copytree
# from shutil import rmtree

"""Guidelines for optimizing auto-multithresh parameters:
Save the masks for each major cycle. To save the intermediate masks,
type the following on the casa command line:"""
os.environ['SAVE_ALL_AUTOMASKS']="true"


def tclean_wrapper_continuum(vis,
                            imagename,
                            line,
                            imsize=None,
                            cellsize=None,
                            robust=0.5,
                            vres_version='v1',
                            uvtaper=[]):
    """
    Master tclean wrapper function to image all lines.

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
    tclean_sv_clean_extensions    = ['.residual', '.psf'] #  .image and .JvM_convolved_model will be added to this
    tclean_rm_clean_extensions    = ['.sumwt']
    JvM_extensions                = ['.image','.JvMcorr.image','.JvMcorr_lowres.image']
    JvM_pbcor_extensions          = [ext+'.pbcor' for ext in JvM_extensions]

    """ Make a dirty image for determining the cleaning threshold """
    imagename +='.dirty'

    # for ext in tclean_rm_extensions:
    #     print('Deleting previous file: ', imagename+ext)
    #     os.system('rm -rf '+ imagename+ext)

    print("Making a dirty image of "+(line)+"...")
    uvtaper = ''

    imsize      = 2048,                  # to image to FWHM of primary beam; FOV 40 arcsec (diameter)
    cell        = '0.02arcsec',          # samples bmin ~9 times
    scales      = [1, 5, 15, 30, 50,100] # x0.02arcsec = 0.02, 0.1, 0.3, 0.6, 1., 2. arcsec
    # Dec 2022: scales  = [0, 5, 15, 25] # x0.04arcsec = 0, 0.2, 0.6, 1. arcsec
             # MAPS scales: [0, 5, 15, 25] x0.02arcsec = 0., 0.1, 0.3, 0.5 arcsec
                                # MAPS Huang on GM Aur: [0, 0.4, 1′′, 2′′]

    start       = dlines.line_dict[line][vres_version+'_start']
    width       = dlines.line_dict[line][vres_version+'_width']
    nchan       = dlines.line_dict[line][vres_version+'_nchan']
    linefreq    = dlines.line_dict[line]['freq']
    print('linefreq = ', linefreq)

    # casatasks.tclean(vis                    = vis,              # msfile to image
    #                  imagename              = imagename,        # file names preceding .image, .residual, etc.
    #                  specmode               = 'cube',
    #                  restfreq               = linefreq,
    #                  start                  = start,
    #                  width                  = width,
    #                  nchan                  = nchan,
    #                  outframe               = 'lsrk',
    #                  veltype                = 'radio',
    #                  deconvolver            = 'multiscale',
    #                  scales                 = scales,
    #                  weighting              = 'briggsbwtaper',  # Ryan's suggestion over 'briggs' (used by MAPS)
    #                  robust                 = robust,
    #                  imsize                 = imsize,
    #                  cell                   = cell,
    #                  spw                    = '',
    #                  niter                  = 0,                # to make a dirty image
    #                  interactive            = False,
    #                  perchanweightdensity   = True,             # Ryan's suggestion over False (used by MAPS)
    #                  restoringbeam          = 'common',
    #                  cycleniter             = 300,
    #                  cyclefactor            = 1,
    #                  uvtaper                = uvtaper,
    #                  savemodel              = 'none'
    #                  )
    #
    # print("Exporting dirty image products of the "+line+" line to FITS...")
    # for ext in tclean_sv_dirty_extensions:
    #     exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)
    for ext in tclean_rm_dirty_extensions:
        print('(Space saving) Deleting this dirty file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)

    print("Dirty image complete, estimating rms noise in line-free channels...")
    mask_for_estimating_rms = imagename+'.mask_for_estimating_rms.image'
    if os.path.isfile(mask_for_estimating_rms):
        print('--> mask_for_estimating_rms.image already exists, and we will use it now...')
        rms = casatasks.imstat(imagename=imagename+'.image', mask='"{}" < 1.0'.format(mask_for_estimating_rms))['rms'][0]
    else:
        print('--> mask_for_estimating_rms.image does not exist yet; we are making it now...')
        rms = make_mask_of_linefree_channels(image           = imagename+'.image',
                                             inc             = ddisk.disk_dict['incl'],
                                             PA              = ddisk.disk_dict['PA_gofish'],
                                             mstar           = ddisk.disk_dict['M_star'], # ideally this would be dynamical, measured with gofish/eddy
                                             dist            = ddisk.disk_dict['distance'],
                                             vlsr            = ddisk.disk_dict['v_sys']*1000., # needs m/s
                                             v_min           = dmask.mask_dict[line+'_estimate_rms']['v_min'], # in m/s
                                             v_max           = dmask.mask_dict[line+'_estimate_rms']['v_max'], # in m/s
                                             r_max           = dmask.mask_dict[line+'_estimate_rms']['r_max'], # in arcsec
                                             restfreqs       = linefreq,
                                             export_FITS     = True,
                                             estimate_rms    = True) # prints and returns rms (Jy/beam) outside mask

    print("Estimated rms noise outside of the masked regions: %.2f mJy/beam"%(rms*1e3))
    threshold   = "%.8f" %(4.*rms*1e3)+'mJy'

    # if Keplerian mask for cleaning...
    # if line=='SO': # Make a keplerian_mask
    #     # Returns rms in Jy/beam
    #     rms = make_keplerian_mask(image           = imagename+'.image',
    #                               inc             = ddisk.disk_dict['incl'],
    #                               PA              = ddisk.disk_dict['PA_gofish'],
    #                               mstar           = ddisk.disk_dict['M_star'], # ideally this would be dynamical, measured with gofish/eddy
    #                               dist            = ddisk.disk_dict['distance'],
    #                               vlsr            = ddisk.disk_dict['v_sys']*1000., # needs m/s
    #                               restfreqs       = linefreq,
    #                               export_FITS     = True,
    #                               estimate_rms    = True, # prints and returns rms (Jy/beam) outside mask
    #                               **dmask.mask_dict[line])
    # else:
    #     print("Code for 12CO, 13CO, C18O not yet written")



    """ Clean down to the cleaning threshold """
    imagename = imagename.replace('.dirty', '.clean')

    for ext in tclean_rm_extensions:
        print('Deleting previous file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)

    print("Starting to clean the "+line+" line down to threshold of "+threshold+"...")
    rec = casatasks.tclean(vis                    = vis,              # msfile to image
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
                           niter                  = 10000,            # to reach the threshold
                           interactive            = False,
                           perchanweightdensity   = True,             # Ryan's suggestion over False (used by MAPS)
                           restoringbeam          = 'common',
                           cycleniter             = 300,
                           cyclefactor            = 1,
                           uvtaper                = uvtaper,
                           savemodel              = 'none',

                           # Keplerian mask:
                           # mask                   = imagename.replace('.clean', '.dirty.mask.image',))

                           # Use broad initial mask (the one used for noise estimation):
                           mask                   = mask_for_estimating_rms,
                           # usemask                = 'user',
                           restart                = True,

                           # Automasking Parameters below this line
                           usemask           = 'auto-multithresh',
                           sidelobethreshold = 2.0,    # changed; Table of Standard values: 12m (long) b75>300m = 3.0
                           noisethreshold    = 4.0,    # changed; Table of Standard values: 12m (long) b75>300m = 5.0
                           lownoisethreshold = 1.5,    #          Table of Standard values: 12m (long) b75>300m = 1.5
                           minbeamfrac       = 0.3,    #          Table of Standard values: 12m (long) b75>300m = 0.3
                           growiterations    = 75,     #          controls the maximum number of iterations that binary dilation performs. A value between 75 and 100 is usually adequate.
                           negativethreshold = 7.0,    #          Table of Standard values: 12m (long) b75>300m = 7.0
                           verbose           = True)
    print("Saving summary log file of tcleaning process...")
    np.save(imagename+'.tclean.summary.npy', rec)

    print("Exporting clean image products of the "+line+" line to FITS...")
    for ext in tclean_sv_clean_extensions:
        exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)
    for ext in tclean_rm_clean_extensions:
        print('(Space saving) Deleting this clean file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)

    # print("Creating imview raster of .clean.image with .clean.mask overlaid...")
    # frame_dir = imagename+'.check.masking/'
    # os.system('mkdir '+frame_dir)
    # for chan_idx in range(nchan):
    #     imview(raster   = {'file': imagename+'.image'},
    #            contour  = {'file': imagename+'.mask'},
    #            zoom     = {'channel': chan_idx},
    #            out      = frame_dir+'%03i.png'%(chan_idx))


    """ Do JvM correction (primary beam correction done concurrently) """

    print('Starting JvM correction of the '+line+' line...')
    JvM_epsilon = do_JvM_correction_and_get_epsilon(root=imagename)
    # Expect: WARN	ImageExprCalculator::compute	image units are not the same: 'Jy/beam' vs ''. Proceed with caution. Output image metadata will be copied from .clean_convolved_model_temp.image

    print("Primary beam correcting the "+line+" line...")
    for ext,pbcor_ext in zip(JvM_extensions,JvM_pbcor_extensions):
        impbcor(imagename=imagename+ext, pbimage=imagename+'.pb',
                outfile=imagename+pbcor_ext, overwrite=True)

    print("Finished primary beam correction, exporting FITS of the "+line+" line...")
    for ext in JvM_extensions + JvM_pbcor_extensions:# + ['.mask']:
        exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)

    print("Saving imaging metrics and image metrics...")
    imaging_metrics = {'.dirty rms (mJy/beam)' : rms*1e3,
                       'tclean threshold' : threshold,
                       'JvM epsilon' : JvM_epsilon,
                       }
    imaging_info = pd.DataFrame(data=imaging_metrics, index=[0])
    imaging_info.to_csv(imagename.replace('.clean', '.imaging_info')+'.csv')
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

    print("Done processing the "+line+" line! (for robust %.1f)"%(robust))




"""
######################################################
################## IMAGE THE LINES ###################
######################################################
"""

molecules       = ['13CO']#'SO', '13CO', '12CO', 'C18O']#, 'SO']
vres_version    = 'v1' # 14-Jan-2023

for line in molecules:
    for robust in [1.]:
        for cont in ['']:#, '_wcont']:
            os.system('mkdir '+ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont)
            vis             = ddata.data_dict[line+cont]
            robust          = robust
            imagename       = ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont+'/ABAur_'+line

            print("################################################")
            print("###### About to start imaging measurement set: ", vis)
            print("###### Creating files whose names will start with: ", imagename)
            print("###### Image cubes will have spectral resolution defined by version: ", vres_version)
            print("###### And Briggs robust weighting: ", robust)
            print("################################################")

            tclean_wrapper_continuum(vis            = vis,
                                     imagename      = imagename,
                                     line           = line,
                                     robust         = robust,
                                     vres_version   = vres_version
                                 )

"""
######################################################
#### Adjustments to this script still to be made #####
######################################################

- adjust rms noise estimation to use first 10 channels instead; mask input is not working
- try to get tclean to start up with user specified mask

- change spectral lims to be more snug
    - change mask_for_estimating_rms to be more snug?

- Save imview frames of images made overlaying automultithresh mask
    - On top of dirty image, residuals
- Confirm broad mask is used by automultithresh as initial mask
- Save image_metrics csv
- Flexibility for different robust parameters and uv tapers
    - Make cell size dependent on those choices
- Flexibility for continuum-subtracted or not
    - Change the masks accordingly

DONE - Save the JvM convolved model?
DONE - Finess/decide on keplerian mask parameters for SO
DONE - Figure out good masks for 12CO, 13CO, C18O, save them
CAN'T - Save the .last files for posterity [.LAST FILES AREN'T CREATED WITH MODULAR CASA]
"""

"""
Spectral extent of line emission, measured in (shallowly) cleaned robust=1.5 image cubes:
12CO: 0.222 km/s --> 12.066 km/s
13CO: 2.919 km/s --> 9.000 km/s
C18O: 3.498 km/s --> 8.538 km/s
SO:   4.422 km/s --> 7.632 km/s
"""


# Notes:

# PSF is unitless
# model is Jy/pixel
# residual is unitless
# pb is unitless

    # 12CO emission starts at ~11.5 km/s or even larger; faint but there. ends at ~0.7 km/s, faint but also there farther out.
    # 13CO emission starts at 9.1 km/s; ends at about 3 km/s.
    # C18 emission starts at 8.1 km/s; ends at 3.7 km/s
    # SO emission starts around 7.5 km/s; ends at ~ 4.17 km/s

# thresholds = ['3.43mJy', '6.69mJy', '1.64mJy', '1.68mJy'] # from TM1 pipeline weblog


# MAPS Huang on GM Aur:
# multiscale scales: [0, 0.4, 1′′, 2′′]
# The mask was initialized with full coverage of the primary beam from 5.2–6.4 km s−1, where the emission is the broadest, because auto-multithresh algorithm does not readily mask diffuse emission.
# After some experimentation, the following auto-multithresh parameters were selected:
# sidelobethreshold = 3.0,
# noisethreshold = 4.0,
# lownoisethreshold = 1.5,
# and minbeamfrac = 0.3.
# The CLEAN threshold was set to 5 mJy, corresponding to ∼3× the rms of line-free channels in the dirty image.

# huang 2022 on DR Tau
# sidelobethreshold=2.0,
# noisethreshold=4.0,
# minbeamfrac=0.3,
# negativethreshold=7.0.



# sys.exit()
