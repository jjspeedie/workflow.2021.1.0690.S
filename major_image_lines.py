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
This script began as a copy of image_lines.py, branched 8-Mar-2023.
Here we adopt a "cautious" cleaning approach, wherein we force frequent
major cycles, and clean with a broad mask.

To run this script, do:
source modularcasa/bin/activate
(modularcasa) python major_image_lines.py
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
# from keplerian_mask import make_mask_for_diffuse_emission
# from keplerian_mask import make_mask_from_model
# from calc_uvtaper import calc_taper
from shutil import copytree
# from shutil import rmtree

# """Guidelines for optimizing auto-multithresh parameters:
# Save the masks for each major cycle. To save the intermediate masks,
# type the following on the casa command line:"""
# os.environ['SAVE_ALL_AUTOMASKS']="true"
#

def tclean_wrapper_line(vis,
                        imagename,
                        # maskname,
                        line,
                        imsize=None,
                        cellsize=None,
                        robust=0.5,
                        vres_version='v8',
                        uvtaper=[]):
    """
    Master tclean wrapper function to image a line.

    Breakdown of the procedure:
        - Make a dirty image cube
        - Estimate the rms noise in the dirty image cube in line free channels
            (makes a mask to do this if one does not already exist)
        - Clean cautiously down to a threshold of 5x the estimated rms noise,
            using frequent major cycles and a broad mask
        - Performs JvM correction and primary beam correction
        - Saves a csv of metrics used for imaging and attempts to save a
            tclean summary dictionary (CASA 6.4 bug?)
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
    scales      = [1, 5, 15, 30, 50] # x0.02arcsec = 0.02, 0.1, 0.3, 0.6, 1. arcsec
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

    print("Exporting dirty image products of the "+line+" line to FITS...")
    for ext in tclean_sv_dirty_extensions:
        exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)
    for ext in tclean_rm_dirty_extensions:
        print('(Space saving) Deleting this dirty file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)


    print("Dirty image complete, estimating rms noise in the dirty image...")
    rms = casatasks.imstat(imagename=imagename+'.image', chans='0~9')['rms'][0] # 'start' params give buffer of 20 channels before emission begins
    print("Estimated rms noise in the full FOV of the first 10 channels: %.2f mJy/beam"%(rms*1e3))
    threshold   = "%.8f" %(5.*rms*1e3)+'mJy'


    # if ((line=='12CO') | (line=='13CO') | (line=='C18O')):
    #     print("The "+line+" line contains diffuse emission at central channels that auto-multithresh struggles to mask. We shall help it...")
    #
    #     print("Have you already made a mask for kickstarting auto-multithresh?")
    #     initial_mask = maskname
    #     if os.path.exists(initial_mask):
    #         print('--> Yes, '+maskname+' already exists, and we will use it now...')
    #     else:
    #         print('--> No, '+maskname+' does not exist yet; so we have to exit.')
    #         sys.exit()
    #         # make_mask_for_diffuse_emission(image           = imagename+'.image',
    #         #                                inc             = ddisk.disk_dict['incl'],
    #         #                                PA              = ddisk.disk_dict['PA_gofish'],
    #         #                                mstar           = ddisk.disk_dict['M_star'], # ideally this would be dynamical, measured with gofish/eddy
    #         #                                dist            = ddisk.disk_dict['distance'],
    #         #                                vlsr            = ddisk.disk_dict['v_sys']*1000., # needs m/s
    #         #                                v_min           = dmask.mask_dict[line+'_diffuse_emission']['v_min'], # in m/s
    #         #                                v_max           = dmask.mask_dict[line+'_diffuse_emission']['v_max'], # in m/s
    #         #                                r_max           = dmask.mask_dict[line+'_diffuse_emission']['r_max'], # in arcsec
    #         #                                restfreqs       = linefreq,
    #         #                                export_FITS     = True, # this could be false
    #         #                                estimate_rms    = False) # this won't be an accurate estimate because there is emission outside the diffuse mask
    #
    # elif line=='SO': # make a keplerian_mask
    #     print("The "+line+" line shall be initialized with a Keplerian mask...")
    #
    #     print("Have you already made a mask for kickstarting auto-multithresh with a keplerian mask?")
    #     initial_mask = imagename+'.initial_mask_keplerian.image'
    #     if os.path.exists(initial_mask):
    #         print('--> Yes, initial_mask_keplerian.image already exists, and we will use it now...')
    #     else:
    #         print('--> No, initial_mask_keplerian.image does not exist yet; we are making it now...')
    #         make_keplerian_mask(image           = imagename+'.image',
    #                             inc             = ddisk.disk_dict['incl'],
    #                             PA              = ddisk.disk_dict['PA_gofish'],
    #                             mstar           = ddisk.disk_dict['M_star'], # ideally this would be dynamical, measured with gofish/eddy
    #                             dist            = ddisk.disk_dict['distance'],
    #                             vlsr            = ddisk.disk_dict['v_sys']*1000., # needs m/s
    #                             restfreqs       = linefreq,
    #                             export_FITS     = True,
    #                             estimate_rms    = True, # prints and returns rms (Jy/beam) outside mask, nice to compare with our estimate
    #                             **dmask.mask_dict[line+'_keplerian'])

    """ Clean down to the cleaning threshold """
    imagename = imagename.replace('.dirty', '.clean')

    for ext in tclean_rm_extensions:
        print('Deleting previous file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)


    # print("Doing initial niter=1 clean to kickstart the initial mask...")
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
    #                  threshold              = threshold,
    #                  perchanweightdensity   = True,             # Ryan's suggestion over False (used by MAPS)
    #                  restoringbeam          = 'common',
    #                  cyclefactor            = 1,
    #                  uvtaper                = uvtaper,
    #                  savemodel              = 'none',
    #
    #                  niter                  = 1,            # these two need to be the same, and we only mean to kickstart
    #                  cycleniter             = 1,            # these two need to be the same, and we only mean to kickstart
    #                  interactive            = False,
    #                  mask                   = initial_mask,
    #                  usemask                = 'user')

    # if os.path.exists(imagename+'.autothresh1'):
    #     print("Making a copy of .clean.autothresh1 for posterity and deleting original so tclean doesn't throw an error...")
    #     copytree(imagename+'.autothresh1', imagename+'.autothresh1_initial')
    #     os.system('rm -rf '+imagename+'.autothresh1')


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
                           niter                  = 1000000,            # to reach the threshold
                           interactive            = False,
                           perchanweightdensity   = True,             # Ryan's suggestion over False (used by MAPS)
                           restoringbeam          = 'common',
                           uvtaper                = uvtaper,
                           savemodel              = 'none',           # we don't wish to alter the ms's model column

                           usemask                = 'pb',             # use a broad mask
                           pbmask                 = 0.2,              # use a broad mask

                           # Cautious/conservative clean:
                           cycleniter             = 20,               # Jess: previously 300, then 100; Maximum number of minor-cycle iterations (per plane) before triggering a major cycle
                           cyclefactor            = 3.0,              # Ryan: 3x max_psf_sidelobe_level as minor cycle threshold (default is 1.0)
                           gain                   = 0.02,             # Ryan: assign clean component peaks to 2% of pixel value (default is 0.1)
                           minpsffraction         = 0.5)              # PHANGS: cycle threshold is never lower than 0.5 times the peak residual (default: 0.05)

                           # fullsummary            = True)           # attempt to access the summary dictionary


    print("Saving summary log file of tcleaning process...")
    np.save(imagename+'.tclean.summary.npy', rec)
    print(rec)

    print("Exporting clean image products of the "+line+" line to FITS...")
    for ext in tclean_sv_clean_extensions:
        exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)
    for ext in tclean_rm_clean_extensions:
        print('(Space saving) Deleting this clean file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)


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

molecules       = ['13CO']# note v7 not implemented in dict_lines
vres_version    = 'v10' # 16-Mar-2023

for line in molecules:
    for robust in [0.5]:
        for cont in ['']:#, '_wcont']:
            os.system('mkdir '+ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont)
            vis             = ddata.data_dict[line+cont]
            robust          = robust
            imagename       = ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont+'/ABAur_'+line
            # maskname        = ddata.data_dict['NRAO_path']+'images_lines/'+line+'/masks/'+'v2'+'_robust1.5_10mean_model.mask'

            print("################################################")
            print("###### About to start imaging measurement set: ", vis)
            print("###### Creating files whose names will start with: ", imagename)
            print("###### Image cubes will have spectral resolution defined by version: ", vres_version)
            print("###### And Briggs robust weighting: ", robust)
            # print("###### Auto-multithresh will be kickstarted with mask: ", maskname)
            print("################################################")

            tclean_wrapper_line(vis            = vis,
                                imagename      = imagename,
                                line           = line,
                                robust         = robust,
                                vres_version   = vres_version,
                                # maskname       = maskname
                                )


"""
######################################################
#### Adjustments to this script still to be made #####
######################################################

- Consider using smallscalebias to bias towards smaller scales?

- Save image_metrics csv
- Flexibility for different robust parameters and uv tapers
    - Make cell size dependent on those choices
- Flexibility for continuum-subtracted or not
    - Change the masks accordingly

DONE - adjust rms noise estimation to use first 10 channels instead; mask input is not working
DONE - try to get tclean to start up with user specified mask
DONE - Save the JvM convolved model?
DONE - Finess/decide on keplerian mask parameters for SO
DONE - Figure out good masks for 12CO, 13CO, C18O, save them
CAN'T - Save the .last files for posterity [.LAST FILES AREN'T CREATED WITH MODULAR CASA]
CAN'T - Save imview frames of images made overlaying automultithresh mask
CAN'T - ADD fullsummary=True (!!) [tclean throws error, "unexpected keyword argument"]
"""



"""
Spectral extent of line emission, measured in (shallowly) cleaned robust=1.5 image cubes:
12CO: 0.222 km/s --> 12.066 km/s
13CO: 2.919 km/s --> 9.000 km/s
C18O: 3.498 km/s --> 8.538 km/s
SO:   4.422 km/s --> 7.632 km/s
"""

# TRACEBACK OF FULLSUMMARY=TRUE THROWING ERROR (not printed to casalog: casa-20230309-020103.log)
# Traceback (most recent call last):
#   File "major_image_lines.py", line 373, in <module>
#     vres_version   = vres_version,
#   File "major_image_lines.py", line 267, in tclean_wrapper_line
#     fullsummary            = True)             # attempt to access the summary dictionary
# TypeError: __call__() got an unexpected keyword argument 'fullsummary'



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
# multiscale scales: [0, 0.4, 1, 2] arcsec
# The mask was initialized with full coverage of the primary beam from 5.2-6.4 km s^1, where the emission is the broadest, because auto-multithresh algorithm does not readily mask diffuse emission.
# After some experimentation, the following auto-multithresh parameters were selected:
# sidelobethreshold = 3.0,
# noisethreshold = 4.0,
# lownoisethreshold = 1.5,
# and minbeamfrac = 0.3.
# The CLEAN threshold was set to 5 mJy, corresponding to ~3x the rms of line-free channels in the dirty image.

# huang 2022 on DR Tau
# sidelobethreshold=2.0,
# noisethreshold=4.0,
# minbeamfrac=0.3,
# negativethreshold=7.0.



# sys.exit()
