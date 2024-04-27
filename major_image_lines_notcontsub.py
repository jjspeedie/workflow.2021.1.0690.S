# Run this script in a virtual environment, made with the following:
# python3.8 -m venv modularcasa24 (2024 version due to following error)
            ###   File "/usr/lib64/python3.6/importlib/__init__.py", line 126, in import_module
            ###     return _bootstrap._gcd_import(name[level:], package, level)
            ### ModuleNotFoundError: No module named '_synthesismaskhandler'
# source modularcasa24/bin/activate
# (modularcasa24) pip install --upgrade pip wheel
# (modularcasa24) pip install casatasks==6.6.3.22
# (modularcasa24) pip install casatools==6.6.3.22
# (modularcasa24) pip install casadata==2024.1.15
# (modularcasa24) pip install numba
# (modularcasa24) pip install astropy
# (modularcasa24) pip install shutils
# (modularcasa24) pip install pandas
'''
ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

Script to image the lines with tclean.
This script began as a copy of image_lines.py, branched 8-Mar-2023.
Here we adopt a "cautious" cleaning approach, wherein we force frequent
major cycles, and clean with a broad mask.

To run this script, do:
source modularcasa/bin/activate
(modularcasa24) python major_image_lines_notcontsub.py
'''
import os
import sys
import numpy as np
import pandas as pd
# import casatools
# tb = casatools.table()
# ia = casatools.image()
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



def tclean_wrapper_line(vis,
                        imagename,
                        # maskname,
                        line,
                        imsize=None,
                        cellsize=None,
                        robust=0.5,
                        vres_version='v11',
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
    if line=='12CO':
        scales.append(100) # cleaning takes too long for 12CO without the 2. arcsec scale

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
                           niter                  = 1000000,            # to reach the threshold
                           interactive            = False,
                           perchanweightdensity   = True,             # Ryan's suggestion over False (used by MAPS)
                           restoringbeam          = 'common',
                           uvtaper                = uvtaper,
                           savemodel              = 'none',           # we don't wish to alter the ms's model column

                           usemask                = 'pb',             # use a broad mask
                           pbmask                 = 0.2,              # use a broad mask

                           # Cautious/conservative clean:
                           cycleniter             = 80,               # Jess: previously 300, then 100, then 20; Maximum number of minor-cycle iterations (per plane) before triggering a major cycle
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

molecules       = ['SO']
vres_version    = 'v11' # 26-Apr-2024

for line in molecules:
    for robust in [0.5, 1.5]:
        for cont in ['_wcont']: #['']:#, '_wcont']:
            os.system('mkdir '+ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+cont+'_robust'+str(robust))
            vis             = ddata.data_dict[line+cont]
            robust          = robust
            imagename       = ddata.data_dict['NRAO_path']+'images_lines/'+line+'/'+vres_version+cont+'_robust'+str(robust)+'/ABAur_'+line+cont
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



# sys.exit()
