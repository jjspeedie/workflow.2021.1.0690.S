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

Script to image the continuum with tclean.

To run this script, do:
source modularcasa/bin/activate
(modularcasa) python image_continuum.py
'''
import os
import sys
import numpy as np
import pandas as pd
import pickle
import casatools
tb = casatools.table()
import casatasks
from casatasks import impbcor
from casatasks import exportfits
import dictionary_data as ddata # contains data_dict
import dictionary_mask as dmask # contains mask_dict
from JvM_correction_casa6 import do_JvM_correction_and_get_epsilon
# from calc_uvtaper import calc_taper

def estimate_rms(imagename, region=''):
    """ Use CASA's imstat task to estimate the rms noise outside the given region (annulus). """
    rms = casatasks.imstat(imagename = imagename, region = region)['rms'][0]*1e6
    print("# rms: %.2f uJy/beam" %(rms))
    return rms

def estimate_disk_flux(imagename, mask=''):
    """ Use CASA's imstat task to estimate the total flux inside the given region (mask). """
    disk_stats  = casatasks.imstat(imagename = imagename, region = mask)
    disk_flux   = disk_stats['flux'][0]*1e3
    return disk_flux

def estimate_peak_intensity(imagename, mask=''):
    """ Use CASA's imstat task to estimate the peak intensity inside the given region (mask). """
    disk_stats      = casatasks.imstat(imagename = imagename, region = mask)
    peak_intensity  = disk_stats['max'][0]*1e3
    return peak_intensity

def estimate_SNR(imagename, mask='', region=''):
    """ Use CASA's imstat task to estimate the SNR as peak intensity inside the
    given region (mask) divided by the rms noise outside the given region (annulus). """
    peak_intensity  = estimate_peak_intensity(imagename=imagename, mask=mask)
    rms             = estimate_rms(imagename=imagename, region=region)
    SNR             = peak_intensity/(rms/1e3) # rms in [uJy/beam], peak in [mJy/beam]
    print("SNR of the image: %.2f"%(SNR))
    return SNR

def tclean_wrapper_continuum(vis, imagename, mask='', region='', imsize=None,
                             cellsize=None, robust=0.5, uvtaper=[]):
    """
    Wrapper for tclean with the parameters we desire for continuum imaging.

    Breakdown of the procedure:
        - Make a dirty image; use that to estimate rms noise for cleaning threshold
        - Clean down to the threshold
        - Do JvM correction, and pbcorrection of JvM images simultaneously
        - Saves csv of key metrics used during the cleaning process
        - Saves csv of key metrics of the created images
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

    imsize      = 2048,                  # to image to FWHM of primary beam; FOV 40 arcsec (diameter)
    cell        = '0.02arcsec',          # samples bmin ~9 times

    print("Making a dirty image of the continuum...")
    casatasks.tclean(vis              = vis, # msfile to image
                     imagename        = imagename, # file names preceding .image, .residual, etc.
                     specmode         = 'mfs', # to make a continuum image
                     deconvolver      = 'hogbom', # better than multiscale, for continuum rings
                     weighting        = 'briggs',
                     robust           = robust,
                     imsize           = imsize,
                     cell             = cell,
                     mask             = mask,
                     spw              = '',
                     niter            = 0, # to make a dirty image
                     interactive      = False,
                     cycleniter       = 300,
                     cyclefactor      = 1,
                     smallscalebias   = 0.6, # CASA's default
                     nterms           = 1, # Number of Taylor coefficients in the spectral model; nterms=1 : Assume flat spectrum source
                     uvtaper          = uvtaper,
                     savemodel        = 'none')

    print("Exporting dirty image products of the continuum to FITS...")
    for ext in tclean_sv_dirty_extensions:
        exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)
    for ext in tclean_rm_dirty_extensions:
        print('(Space saving) Deleting this dirty file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)

    rms         = estimate_rms(imagename=imagename+'.image', region=region)
    threshold   = "%.8f" %(2.*rms/1e3)+'mJy'

    """ Clean down to the cleaning threshold """
    imagename = imagename.replace('.dirty', '.clean')

    for ext in tclean_rm_extensions:
        print('Deleting previous file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)

    print("Starting to clean the continuum down to threshold of "+threshold+"...")
    rec =   casatasks.tclean(vis              = vis, # msfile to image
                             imagename        = imagename, # file names preceding .image, .residual, etc.
                             specmode         = 'mfs', # to make a continuum image
                             deconvolver      = 'hogbom', # better than multiscale, for continuum rings
                             weighting        = 'briggs',
                             robust           = robust,
                             imsize           = imsize,
                             cell             = cell,
                             mask             = mask,
                             spw              = '',
                             niter            = 100000,
                             threshold        = threshold,
                             interactive      = False,
                             cycleniter       = 300,
                             cyclefactor      = 1,
                             smallscalebias   = 0.6, # CASA's default
                             gain             = 0.2, # slightly low
                             nterms           = 1, # Number of Taylor coefficients in the spectral model; nterms=1 : Assume flat spectrum source
                             uvtaper          = uvtaper,
                             savemodel        = 'none')#,
                             # fullsummary      = True) # This is an unexpected keyword argument...?!

    print("Saving summary log file of tcleaning process...")
    np.save(imagename+'.tclean.summary.npy', rec)
    print(rec)

    print("Exporting clean image products of the continuum to FITS...")
    for ext in tclean_sv_clean_extensions:
        exportfits(imagename+ext, imagename+ext+'.fits', dropstokes=True, overwrite=True)
    for ext in tclean_rm_clean_extensions:
        print('(Space saving) Deleting this clean file: ', imagename+ext)
        os.system('rm -rf '+ imagename+ext)


    """ Do JvM correction (primary beam correction done concurrently) """

    print('Starting JvM correction of the continuum...')
    JvM_epsilon = do_JvM_correction_and_get_epsilon(root=imagename)
    # Expect: WARN	ImageExprCalculator::compute	image units are not the same: 'Jy/beam' vs ''. Proceed with caution. Output image metadata will be copied from .clean_convolved_model_temp.image

    print("Primary beam correcting the continuum...")
    for ext,pbcor_ext in zip(JvM_extensions,JvM_pbcor_extensions):
        impbcor(imagename=imagename+ext, pbimage=imagename+'.pb',
                outfile=imagename+pbcor_ext, overwrite=True)

    print("Finished primary beam correction, exporting FITS of the continuum...")
    for ext in JvM_extensions + JvM_pbcor_extensions + ['.mask']:
        export_cube = imagename+ext
        exportfits(imagename=export_cube, fitsimage=export_cube+'.fits',
                   dropstokes=True, overwrite=True)

    print("Saving imaging metrics and image metrics...")
    imaging_metrics = {'.dirty rms (uJy/beam)' : rms,
                       'tclean threshold' : threshold,
                       'JvM epsilon' : JvM_epsilon,
                       '.clean.residual peak (unitless)' : estimate_peak_intensity(imagename=imagename+'.residual', mask=mask),
                       }
    imaging_info = pd.DataFrame(data=imaging_metrics, index=[0])
    imaging_info.to_csv(imagename.replace('.clean', '.imaging_info')+'.csv')
    print("Saved imaging info csv.")

    rows = []
    for ext in JvM_extensions + JvM_pbcor_extensions:
        image_metrics = {}
        image_metrics['peak intensity (mJy/beam)'] = estimate_peak_intensity(imagename=imagename+ext, mask=mask)
        image_metrics['disk flux (mJy)'] = estimate_disk_flux(imagename=imagename+ext, mask=mask)
        image_metrics['rms noise (uJy/beam)'] = estimate_rms(imagename=imagename+ext, region=region)
        image_metrics['SNR'] = estimate_SNR(imagename=imagename+ext, mask=mask, region=region)
        rows.append(pd.Series(image_metrics, name=ext))

    image_info = pd.concat(rows, axis=1)
    image_info.transpose().to_csv(imagename.replace('.clean', '.image_info')+'.csv')
    print("Saved image info csv.")

    print("Done processing the continuum! (for robust %.1f)"%(robust))






"""
######################################################
################ IMAGE THE CONTINUUM #################
######################################################
"""

for robust in [0.5, 1.0, 1.5]:
    vis             = ddata.data_dict['continuum']
    mask            = dmask.mask_dict['continuum']['circle mask']
    region          = dmask.mask_dict['continuum']['noise annulus']
    robust          = robust
    imagename       = ddata.data_dict['NRAO_path']+'images_continuum/'+'robust'+str(robust)+'/ABAur_continuum'

    tclean_wrapper_continuum(vis        = vis,
                             imagename  = imagename,
                             mask       = mask,
                             region     = region,
                             robust     = robust
                             )



# SAVE .LASTS?
# MAKE CELL SIZE ROBUST PARAM. DEPENDENT; BEAM NOT SUFFICIENTLY SAMPLED IN ROBUST 0 IMAGE

#arguments to set every time:
# robust value
# taper or no taper

# all other arguments same

#PSF is unitless
#model is Jy/pixel
#residual is unitless
#pb is unitless

# dirty image (get rms)
# clean to to different levels of rms
# do all robust values
# do some uv tapering?


















# sys.exit()
