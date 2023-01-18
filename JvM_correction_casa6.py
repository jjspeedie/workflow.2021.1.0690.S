import numpy as np
import shutil
import casatasks
import casatools
import os

from casatasks import exportfits # Jess added this line

def gaussian_eval(params, data, center):
    """Returns a gaussian with the given parameters"""
    width_x, width_y, rotation = params
    rotation = 90-rotation
    width_x = float(width_x)
    width_y = float(width_y)

    rotation = np.deg2rad(rotation)
    x, y = np.indices(data.shape) - center

    xp = x * np.cos(rotation) - y * np.sin(rotation)
    yp = x * np.sin(rotation) + y * np.cos(rotation)
    g = 1.*np.exp(
        -(((xp)/width_x)**2+
          ((yp)/width_y)**2)/2.)
    return g

def do_JvM_correction_and_get_epsilon(root, taper_match=None):
    # Get the psf file to fit
    psf_file = root + '.psf'
    model_file = root + '.model'
    residual_file = root + '.residual'
    npix_window = 201

    # Open psf and read off the metadata
    ia = casatools.image()
    ia.open(psf_file)
    psf_data_raw = ia.getregion()
    hdr = ia.summary(list=False)
    ia.close()

    delta = np.abs(hdr['incr'][0]*206265)
    try:
        rb = hdr['restoringbeam']
        major = rb['major']['value']
        minor = rb['minor']['value']
        phi = rb['positionangle']['value']
    except:
        major = hdr['perplanebeams']['beams']['*0']['*0']['major']['value']
        minor = hdr['perplanebeams']['beams']['*0']['*0']['minor']['value']
        phi = hdr['perplanebeams']['beams']['*0']['*0']['positionangle']['value']

    print("The CASA fitted beam is " + str(major) + "x" + str(minor) + '" at ' + str(phi) + "deg")

    npix = psf_data_raw.shape[0]         # Assume image is square

    # Check if image cube, or just single psf; this example doesn't handle the full
    # polarization case - implicitly assumes we can drop Stokes
    # If single psf, add an axis so we can use a single loop
    psf_data = np.squeeze(psf_data_raw)
    if len(psf_data.shape) == 2:
        psf_data = np.expand_dims(psf_data, axis=2)

    # Roll the axes to make looping more straightforward
    psf_rolled = np.rollaxis(psf_data,2)

    # Window out the region we want to consider
    i_min = int(npix/2-(npix_window-1)/2)
    i_max = int(npix/2+(npix_window-1)/2 + 1)
    for j in range(10):
        psf_windowed = psf_rolled[12*j][i_min:i_max,i_min:i_max]
        if np.sum(psf_windowed) > 0:
            break

    # Mask out anything beyond the first null
    psf_windowed[psf_windowed<0.] = -1.
    psf_windowed = np.fft.fftshift(psf_windowed)

    for i in range(psf_windowed.shape[0]):
        left_edge = np.argmax(psf_windowed[i] < 0.)
        right_edge = 201-np.argmax(psf_windowed[i][::-1] < 0.)
        psf_windowed[i][left_edge:right_edge] = 0.

    psf_windowed = np.fft.fftshift(psf_windowed)

    # Create a clean beam to evaluate against
    clean_beam = gaussian_eval([major/2.355/delta, minor/2.355/delta, phi], psf_windowed,
                               (npix_window-1)/2)

    # Calculate epsilon
    epsilon = np.sum(clean_beam)/np.sum(psf_windowed)
    print("Epsilon = " + str(epsilon))


    if taper_match:
        # create the convolved model
        convolved_temp_image = '{:s}_convolved_model_temp.image'.format(root)
        casatasks.imsmooth(imagename=model_file, major=str(major)+'arcsec', minor=str(minor)+'arcsec',
                 pa=str(phi)+'deg', targetres=True, outfile=convolved_temp_image)

        # doing the correction
        try:
            shutil.rmtree(root+".JvMcorr.temp.image")
        except:
            pass
        casatasks.immath(imagename=[convolved_temp_image, residual_file],
               expr='IM0 + ' + str(epsilon) + '*IM1', outfile=root+".JvMcorr.temp.image")


        # now smooth to taper_match
        try:
            shutil.rmtree(root+".JvMcorr.image")
        except:
            pass

        casatasks.imsmooth(imagename=root+".JvMcorr.temp.image", major=str(taper_match)+'arcsec', minor=str(taper_match)+'arcsec',
             pa=str(phi)+'deg', targetres=True, outfile=root+".JvMcorr.image")

        if not os.path.exists(root+".JvMcorr.image"):
            casatasks.imsmooth(imagename=root+".JvMcorr.temp.image", major=str(taper_match*1.05)+'arcsec', minor=str(taper_match*1.05)+'arcsec',
                 pa=str(phi)+'deg', targetres=True, outfile=root+".JvMcorr.image")
            if not os.path.exists(root+".JvMcorr.image"):
                casatasks.imsmooth(imagename=root+".JvMcorr.temp.image", major=str(taper_match*1.1)+'arcsec', minor=str(taper_match*1.1)+'arcsec',
                     pa=str(phi)+'deg', targetres=True, outfile=root+".JvMcorr.image")

        print("Wrote " + root + ".JvMcorr.image")

        # clean up
        shutil.rmtree(convolved_temp_image)
        shutil.rmtree(root+".JvMcorr.temp.image")

    else:
        # Regular JvM correction
        # create the convolved model
        convolved_temp_image = '{:s}_convolved_model_temp.image'.format(root)
        casatasks.imsmooth(imagename=model_file, major=str(major)+'arcsec', minor=str(minor)+'arcsec',
                 pa=str(phi)+'deg', targetres=True, outfile=convolved_temp_image)

        # doing the correction
        try:
            shutil.rmtree(root+".JvMcorr.image")
        except:
            pass
        casatasks.immath(imagename=[convolved_temp_image, residual_file], expr='IM0 + ' + str(epsilon) + '*IM1',
                outfile=root+".JvMcorr.image", imagemd=convolved_temp_image) # Jess specified imagemd
        print("Wrote " + root + ".JvMcorr.image")

        # clean up
        exportfits(convolved_temp_image, convolved_temp_image+'.fits', dropstokes=True, overwrite=True) # Jess added this line
        shutil.rmtree(convolved_temp_image)

        # 'lowres' JvM correction
        # create the convolved model
        casatasks.imsmooth(imagename=model_file, major=str(np.sqrt(1./epsilon)*major)+'arcsec',
                 minor=str(np.sqrt(1./epsilon)*minor)+'arcsec', pa=str(phi)+'deg',
                 targetres=True, outfile=convolved_temp_image)

        # doing the correction
        try:
            shutil.rmtree(root+".JvMcorr_lowres.image")
        except:
            pass
        casatasks.immath(imagename=[convolved_temp_image, residual_file], expr='IM0 + IM1',
               outfile=root+".JvMcorr_lowres.image", imagemd=convolved_temp_image) # Jess specified imagemd
        print("Wrote " + root + ".JvMcorr_lowres.image")

        # clean up
        shutil.rmtree(convolved_temp_image)

    return epsilon
