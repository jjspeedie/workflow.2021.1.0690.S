# Run this script in a virtual environment, made with the following:
# python3 -m venv modularcasa
# source modularcasa/bin/activate
# (modularcasa) pip install bettermoments
'''
ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

Script to make moment maps using Rich Teague's bettermoments.

To run this script, do:
source modularcasa/bin/activate
(modularcasa) python make_momentmaps_bettermoments_masked.py
'''

import sys, os
import numpy as np
import bettermoments as bm
import dictionary_data as ddata # contains data_dict

extensions = ['.JvMcorr.image', '.JvMcorr.image.pbcor']#'.image', '.image.pbcor', '.JvMcorr_lowres.image', '.JvMcorr_lowres.image.pbcor' ]
molecules  = ['12CO']
robusts    = ['robust0.0']#, 'robust0.5']
vres_version = 'v12'
masks = ['keplerian_mask', 'anti_keplerian_mask']

smooth = 0      # doesn't really seem to make a difference when vres is 41 or 82 m/s
polyorder = 0   # doesn't really seem to make a difference when vres is 41 or 82 m/s

for robust in robusts:
    for ext in extensions:
        for molecule in molecules:
            for mask in masks:
                filename = ddata.data_dict['NRAO_path']+'images_lines/'+molecule+'/'+vres_version+'_'+robust+'/ABAur_'+molecule+'.clean'+ext
                maskname = ddata.data_dict['NRAO_path']+'images_lines/'+molecule+'/'+vres_version+'_'+robust+'/ABAur_'+molecule+'.clean.'+mask+'.fits'

                savedir  = ddata.data_dict['NRAO_path']+'images_lines/'+molecule+'/'+vres_version+'_'+robust+'/moments_bm_'+mask+'/'
                os.system('mkdir '+savedir)

                print('Working on '+robust+', '+ext+', '+molecule+'...')

                data, velax = bm.load_cube(filename+'.fits')

                # Noise estimation on the original data (not spectrally smoothed - this way the noise will be higher)
                rms = bm.estimate_RMS(data=data, N=9) # cube channels were selected to have 10 line-free channels at beginning and end

                # Spectral smoothing: convolution with top-hat or a Savitzky-Golay filter.
                smoothed_data = bm.smooth_data(data=data, smooth=smooth, polyorder=polyorder) # no smoothing

                # Optional: user-defined mask; IN THIS SCRIPT WE DO USE A MASK
                user_mask = bm.get_user_mask(data=data, user_mask_path=maskname)

                # Mask in the spatial plane: "sigma-clip"
                if molecule=='SO':
                    Xsigma = 3.0 # after experimentation, 3 sigma determined to be best for SO
                else:
                    Xsigma = 5.0 # with 12CO, 13CO, C18O signal is strong enough for 5sigma clip
                threshold_mask = bm.get_threshold_mask(data=data,
                                                       clip=(-np.inf, Xsigma), # to remove all pixels below Xsigma, including high significance but negative pixels
                                                       smooth_threshold_mask=3.0) # the threshold mask is smoothed/conservative

                # Mask along the velocity axis
                channel_mask = bm.get_channel_mask(data=data,
                                                  firstchannel=0, # we don't include first 10 channels
                                                  lastchannel=-1) # we don't include last 10 channels

                # Mask combination
                mask = bm.get_combined_mask(user_mask=user_mask,
                                            threshold_mask=threshold_mask,
                                            channel_mask=channel_mask,
                                            combine='and')
                masked_data = smoothed_data * mask

                # Collapse the data
                print('Making moment 0...')
                moments = bm.collapse_zeroth(velax=velax, data=masked_data, rms=rms)
                bm.save_to_FITS(moments=moments, method='zeroth', path=filename+'.fits')
                # makes M0, dM0
                print('Making moment 1...')
                moments = bm.collapse_first(velax=velax, data=masked_data, rms=rms)
                bm.save_to_FITS(moments=moments, method='first', path=filename+'.fits')
                # makes M1, dM1
                print('Making moment 2...')
                moments = bm.collapse_second(velax=velax, data=masked_data, rms=rms)
                bm.save_to_FITS(moments=moments, method='second', path=filename+'.fits')
                # makes M2, dM2
                print('Making moment 9...')
                moments = bm.collapse_ninth(velax=velax, data=masked_data, rms=rms)
                bm.save_to_FITS(moments=moments, method='ninth', path=filename+'.fits')
                # makes M9, dM9
                print('Making quadratic moment...')
                moments = bm.collapse_quadratic(velax=velax, data=masked_data, rms=rms)
                bm.save_to_FITS(moments=moments, method='quadratic', path=filename+'.fits')
                # makes Fnu, dFnu, v0, and dv0

                bm_exts = ['_M0', '_dM0', '_M1', '_dM1', '_M2', '_dM2', '_M9', '_dM9', '_Fnu', '_dFnu', '_v0', '_dv0']
                for bm_ext in bm_exts:
                    os.system('mv '+filename+bm_ext+'.fits '+savedir)




# sys.exit()
