"""
ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

Dictionary of masking parameters for all lines and the continuum.
"""


mask_dict = {}

mask_dict['continuum'] = {'circle mask' : 'circle[[%s, %s], %.1farcsec]' %('04h55m45.8549s', '+30.33.03.733', 2.), # J2000 common/aligned phase center
                          'noise annulus' : "annulus[[%s, %s],['%.2farcsec', '10.0arcsec']]" % ('04h55m45.8549s', '+30.33.03.733', 4.) # J2000 common/aligned phase center
                         }

# The following is for cleaning, using with Rich Teague's keplerian_mask, based on experimentation:
mask_dict['SO']   = {'r_max': 3.0,          # Maximum radius in [arcsec] of the mask.
                     'dV0': 400.0,          # The Doppler width of the line in [m/s] at 1 arcsec.
                     'dVq': -0.5,           # The exponent of the power law describing the Doppler width as a function of radius.
                     'zr': 0.1,             # For elevated emission, the z/r value.
                     'target_res': 0.5      # Instead of scaling the CLEAN beam for the convolution kernel, specify the FWHM of the convolution kernel directly.
                     }

# The following is for  kickstarting auto-multithresh to help it capture diffuse emission using make_mask_for_diffuse_emission() from keplerian_mask, which is a hack by Jess:
# mask_dict['SO_diffuse_emission']   = {'r_max': 15.0,         # Maximum radius in [arcsec] of the mask.
#                                       'v_min': 3000.,        # Maximum velocity in [m/s] of the mask.
#                                       'v_max': 8000.,        # Maximum velocity in [m/s] of the mask.
#                                       }
# mask_dict['C18O_diffuse_emission']   = {'r_max': 15.0,         # Maximum radius in [arcsec] of the mask.
#                                         'v_min': 2000.,        # Maximum velocity in [m/s] of the mask.
#                                         'v_max': 9000.,        # Maximum velocity in [m/s] of the mask.
#                                         }
mask_dict['13CO_diffuse_emission']   = {'r_max': 15.0,         # Maximum radius in [arcsec] of the mask.
                                        'v_min': 4926.,        # Maximum velocity in [m/s] of the mask.
                                        'v_max': 6900.,        # Maximum velocity in [m/s] of the mask.
                                        }
mask_dict['12CO_diffuse_emission']   = {'r_max': 15.0,         # Maximum radius in [arcsec] of the mask.
                                        'v_min': -4212.,        # Maximum velocity in [m/s] of the mask.
                                        'v_max': 7362.,        # Maximum velocity in [m/s] of the mask.
                                        }



# # The following is for estimating noise, using with make_mask_of_linefree_channels() from keplerian_mask, which is a hack by Jess:
# mask_dict['SO_estimate_rms']   = {'r_max': 15.0,         # Maximum radius in [arcsec] of the mask.
#                                   'v_min': 3000.,        # Maximum velocity in [m/s] of the mask.
#                                   'v_max': 8000.,        # Maximum velocity in [m/s] of the mask.
#                                   }
# mask_dict['C18O_estimate_rms']   = {'r_max': 15.0,         # Maximum radius in [arcsec] of the mask.
#                                     'v_min': 2000.,        # Maximum velocity in [m/s] of the mask.
#                                     'v_max': 9000.,        # Maximum velocity in [m/s] of the mask.
#                                     }
# mask_dict['13CO_estimate_rms']   = {'r_max': 15.0,         # Maximum radius in [arcsec] of the mask.
#                                     'v_min': 1500.,        # Maximum velocity in [m/s] of the mask.
#                                     'v_max': 10000.,        # Maximum velocity in [m/s] of the mask.
#                                     }
# mask_dict['12CO_estimate_rms']   = {'r_max': 15.0,         # Maximum radius in [arcsec] of the mask.
#                                     'v_min': -500.,        # Maximum velocity in [m/s] of the mask.
#                                     'v_max': 12000.,        # Maximum velocity in [m/s] of the mask.
#                                     }


# Keplerian mask not used for the following molecules:
# mask_dict['C18O'] = {'r_max': 4.0,          # Maximum radius in [arcsec] of the mask.
#                      'dV0': 400.0,          # The Doppler width of the line in [m/s] at 1 arcsec.
#                      'dVq': -0.5,           # The exponent of the power law describing the Doppler width as a function of radius.
#                      'zr': 0.1,             # For elevated emission, the z/r value.
#                      'target_res': 0.5      # Instead of scaling the CLEAN beam for the convolution kernel, specify the FWHM of the convolution kernel directly.
#                      }
# mask_dict['13CO'] = {'r_max': 5.0,          # Maximum radius in [arcsec] of the mask.
#                      'dV0': 400.0,          # The Doppler width of the line in [m/s] at 1 arcsec.
#                      'dVq': -0.5,           # The exponent of the power law describing the Doppler width as a function of radius.
#                      'zr': 0.2,             # For elevated emission, the z/r value.
#                      'target_res': 1.0      # Instead of scaling the CLEAN beam for the convolution kernel, specify the FWHM of the convolution kernel directly.
#                      }
# mask_dict['12CO'] = {'r_max': 6.0,          # Maximum radius in [arcsec] of the mask.
#                      'dV0': 400.0,          # The Doppler width of the line in [m/s] at 1 arcsec.
#                      'dVq': -0.5,           # The exponent of the power law describing the Doppler width as a function of radius.
#                      'zr': 0.3,             # For elevated emission, the z/r value.
#                      'target_res': 2.0      # Instead of scaling the CLEAN beam for the convolution kernel, specify the FWHM of the convolution kernel directly.
#                      }















# MAPS maskdictionary.py
#
# """
# Dictionary of masking parameters for all lines other than 12CO. This should be
# used with keplerian_mask.py found at:
# > https://github.com/richteague/keplerian_mask
# Note that this has slightly different argument names than in
# disk_dictionary.py.
# """
#
#
# mask_dictionary = {}
#
# mask_dictionary['AS_209'] = {'r_max': 1.8,
#                              'dV0': 400.0,
#                              'dVq': -0.5,
#                              'zr': 0.1,
#                              'target_res': 0.5}
#
# mask_dictionary['HD_163296'] = {'r_max': 4.2,
#                                 'dV0': 500.0,
#                                 'dVq': -1.0,
#                                 'zr': 0.3,
#                                 'target_res': 0.5}
#
# mask_dictionary['MWC_480'] = {'r_max': 3.0,
#                               'dV0': 300.0,
#                               'dVq': -0.5,
#                               'zr': 0.3,
#                               'target_res': 0.5}
#
# mask_dictionary['IM_Lup'] = {'r_max': 4.0,
#                              'dV0': 400.0,
#                              'dVq': -1.0,
#                              'zr': 0.3,
#                              'target_res': 0.5}
#
# mask_dictionary['GM_Aur'] = {'r_max': 3.0,
#                              'dV0': 500.0,
#                              'dVq': -0.5,
#                              'zr': 0.2,
#                              'target_res': 0.5}
