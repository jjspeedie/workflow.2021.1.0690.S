"""
ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

Dictionary of physical properties of the AB Aur system for many different uses.
When applicable/possible, I've tried to cite the *original* paper who found each value.
"""

disk_dict = {'name': 'AB Aur',
             'distance': 155.9403,                      # +/- 0.9045; source distance in pc, from Gaia DR3 archive calculated as 1000/(parallax), error as 1000*(error in p)/(p^2)
             'incl': 23.2,                              # +/- ?; inclination in degrees, Huang et al. (2020) (from the dust)
             'PA': 54.,                                 # +/- ?; position angle in degrees, Tang et al. (2017)
             'PA_gofish': 54. + 180.0,                  # +/- ?; position angle in degrees, corrected for gofish
             'M_star': 2.4,                             # +/- 0.2; stellar mass in solar masses, DeWarf et al. 2003, Herbertz et al. 1991

             'L_star_23': 59, 				# +/- 5; stellar luminosity in solar luminosities, estimated from radiative transfer modelling, Currie et al. 2023
             'L_star_03': 38, 				# +/- ?; stellar luminosity in solar luminosities, deWarf et al. 2003. Note, given as ~38
             'L_star_98': 47, 				# +17/-10; stellar luminosity in solar luminosities, van den Ancker et al. 1998
             'L_star_92': 57, 				# +/- ?; stellar luminosity in solar luminosities, estimated from Hillenbrand et al. 1992

             'T_eff_23': 9770, 				# +/- 100; effective temperature in Kelvin, estimated from spectral type, Currie et al. 2023
             'T_eff_01': 9500, 				# +/- ?; effective temperature in Kelvin, estimated from SED modeling, Natta et al. 2001
             'T_eff_92': 9772, 				# +/- ?; effective temperature in Kelvin, first reported by Hillenbrand 1992

             'v_sys': 5.85,                             # +/- 0.02; LSR systemic velocity [km/s], from Tang et al. (2012) outer regions (note 5.73 km/s for inner regions)
             'RA_phase_center': '04:55:45.854900',      # J2000 coordinates of LB1 phase center, to which all other EBs were aligned during reduction
             'Dec_phase_center': '+30.33.03.73320'      # J2000 coordinates of LB1 phase center, to which all other EBs were aligned during reduction
}


# Other properties that MAPS included in their disk dictionaries (GM AUR):
             # 'PA_diskprojection': 57.17 + 180.0, #position angle in degrees, corrected for diskprojection
             # 'M_star_dynamical': 1.1, # from Rich Teague's gofish
             # 'logMdot': -8.1, #log of stellar accretion rate in solar masses/yr (GM Aur is variable!), Ingleby et al. 2013
             # '12CO_extent' : 14.0, # km/s, measured in v0 images


