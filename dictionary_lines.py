"""
ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

Dictionary of emission line properties and spectral parameters for imaging.
"""

line_dict = {}

line_dict['12CO'] = {'qn'    : '2-1',                # Transition
                     'freq'  : '230.538000GHz',      # Rest frequency
                     'nchan' : 1920,                 # Native number of channels
                     'width' : '0.039684km/s' ,      # Native spectral resolution

                     # Spectral resolution for imaging, version 1 (14-Jan-2023):
                     # 'v1_nchan' : 480,
                     # 'v1_start' : '-4.23km/s',
                     # 'v1_width' : '0.042km/s',
                     # Spectral resolution for imaging, version 2 (17-Jan-2023):
                     # 'v2_nchan' : 302,
                     # 'v2_start' : '-0.198km/s',
                     # 'v2_width' : '0.042km/s',
                     # Spectral resolution for imaging, version 3 (26-Jan-2023):
                     'v3_nchan' : 302,
                     'v3_start' : '-0.198km/s',
                     'v3_width' : '0.042km/s'
                    }

line_dict['13CO'] = {'qn'    : '2-1',                # Transition
                     'freq'  : '220.398684GHz',      # Rest frequency
                     'nchan' : 1920,                 # Native number of channels
                     'width' : '0.041510km/s' ,      # Native spectral resolution

                     # Spectral resolution for imaging, version 1 (14-Jan-2023):
                     # 'v1_nchan' : 480,
                     # 'v1_start' : '-4.23km/s',
                     # 'v1_width' : '0.042km/s',
                     # Spectral resolution for imaging, version 2 (17-Jan-2023):
                     # 'v2_nchan' : 186,
                     # 'v2_start' : '2.070km/s',
                     # 'v2_width' : '0.042km/s',
                     # Spectral resolution for imaging, version 3 (26-Jan-2023):
                     # 'v3_nchan' : 186,
                     # 'v3_start' : '2.070km/s',
                     # 'v3_width' : '0.042km/s'
                     # Spectral resolution for imaging, version 4 (27-Feb-2023):
                     # 'v4_nchan' : 186,
                     # 'v4_start' : '2.070km/s',
                     # 'v4_width' : '0.042km/s'
                     # Spectral grid offset for channel dropping diagnostics, version 5 (7-Mar-2023):
                     # 'v5_nchan' : 181,
                     # 'v5_start' : '2.238km/s', # offset by 5 channels (clean goes in chunks of 10, and rms noise estimation for clean threshold needs 10 channel buffer)
                     # 'v5_width' : '0.042km/s'
                     # Spectral resolution for imaging (forcing major cycles), version 6 (8-Mar-2023):
                     # 'v6_nchan' : 186,
                     # 'v6_start' : '2.070km/s',
                     # 'v6_width' : '0.042km/s'
                     # Spectral resolution for imaging (forcing major cycles), version 7 (10-Mar-2023):
                     'v7_nchan' : 186,
                     'v7_start' : '2.070km/s',
                     'v7_width' : '0.042km/s'
                    }

line_dict['C18O'] = {'qn'    : '2-1',                 # Transition
                     'freq'  : '219.560358GHz',       # Rest frequency
                     'nchan' : 960,                   # Native number of channels
                     'width' : '0.083336km/s',        # Native spectral resolution

                     # Spectral resolution for imaging, version 1 (14-Jan-2023):
                     # 'v1_nchan' : 120,
                     # 'v1_start' : '0.81km/s',
                     # 'v1_width' : '0.084km/s',
                     # Spectral resolution for imaging, version 2 (17-Jan-2023):
                     # 'v2_nchan' : 82,
                     # 'v2_start' : '2.658km/s',
                     # 'v2_width' : '0.084km/s',
                     # Spectral resolution for imaging, version 3 (26-Jan-2023):
                     'v3_nchan' : 82,
                     'v3_start' : '2.658km/s',
                     'v3_width' : '0.084km/s'
                    }

line_dict['SO']   = {'qn'    : '6(5)-5(4)',           # Transition
                     'freq'  : '219.949442GHz',       # Rest frequency
                     'nchan' : 960,                   # Native number of channels
                     'width' : '0.0831883km/s',       # Native spectral resolution

                     # Spectral resolution for imaging, version 1 (14-Jan-2023):
                     # 'v1_nchan' : 120,
                     # 'v1_start' : '0.81km/s',
                     # 'v1_width' : '0.084km/s',
                     # Spectral resolution for imaging, version 2 (17-Jan-2023):
                     # 'v2_nchan' : 60,
                     # 'v2_start' : '3.582km/s',
                     # 'v2_width' : '0.084km/s',
                     # Spectral resolution for imaging, version 3 (26-Jan-2023):
                     'v3_nchan' : 60,
                     'v3_start' : '3.582km/s',
                     'v3_width' : '0.084km/s'
                    }
