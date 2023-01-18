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
                     # 'v1_width' : '0.042km/s'
                     'v1_nchan' : 320,
                     'v1_start' : '-0.198km/s',
                     'v1_width' : '0.042km/s'
                    }

line_dict['13CO'] = {'qn'    : '2-1',                # Transition
                     'freq'  : '220.398684GHz',      # Rest frequency
                     'nchan' : 1920,                 # Native number of channels
                     'width' : '0.041510km/s' ,      # Native spectral resolution

                     # Spectral resolution for imaging, version 1 (14-Jan-2023):
                     # 'v1_nchan' : 480,
                     # 'v1_start' : '-4.23km/s',
                     # 'v1_width' : '0.042km/s'
                     'v1_nchan' : 186,
                     'v1_start' : '2.070km/s',
                     'v1_width' : '0.042km/s'
                    }

line_dict['C18O'] = {'qn'    : '2-1',                 # Transition
                     'freq'  : '219.560358GHz',       # Rest frequency
                     'nchan' : 960,                   # Native number of channels
                     'width' : '0.083336km/s',        # Native spectral resolution

                     # Spectral resolution for imaging, version 1 (14-Jan-2023):
                     # 'v1_nchan' : 120,
                     # 'v1_start' : '0.81km/s',
                     # 'v1_width' : '0.084km/s'
                     'v1_nchan' : 80,
                     'v1_start' : '2.658km/s',
                     'v1_width' : '0.084km/s'
                    }

line_dict['SO']   = {'qn'    : '6(5)-5(4)',           # Transition
                     'freq'  : '219.949442GHz',       # Rest frequency
                     'nchan' : 960,                   # Native number of channels
                     'width' : '0.0831883km/s',       # Native spectral resolution

                     # Spectral resolution for imaging, version 1 (14-Jan-2023):
                     # 'v1_nchan' : 120,
                     # 'v1_start' : '0.81km/s',
                     # 'v1_width' : '0.084km/s'
                     'v1_nchan' : 60,
                     'v1_start' : '3.582km/s',
                     'v1_width' : '0.084km/s'
                    }
