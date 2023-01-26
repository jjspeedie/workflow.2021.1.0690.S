'''
ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

Script to get imview to create pngs of all channels in a cube. Useful for
investigating the results of a clean run.

To run this script, do:
casa -r 6.2.1-7-pipeline-2021.2.0.128
#### <CASA> execfile('dictionary_disk.py') # loads disk_dict
#### <CASA> execfile('dictionary_data.py') # loads data_dict
[instead, be in /data directory, which contains images_lines/]
<CASA> execfile('dictionary_lines.py') # loads line_dict
<CASA> execfile('imview_channels.py')

Note you might need to open the viewer interactively first, close it again, and
then this script will run properly. While the script runs, you should see the
viewer opening and closing automatically.
'''

import os
# execfile('dictionary_disk.py')  # loads disk_dict
# execfile('dictionary_data.py') # loads data_dict
execfile('dictionary_lines.py') # loads line_dict

molecules       = ['C18O', '13CO', '12CO', 'SO']
vres_version    = 'v2' # 17-Jan-2023

for line in molecules:
    for robust in [1.5]:
        for cont in ['']:#, '_wcont']:

            # imagename       = 'images_lines/'+line+'/'+vres_version+'_robust'+str(robust)+cont+'/ABAur_'+line
            imagename       = 'images_lines/'+line+'/'+vres_version+'_initial_mask'+cont+'/ABAur_'+line

            contour_im      = imagename+'.initial.model.copy.initial.mask'

            for ext in ['.initial.image', '.initial.residual', '.initial.model']:

                im          = imagename+ext
                nchan       = line_dict[line][vres_version+'_nchan']

                os.system('mkdir '+im+'.imview_channels')

                # for chan_idx in [0+i for i in range(nchan)]:
                #     imview(raster = {'file':im},
                #           contour = {'file':contour_im},
                #           zoom    = {'file':im, 'channel': chan_idx},
                #           out     = im+'.imview_channels/img_%04i.png'%(chan_idx) )




sys.exit()
