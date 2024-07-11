# python3 -m venv modularcasa
# source modularcasa/bin/activate
# (modularcasa) pip install --upgrade pip wheel
# (modularcasa) pip install casatasks==6.2.1.7 # actually it's this: 6.4.3.27
# (modularcasa) pip install casatools==6.2.1.7 # actually it's this: 6.4.3.27
# (modularcasa) pip install casadata==2022.9.5
# (modularcasa) pip install numba
# (modularcasa) pip install astropy
# (modularcasa) pip install shutils
#
# (modularcasa) python step4_detour.py

######
# ALIGN DATA (go from *initlines.ms to *initlines_shift.ms)
######

import alignment
import sys, os
import casatasks
import dictionary_data as ddata # contains data_dict

"""
######################################################
###### Align LB EBs to a chosen LB reference EB ######
######################################################
"""

"""Select the LB EB to act as the reference (usually the best SNR one)"""
reference_for_LB_alignment  = ddata.data_dict['LB_EB1']['_initlines.ms']
assert 'LB' in reference_for_LB_alignment, 'you need to choose an LB EB for alignment of LBs'
print('Selected reference long baseline execution block: ', reference_for_LB_alignment)

"""All the other EBs will be aligned to the reference EB.
We also include the reference EB itself to make sure the coordinate changes are copied over."""
offset_LB_EBs       = [ddata.data_dict[EB]['_initlines.ms'] for EB in ddata.data_dict['LB_EBs']]
print('List of long baseline execution blocks to be aligned: ', offset_LB_EBs)

"""Update the phase centers for all offset_EBs using the offset values determined
for the continuum, from script step2_phase_alignment.py."""
alignment_offsets   = []
#insert offsets from the alignment output
alignment_offsets.append([0,0])
alignment_offsets.append([-0.019783,0.0046808])
alignment_offsets.append([-0.012707,0.001336])
alignment_offsets.append([0.023685,-0.021869])
alignment_offsets.append([-0.020729,-0.01054])
alignment_offsets.append([-0.0045512,-0.0276])

print('Aligning the long baseline execution blocks with the following offsets: ', alignment_offsets)
alignment.align_measurement_sets(reference_ms       = reference_for_LB_alignment,
                                 align_ms           = offset_LB_EBs,
                                 align_offsets      = alignment_offsets)

"""Check that the phase center of each EB is: 04:55:45.854900 +30.33.03.73320 J2000"""
for EB in ddata.data_dict['LB_EBs']:
    vis          = ddata.data_dict['NRAO_path']+ddata.data_dict[EB]['_initlines_shift.ms']
    casatasks.listobs(vis=vis, listfile=vis+'.listobs.txt')



"""
######################################################
###### Align SB EBs to a chosen LB reference EB ######
######################################################
"""

"""Select the LB EB to act as the reference (usually the best SNR one)"""
"""Previously this was the concatenation of LB EBs, but now not necessary"""
reference_for_SB_alignment  = ddata.data_dict['LB_EB1']['_initlines.ms'] # same as before
print('Selected reference long baseline execution block: ', reference_for_LB_alignment)

offset_SB_EBs = [ddata.data_dict[EB]['_initlines.ms'] for EB in ddata.data_dict['SB_EBs']]
print('List of short baseline execution blocks to be aligned: ', offset_SB_EBs)

"""Update the phase centers for all offset_EBs using the offset values determined
for the continuum, from script step2_phase_alignment.py."""
alignment_offsets   = []
#insert offsets from the alignment output
alignment_offsets.append([-0.013133,0.03949])
alignment_offsets.append([0.11922,-0.19222])

print('Aligning the long baseline execution blocks with the following offsets: ', alignment_offsets)
alignment.align_measurement_sets(reference_ms       = reference_for_SB_alignment,
                                 align_ms           = offset_SB_EBs,
                                 align_offsets      = alignment_offsets)

"""Check that the phase center of each EB is: 04:55:45.854900 +30.33.03.73320 J2000"""
for EB in ddata.data_dict['SB_EBs']:
    vis          = ddata.data_dict['NRAO_path']+ddata.data_dict[EB]['_initlines_shift.ms']
    casatasks.listobs(vis=vis, listfile=vis+'.listobs.txt')





sys.exit()
