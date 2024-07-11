# python3 -m venv modularcasa
# source modularcasa/bin/activate
# (modularcasa) pip install --upgrade pip wheel
# (modularcasa) pip install casatasks==6.2.1.7 # 6.4.3.27
# (modularcasa) pip install casatools==6.2.1.7 # 6.4.3.27
# (modularcasa) pip install casadata==2022.9.5
# (modularcasa) pip install numba
# (modularcasa) pip install astropy
# (modularcasa) pip install shutils
#
# (modularcasa) python step2_phase_alignment.py

######
# ALIGN DATA (go from *initcont_selfcal.ms to *initcont_shift.ms)
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
reference_for_LB_alignment  = ddata.data_dict['LB_EB1']['_initcont_selfcal.ms']
assert 'LB' in reference_for_LB_alignment, 'you need to choose an LB EB for alignment of LBs'
print('Selected reference long baseline execution block: ', reference_for_LB_alignment)

"""All the other EBs will be aligned to the reference EB.
We also include the reference EB itself to make sure the coordinate changes are copied over."""
offset_LB_EBs       = [ddata.data_dict[EB]['_initcont_selfcal.ms'] for EB in ddata.data_dict['LB_EBs']]
print('List of long baseline execution blocks to be aligned: ', offset_LB_EBs)

"""Select the continuum spw with the large bandwidth"""
continuum_spw_id    = 0
print('Continuum spectral window in all (ALL) execution blocks: ', continuum_spw_id)


"""Find the relative offsets and update the phase centers for all offset_EBs."""
"""cell_size defines the size of the uv grid"""

alignment_npix = {'LB':2000,'SB':500}
alignment_cell_size = {'LB':0.01,'SB':0.04}
# alignment_npix = {'LB':1024,'SB':102} # exoALMA
# alignment_cell_size = {'LB':0.01,'SB':0.1} # exoALMA
alignment_plot_file_template = 'alignment_uv_grid.png' # this will become a suffix
alignment.align_measurement_sets(reference_ms       = reference_for_LB_alignment,
                                 align_ms           = offset_LB_EBs,
                                 npix               = alignment_npix['LB'],
                                 cell_size          = alignment_cell_size['LB'],
                                 spwid              = continuum_spw_id,
                                 plot_uv_grid       = True,
                                 plot_file_template = alignment_plot_file_template)


# This is without per-EB self-cal:
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_LB_EB1_initcont.ms
#no shift, reference MS.
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_LB_EB2_initcont.ms
#requires a shift of [-0.0013834,-0.0032833]
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_LB_EB3_initcont.ms
#requires a shift of [-0.0073534,8.5593e-05]
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_LB_EB4_initcont.ms
#requires a shift of [0.024526,-0.023529]
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_LB_EB5_initcont.ms
#requires a shift of [-0.01779,-0.0037046]
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_LB_EB6_initcont.ms
#requires a shift of [0.0024483,-0.03545]
# alignment_offsets['LB_EB1'] = [0,0]
# alignment_offsets['LB_EB2'] = [-0.0013834,-0.0032833]
# alignment_offsets['LB_EB3'] = [-0.0073534,8.5593e-05]
# alignment_offsets['LB_EB4'] = [0.024526,-0.023529]
# alignment_offsets['LB_EB5'] = [-0.01779,-0.0037046]
# alignment_offsets['LB_EB6'] = [0.0024483,-0.03545]

# This is with per-EB self-cal:
#New coordinates for workflow/step1/ABAur_LB_EB1_initcont_selfcal.ms
#no shift, reference MS.
#New coordinates for workflow/step1/ABAur_LB_EB2_initcont_selfcal.ms
#requires a shift of [-0.019783,0.0046808]
#New coordinates for workflow/step1/ABAur_LB_EB3_initcont_selfcal.ms
#requires a shift of [-0.012707,0.001336]
#New coordinates for workflow/step1/ABAur_LB_EB4_initcont_selfcal.ms
#requires a shift of [0.023685,-0.021869]
#New coordinates for workflow/step1/ABAur_LB_EB5_initcont_selfcal.ms
#requires a shift of [-0.020729,-0.01054]
#New coordinates for workflow/step1/ABAur_LB_EB6_initcont_selfcal.ms
#requires a shift of [-0.0045512,-0.0276]

for EB in ddata.data_dict['LB_EBs']:
    initcont_selfcal          = ddata.data_dict[EB]['_initcont_selfcal.ms']
    initcont_selfcal_shift    = initcont_selfcal.replace('.ms', '_shift.ms')
    os.system('mv '+initcont_selfcal_shift+' ./workflow/step2')
os.system('mv ./workflow/step1/*'+alignment_plot_file_template+'  ./workflow/step2')


alignment_offsets   = {}
#insert offsets from the alignment output
alignment_offsets['LB_EB1'] = [0,0]
alignment_offsets['LB_EB2'] = [-0.019783,0.0046808]
alignment_offsets['LB_EB3'] = [-0.012707,0.001336]
alignment_offsets['LB_EB4'] = [0.023685,-0.021869]
alignment_offsets['LB_EB5'] = [-0.020729,-0.01054]
alignment_offsets['LB_EB6'] = [-0.0045512,-0.0276]

shifted_LB_EBs = [ddata.data_dict[EB]['_initcont_selfcal_shift.ms'] for EB in ddata.data_dict['LB_EBs']]
print('List of long baseline execution blocks that have been aligned: ', shifted_LB_EBs)

"""to check if alignment worked, calculate shift again and verify that shifts are small (i.e.
a fraction of the cell size):"""
for i,shifted_ms in enumerate(shifted_LB_EBs):
    if i==0: #
        #for some reason the fitter fails when computing the offset of an EB to itself,
        #so we skip the ref EB
        continue
    offset = alignment.find_offset(reference_ms     = reference_for_LB_alignment,
                                   offset_ms        = shifted_ms,
                                   npix             = alignment_npix['LB'],
                                   cell_size        = alignment_cell_size['LB'],
                                   spwid            = continuum_spw_id)
    print(f'#offset for {shifted_ms}: ',offset)

# This is without per-EB self-cal:
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_LB_EB1_initcont_shift.ms:  [-4.51798589e-10 -3.36817461e-09]
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_LB_EB2_initcont_shift.ms:  [-1.42397170e-06 -3.31838587e-05]
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_LB_EB3_initcont_shift.ms:  [-3.37030700e-05 -2.14097102e-05]
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_LB_EB4_initcont_shift.ms:  [ 0.00011575 -0.0001446 ]
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_LB_EB5_initcont_shift.ms:  [-4.58253229e-05  1.04505131e-04]
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_LB_EB6_initcont_shift.ms:  [-0.00015515 -0.00030391]

# This is with per-EB self-cal:
#offset for workflow/step2/ABAur_LB_EB2_initcont_selfcal_shift.ms:  [-0.00053475 -0.00016906]
#offset for workflow/step2/ABAur_LB_EB3_initcont_selfcal_shift.ms:  [1.71563046e-04 2.47059087e-05]
#offset for workflow/step2/ABAur_LB_EB4_initcont_selfcal_shift.ms:  [ 5.09346597e-04 -1.23296973e-05]
#offset for workflow/step2/ABAur_LB_EB5_initcont_selfcal_shift.ms:  [-7.52639036e-05 -3.63954304e-04]
#offset for workflow/step2/ABAur_LB_EB6_initcont_selfcal_shift.ms:  [ 0.00052468 -0.00087966]


"""Merge shifted LB EBs for aligning SB EBs"""
# print('Concatenating shifted LB EBs for aligning SB EBs...')
path, _ = os.path.split(shifted_LB_EBs[0])
LB_concat_shifted = path + '/ABAur_LB_concat_shifted.ms'
os.system(f'rm -rf {LB_concat_shifted}')
casatasks.concat(vis            = shifted_LB_EBs,
                 concatvis      = LB_concat_shifted,
                 dirtol         = '0.1arcsec',
                 copypointing   = False)


"""
######################################################
######## Align SB EBs to the LB concatenation ########
######################################################
"""

"""Align SB EBs to concat shifted LB EBs"""
reference_for_SB_alignment = LB_concat_shifted

offset_SB_EBs = [ddata.data_dict[EB]['_initcont_selfcal.ms'] for EB in ddata.data_dict['SB_EBs']]
print('List of short baseline execution blocks to be aligned: ', offset_SB_EBs)

alignment.align_measurement_sets(reference_ms       = reference_for_SB_alignment,
                                 align_ms           = offset_SB_EBs,
                                 npix               = alignment_npix['SB'],
                                 cell_size          = alignment_cell_size['SB'],
                                 spwid              = continuum_spw_id,
                                 plot_uv_grid       = True,
                                 plot_file_template = alignment_plot_file_template)

# This is without per-EB self-cal:
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_SB_EB1_initcont.ms
#requires a shift of [-0.012473,0.064629]
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_SB_EB2_initcont.ms
#requires a shift of [0.11748,-0.16897]
# results when aligning to reference LB alone, and not the concatenation: (not different at all...!)
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_SB_EB1_initcont.ms
#requires a shift of [-0.012473,0.064629]
#New coordinates for /arc/projects/abaur/workflow/step1/ABAur_SB_EB2_initcont.ms
#requires a shift of [0.11748,-0.16897]
# alignment_offsets['SB_EB1'] = [-0.012473,0.064629]
# alignment_offsets['SB_EB2'] = [0.11748,-0.16897]

# This is with per-EB self-cal:
#New coordinates for workflow/step1/ABAur_SB_EB1_initcont_selfcal.ms
#requires a shift of [-0.013133,0.03949]
#New coordinates for workflow/step1/ABAur_SB_EB2_initcont_selfcal.ms
#requires a shift of [0.11922,-0.19222]

for EB in ddata.data_dict['SB_EBs']:
    initcont_selfcal          = ddata.data_dict[EB]['_initcont_selfcal.ms']
    initcont_selfcal_shift    = initcont_selfcal.replace('.ms', '_shift.ms')
    os.system('mv '+initcont_selfcal_shift+' ./workflow/step2')
os.system('mv ./workflow/step1/*'+alignment_plot_file_template+'  ./workflow/step2')


#insert offsets from the alignment output
alignment_offsets['SB_EB1'] = [-0.013133,0.03949]
alignment_offsets['SB_EB2'] = [0.11922,-0.19222]

shifted_SB_EBs = [ddata.data_dict[EB]['_initcont_selfcal_shift.ms'] for EB in ddata.data_dict['SB_EBs']]
print('List of short baseline execution blocks that have been aligned: ', shifted_LB_EBs)

"""Again: to check if alignment worked, calculate shift again and verify that shifts are small (i.e.
a fraction of the cell size):"""
for shifted_ms in shifted_SB_EBs:
    offset = alignment.find_offset(reference_ms     = reference_for_SB_alignment,
                                   offset_ms        = shifted_ms,
                                   npix             = alignment_npix['SB'],
                                   cell_size        = alignment_cell_size['SB'],
                                   spwid            = continuum_spw_id)
    print(f'#offset for {shifted_ms}: ',offset)

# This is without per-EB self-cal:
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_SB_EB1_initcont_shift.ms:  [0.00044042 0.00118166]
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_SB_EB2_initcont_shift.ms:  [ 0.00026603 -0.00286368]

# This is with per-EB self-cal:
#offset for workflow/step2/ABAur_SB_EB1_initcont_selfcal_shift.ms:  [0.00011743 0.00072109]
#offset for workflow/step2/ABAur_SB_EB2_initcont_selfcal_shift.ms:  [-0.00033063 -0.00405134]

# Jess: I guess this is here in the exoALMA script just out of interest
# for shifted_ms in shifted_SB_EBs[1:]:
#     offset = alignment.find_offset(reference_ms=shifted_SB_EBs[0],
#                                    offset_ms=shifted_ms,npix=alignment_npix['SB'],
#                                    cell_size=alignment_cell_size['SB'],
#                                    spwid=continuum_spw_id)
#     print(f'#offset for {shifted_ms} to SB EB0: ',offset)
# This is without per-EB self-cal:
#offset for /arc/projects/abaur/workflow/step2_noselfcal/ABAur_SB_EB2_initcont_shift.ms to SB EB0:  [-0.01041825 -0.02799446]
# This is with per-EB self-cal:
#offset for workflow/step2/ABAur_SB_EB2_initcont_selfcal_shift.ms to SB EB0:  [-0.01378587 -0.01481389]

"""Check that the phase center of each EB is: 04:55:45.854900 +30.33.03.73320 J2000"""
for EB in ddata.data_dict['EBs']:
    vis          = ddata.data_dict['NRAO_path']+ddata.data_dict[EB]['_initcont_selfcal_shift.ms']
    casatasks.listobs(vis=vis, listfile=vis+'.listobs.txt')

"""Merge shifted SB EBs for continuing into step 3: self calibration"""
print('Concatenating shifted SB EBs...')
path, _ = os.path.split(shifted_SB_EBs[0])
SB_concat_shifted = path + '/ABAur_SB_concat_shifted_contp0.ms'
os.system(f'rm -rf {SB_concat_shifted}')
casatasks.concat(vis            = shifted_SB_EBs,
                 concatvis      = SB_concat_shifted,
                 dirtol         = '0.1arcsec',
                 copypointing   = False)

"""Move concatenated shifted LB and SB EBs to step 3 folder"""
os.system('mv '+SB_concat_shifted+' '+ddata.data_dict['NRAO_path']+'workflow/step3')
os.system('mv '+LB_concat_shifted+' '+ddata.data_dict['NRAO_path']+'workflow/step3')

sys.exit()

# exoALMA: Remove the '_selfcal' part of the names to match the naming convention below.
# for shifted_EB in shifted_LB_EBs+shifted_SB_EBs:
#     os.system('mv {} {}'.format(shifted_EB, shifted_EB.replace('_selfcal', '')))
