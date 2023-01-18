######
# ALIGN DATA (go from *initcont_selfcal.ms to *initcont_shift.ms)
######

# Select the LB EB to act as the reference (usually the best SNR one).

reference_for_LB_alignment = f'{prefix}_LB_EB3_initcont_selfcal.ms'
assert 'LB' in reference_for_LB_alignment,\
            'you need to choose an LB EB for alignment of LB'

alignment_offsets = {}

# All the other EBs will be aligned to the reference EB
# We also include the reference EB itself to make sure the coordinate changes are
# copied over.

offset_LB_EBs = ['{}_{}_initcont_selfcal.ms'.format(prefix, params['name'])
                 for params in data_params_LB.values()]

#select the continuum spw with the large bandwidth
continuum_spw_id = 1
# Find the relative offsets and update the phase centers for all offset_EBs.
# cell_size defines the size of the uv grid
alignment_npix = {'LB':1024,'SB':102}
alignment_cell_size = {'LB':0.01,'SB':0.1}
alignment_plot_file_template = os.path.join(EB_selfcal_shift_folder,
                                            'alignment_uv_grid.png')
alignment.align_measurement_sets(reference_ms=reference_for_LB_alignment,
                                 align_ms=offset_LB_EBs,npix=alignment_npix['LB'],
                                 cell_size=alignment_cell_size['LB'],
                                 spwid=continuum_spw_id,plot_uv_grid=True,
                                 plot_file_template=alignment_plot_file_template)

#New coordinates for CQ_Tau_LB_EB0_initcont_selfcal.ms
#requires a shift of [0.0078172,-0.0037858].

#New coordinates for CQ_Tau_LB_EB1_initcont_selfcal.ms
#requires a shift of [0.0088372,-0.0068649].

#New coordinates for CQ_Tau_LB_EB2_initcont_selfcal.ms
#requires a shift of [0.0055495,0.0025618].

#New coordinates for CQ_Tau_LB_EB3_initcont_selfcal.ms
#no shift, reference MS.


#insert offsets from the alignment output
alignment_offsets['LB_EB0'] = [0.0078172,-0.0037858]
alignment_offsets['LB_EB1'] = [0.0088372,-0.0068649]
alignment_offsets['LB_EB2'] = [0.0055495,0.0025618]
alignment_offsets['LB_EB3'] = [0,0]

shifted_LB_EBs = [EB.replace('.ms','_shift.ms') for EB in offset_LB_EBs]

#to check if alignment worked, calculate shift again and verify that shifts are small (i.e.
#a fraction of the cell size):

for shifted_ms in shifted_LB_EBs:
    if shifted_ms == reference_for_LB_alignment.replace('.ms','_shift.ms'):
        #for some reason the fitter fails when computing the offset of an EB to itself,
        #so we skip the ref EB
        continue
    offset = alignment.find_offset(reference_ms=reference_for_LB_alignment,
                                   offset_ms=shifted_ms,npix=alignment_npix['LB'],
                                   cell_size=alignment_cell_size['LB'],
                                   spwid=continuum_spw_id)
    print(f'#offset for {shifted_ms}: ',offset)

# offset for CQ_Tau_LB_EB0_initcont_selfcal_shift.ms:  [ 1.02707611e-03 -1.27832463e-05]
# offset for CQ_Tau_LB_EB1_initcont_selfcal_shift.ms:  [7.37408140e-04 3.02619741e-06]
# offset for CQ_Tau_LB_EB2_initcont_selfcal_shift.ms:  [1.26997133e-04 4.28586053e-05]


# Merge shifted LB EBs for aligning SB EBs
LB_concat_shifted = f'{prefix}_LB_concat_shifted.ms'
os.system(f'rm -rf {LB_concat_shifted}')
concat(vis=shifted_LB_EBs,concatvis=LB_concat_shifted,dirtol='0.1arcsec',
       copypointing=False)


# Align SB EBs to concat shifted LB EBs
reference_for_SB_alignment = LB_concat_shifted

offset_SB_EBs = ['{}_{}_initcont_selfcal.ms'.format(prefix, params['name'])
                 for params in data_params_SB.values()]

alignment.align_measurement_sets(reference_ms=reference_for_SB_alignment,
                                 align_ms=offset_SB_EBs,npix=alignment_npix['SB'],
                                 cell_size=alignment_cell_size['SB'],
                                 spwid=continuum_spw_id,plot_uv_grid=True,
                                 plot_file_template=alignment_plot_file_template)

#New coordinates for CQ_Tau_SB_EB0_initcont_selfcal.ms
#requires a shift of [0.034742,-0.059789].
#New coordinates for CQ_Tau_SB_EB1_initcont_selfcal.ms
#requires a shift of [0.067768,-0.099577].


alignment_offsets['SB_EB0'] = [0.034742,-0.059789]
alignment_offsets['SB_EB1'] = [0.067768,-0.099577]

shifted_SB_EBs = [EB.replace('.ms','_shift.ms') for EB in offset_SB_EBs]

#check by calculating offset again
for shifted_ms in shifted_SB_EBs:
    offset = alignment.find_offset(reference_ms=reference_for_SB_alignment,
                                   offset_ms=shifted_ms,npix=alignment_npix['SB'],
                                   cell_size=alignment_cell_size['SB'],
                                   spwid=continuum_spw_id)
    print(f'#offset for {shifted_ms}: ',offset)
#offset for CQ_Tau_SB_EB0_initcont_selfcal_shift.ms:  [ 0.00030234 -0.00045335]
#offset for CQ_Tau_SB_EB1_initcont_selfcal_shift.ms:  [ 0.00116181 -0.00196541]
for shifted_ms in shifted_SB_EBs[1:]:
    offset = alignment.find_offset(reference_ms=shifted_SB_EBs[0],
                                   offset_ms=shifted_ms,npix=alignment_npix['SB'],
                                   cell_size=alignment_cell_size['SB'],
                                   spwid=continuum_spw_id)
    print(f'#offset for {shifted_ms} to SB EB0: ',offset)


# Remove the '_selfcal' part of the names to match the naming convention below.
for shifted_EB in shifted_LB_EBs+shifted_SB_EBs:
    os.system('mv {} {}'.format(shifted_EB, shifted_EB.replace('_selfcal', '')))
