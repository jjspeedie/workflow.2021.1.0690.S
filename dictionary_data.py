"""
ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

Dictionary of data file names, paths, and properties used during reduction.
"""

import numpy as np
# execfile('dictionary_disk.py') # either this or the below line, depending if working in casa or modularcasa
import dictionary_disk as ddisk # and replace disk_dict with ddisk.disk_dict

"""
######################################################
########### Dictionary for data properties ###########
######################################################
"""

"""######### Properties common to all EBs #########"""

CANFAR_path     = '/arc/projects/abaur/'
# NRAO_path       = '/lustre/cv/observers/cv-13687/data/'   # short term account (~Nov-Dec 2022)
NRAO_path       = '/lustre/cv/users/jspeedie/data/'         # long term account (1 year starting Jan. 2023)
LB_raw_path     = '2021.1.00690.S/science_goal.uid___A001_X15a2_Xb69/group.uid___A001_X15a2_Xb6a/member.uid___A001_X15a2_Xb6b/calibrated/'
SB_raw_path     = '2021.1.00690.S/science_goal.uid___A001_X15a2_Xb69/group.uid___A001_X15a2_Xb6a/member.uid___A001_X15a2_Xb6d/calibrated/'
step1_path      = 'workflow/step1/'
step2_path      = 'workflow/step2/'
step3_path      = 'workflow/step3/'
step4_path      = 'workflow/step4/'
step5_path      = 'workflow/step5/'
prefix          = 'ABAur'

continuum               = NRAO_path+'ABAur_continuum.bin30s.ms'
line_12CO               = NRAO_path+'ABAur_12CO.bin30s.ms.contsub'
line_13CO               = NRAO_path+'ABAur_13CO.bin30s.ms.contsub'
line_C18O               = NRAO_path+'ABAur_C18O.bin30s.ms.contsub'
line_SO                 = NRAO_path+'ABAur_SO.bin30s.ms.contsub'
line_12CO_wcont         = NRAO_path+'ABAur_12CO.bin30s.ms'
line_13CO_wcont         = NRAO_path+'ABAur_13CO.bin30s.ms'
line_C18O_wcont         = NRAO_path+'ABAur_C18O.bin30s.ms'
line_SO_wcont           = NRAO_path+'ABAur_SO.bin30s.ms'

SB_EB1_msplitcalsource = SB_raw_path+'uid___A002_Xf7ad58_Xd406.ms.split.cal.source'  # 34.0 GB
SB_EB2_msplitcalsource = SB_raw_path+'uid___A002_Xf8f6a9_X15c79.ms.split.cal.source' # 32.3 GB
LB_EB1_msplitcalsource = LB_raw_path+'uid___A002_Xfb8480_X6a5b.ms.split.cal.source'  # 39.5 GB
LB_EB2_msplitcalsource = LB_raw_path+'uid___A002_Xfb8480_X178cb.ms.split.cal.source' # 41.5 GB
LB_EB3_msplitcalsource = LB_raw_path+'uid___A002_Xfbb255_X2286.ms.split.cal.source'  # 30.5 GB
LB_EB4_msplitcalsource = LB_raw_path+'uid___A002_Xfbb9c7_X6731.ms.split.cal.source'  # 41.5 GB
LB_EB5_msplitcalsource = LB_raw_path+'uid___A002_Xfbb9c7_Xd621.ms.split.cal.source'  # 45.4 GB
LB_EB6_msplitcalsource = LB_raw_path+'uid___A002_Xfbb9c7_Xe05e.ms.split.cal.source'  # 39.6 GB

SB_EB1_initcont        = step1_path+prefix+'_SB_EB1_initcont.ms'
SB_EB2_initcont        = step1_path+prefix+'_SB_EB2_initcont.ms'
LB_EB1_initcont        = step1_path+prefix+'_LB_EB1_initcont.ms'
LB_EB2_initcont        = step1_path+prefix+'_LB_EB2_initcont.ms'
LB_EB3_initcont        = step1_path+prefix+'_LB_EB3_initcont.ms'
LB_EB4_initcont        = step1_path+prefix+'_LB_EB4_initcont.ms'
LB_EB5_initcont        = step1_path+prefix+'_LB_EB5_initcont.ms'
LB_EB6_initcont        = step1_path+prefix+'_LB_EB6_initcont.ms'

SB_EB1_initcont_selfcal        = step1_path+prefix+'_SB_EB1_initcont_selfcal.ms'
SB_EB2_initcont_selfcal        = step1_path+prefix+'_SB_EB2_initcont_selfcal.ms'
LB_EB1_initcont_selfcal        = step1_path+prefix+'_LB_EB1_initcont_selfcal.ms'
LB_EB2_initcont_selfcal        = step1_path+prefix+'_LB_EB2_initcont_selfcal.ms'
LB_EB3_initcont_selfcal        = step1_path+prefix+'_LB_EB3_initcont_selfcal.ms'
LB_EB4_initcont_selfcal        = step1_path+prefix+'_LB_EB4_initcont_selfcal.ms'
LB_EB5_initcont_selfcal        = step1_path+prefix+'_LB_EB5_initcont_selfcal.ms'
LB_EB6_initcont_selfcal        = step1_path+prefix+'_LB_EB6_initcont_selfcal.ms'

SB_EB1_initcont_selfcal_shift        = step2_path+prefix+'_SB_EB1_initcont_selfcal_shift.ms'
SB_EB2_initcont_selfcal_shift        = step2_path+prefix+'_SB_EB2_initcont_selfcal_shift.ms'
LB_EB1_initcont_selfcal_shift        = step2_path+prefix+'_LB_EB1_initcont_selfcal_shift.ms'
LB_EB2_initcont_selfcal_shift        = step2_path+prefix+'_LB_EB2_initcont_selfcal_shift.ms'
LB_EB3_initcont_selfcal_shift        = step2_path+prefix+'_LB_EB3_initcont_selfcal_shift.ms'
LB_EB4_initcont_selfcal_shift        = step2_path+prefix+'_LB_EB4_initcont_selfcal_shift.ms'
LB_EB5_initcont_selfcal_shift        = step2_path+prefix+'_LB_EB5_initcont_selfcal_shift.ms'
LB_EB6_initcont_selfcal_shift        = step2_path+prefix+'_LB_EB6_initcont_selfcal_shift.ms'

LB_concat_shifted               = step3_path+prefix+'_LB_concat_shifted.ms'
SB_concat_shifted_contp0        = step3_path+prefix+'_SB_concat_shifted_contp0.ms'
BB_concat_shifted_contp0        = step3_path+prefix+'_BB_concat_shifted_contp0.ms' # BB stands for both baselines

SB_contp1_cal        = step3_path+prefix+'_SB_concat_shifted_contp1.cal'
SB_contp2_cal        = step3_path+prefix+'_SB_concat_shifted_contp2.cal'
SB_contp3_cal        = step3_path+prefix+'_SB_concat_shifted_contp3.cal'
SB_contp4_cal        = step3_path+prefix+'_SB_concat_shifted_contp4.cal'
SB_contp5_cal        = step3_path+prefix+'_SB_concat_shifted_contp5.cal'
SB_contp6_cal        = step3_path+prefix+'_SB_concat_shifted_contp6.cal'
SB_contap_cal        = step3_path+prefix+'_SB_concat_shifted_contap.cal'
BB_contp1_cal        = step3_path+prefix+'_BB_concat_shifted_contp1.cal'
BB_contp2_cal        = step3_path+prefix+'_BB_concat_shifted_contp2.cal'
BB_contp3_cal        = step3_path+prefix+'_BB_concat_shifted_contp3.cal'
BB_contp4_cal        = step3_path+prefix+'_BB_concat_shifted_contp4.cal'
BB_contap_cal        = step3_path+prefix+'_BB_concat_shifted_contap.cal'

final_continuum                 = step3_path+prefix+'_continuum.ms'

line_spws       = np.array([1, 2, 3, 4]) # spws containing emission lines (SO, C18O, 13CO, 12CO)
line_rest_freqs = np.array([2.19949442e11, 2.19560358e11, 2.20398684e11, 2.30538000e11]) # rest frequencies of the emission lines (SO, C18O, 13CO, 12CO)

SB_EB1_initlines        = step4_path+prefix+'_SB_EB1_initlines.ms'
SB_EB2_initlines        = step4_path+prefix+'_SB_EB2_initlines.ms'
LB_EB1_initlines        = step4_path+prefix+'_LB_EB1_initlines.ms'
LB_EB2_initlines        = step4_path+prefix+'_LB_EB2_initlines.ms'
LB_EB3_initlines        = step4_path+prefix+'_LB_EB3_initlines.ms'
LB_EB4_initlines        = step4_path+prefix+'_LB_EB4_initlines.ms'
LB_EB5_initlines        = step4_path+prefix+'_LB_EB5_initlines.ms'
LB_EB6_initlines        = step4_path+prefix+'_LB_EB6_initlines.ms'

SB_EB1_initlines_shift        = step4_path+prefix+'_SB_EB1_initlines_shift.ms'
SB_EB2_initlines_shift        = step4_path+prefix+'_SB_EB2_initlines_shift.ms'
LB_EB1_initlines_shift        = step4_path+prefix+'_LB_EB1_initlines_shift.ms'
LB_EB2_initlines_shift        = step4_path+prefix+'_LB_EB2_initlines_shift.ms'
LB_EB3_initlines_shift        = step4_path+prefix+'_LB_EB3_initlines_shift.ms'
LB_EB4_initlines_shift        = step4_path+prefix+'_LB_EB4_initlines_shift.ms'
LB_EB5_initlines_shift        = step4_path+prefix+'_LB_EB5_initlines_shift.ms'
LB_EB6_initlines_shift        = step4_path+prefix+'_LB_EB6_initlines_shift.ms'

SB_EB1_initlines_SBselfcal        = step4_path+prefix+'_SB_EB1_initlines_SBselfcal.ms'
SB_EB2_initlines_SBselfcal        = step4_path+prefix+'_SB_EB2_initlines_SBselfcal.ms'
SB_EB1_initlines_selfcal        = step4_path+prefix+'_SB_EB1_initlines_selfcal.ms'
SB_EB2_initlines_selfcal        = step4_path+prefix+'_SB_EB2_initlines_selfcal.ms'
LB_EB1_initlines_selfcal        = step4_path+prefix+'_LB_EB1_initlines_selfcal.ms'
LB_EB2_initlines_selfcal        = step4_path+prefix+'_LB_EB2_initlines_selfcal.ms'
LB_EB3_initlines_selfcal        = step4_path+prefix+'_LB_EB3_initlines_selfcal.ms'
LB_EB4_initlines_selfcal        = step4_path+prefix+'_LB_EB4_initlines_selfcal.ms'
LB_EB5_initlines_selfcal        = step4_path+prefix+'_LB_EB5_initlines_selfcal.ms'
LB_EB6_initlines_selfcal        = step4_path+prefix+'_LB_EB6_initlines_selfcal.ms'

SB_EB1_initlines_selfcal_contsub        = step4_path+prefix+'_SB_EB1_initlines_selfcal.ms.contsub'
SB_EB2_initlines_selfcal_contsub        = step4_path+prefix+'_SB_EB2_initlines_selfcal.ms.contsub'
LB_EB1_initlines_selfcal_contsub        = step4_path+prefix+'_LB_EB1_initlines_selfcal.ms.contsub'
LB_EB2_initlines_selfcal_contsub        = step4_path+prefix+'_LB_EB2_initlines_selfcal.ms.contsub'
LB_EB3_initlines_selfcal_contsub        = step4_path+prefix+'_LB_EB3_initlines_selfcal.ms.contsub'
LB_EB4_initlines_selfcal_contsub        = step4_path+prefix+'_LB_EB4_initlines_selfcal.ms.contsub'
LB_EB5_initlines_selfcal_contsub        = step4_path+prefix+'_LB_EB5_initlines_selfcal.ms.contsub'
LB_EB6_initlines_selfcal_contsub        = step4_path+prefix+'_LB_EB6_initlines_selfcal.ms.contsub'

"""######### Parameter choices for Step 1 #########"""

'''For spectral averaging: Spectral windows to be included. spw 0 is the only
pure-continuum spw though'''
cont_spws       = '0, 1, 2, 3, 4'

'''For spectral averaging: Velocity ranges to flag out, about the line'''
velocity_ranges = np.array([np.array([ddisk.disk_dict['v_sys']-3., ddisk.disk_dict['v_sys']+3.]),
                            np.array([ddisk.disk_dict['v_sys']-4., ddisk.disk_dict['v_sys']+4.]),
                            np.array([ddisk.disk_dict['v_sys']-6., ddisk.disk_dict['v_sys']+6.]),
                            np.array([ddisk.disk_dict['v_sys']-10., ddisk.disk_dict['v_sys']+10.])])

'''For spectral averaging: The non-continuum spectral windows have total
bandwidth less than 250 MHz (58.6 MHz), but we include them regardless using
width = number of channels (so, all channels are averaged into one bin)'''
width_array = [8, 960, 960, 1920, 1920] # number of channels to average together

'''For initial per-EB continuum images: Get some estimation of how deeply to clean,
from the weblog. These are overestimations because the weblog combined EBs, but
it works well enough'''
SB_pipeline_cont_cleanthresh           = 0.18, #mJy; from weblog hif_makeimages (aggregate)
SB_pipeline_cont_cleanthresh_perspw    = [0.18, 0.66, 0.57, 0.78, 0.7], #mJy; from weblog hif_makeimages (mfs)
LB_pipeline_cont_cleanthresh           = 0.060, #mJy; from weblog hif_makeimages (aggregate)
LB_pipeline_cont_cleanthresh_perspw    = [0.062, 0.26, 0.25, 0.2, 0.11], #mJy; from weblog hif_makeimages (mfs)

"""######### Parameter choices for Step 4 #########"""

'''For extracting the spectral lines: We average the continuum spw down, and do
no averaging on the line spectral windows'''
line_width_array = [128, 1, 1, 1, 1] # number of channels to average together


"""############### Data dictionary ################"""

data_dict   = {'CANFAR_path' : CANFAR_path,
               'NRAO_path' : NRAO_path,
               'SB_raw_path' : SB_raw_path,
               'LB_raw_path' : LB_raw_path,
               'prefix' : prefix,
               'EBs' : ['SB_EB1', 'SB_EB2', 'LB_EB1', 'LB_EB2', 'LB_EB3', 'LB_EB4', 'LB_EB5', 'LB_EB6'],
               'LB_EBs' : ['LB_EB1', 'LB_EB2', 'LB_EB3', 'LB_EB4', 'LB_EB5', 'LB_EB6'],
               'SB_EBs' : ['SB_EB1', 'SB_EB2'],

               'continuum' : continuum,
               '12CO' : line_12CO,
               '13CO' : line_13CO,
               'C18O' : line_C18O,
               'SO' : line_SO,
               '12CO_wcont' : line_12CO_wcont,
               '13CO_wcont' : line_13CO_wcont,
               'C18O_wcont' : line_C18O_wcont,
               'SO_wcont' : line_SO_wcont,

               'all_selfcaltables' : [SB_contp1_cal, SB_contp2_cal, SB_contp3_cal, SB_contp4_cal, SB_contp5_cal, SB_contp6_cal, SB_contap_cal, BB_contp1_cal, BB_contp2_cal, BB_contp3_cal, BB_contp4_cal, BB_contap_cal],
               'SB_selfcaltables' : [SB_contp1_cal, SB_contp2_cal, SB_contp3_cal, SB_contp4_cal, SB_contp5_cal, SB_contp6_cal, SB_contap_cal],
               'BB_selfcaltables' : [BB_contp1_cal, BB_contp2_cal, BB_contp3_cal, BB_contp4_cal, BB_contap_cal],

               'SB_concat': {'contp0' : SB_concat_shifted_contp0,
                            'refants_list' : 'DV08, DA57, DA52', # first 3 in common to SB_EB1 and SB_EB2!
                            'observations' : ['0', '1'] # this ID's the execution blocks
                            },
               'LB_concat': {'contp0' : LB_concat_shifted,
                            'refants_list' : 'DA58, DA42, DA55, DV03', # scores: 17, 19, 18, 20
                            'observations' : ['0', '1', '2', '3', '4', '5'] # this ID's the execution blocks
                            },
               'BB_concat': {'contp0' : BB_concat_shifted_contp0, # BB stands for both baselines
                            'refants_list' : 'DA58, DA42, DA55, DV03', # scores: 17, 19, 18, 20
                            'observations' : ['0', '1', '2', '3', '4', '5', '6', '7'] # this ID's the execution blocks
                            },
               'continuum' : final_continuum,

               'SB_EB1': {'.ms.split.cal.source' : SB_EB1_msplitcalsource,
                           '_initcont.ms' : SB_EB1_initcont,
                           '_initlines.ms' : SB_EB1_initlines,
                           '_initlines_shift.ms' : SB_EB1_initlines_shift,
                           '_initlines_SBselfcal.ms' : SB_EB1_initlines_SBselfcal,
                           '_initlines_selfcal.ms' : SB_EB1_initlines_selfcal,
                           '_initlines_selfcal.ms.contsub' : SB_EB1_initlines_selfcal_contsub,
                           '_initcont_selfcal.ms' : SB_EB1_initcont_selfcal,
                           '_initcont_selfcal_shift.ms' : SB_EB1_initcont_selfcal_shift,
                           'name' : 'SB_EB1',
                           'field': 'AB_Aur',
                           'line_spws': line_spws,
                           'line_freqs':line_rest_freqs,
                           'velocity_ranges':velocity_ranges,
                           'cont_spws': cont_spws,
                           'width_array': width_array,
                           'line_width_array': line_width_array,
                           'pipeline_cont_cleanthresh': SB_pipeline_cont_cleanthresh,
                           'pipeline_cont_cleanthresh_perspw': SB_pipeline_cont_cleanthresh_perspw,
                           'refants_list':'DV08, DA57, DA52, DV14, DA60',
                           'SB_selfcaltables_spwmap' : [[0,1,2,3,4],  [0,1,2,3,4],   [0,0,0,0,0],   [0,0,0,0,0],   [0,0,0,0,0],    [0,0,0,0,0],   [0,1,2,3,4]],
                           'BB_selfcaltables_spwmap' : [[0,0,0,0,0],  [0,0,0,0,0],   [0,0,0,0,0],   [0,0,0,0,0],   [0,0,0,0,0]]
                          },
               'SB_EB2': {'.ms.split.cal.source' : SB_EB2_msplitcalsource,
                           '_initcont.ms' : SB_EB2_initcont,
                           '_initlines.ms' : SB_EB2_initlines,
                           '_initlines_shift.ms' : SB_EB2_initlines_shift,
                           '_initlines_SBselfcal.ms' : SB_EB2_initlines_SBselfcal,
                           '_initlines_selfcal.ms' : SB_EB2_initlines_selfcal,
                           '_initlines_selfcal.ms.contsub' : SB_EB2_initlines_selfcal_contsub,
                           '_initlines_selfcal.ms.contsub' : SB_EB2_initlines_selfcal_contsub,
                           '_initcont_selfcal.ms' : SB_EB2_initcont_selfcal,
                           '_initcont_selfcal_shift.ms' : SB_EB2_initcont_selfcal_shift,
                           'name' : 'SB_EB2',
                           'field': 'AB_Aur',
                           'line_spws': line_spws,
                           'line_freqs':line_rest_freqs,
                           'velocity_ranges':velocity_ranges,
                           'cont_spws': cont_spws,
                           'width_array': width_array,
                           'line_width_array': line_width_array,
                           'pipeline_cont_cleanthresh': SB_pipeline_cont_cleanthresh,
                           'pipeline_cont_cleanthresh_perspw': SB_pipeline_cont_cleanthresh_perspw,
                           'refants_list':'DV08, DA57, DA52, DA60, DV05',
                           'SB_selfcaltables_spwmap' : [[5,6,7,8,9],  [5,6,7,8,9],   [5,5,5,5,5],   [5,5,5,5,5],   [5,5,5,5,5],    [5,5,5,5,5],   [5,6,7,8,9]],
                           'BB_selfcaltables_spwmap' : [[5,5,5,5,5],  [5,5,5,5,5],   [5,5,5,5,5],   [5,5,5,5,5],   [5,5,5,5,5]]
                          },
               'LB_EB1': {'.ms.split.cal.source' : LB_EB1_msplitcalsource,
                           '_initcont.ms' : LB_EB1_initcont,
                           '_initlines.ms' : LB_EB1_initlines,
                           '_initlines_shift.ms' : LB_EB1_initlines_shift,
                           '_initlines_selfcal.ms' : LB_EB1_initlines_selfcal,
                           '_initlines_selfcal.ms.contsub' : LB_EB1_initlines_selfcal_contsub,
                           '_initcont_selfcal.ms' : LB_EB1_initcont_selfcal,
                           '_initcont_selfcal_shift.ms' : LB_EB1_initcont_selfcal_shift,
                           'name' : 'LB_EB1',
                           'field' : 'AB_Aur',
                           'line_spws': line_spws,
                           'line_freqs':line_rest_freqs,
                           'velocity_ranges':velocity_ranges,
                           'cont_spws': cont_spws,
                           'width_array': width_array,
                           'line_width_array': line_width_array,
                           'pipeline_cont_cleanthresh': LB_pipeline_cont_cleanthresh,
                           'pipeline_cont_cleanthresh_perspw': LB_pipeline_cont_cleanthresh_perspw,
                           'refants_list':'DV03, DA58, DV01, DA42, DA55',
                           'BB_selfcaltables_spwmap' : [[10,10,10,10,10],[10,10,10,10,10],[10,10,10,10,10],[10,10,10,10,10],[10,10,10,10,10]]
                          },
               'LB_EB2': {'.ms.split.cal.source' : LB_EB2_msplitcalsource,
                           '_initcont.ms' : LB_EB2_initcont,
                           '_initlines.ms' : LB_EB2_initlines,
                           '_initlines_shift.ms' : LB_EB2_initlines_shift,
                           '_initlines_selfcal.ms' : LB_EB2_initlines_selfcal,
                           '_initlines_selfcal.ms.contsub' : LB_EB2_initlines_selfcal_contsub,
                           '_initcont_selfcal.ms' : LB_EB2_initcont_selfcal,
                           '_initcont_selfcal_shift.ms' : LB_EB2_initcont_selfcal_shift,
                           'name' : 'LB_EB2',
                           'field' : 'AB_Aur',
                           'line_spws': line_spws,
                           'line_freqs':line_rest_freqs,
                           'velocity_ranges':velocity_ranges,
                           'cont_spws': cont_spws,
                           'width_array': width_array,
                           'line_width_array': line_width_array,
                           'pipeline_cont_cleanthresh': LB_pipeline_cont_cleanthresh,
                           'pipeline_cont_cleanthresh_perspw': LB_pipeline_cont_cleanthresh_perspw,
                           'refants_list':'DV04, DA55, DV03, DA58, DA42',
                           'BB_selfcaltables_spwmap' : [[15,15,15,15,15],[15,15,15,15,15],[15,15,15,15,15],[15,15,15,15,15],[15,15,15,15,15]]
                          },
               'LB_EB3': {'.ms.split.cal.source' : LB_EB3_msplitcalsource,
                           '_initcont.ms' : LB_EB3_initcont,
                           '_initlines.ms' : LB_EB3_initlines,
                           '_initlines_shift.ms' : LB_EB3_initlines_shift,
                           '_initlines_selfcal.ms' : LB_EB3_initlines_selfcal,
                           '_initlines_selfcal.ms.contsub' : LB_EB3_initlines_selfcal_contsub,
                           '_initcont_selfcal.ms' : LB_EB3_initcont_selfcal,
                           '_initcont_selfcal_shift.ms' : LB_EB3_initcont_selfcal_shift,
                           'name' : 'LB_EB3',
                           'field' : 'AB_Aur',
                           'line_spws': line_spws,
                           'line_freqs':line_rest_freqs,
                           'velocity_ranges':velocity_ranges,
                           'cont_spws': cont_spws,
                           'width_array': width_array,
                           'line_width_array': line_width_array,
                           'pipeline_cont_cleanthresh': LB_pipeline_cont_cleanthresh,
                           'pipeline_cont_cleanthresh_perspw': LB_pipeline_cont_cleanthresh_perspw,
                           'refants_list':'DA55, DA42, DA58, DV03, DV04',
                           'BB_selfcaltables_spwmap' : [[20,20,20,20,20],[20,20,20,20,20],[20,20,20,20,20],[20,20,20,20,20],[20,20,20,20,20]]
                          },
               'LB_EB4': {'.ms.split.cal.source' : LB_EB4_msplitcalsource,
                           '_initcont.ms' : LB_EB4_initcont,
                           '_initlines.ms' : LB_EB4_initlines,
                           '_initlines_shift.ms' : LB_EB4_initlines_shift,
                           '_initlines_selfcal.ms' : LB_EB4_initlines_selfcal,
                           '_initlines_selfcal.ms.contsub' : LB_EB4_initlines_selfcal_contsub,
                           '_initcont_selfcal.ms' : LB_EB4_initcont_selfcal,
                           '_initcont_selfcal_shift.ms' : LB_EB4_initcont_selfcal_shift,
                           'name' : 'LB_EB4',
                           'field' : 'AB_Aur',
                           'line_spws': line_spws,
                           'line_freqs':line_rest_freqs,
                           'velocity_ranges':velocity_ranges,
                           'cont_spws': cont_spws,
                           'width_array': width_array,
                           'line_width_array': line_width_array,
                           'pipeline_cont_cleanthresh': LB_pipeline_cont_cleanthresh,
                           'pipeline_cont_cleanthresh_perspw': LB_pipeline_cont_cleanthresh_perspw,
                           'refants_list':'DV03, DA58, DV01, DA42, DA55',
                           'BB_selfcaltables_spwmap' : [[25,25,25,25,25],[25,25,25,25,25],[25,25,25,25,25],[25,25,25,25,25],[25,25,25,25,25]]
                          },
               'LB_EB5': {'.ms.split.cal.source' : LB_EB5_msplitcalsource,
                           '_initcont.ms' : LB_EB5_initcont,
                           '_initlines.ms' : LB_EB5_initlines,
                           '_initlines_shift.ms' : LB_EB5_initlines_shift,
                           '_initlines_selfcal.ms' : LB_EB5_initlines_selfcal,
                           '_initlines_selfcal.ms.contsub' : LB_EB5_initlines_selfcal_contsub,
                           '_initcont_selfcal.ms' : LB_EB5_initcont_selfcal,
                           '_initcont_selfcal_shift.ms' : LB_EB5_initcont_selfcal_shift,
                           'name' : 'LB_EB5',
                           'field' : 'AB_Aur',
                           'line_spws': line_spws,
                           'line_freqs':line_rest_freqs,
                           'velocity_ranges':velocity_ranges,
                           'cont_spws': cont_spws,
                           'width_array': width_array,
                           'line_width_array': line_width_array,
                           'pipeline_cont_cleanthresh': LB_pipeline_cont_cleanthresh,
                           'pipeline_cont_cleanthresh_perspw': LB_pipeline_cont_cleanthresh_perspw,
                           'refants_list':'DA60, DA42, DA58, DA55, DV03',
                           'BB_selfcaltables_spwmap' : [[30,30,30,30,30],[30,30,30,30,30],[30,30,30,30,30],[30,30,30,30,30],[30,30,30,30,30]]
                          },
               'LB_EB6': {'.ms.split.cal.source' : LB_EB6_msplitcalsource,
                           '_initcont.ms' : LB_EB6_initcont,
                           '_initlines.ms' : LB_EB6_initlines,
                           '_initlines_shift.ms' : LB_EB6_initlines_shift,
                           '_initlines_selfcal.ms' : LB_EB6_initlines_selfcal,
                           '_initlines_selfcal.ms.contsub' : LB_EB6_initlines_selfcal_contsub,
                           '_initcont_selfcal.ms' : LB_EB6_initcont_selfcal,
                           '_initcont_selfcal_shift.ms' : LB_EB6_initcont_selfcal_shift,
                           'name' : 'LB_EB6',
                           'field' : 'AB_Aur',
                           'line_spws': line_spws,
                           'line_freqs':line_rest_freqs,
                           'velocity_ranges':velocity_ranges,
                           'cont_spws': cont_spws,
                           'width_array': width_array,
                           'line_width_array': line_width_array,
                           'pipeline_cont_cleanthresh': LB_pipeline_cont_cleanthresh,
                           'pipeline_cont_cleanthresh_perspw': LB_pipeline_cont_cleanthresh_perspw,
                           'refants_list':'DA55, DA42, DV04, DV25, DA58',
                           'BB_selfcaltables_spwmap' : [[35,35,35,35,35],[35,35,35,35,35],[35,35,35,35,35],[35,35,35,35,35],[35,35,35,35,35]]
                          }
               }
