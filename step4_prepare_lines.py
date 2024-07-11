'''
This script was written for CASA 6.4.1.12
(on NRAO computers: casa -r 6.4.1-12-pipeline-2022.2.0.64)

ALMA Program ID: 2021.1.00690.S (PI: R. Dong)

Datasets (in order of date observed):
AB_Aur_a_06_TM2: member.uid___A001_X15a2_Xb6d (short baseline, 2 EB's)
    SB1: uid___A002_Xf7ad58_Xd406
         Observed 2022-04-19 7:25:22 - 2022-04-19  8:20:05
            04:55:45.852430 +30.33.03.75455 # RA and Dec in ICRS from Field section of listobs
         [uid___A002_Xf8d822_Xa976 (not used, QA fail)]
    SB2: uid___A002_Xf8f6a9_X15c79
         Observed 2022-05-15 18:05:13 - 2022-05-15 18:49:47
            04:55:45.852624 +30.33.03.75287
AB_Aur_a_06_TM1: member.uid___A001_X15a2_Xb6b (long baseline, 6 EB's)
    LB1: uid___A002_Xfb8480_X6a5b
         Observed 2022-07-17 14:07:17 - 2022-07-17 15:25:26
            04:55:45.853150 +30.33.03.74953
    LB2: uid___A002_Xfb8480_X178cb
         Observed 2022-07-19 12:02:10 - 2022-07-19 13:21:27
            04:55:45.853163 +30.33.03.74944
    LB3: uid___A002_Xfbb255_X2286
         Observed 2022-07-20 14:03:55 - 2022-07-20 14:59:55
            04:55:45.853170 +30.33.03.74939
    LB4: uid___A002_Xfbb9c7_X6731
         Observed 2022-07-21 14:01:35 - 2022-07-21 15:20:52
            04:55:45.853177 +30.33.03.74934
    LB5: uid___A002_Xfbb9c7_Xd621
         Observed 2022-07-22 11:30:28 - 2022-07-22 12:49:58
            04:55:45.853183 +30.33.03.74930
    LB6: uid___A002_Xfbb9c7_Xe05e
         Observed 2022-07-22 12:59:16 - 2022-07-22 14:18:16
            04:55:45.853183 +30.33.03.74929

All 8 EB's have the same spectral setup:
spw23	continuum
spw25	SO_3__v=0_6(5)-5(4)
spw27	C18O_2-1
spw29	13CO_v=0_2-1 **science goal window
spw31	CO_v=0_2-1

reducer: J. Speedie
'''

""" Starting matter """
import os
execfile('dictionary_disk.py') # loads disk_dict
execfile('dictionary_data.py') # loads data_dict
execfile('step1_utils.py') # loads multiple functions
import sys
sys.path.append(data_dict['NRAO_path']+'analysis_scripts') # path to analysis_scripts/
import analysisUtils as au

'''
Note on where we begin:
We downloaded the entire data set from the ALMA Archive and restored the
pipeline calibration using scriptForPI.py to create .ms files for each execution.
We then split out the science spectral windows, using the code:
    # in CASA, in e.g., member.uid___A001_X15a2_Xb6d/calibrated/
    import glob
    vislist = glob.glob('*[!_t].ms')  # should be ['uid___A002_Xf7ad58_Xd406.ms', 'uid___A002_Xf8f6a9_X15c79.ms']
    for myvis in vislist:
        msmd.open(myvis)
        targetspws = msmd.spwsforintent('OBSERVE_TARGET*')
        sciencespws = []
        for myspw in targetspws:
            if msmd.nchan(myspw)>4:
                sciencespws.append(myspw)
        sciencespws = ','.join(map(str,sciencespws))
        msmd.close()
        split(vis=myvis,outputvis=myvis+'.split.cal',spw=sciencespws)
So, the .ms.split.cal measurement sets contain the science target AND the calibrators.
We then saved the original flags, using the code:
    # in CASA, in e.g., member.uid___A001_X15a2_Xb6d/calibrated/
    import glob
    vislist=glob.glob('*.ms.split.cal')
    for vis in vislist:
        flagmanager(vis=vis,
                    mode='save',
                    versionname='original_flags')
After that, we split off the science target, and saved those flags, using the code:
    for vis in vislist:
        outputvis = vis+'.source'
        split(vis=vis,
          intent='*TARGET*', # split off the target sources by intent
          outputvis=outputvis,
          datacolumn='data')
        os.system('rm -rf ' + outputvis + '.flagversions')
        flagmanager(vis=outputvis,
                    mode='save',
                    versionname='original_flags')

In summary, this script starts with 8 separate measurement sets (.ms.split.cal.target)
which contain only the science spectral windows for only the science target.
'''

"""
######################################################
############# Additional manual flagging #############
######################################################
"""
'''
Possible things to flag have been identified, but opting for now to rely on
self calibration to fix them.
'''

"""
######################################################
######## SPLIT OUT THE LINE SPECTRAL WINDOWS #########
######################################################
"""

print('\nSplitting out the line spectral windows, with the continuum spw averaged down...')
for EB in data_dict['EBs']:
    inputvis            = data_dict['NRAO_path']+data_dict[EB]['.ms.split.cal.source']
    outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines.ms']
    print('inputvis = ', inputvis)
    print('outputvis = ', outputvis)

    os.system('rm -rf '+outputvis)
    split(vis         = inputvis,
          field       = 'AB_Aur',
          spw         = data_dict[EB]['cont_spws'], # in hindsight, poorly named; this is all spws
          outputvis   = outputvis,
          width       = data_dict[EB]['line_width_array'], # averages the cont spw (0) down to 1 channel
          datacolumn  = 'data',
          keepflags   = False)
    print("Spectrally averaged continuum dataset saved to %s" % outputvis)
    listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')


"""
######################################################
###### ALIGN EB'S TO THE REFERENCE EB (LB EB1) #######
######################################################
"""
'''
Here, we have to take a detour into modular casa.
Please go run the script called step4_detour.py using modularcasa venv.
Then come back and continue with the initlines_shift.ms's that that script creates.
'''


"""
####################################################################################
###### APPLY SAME SELF CALIBRATION SOLUTIONS AS DETERMINED FOR THE CONTINUUM #######
####################################################################################
"""

'''
First, we apply the SB self-cal solutions to the SB execution blocks.

SB_contp1.cal was gained on SB 1 and SB2 with combine='', so it contains spws [0,1,2,3,4,5,6,7,8,9]
If we want to apply it to SB 1 (lines) only, spwmap=[0,1,2,3,4]
                   and to SB 2 (lines) only, spwmap=[5,6,7,8,9]

SB_contp2.cal was gained on SB 1 and SB2 with combine='', so it contains spws [0,1,2,3,4,5,6,7,8,9]
If we want to apply it to SB 1 (lines) only, spwmap=[0,1,2,3,4]
                   and to SB 2 (lines) only, spwmap=[5,6,7,8,9]

SB_contp3.cal was gained on SB 1 and SB2 with combine=spw, so it contains spws [0,0,0,0,0,5,5,5,5,5]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]

SB_contp4.cal was gained on SB 1 and SB2 with combine=spw, so it contains spws [0,0,0,0,0,5,5,5,5,5]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]

SB_contp5.cal was gained on SB 1 and SB2 with combine=spw, so it contains spws [0,0,0,0,0,5,5,5,5,5]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]

SB_contp6.cal was gained on SB 1 and SB2 with combine=spw, so it contains spws [0,0,0,0,0,5,5,5,5,5]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]

SB_contap.cal was gained on SB 1 and SB2 with combine='', so it contains spws [0,1,2,3,4,5,6,7,8,9]
If we want to apply it to SB 1 (lines) only, spwmap=[0,1,2,3,4]
                   and to SB 2 (lines) only, spwmap=[5,6,7,8,9]


To apply all SB caltables to SB1, gaintable=[SB_contp1.cal, SB_contp2.cal, SB_contp3.cal, SB_contp4.cal, SB_contp5.cal, SB_contp6.cal, SB_contap.cal]
                                  and spwmap=[[0,1,2,3,4],  [0,1,2,3,4],   [0,0,0,0,0],   [0,0,0,0,0],   [0,0,0,0,0],    [0,0,0,0,0],   [0,1,2,3,4]]
To apply all SB caltables to SB2, gaintable=[SB_contp1.cal, SB_contp2.cal, SB_contp3.cal, SB_contp4.cal, SB_contp5.cal, SB_contp6.cal, SB_contap.cal]
                                  and spwmap=[[5,6,7,8,9],  [5,6,7,8,9],   [5,5,5,5,5],   [5,5,5,5,5],   [5,5,5,5,5],    [5,5,5,5,5],   [5,6,7,8,9]]



BB_contp1.cal was gained on SB1,SB2,LB1,LB2,LB3,LB4,LB5,LB6 with combine='spw', so it contains spws [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]
                   and to LB 1 (lines) only, spwmap=[10,10,10,10,10]
                   and to LB 2 (lines) only, spwmap=[15,15,15,15,15]
                   and to LB 3 (lines) only, spwmap=[20,20,20,20,20]
                   and to LB 4 (lines) only, spwmap=[25,25,25,25,25]
                   and to LB 5 (lines) only, spwmap=[30,30,30,30,30]
                   and to LB 6 (lines) only, spwmap=[35,35,35,35,35]

BB_contp2.cal was gained on SB1,SB2,LB1,LB2,LB3,LB4,LB5,LB6 with combine='spw', so it contains spws [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]
                   and to LB 1 (lines) only, spwmap=[10,10,10,10,10]
                   and to LB 2 (lines) only, spwmap=[15,15,15,15,15]
                   and to LB 3 (lines) only, spwmap=[20,20,20,20,20]
                   and to LB 4 (lines) only, spwmap=[25,25,25,25,25]
                   and to LB 5 (lines) only, spwmap=[30,30,30,30,30]
                   and to LB 6 (lines) only, spwmap=[35,35,35,35,35]

BB_contp3.cal was gained on SB1,SB2,LB1,LB2,LB3,LB4,LB5,LB6 with combine='spw', so it contains spws [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]
                   and to LB 1 (lines) only, spwmap=[10,10,10,10,10]
                   and to LB 2 (lines) only, spwmap=[15,15,15,15,15]
                   and to LB 3 (lines) only, spwmap=[20,20,20,20,20]
                   and to LB 4 (lines) only, spwmap=[25,25,25,25,25]
                   and to LB 5 (lines) only, spwmap=[30,30,30,30,30]
                   and to LB 6 (lines) only, spwmap=[35,35,35,35,35]

BB_contp4.cal was gained on SB1,SB2,LB1,LB2,LB3,LB4,LB5,LB6 with combine='spw', so it contains spws [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]
                   and to LB 1 (lines) only, spwmap=[10,10,10,10,10]
                   and to LB 2 (lines) only, spwmap=[15,15,15,15,15]
                   and to LB 3 (lines) only, spwmap=[20,20,20,20,20]
                   and to LB 4 (lines) only, spwmap=[25,25,25,25,25]
                   and to LB 5 (lines) only, spwmap=[30,30,30,30,30]
                   and to LB 6 (lines) only, spwmap=[35,35,35,35,35]

BB_contap.cal was gained on SB1,SB2,LB1,LB2,LB3,LB4,LB5,LB6 with combine='spw', so it contains spws [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35]
If we want to apply it to SB 1 (lines) only, spwmap=[0,0,0,0,0]
                   and to SB 2 (lines) only, spwmap=[5,5,5,5,5]
                   and to LB 1 (lines) only, spwmap=[10,10,10,10,10]
                   and to LB 2 (lines) only, spwmap=[15,15,15,15,15]
                   and to LB 3 (lines) only, spwmap=[20,20,20,20,20]
                   and to LB 4 (lines) only, spwmap=[25,25,25,25,25]
                   and to LB 5 (lines) only, spwmap=[30,30,30,30,30]
                   and to LB 6 (lines) only, spwmap=[35,35,35,35,35]


To apply all BB caltables to SB1, gaintable=[BB_contp1.cal, BB_contp2.cal, BB_contp3.cal, BB_contp4.cal, BB_contap.cal]
                                  and spwmap=[[0,0,0,0,0],  [0,0,0,0,0],   [0,0,0,0,0],   [0,0,0,0,0],   [0,0,0,0,0]]
To apply all BB caltables to SB2, gaintable=[BB_contp1.cal, BB_contp2.cal, BB_contp3.cal, BB_contp4.cal, BB_contap.cal]
                                  and spwmap=[[5,5,5,5,5],  [5,5,5,5,5],   [5,5,5,5,5],   [5,5,5,5,5],   [5,5,5,5,5]]
To apply all BB caltables to LB1, gaintable=[BB_contp1.cal,    BB_contp2.cal,   BB_contp3.cal,   BB_contp4.cal,   BB_contap.cal]
                                  and spwmap=[[10,10,10,10,10],[10,10,10,10,10],[10,10,10,10,10],[10,10,10,10,10],[10,10,10,10,10]]
To apply all BB caltables to LB2, gaintable=[BB_contp1.cal,    BB_contp2.cal,   BB_contp3.cal,   BB_contp4.cal,   BB_contap.cal]
                                  and spwmap=[[15,15,15,15,15],[15,15,15,15,15],[15,15,15,15,15],[15,15,15,15,15],[15,15,15,15,15]]
To apply all BB caltables to LB3, gaintable=[BB_contp1.cal,    BB_contp2.cal,   BB_contp3.cal,   BB_contp4.cal,   BB_contap.cal]
                                  and spwmap=[[20,20,20,20,20],[20,20,20,20,20],[20,20,20,20,20],[20,20,20,20,20],[20,20,20,20,20]]
To apply all BB caltables to LB4, gaintable=[BB_contp1.cal,    BB_contp2.cal,   BB_contp3.cal,   BB_contp4.cal,   BB_contap.cal]
                                  and spwmap=[[25,25,25,25,25],[25,25,25,25,25],[25,25,25,25,25],[25,25,25,25,25],[25,25,25,25,25]]
To apply all BB caltables to LB5, gaintable=[BB_contp1.cal,    BB_contp2.cal,   BB_contp3.cal,   BB_contp4.cal,   BB_contap.cal]
                                  and spwmap=[[30,30,30,30,30],[30,30,30,30,30],[30,30,30,30,30],[30,30,30,30,30],[30,30,30,30,30]]
To apply all BB caltables to LB6, gaintable=[BB_contp1.cal,    BB_contp2.cal,   BB_contp3.cal,   BB_contp4.cal,   BB_contap.cal]
                                  and spwmap=[[35,35,35,35,35],[35,35,35,35,35],[35,35,35,35,35],[35,35,35,35,35],[35,35,35,35,35]]
'''


"""
####################################################
###### APPLY THE SB CALTABLES TO THE SB DATA #######
####################################################
"""
# this takes ~1hr per EB (so 2hr total)
for EB in data_dict['SB_EBs']:
    inputvis            = data_dict['NRAO_path']+data_dict[EB]['_initlines_shift.ms']
    os.system('cp -r '+inputvis+' '+inputvis.replace('.ms', '.keepsafe.ms'))
    outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_SBselfcal.ms'] # half selfcal
    dir4lasts           = outputvis.replace('.ms', '/')
    os.system('mkdir '+dir4lasts)

    vis         = inputvis
    gaintable   = data_dict['SB_selfcaltables']
    spw         = ''
    spwmap      = data_dict[EB]['SB_selfcaltables_spwmap']
    applymode   = 'calflag'
    interp      = ['linearPD' for gaintab in gaintable]
    calwt       = True

    applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, applymode=applymode, calwt=calwt, spwmap=spwmap)

    os.system('rm -rf ' + outputvis)
    split(vis=vis, outputvis=outputvis, datacolumn='corrected', keepflags=False)
    listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')

    for spwi in data_dict[EB]['line_spws']:
        plotms(vis=inputvis, yaxis='amp', xaxis='UVdist', avgtime='1e8', spw=str(spwi), plotfile=inputvis+'_amp-vs-UVdist_spw'+str(spwi)+'.png', showgui=False, ydatacolumn='data', overwrite=True)
        plotms(vis=outputvis, yaxis='amp', xaxis='UVdist', avgtime='1e8', spw=str(spwi), plotfile=outputvis+'_amp-vs-UVdist_spw'+str(spwi)+'.png', showgui=False, ydatacolumn='data', overwrite=True)

    os.system('mv *.last '+dir4lasts)
    os.system('rm -rf ' + inputvis)


"""
####################################################
###### APPLY THE BB CALTABLES TO THE SB DATA #######
####################################################
"""

# this takes ~40min per EB (so 1hr20min total)
for EB in data_dict['SB_EBs']:
    inputvis            = data_dict['NRAO_path']+data_dict[EB]['_initlines_SBselfcal.ms'] # start with the "half selfcal" (SBcaltable-selfcal'ed SB data)
    os.system('cp -r '+inputvis+' '+inputvis.replace('.ms', '.keepsafe.ms'))
    outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms'] # full selfcal
    dir4lasts           = outputvis.replace('.ms', '/')
    os.system('mkdir '+dir4lasts)

    vis         = inputvis
    gaintable   = data_dict['BB_selfcaltables']
    spw         = ''
    spwmap      = data_dict[EB]['BB_selfcaltables_spwmap']
    applymode   = 'calonly' # 0% of data flagged
    interp      = ['linearPD' for gaintab in gaintable]
    calwt       = True

    applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, applymode=applymode, calwt=calwt, spwmap=spwmap)

    os.system('rm -rf ' + outputvis)
    split(vis=vis, outputvis=outputvis, datacolumn='corrected', keepflags=False, timebin='30s') # time average now to save space
    listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')

    for spwi in data_dict[EB]['line_spws']:
        plotms(vis=inputvis, yaxis='amp', xaxis='UVdist', avgtime='1e8', spw=str(spwi), plotfile=inputvis+'_amp-vs-UVdist_spw'+str(spwi)+'.png', showgui=False, ydatacolumn='data', overwrite=True)
        plotms(vis=outputvis, yaxis='amp', xaxis='UVdist', avgtime='1e8', spw=str(spwi), plotfile=outputvis+'_amp-vs-UVdist_spw'+str(spwi)+'.png', showgui=False, ydatacolumn='data', overwrite=True)

    os.system('mv *.last '+dir4lasts)
    os.system('rm -rf ' + inputvis)




"""
####################################################
###### APPLY THE BB CALTABLES TO THE LB DATA #######
####################################################
"""
# 1hr 10min per EB
for EB in data_dict['LB_EBs']:
    print('/n--> Working on ', EB)
    inputvis            = data_dict['NRAO_path']+data_dict[EB]['_initlines_shift.ms']
    # os.system('cp -r '+inputvis+' '+inputvis.replace('.ms', '.keepsafe.ms')) # too much space
    outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms']
    dir4lasts           = outputvis.replace('.ms', '/')
    os.system('mkdir '+dir4lasts)

    vis         = inputvis
    gaintable   = data_dict['BB_selfcaltables']
    spw         = ''
    spwmap      = data_dict[EB]['BB_selfcaltables_spwmap']
    applymode   = 'calonly' # 0% of data will be flagged
    interp      = ['linearPD' for gaintab in gaintable]
    calwt       = True

    applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, applymode=applymode, calwt=calwt, spwmap=spwmap)

    os.system('rm -rf ' + outputvis)
    split(vis=vis, outputvis=outputvis, datacolumn='corrected', keepflags=False, timebin='30s') # time average now to save space
    listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')

    for spwi in data_dict[EB]['line_spws']:
        plotms(vis=inputvis, yaxis='amp', xaxis='UVdist', avgtime='1e8', spw=str(spwi), plotfile=inputvis+'_amp-vs-UVdist_spw'+str(spwi)+'.png', showgui=False, ydatacolumn='data', overwrite=True)
        plotms(vis=outputvis, yaxis='amp', xaxis='UVdist', avgtime='1e8', spw=str(spwi), plotfile=outputvis+'_amp-vs-UVdist_spw'+str(spwi)+'.png', showgui=False, ydatacolumn='data', overwrite=True)

    os.system('mv *.last '+dir4lasts)
    os.system('rm -rf ' + inputvis)


"""
####################################################
########## PERFORM CONTINUUM SUBTRACTION ###########
####################################################
"""

for EB in data_dict['EBs']: # this takes 1hr10min per LB EB
    inputvis            = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms']
    outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms.contsub']
    spw                 = '1, 2, 3, 4' # only line spws; they will be renumbered 0,1,2,3
    combine             = '' # break at scan, field, and spw
    fitorder            = 1 # order of polynomial fit; default is 0; DSHARP uses 1 so we copy them
    solint              = 'int' # Timescale for per-baseline fit. default (recommended): ‘int’, i.e. no time averaging, do a fit for each integration and let the noisy fits average out in the image.continuum
    excludechans        = False # the channels in fitspw will be used to fit the continuum
    want_cont           = False # we don't need another ms to hold the continuum estimate

    flagchannels_string = get_flagchannels(ms_file          = inputvis,
                                           ms_dict          = data_dict[EB],
                                           velocity_range   = data_dict[EB]['velocity_ranges'])
    fitspw = au.invertChannelRanges(flagchannels_string, vis=inputvis)
    # Flagchannels input string for SB_EB1: '1:517~589,       2:504~600,       3:961~1251,        4:556~1059'
    #                               fitspw:  1:0~516;590~959, 2:0~503;601~959, 3:0~960;1252~1919, 4:0~555;1060~1919
    # Flagchannels input string for SB_EB2: '1:517~589,       2:504~600,       3:962~1251,        4:554~1058'
    #                               fitspw:  1:0~516;590~959, 2:0~503;601~959, 3:0~961;1252~1919, 4:0~553;1059~1919
    # Flagchannels input string for LB_EB1: '1:516~589,       2:504~600,       3:960~1249,        4:553~1057'
    #                               fitspw:  1:0~515;590~959, 2:0~503;601~959, 3:0~959;1250~1919, 4:0~552;1058~1919
    # Flagchannels input string for LB_EB2: '1:517~589,       2:504~600,       3:960~1249,        4:553~1057'
    #                               fitspw:  1:0~516;590~959, 2:0~503;601~959, 3:0~959;1250~1919, 4:0~552;1058~1919
    # Flagchannels input string for LB_EB3: '1:517~589,       2:505~601,       3:960~1249,        4:553~1057'
    #                               fitspw:  1:0~516;590~959, 2:0~504;602~959, 3:0~959;1250~1919, 4:0~552;1058~1919
    # Flagchannels input string for LB_EB4: '1:516~589,       2:504~600,       3:960~1249,        4:553~1057'
    #                               fitspw:  1:0~515;590~959, 2:0~503;601~959, 3:0~959;1250~1919, 4:0~552;1058~1919
    # Flagchannels input string for LB_EB5: '1:516~588,       2:504~600,       3:960~1250,        4:553~1057'
    #                               fitspw:  1:0~515;589~959, 2:0~503;601~959, 3:0~959;1251~1919, 4:0~552;1058~1919
    # Flagchannels input string for LB_EB6: '1:517~589,       2:505~601,       3:960~1249,        4:553~1057'
    #                               fitspw:  1:0~516;590~959, 2:0~504;602~959, 3:0~959;1250~1919, 4:0~552;1058~1919

    os.system('rm -rf ' + outputvis)
    uvcontsub(vis=inputvis, spw=spw, combine=combine, fitorder=fitorder,
                       solint=solint, excludechans=excludechans, want_cont=want_cont,
                       fitspw=fitspw)
    listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')
    for spwi in data_dict[EB]['line_spws']:
        plotms(vis=inputvis, yaxis='amp', xaxis='channel', avgtime='1e8', avgspw=False, avgscan=True, avgbaseline=True, spw=str(spwi), plotfile=inputvis+'_amp-vs-channel_spw'+str(spwi)+'.png', showgui=False, overwrite=True)
        plotms(vis=outputvis, yaxis='amp', xaxis='channel', avgtime='1e8', avgspw=False, avgscan=True, avgbaseline=True, spw=str(spwi-1), plotfile=outputvis+'_amp-vs-channel_spw'+str(spwi)+'.png', showgui=False, overwrite=True)




"""
#############################################################
###################### SPLIT OUT SPWS #######################
#############################################################
"""

""" Note we originally performed velocity regridding with cvel2 at this time.
But as we later would like to image on an arbitrary velocity regrid, we will not
do any regridding. """

""" First do non-continuum-subtracted ms's (contains spws 0,1,2,3,4) """
for EB in data_dict['EBs']:
    for i,spwi in enumerate(['1', '2', '3', '4']):
        inputvis            = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms']
        outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms'].replace('.ms', '_spw'+str(i)+'.ms') # want i instead of spwi to align with contsub indexing

        os.system('rm -rf ' + outputvis)
        split(vis=inputvis, outputvis=outputvis, spw=spwi, datacolumn='data') # here, spws are indexed 1,2,3,4
        listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')

""" Second do continuum-subtracted ms's (contains spws 0,1,2,3) """
for EB in data_dict['EBs']:
    for i,spwi in enumerate(['0', '1', '2', '3']):
        inputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms.contsub']
        outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms.contsub'].replace('.ms', '_spw'+spwi+'.ms') # want spwi now

        os.system('rm -rf ' + outputvis)
        split(vis=inputvis, outputvis=outputvis, spw=spwi, datacolumn='data') # here, spws are indexed 0,1,2,3
        listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')


"""
####################################################################################
########### COMBINE FINAL MEASUREMENT SETS, ONE FOR EACH LINE, AND SAVE ############
####################################################################################
"""
molecules = ['SO', 'C18O', '13CO', '12CO']

""" First do non-continuum-subtracted ms's """
for spwi,molecule in enumerate(molecules):
    ms_list_to_concatenate = []
    for EB in data_dict['EBs']:
        ms_list_to_concatenate.append(data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms'].replace('.ms', '_spw'+str(spwi)+'.ms'))

    print('For molecule: ', molecule)
    print('Combining spectral window '+str(spwi)+' across these measurement sets:', ms_list_to_concatenate)

    final_line_ms = 'ABAur_'+molecule+'.bin30s.ms'

    os.system('rm -rf %s*' % final_line_ms)
    concat(vis          = ms_list_to_concatenate,
           concatvis    = final_line_ms,
           dirtol       = '0.1arcsec',
           copypointing = False)
    listobs(vis=final_line_ms, listfile=final_line_ms+'.listobs.txt')
    os.system('tar cvzf backups/' + final_line_ms+'.tgz ' + final_line_ms)


""" Second do continuum-subtracted ms's (contains spws 0,1,2,3) """
for spwi,molecule in enumerate(molecules):
    ms_list_to_concatenate = []
    for EB in data_dict['EBs']:
        ms_list_to_concatenate.append(data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms'].replace('.ms', '_spw'+str(spwi)+'.ms')+'.contsub')

    print('For molecule: ', molecule)
    print('Combining spectral window '+str(spwi)+' across these measurement sets:', ms_list_to_concatenate)

    final_line_ms = 'ABAur_'+molecule+'.bin30s.ms.contsub'

    os.system('rm -rf %s*' % final_line_ms)
    concat(vis          = ms_list_to_concatenate,
           concatvis    = final_line_ms,
           dirtol       = '0.1arcsec',
           copypointing = False)
    listobs(vis=final_line_ms, listfile=final_line_ms+'.listobs.txt')
    os.system('tar cvzf backups/' + final_line_ms+'.tgz ' + final_line_ms)

os.system('mv ABAur*.tgz backups/')




# sys.exit()



# BELOW THIS IS A REPEAT OF THE LAST TWO STEPS, BUT WITH CVEL2 APPLIED
# ORIGINALLY PERFORMED ON DECEMBER 19 2022
# REALIZED WE DO NOT WANT TO CVEL2 THE MS'S IN JANUARY 2023 (THANKS RICH!)
























































"""
######################################################################################
########### PERFORM VELOCITY REGRIDDING AND SPLIT OUT SPWS SIMULTANEOUSLY ############
######################################################################################
"""

""" Extremely strangely, trying this on all spectral windows gave the cvel2 error:
'Task cvel2 raised an exception of class RuntimeError with the following message:
Currently the option 'combinespws' is only supported when the number of channels
is the same for all the spw's selected. One of the SPWs selected has 1 channels,
but another selected SPW has 960 channels.'
It's strange because there is no combinespws parameter. But I figure we can go
one spw at a time, and split in the process. """

""" First do non-continuum-subtracted ms's (contains spws 0,1,2,3,4) """
# for EB in data_dict['EBs']:
#     for i,spwi in enumerate(['1', '2', '3', '4']):
#         inputvis            = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms']
#         outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms'].replace('.ms', '_spw'+str(i)+'.ms') # want i instead of spwi to align with contsub indexing
#
#         os.system('rm -rf ' + outputvis)
#         split(vis=inputvis, outputvis=outputvis, spw=spwi, datacolumn='data') # here, spws are indexed 1,2,3,4
#         listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')
#
#         tb.open(outputvis+'/SPECTRAL_WINDOW')
#         chanfreqs = tb.getcol('CHAN_FREQ', startrow = 0, nrow = 1) # split spw id will always just be 0, contains 1 spw now
#         tb.close()
#         print('--> Before cvel2: '+str(chanfreqs[0][0]/1e9)+' GHz ... '+str(chanfreqs[-1][0]/1e9)+' GHz - TOPO frame (total of '+str(len(chanfreqs))+' channels)')
#
#         os.system('rm -rf ' + outputvis+'.cvel2')
#         cvel2(vis=outputvis, outputvis=outputvis+'.cvel2', spw='0', outframe='LSRK') # split spw id will always just be 0, contains 1 spw now
#         listobs(vis=outputvis+'.cvel2', listfile=outputvis+'.cvel2'+'.listobs.txt', overwrite=True)
#         os.system('rm -rf ' + outputvis) # don't need this now
#
#         tb.open(outputvis+'.cvel2'+'/SPECTRAL_WINDOW')
#         chanfreqs = tb.getcol('CHAN_FREQ', startrow = 0, nrow = 1) # split spw id will always just be 0, contains 1 spw now
#         tb.close()
#         print('--> After cvel2: '+str(chanfreqs[0][0]/1e9)+' GHz ... '+str(chanfreqs[-1][0]/1e9)+' GHz - LSRK frame (total of '+str(len(chanfreqs))+' channels)')


""" Second do continuum-subtracted ms's (contains spws 0,1,2,3) """
# for EB in data_dict['EBs']:
#     for i,spwi in enumerate(['0', '1', '2', '3']):
#         inputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms.contsub']
#         outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms.contsub'].replace('.ms', '_spw'+spwi+'.ms') # want spwi now
#
#         os.system('rm -rf ' + outputvis)
#         split(vis=inputvis, outputvis=outputvis, spw=spwi, datacolumn='data') # here, spws are indexed 0,1,2,3
#         listobs(vis=outputvis, listfile=outputvis+'.listobs.txt')
#
#         tb.open(outputvis+'/SPECTRAL_WINDOW')
#         chanfreqs = tb.getcol('CHAN_FREQ', startrow = 0, nrow = 1) # split spw id will always just be 0, contains 1 spw now
#         tb.close()
#         print('--> Before cvel2: '+str(chanfreqs[0][0]/1e9)+' GHz ... '+str(chanfreqs[-1][0]/1e9)+' GHz - TOPO frame (total of '+str(len(chanfreqs))+' channels)')
#
#         os.system('rm -rf ' + outputvis+'.cvel2')
#         cvel2(vis=outputvis, outputvis=outputvis+'.cvel2', spw='0', outframe='LSRK') # split spw id will always just be 0, contains 1 spw now
#         listobs(vis=outputvis+'.cvel2', listfile=outputvis+'.cvel2'+'.listobs.txt')
#         os.system('rm -rf ' + outputvis) # don't need this now
#
#         tb.open(outputvis+'.cvel2'+'/SPECTRAL_WINDOW')
#         chanfreqs = tb.getcol('CHAN_FREQ', startrow = 0, nrow = 1) # split spw id will always just be 0, contains 1 spw now
#         tb.close()
#         print('--> After cvel2: '+str(chanfreqs[0][0]/1e9)+' GHz ... '+str(chanfreqs[-1][0]/1e9)+' GHz - LSRK frame (total of '+str(len(chanfreqs))+' channels)')


# """ Point of interest only: Save the channel frequencies before and after cvel2 """
# import numpy as np
#
# for EB in data_dict['EBs']:
#     for i,spwi in enumerate(['1', '2', '3', '4']):
#         before_cvel2           = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms']
#         tb.open(before_cvel2+'/SPECTRAL_WINDOW')
#         before_cvel2_chanfreqs = tb.getcol('CHAN_FREQ', startrow = i+1, nrow = 1) # 1,2,3,4
#         tb.close()
#         np.savetxt(data_dict['NRAO_path']+'workflow/step4/'+EB+'_spw'+str(i)+'_chanfreqs_beforecvel2.txt', before_cvel2_chanfreqs)
#
#         after_cvel2            = data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms'].replace('.ms', '_spw'+str(i)+'.ms.cvel2') # 0,1,2,3
#         tb.open(after_cvel2+'/SPECTRAL_WINDOW')
#         after_cvel2_chanfreqs = tb.getcol('CHAN_FREQ', startrow = 0, nrow = 1) # will always be 0
#         tb.close()
#         np.savetxt(data_dict['NRAO_path']+'workflow/step4/'+EB+'_spw'+str(i)+'_chanfreqs_aftercvel2.txt', after_cvel2_chanfreqs)



"""
####################################################################################
########### COMBINE FINAL MEASUREMENT SETS, ONE FOR EACH LINE, AND SAVE ############
####################################################################################
"""
molecules = ['SO', 'C18O', '13CO', '12CO']

""" First do non-continuum-subtracted ms's """
# for spwi,molecule in enumerate(molecules):
#     ms_list_to_concatenate = []
#     for EB in data_dict['EBs']:
#         ms_list_to_concatenate.append(data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms'].replace('.ms', '_spw'+str(spwi)+'.ms')+'.cvel2')
#
#     print('For molecule: ', molecule)
#     print('Combining spectral window '+str(spwi)+' across these measurement sets:', ms_list_to_concatenate)
#
#     final_line_ms = 'ABAur_'+molecule+'.bin30s.ms.cvel2'
#
#     # os.system('rm -rf %s*' % final_line_ms)
#     # concat(vis          = ms_list_to_concatenate,
#     #        concatvis    = final_line_ms,
#     #        dirtol       = '0.1arcsec',
#     #        copypointing = False)
#     listobs(vis=final_line_ms, listfile=final_line_ms+'.listobs.txt')
#     os.system('tar cvzf backups/' + final_line_ms+'.tgz ' + final_line_ms)


""" Second do continuum-subtracted ms's (contains spws 0,1,2,3) """
# for spwi,molecule in enumerate(molecules):
#     ms_list_to_concatenate = []
#     for EB in data_dict['EBs']:
#         ms_list_to_concatenate.append(data_dict['NRAO_path']+data_dict[EB]['_initlines_selfcal.ms'].replace('.ms', '_spw'+str(spwi)+'.ms')+'.contsub.cvel2')
#
#     print('For molecule: ', molecule)
#     print('Combining spectral window '+str(spwi)+' across these measurement sets:', ms_list_to_concatenate)
#
#     final_line_ms = 'ABAur_'+molecule+'.bin30s.ms.contsub.cvel2'
#
#     os.system('rm -rf %s*' % final_line_ms)
#     concat(vis          = ms_list_to_concatenate,
#            concatvis    = final_line_ms,
#            dirtol       = '0.1arcsec',
#            copypointing = False)
#     listobs(vis=final_line_ms, listfile=final_line_ms+'.listobs.txt')
#     os.system('tar cvzf backups/' + final_line_ms+'.tgz ' + final_line_ms)
#
# os.system('mv ABAur*.tgz backups')








# sys.exit()
