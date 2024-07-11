'''
This script was written for CASA 6.2.1.7
(on NRAO computers: casa -r 6.2.1-7-pipeline-2021.2.0.128)
(though we switch halfway through to casa -r 6.4.1-12-pipeline-2022.2.0.64)

ALMA Program ID: 2021.1.00690.S (PI: R. Dong)
reducer: J. Speedie

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
'''

""" Starting matter """
import os
execfile('dictionary_data.py') # loads data_dict
execfile('selfcal_utils.py') # necessary for an initial round of selfcal at the end

"""
######################################################
###### SELF-CAL OF ALIGNED SHORT-BASELINE EB'S #######
######################################################

For reference in choosing solint, here's the number of seconds in each scan, in
each execution block:

SB_EB1: 486 486 486 486 486 30
SB_EB2: 486 486 486 486

The "average interval (seconds)" of every scan in every execution block is 6.05
seconds. (Representing smallest possible solint).
"""

"""Start ourselves off with an orientation"""
listobs(vis=data_dict['NRAO_path']+data_dict['SB_concat']['contp0'], listfile=data_dict['NRAO_path']+data_dict['SB_concat']['contp0']+'.listobs.txt')
listobs(vis=data_dict['NRAO_path']+data_dict['LB_concat']['contp0'], listfile=data_dict['NRAO_path']+data_dict['LB_concat']['contp0']+'.listobs.txt')


"""Set image plane cell and image size:
imsize: needs to be big enough to encompass the FOV.
        FOV (FWHM of primary beam) = 1.13 * lambda_obs / (D=12 m) = 0.00012274 rad = 25.316 arcsec
cellsize: needs to be small enough to well sample the beam (5-8 times)
        beamsize (arcsec) = (lambda_obs / maximum baseline) * 206265
            SB1, SB2 maximum baseline: 500.2 m, 680.0 m
            SB1, SB2 beam size: 0.53771036 arcsec, 0.39537526 arcsec
                # pipeline SB cell: 0.12 arsec
            LB1-6 maximum baseline: 2.6 km
            LB1-6 beam size: 0.10340584 arcsec
                # pipeline LB cell: 0.031 arcsec"""
SB_cellsize     = '0.04arcsec' # samples beam 9-12 times
LB_cellsize     = '0.01arcsec' # samples beam ~10 times
SB_imsize       = 500 # to image FOV of 20 arcsec (radius)
LB_imsize       = 2000 # to image FOV of 20 arcsec (radius)

""" Define simple masks and clean scales for imaging """
mask_diam   = 2. 	# diameter of mask in arcsec; for LB data (SB data needs to be bigger)
mask_ra     = '04h55m45.8524s' # roughly; estimated in CARTA/DS9
mask_dec    = '+30.33.03.755'  # roughly; estimated in CARTA/DS9
SB_mask     = 'circle[[%s, %s], %.1farcsec]' %(mask_ra, mask_dec, mask_diam+0.8)
LB_mask     = 'circle[[%s, %s], %.1farcsec]' %(mask_ra, mask_dec, mask_diam)

""" Define a noise annulus, measure the peak SNR in map """
noise_annulus = "annulus[[%s, %s],['%.2farcsec', '9.00arcsec']]" % (mask_ra, mask_dec, 1.1*mask_diam) # outer annulus radius (9arcsec) almost = FOV


"""
################ IMAGE THE STARTING MS FOR MODEL-GENERATION AND COMPARISON #################
"""

# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'],
              imagename     = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = SB_imsize,
              cellsize      = SB_cellsize,
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
number (interactive) iterations performed: 2

imagename         = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
dir4lasts         = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('.ms', '/')
os.system('mkdir '+dir4lasts)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'SB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict['contp0.ms']:  {'SB_concat': {'beammajor': 0.90156871080396, 'beamminor': 0.6002742648124799, 'beampa': -7.856073856354, 'disk_flux': 96.75444235717598, 'peak_intensity': 13.460222631692886, 'rms': 68.33170657680805, 'SNR': 196.9835572094627}}
# selfcal_dict:  {'SB_concat': {'beammajor': 0.90191823244092, 'beamminor': 0.60061353445044, 'beampa': -7.856902599335, 'disk_flux': 96.24543315326115, 'peak_intensity': 13.490645214915276, 'rms': 77.00131321364115, 'SNR': 175.20019661853433}}
os.system('mv *.last '+dir4lasts)



"""
################ PERFORM 1st ROUND OF SELFCAL #################
"""
"""We need to not combine spectral windows at least once, to remove any potential per-spw phase offsets"""

inputvis          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0']
caltable          = inputvis.replace('p0.ms', 'p1.cal')
outputvis         = inputvis.replace('p0.ms', 'p1.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['SB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = '' #
solint      = 'inf' #
minsnr      = 2.5 # don't want this to be too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['SB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=1 at 2022/04/19/20:04:39.4
# 1 of 26 solutions flagged due to SNR < 2.5 in spw=1 at 2022/04/19/20:19:10.7
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=6 at 2022/05/15/18:17:09.1
# 1 of 33 solutions flagged due to SNR < 2.5 in spw=8 at 2022/05/15/18:44:40.3

""" Apply the solutions """
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp,  calwt=calwt)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')

""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = SB_imsize,
              cellsize      = SB_cellsize,
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 3

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'SB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict['contp1.ms']:  {'SB_concat': {'beammajor': 0.90191823244092, 'beamminor': 0.60061353445044, 'beampa': -7.856902599335, 'disk_flux': 96.75365587533217, 'peak_intensity': 13.484212569892406, 'rms': 68.46133555154267, 'SNR': 196.96099208786998}}
# selfcal_dict:  {'SB_concat': {'beammajor': 0.9019809961318801, 'beamminor': 0.60065919160848, 'beampa': -7.844259262085, 'disk_flux': 96.75450787390174, 'peak_intensity': 13.490589335560799, 'rms': 68.51391554740368, 'SNR': 196.9029098362781}}

os.system('mv *.last '+dir4lasts)


"""
################ PERFORM 2nd ROUND OF SELFCAL #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('p0.ms', 'p1.ms')
caltable          = inputvis.replace('p1.ms', 'p2.cal')
outputvis         = inputvis.replace('p1.ms', 'p2.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['SB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = '' #
solint      = '243s' # half of a scan
minsnr      = 2.5 # don't want this to be too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['SB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)

# looks like now or, latest, next round, is the time to combine spws
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=1 at 2022/04/19/19:35:26.7
# 2 of 36 solutions flagged due to SNR < 2.5 in spw=2 at 2022/04/19/19:35:26.7
# 2 of 35 solutions flagged due to SNR < 2.5 in spw=3 at 2022/04/19/19:35:26.7
# 2 of 35 solutions flagged due to SNR < 2.5 in spw=3 at 2022/04/19/19:39:30.2
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=4 at 2022/04/19/19:35:26.7
# 2 of 36 solutions flagged due to SNR < 2.5 in spw=4 at 2022/04/19/19:39:30.2
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=1 at 2022/04/19/19:44:24.1
# 2 of 36 solutions flagged due to SNR < 2.5 in spw=2 at 2022/04/19/19:48:27.7
# 2 of 35 solutions flagged due to SNR < 2.5 in spw=3 at 2022/04/19/19:48:27.7
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=4 at 2022/04/19/19:44:24.1
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=4 at 2022/04/19/19:48:27.6
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=1 at 2022/04/19/19:53:43.1
# 1 of 35 solutions flagged due to SNR < 2.5 in spw=1 at 2022/04/19/19:57:46.8
# 1 of 35 solutions flagged due to SNR < 2.5 in spw=3 at 2022/04/19/19:53:43.3
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=4 at 2022/04/19/19:57:46.9
# 1 of 35 solutions flagged due to SNR < 2.5 in spw=1 at 2022/04/19/20:02:42.2
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=2 at 2022/04/19/20:06:45.8
# 1 of 34 solutions flagged due to SNR < 2.5 in spw=3 at 2022/04/19/20:06:45.7
# 1 of 35 solutions flagged due to SNR < 2.5 in spw=4 at 2022/04/19/20:06:45.7
# 2 of 23 solutions flagged due to SNR < 2.5 in spw=3 at 2022/04/19/20:16:04.4
# 1 of 24 solutions flagged due to SNR < 2.5 in spw=4 at 2022/04/19/20:16:04.3
# 1 of 23 solutions flagged due to SNR < 2.5 in spw=3 at 2022/04/19/20:19:10.7
# 1 of 33 solutions flagged due to SNR < 2.5 in spw=6 at 2022/05/15/18:19:22.3
# 2 of 34 solutions flagged due to SNR < 2.5 in spw=8 at 2022/05/15/18:15:18.8
# 2 of 33 solutions flagged due to SNR < 2.5 in spw=8 at 2022/05/15/18:19:22.3
# 2 of 33 solutions flagged due to SNR < 2.5 in spw=9 at 2022/05/15/18:15:18.5
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:28:21.1
# 1 of 33 solutions flagged due to SNR < 2.5 in spw=6 at 2022/05/15/18:24:17.1
# 4 of 33 solutions flagged due to SNR < 2.5 in spw=6 at 2022/05/15/18:28:20.7
# 2 of 33 solutions flagged due to SNR < 2.5 in spw=7 at 2022/05/15/18:28:20.1
# 2 of 33 solutions flagged due to SNR < 2.5 in spw=8 at 2022/05/15/18:24:16.6
# 4 of 32 solutions flagged due to SNR < 2.5 in spw=8 at 2022/05/15/18:28:20.5
# 2 of 31 solutions flagged due to SNR < 2.5 in spw=9 at 2022/05/15/18:24:17.2
# 1 of 32 solutions flagged due to SNR < 2.5 in spw=9 at 2022/05/15/18:28:20.6
# 1 of 33 solutions flagged due to SNR < 2.5 in spw=6 at 2022/05/15/18:33:38.2
# 3 of 32 solutions flagged due to SNR < 2.5 in spw=6 at 2022/05/15/18:37:42.1
# 1 of 34 solutions flagged due to SNR < 2.5 in spw=7 at 2022/05/15/18:33:42.1
# 3 of 32 solutions flagged due to SNR < 2.5 in spw=8 at 2022/05/15/18:33:38.2
# 3 of 31 solutions flagged due to SNR < 2.5 in spw=8 at 2022/05/15/18:37:42.2
# 1 of 32 solutions flagged due to SNR < 2.5 in spw=9 at 2022/05/15/18:33:38.3
# 1 of 33 solutions flagged due to SNR < 2.5 in spw=6 at 2022/05/15/18:42:40.1
# 3 of 33 solutions flagged due to SNR < 2.5 in spw=6 at 2022/05/15/18:46:39.9
# 2 of 34 solutions flagged due to SNR < 2.5 in spw=7 at 2022/05/15/18:46:39.8
# 1 of 32 solutions flagged due to SNR < 2.5 in spw=8 at 2022/05/15/18:46:40.0
# 2 of 33 solutions flagged due to SNR < 2.5 in spw=9 at 2022/05/15/18:42:41.4
# 1 of 33 solutions flagged due to SNR < 2.5 in spw=9 at 2022/05/15/18:46:39.9

""" Apply the solutions """
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp,  calwt=calwt)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')

""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = SB_imsize,
              cellsize      = SB_cellsize,
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 3

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'SB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict['contp2.ms']:  {'SB_concat': {'beammajor': 0.91481429338452, 'beamminor': 0.6015850901604, 'beampa': -8.591674804688, 'disk_flux': 96.83035098664075, 'peak_intensity': 13.805755414068699, 'rms': 67.89951261198962, 'SNR': 203.32628148542696}}
# selfcal_dict:  {'SB_concat': {'beammajor': 0.9148665666582, 'beamminor': 0.6016330122947999, 'beampa': -8.578067779541, 'disk_flux': 96.83045548541152, 'peak_intensity': 13.81117943674326, 'rms': 67.93956584327374, 'SNR': 203.28624808412158}}

os.system('mv *.last '+dir4lasts)


"""
################ PERFORM 3rd ROUND OF SELFCAL #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('p0.ms', 'p2.ms')
caltable          = inputvis.replace('p2.ms', 'p3.cal')
outputvis         = inputvis.replace('p2.ms', 'p3.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['SB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = 'spw' # now reluctantly combining spectral windows
solint      = '120s' # quarter of a scan
spwmap      = [0,0,0,0,0,5,5,5,5,5] # since we did combine=spw, now we need a spwmap.
minsnr      = 2.5 # don't want this to be too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['SB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)

# these failed solutions occur in SB EB2 (to be expected)
# 1 of 38 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:16:15.6
# 1 of 38 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:23:14.0
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:29:19.0
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:38:39.6


""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')

""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = SB_imsize,
              cellsize      = SB_cellsize,
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 3


imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'SB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict['contp3.ms']:  {'SB_concat': {'beammajor': 0.9178013801574, 'beamminor': 0.60588830709468, 'beampa': -8.451219558716, 'disk_flux': 96.87789973656545, 'peak_intensity': 14.024244621396065, 'rms': 67.72439317210116, 'SNR': 207.078187998786}}
# selfcal_dict:  {'SB_concat': {'beammajor': 0.9178589582441999, 'beamminor': 0.6059303283690001, 'beampa': -8.436779975891, 'disk_flux': 96.87768215597033, 'peak_intensity': 14.026996679604053, 'rms': 67.73446447962476, 'SNR': 207.0880280425561}}

os.system('mv *.last '+dir4lasts)


"""
################ PERFORM 4th ROUND OF SELFCAL #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('p0.ms', 'p3.ms')
caltable          = inputvis.replace('p3.ms', 'p4.cal')
outputvis         = inputvis.replace('p3.ms', 'p4.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['SB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = 'spw' # still reluctantly combining spectral windows
solint      = '60s' # an eighth of a scan
spwmap      = [0,0,0,0,0,5,5,5,5,5] # since we did combine=spw, still we need a spwmap.
minsnr      = 2.5 # don't want this to be too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['SB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)

# these failed solutions occur in SB EB2 (to be expected)
# 1 of 38 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:13:44.8
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:26:45.7
# 1 of 35 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:29:49.3
# 2 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:37:08.9
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:39:09.8
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:42:03.8
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:43:04.3
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:46:07.9
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:47:08.3
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:48:08.8

""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')


""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = SB_imsize,
              cellsize      = SB_cellsize,
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 3

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'SB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict['contp4.ms']:  {'SB_concat': {'beammajor': 0.9203351736067201, 'beamminor': 0.61289548873884, 'beampa': -7.853649616241, 'disk_flux': 97.01265060549277, 'peak_intensity': 14.433663338422775, 'rms': 65.79302278945445, 'SNR': 219.37984799105863}}
# selfcal_dict:  {'SB_concat': {'beammajor': 0.9203972220421199, 'beamminor': 0.6129397153853999, 'beampa': -7.838306427002, 'disk_flux': 97.01395422578855, 'peak_intensity': 14.434458687901497, 'rms': 65.82876463800486, 'SNR': 219.27281739642527}}

os.system('mv *.last '+dir4lasts)


"""
################ PERFORM 5th ROUND OF SELFCAL #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('p0.ms', 'p4.ms')
caltable          = inputvis.replace('p4.ms', 'p5.cal')
outputvis         = inputvis.replace('p4.ms', 'p5.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['SB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = 'spw,scan' # still reluctantly combining spectral windows; now scans as well
solint      = '30s' # a 16th of a scan
spwmap      = [0,0,0,0,0,5,5,5,5,5] # since we did combine=spw, still we need a spwmap.
minsnr      = 2.5 # don't want this to be too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['SB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)

# we're getting the same flags over and over; maybe a cloud passed between 6:15-6:45pm that day
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:16:43.2
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:17:43.7
# 2 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:21:08.4
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:22:58.9
# 1 of 35 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:26:00.4
# 1 of 35 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:28:33.7
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:35:20.9
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:44:20.0
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:47:23.5


""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')

""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = SB_imsize,
              cellsize      = SB_cellsize,
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 4

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'SB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'SB_concat': {'beammajor': 0.9267931580544, 'beamminor': 0.6219916939735199, 'beampa': -7.578899383545, 'disk_flux': 97.4270423236754, 'peak_intensity': 14.997098594903946, 'rms': 63.08917232357962, 'SNR': 237.712717453083}}

os.system('mv *.last '+dir4lasts)
os.system('cp -r '+outputvis+' '+outputvis.replace('.ms', '.keepsafe.ms'))


"""
################ PERFORM 6th ROUND OF SELFCAL #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('p0.ms', 'p5.ms')
caltable          = inputvis.replace('p5.ms', 'p6.cal')
outputvis         = inputvis.replace('p5.ms', 'p6.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['SB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = 'spw,scan' # still reluctantly combining spectral windows; now scans as well
solint      = '18s' # a 27th of a scan
spwmap      = [0,0,0,0,0,5,5,5,5,5] # since we did combine=spw, still we need a spwmap.
minsnr      = 2.5 # don't want this to be too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['SB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)

# 1 of 36 solutions flagged due to SNR < 2.5 in spw=0 at 2022/04/19/20:19:19.7
# 1 of 37 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:14:18.0
# 1 of 34 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:29:04.0
# 1 of 34 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:29:22.1
# 1 of 34 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:29:40.2
# 1 of 34 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:36:53.7
# 1 of 34 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:37:11.9
# 1 of 35 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:38:24.5
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:47:59.8
# 3 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:48:17.9
# 1 of 36 solutions flagged due to SNR < 2.5 in spw=5 at 2022/05/15/18:48:33.1

""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')

""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = SB_imsize,
              cellsize      = SB_cellsize,
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 4

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'SB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'SB_concat': {'beammajor': 0.9292359352111199, 'beamminor': 0.62644398212448, 'beampa': -7.268037319183, 'disk_flux': 97.56292020763868, 'peak_intensity': 15.533193945884705, 'rms': 60.367925476803954, 'SNR': 257.3087251748157}}

os.system('mv *.last '+dir4lasts)
os.system('cp -r '+outputvis+' '+outputvis.replace('.ms', '.keepsafe.ms'))




"""
################ PERFORM 7th ROUND OF SELFCAL (AMPLITUDE+PHASE) #################
"""


inputvis          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('p0.ms', 'p6.ms')
caltable          = inputvis.replace('p6.ms', 'ap.cal')
outputvis         = inputvis.replace('p6.ms', 'ap.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['SB_concat']['refants_list']
calmode     = 'ap' # amplitude+phase
combine     = '' # in none of the DSHARP continuum scripts do they combine scan for SB ap selfcal; for HD 163296 they combine 'spw'; but in-text they say they don't.
spwmap      = [0,0,0,0,0,5,5,5,5,5] # since we did combine=spw, still we need a spwmap.
solint      = 'inf' # will go up to scan length
minsnr      = 3. # this gets raised
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True
solnorm     = False # I don't know why we're doing this but we are because DSHARP does; default is False anyway

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant, solnorm=solnorm)
for observation in data_dict['SB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)
# combine='spw'
# 1 of 36 solutions flagged due to SNR < 3 in spw=0 at 2022/04/19/20:19:10.5
# 1 of 37 solutions flagged due to SNR < 3 in spw=5 at 2022/05/15/18:17:16.9
# 2 of 36 solutions flagged due to SNR < 3 in spw=5 at 2022/05/15/18:26:06.1
# 3 of 37 solutions flagged due to SNR < 3 in spw=5 at 2022/05/15/18:35:39.5
# 1 of 36 solutions flagged due to SNR < 3 in spw=5 at 2022/05/15/18:44:28.8

#combine=''
# 1 of 36 solutions flagged due to SNR < 3 in spw=2 at 2022/04/19/19:37:26.9
# 4 of 35 solutions flagged due to SNR < 3 in spw=3 at 2022/04/19/19:37:26.4
# 2 of 35 solutions flagged due to SNR < 3 in spw=4 at 2022/04/19/19:37:25.3
# 1 of 34 solutions flagged due to SNR < 3 in spw=3 at 2022/04/19/19:46:22.4
# 1 of 35 solutions flagged due to SNR < 3 in spw=1 at 2022/04/19/19:55:40.1
# 1 of 35 solutions flagged due to SNR < 3 in spw=3 at 2022/04/19/19:55:41.6
# 1 of 35 solutions flagged due to SNR < 3 in spw=1 at 2022/04/19/20:04:42.1
# 1 of 36 solutions flagged due to SNR < 3 in spw=2 at 2022/04/19/20:04:41.3
# 2 of 35 solutions flagged due to SNR < 3 in spw=3 at 2022/04/19/20:04:40.3
# 1 of 35 solutions flagged due to SNR < 3 in spw=4 at 2022/04/19/20:04:40.8
# 3 of 35 solutions flagged due to SNR < 3 in spw=1 at 2022/04/19/20:13:36.8
# 2 of 35 solutions flagged due to SNR < 3 in spw=2 at 2022/04/19/20:13:52.7
# 8 of 33 solutions flagged due to SNR < 3 in spw=3 at 2022/04/19/20:13:38.8
# 3 of 34 solutions flagged due to SNR < 3 in spw=4 at 2022/04/19/20:13:44.5
# 10 of 25 solutions flagged due to SNR < 3 in spw=1 at 2022/04/19/20:19:10.7
# 9 of 30 solutions flagged due to SNR < 3 in spw=2 at 2022/04/19/20:19:10.7
# 12 of 22 solutions flagged due to SNR < 3 in spw=3 at 2022/04/19/20:19:10.8
# 10 of 25 solutions flagged due to SNR < 3 in spw=4 at 2022/04/19/20:19:10.8
# 2 of 34 solutions flagged due to SNR < 3 in spw=6 at 2022/05/15/18:17:17.4
# 1 of 34 solutions flagged due to SNR < 3 in spw=7 at 2022/05/15/18:17:18.6
# 2 of 32 solutions flagged due to SNR < 3 in spw=8 at 2022/05/15/18:17:17.8
# 3 of 33 solutions flagged due to SNR < 3 in spw=9 at 2022/05/15/18:17:16.7
# 2 of 32 solutions flagged due to SNR < 3 in spw=6 at 2022/05/15/18:26:14.2
# 1 of 33 solutions flagged due to SNR < 3 in spw=7 at 2022/05/15/18:26:14.3
# 3 of 31 solutions flagged due to SNR < 3 in spw=8 at 2022/05/15/18:26:14.9
# 2 of 31 solutions flagged due to SNR < 3 in spw=9 at 2022/05/15/18:26:17.6
# 1 of 32 solutions flagged due to SNR < 3 in spw=9 at 2022/05/15/18:35:34.3
# 1 of 31 solutions flagged due to SNR < 3 in spw=6 at 2022/05/15/18:44:32.4
# 1 of 33 solutions flagged due to SNR < 3 in spw=7 at 2022/05/15/18:44:30.9
# 1 of 31 solutions flagged due to SNR < 3 in spw=9 at 2022/05/15/18:44:32.6


""" Apply the solutions"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')


""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = SB_imsize,
              cellsize      = SB_cellsize,
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 4

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'SB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'SB_concat': {'beammajor': 0.9361717700958001, 'beamminor': 0.6296230554580801, 'beampa': -7.48325920105, 'disk_flux': 93.63928425306288, 'peak_intensity': 15.294084325432777, 'rms': 53.25272680004415, 'SNR': 287.1981444792438}}

os.system('mv *.last '+dir4lasts)



"""
###############################################################################
###### SELF-CAL OF COMBINED ALIGNED SHORT-BASELINE + LONG-BASELINE EB'S #######
###############################################################################

For reference in choosing solint, here's the number of seconds in each scan, in
each execution block:

SB_EB1: 486 486 486 486 486 30
SB_EB2: 486 486 486 486

LB_EB1: 121 181 181 181 60  121 181 181 181 60  121 181 181 181 60  121 181 181 181
LB_EB2: 121 181 181 181 30  121 181 181 181 60  121 181 181 181 30  121 181 181 181 60
LB_EB3: 121 181 181 181 30  121 181 181 181 60  121 181 181
LB_EB4: 121 181 181 181 30  121 181 181 181 60  121 181 181 181 30  121 181 181 181 60
LB_EB5: 121 181 181 181 30  121 181 181 181 60  121 181 181 181 30  121 181 181 181 60
LB_EB6: 121 181 181 181 30  121 181 181 181 60  121 181 181 181 30  121 181 181 181 60

The "average interval (seconds)" of every scan in every execution block is 6.05
seconds. (Representing smallest possible solint).
"""

"""Merge the aligned SB+LB executions into a single MS"""

SB_concat_shifted_selfcal  = data_dict['NRAO_path']+data_dict['SB_concat']['contp0'].replace('p0.ms', 'ap.ms')
LB_concat_shifted          = data_dict['NRAO_path']+data_dict['LB_concat']['contp0']
BB_concat_shifted          = data_dict['NRAO_path']+data_dict['BB_concat']['contp0']

os.system('rm -rf %s*' % BB_concat_shifted)
concat(vis          = [SB_concat_shifted_selfcal, LB_concat_shifted],
       concatvis    = BB_concat_shifted,
       dirtol       = '0.1arcsec',
       copypointing = False)

"""Start ourselves off with an orientation"""
listobs(vis=BB_concat_shifted, listfile=BB_concat_shifted+'.listobs.txt')


"""Set image plane cell and image size:
imsize: needs to be big enough to encompass the FOV.
        FOV (FWHM of primary beam) = 1.13 * lambda_obs / (D=12 m) = 0.00012274 rad = 25.316 arcsec
cellsize: needs to be small enough to well sample the beam (5-8 times)
        beamsize (arcsec) = (lambda_obs / maximum baseline) * 206265
            SB1, SB2 maximum baseline: 500.2 m, 680.0 m
            SB1, SB2 beam size: 0.53771036 arcsec, 0.39537526 arcsec
                # pipeline SB cell: 0.12 arsec
            LB1-6 maximum baseline: 2.6 km
            LB1-6 beam size: 0.10340584 arcsec
                # pipeline LB cell: 0.031 arcsec"""
SB_cellsize     = '0.04arcsec' # samples beam 9-12 times
LB_cellsize     = '0.01arcsec' # samples beam ~10 times
SB_imsize       = 500 # to image FOV of 20 arcsec (radius)
LB_imsize       = 2000 # to image FOV of 20 arcsec (radius)

""" Define simple masks and clean scales for imaging """
mask_diam   = 2. 	# diameter of mask in arcsec; for LB data (SB data needs to be bigger)
mask_ra     = '04h55m45.8524s' # roughly; estimated in CARTA/DS9
mask_dec    = '+30.33.03.755'  # roughly; estimated in CARTA/DS9
SB_mask     = 'circle[[%s, %s], %.1farcsec]' %(mask_ra, mask_dec, mask_diam+0.8)
LB_mask     = 'circle[[%s, %s], %.1farcsec]' %(mask_ra, mask_dec, mask_diam)

""" Define a noise annulus, measure the peak SNR in map """
noise_annulus = "annulus[[%s, %s],['%.2farcsec', '9.00arcsec']]" % (mask_ra, mask_dec, 1.1*mask_diam) # outer annulus radius (9arcsec) almost = FOV


"""
################ IMAGE THE STARTING MS FOR MODEL GENERATION AND COMPARISON #################
"""

# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = BB_concat_shifted,
              imagename     = BB_concat_shifted.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask, # it's good
              imsize        = LB_imsize,
              cellsize      = LB_cellsize,
              robust        = 1.0, # raise this from 0.5 to increase SNR; we can image with lower robust later
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 5
# 2022-12-09 21:15:05	WARN	tclean::::casa	Warning! Non-zero values at the edge of the .pb image can cause unexpected aliasing effects! (found value 0.6844037175178528 at index [1996, 983, 0, 0])

imagename         = BB_concat_shifted.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask,
                                           noise_mask     = noise_annulus)
dir4lasts         = BB_concat_shifted.replace('.ms', '/')
os.system('mkdir '+dir4lasts)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'BB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'BB_concat': {'beammajor': 0.37748101353659996, 'beamminor': 0.26945322751998, 'beampa': -9.022125244141, 'disk_flux': 113.30870792425621, 'peak_intensity': 4.4128429144620895, 'rms': 31.834778669739176, 'SNR': 138.61704396446004}}

os.system('mv *.last '+dir4lasts)


"""
To avoid gaincal error:
# For spw=35:
# Current CalTable nchan= 1
# Current CalTable freq = [224698702714]
# Current Solution nchan= 1
# Current Solution freq = [224698702714]
# Diff = [6.103515625e-05]
# 2022-12-08 21:08:52	SEVERE	Calibrater::solve	Caught exception: Mismatch between Solution frequencies and existing CalTable frequencies for spw=35
# 2022-12-08 21:08:52	SEVERE		Exception Reported: Error in Calibrater::solve.
# 2022-12-08 21:08:53	SEVERE	gaincal::::casa	Task gaincal raised an exception of class RuntimeError with the following message: Error in Calibrater::solve.
# 2022-12-08 21:08:53	SEVERE	gaincal::::casa	Exception Reported: Error in gaincal: Error in Calibrater::solve.

We switch to casa 6.4.1-12-pipeline-2022.2.0.64
"""

"""
################ PERFORM 1st ROUND OF SELFCAL #################
"""
"""We want to not combine spectral windows at least once, but it's just too low SNR unfortunately"""

inputvis          = data_dict['NRAO_path']+data_dict['BB_concat']['contp0']
caltable          = inputvis.replace('p0.ms', 'p1.cal')
outputvis         = inputvis.replace('p0.ms', 'p1.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['BB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = 'scan, spw' # too many failed solns with combine='' or 'scan'
solint      = '900s' # a little over half an entire EB
spwmap      = [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35] # since we did combine=spw, still we need a spwmap.
applymode   = 'calonly' # default is 'calflag', which applies solns and flags those with SNR<minSNR. Here we follow DSHARP and don't flag those with SNR<mnSNR; use with extreme caution
minsnr      = 2.5 # DSHARP uses 1.5 for combined, but that seems too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['BB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:49:53.6
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:08:11.7
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:12:23.4
# 9 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:20:20.5
# 10 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:30:07.6
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:38:35.5
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:54:36.3
# 14 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:02:34.8
# 4 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:11:36.0
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:20:19.3
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:14:08.6
# 1 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:22:03.3
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:31:49.8
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:39:42.3
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:53:58.2
# 2 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:11:48.8
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:29:30.0
# 1 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:37:26.0
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:53:44.0
# 2 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:09:54.5
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:17:28.4
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:19:41.1
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:40:41.5
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:49:49.8
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:58:28.9
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:39:06.7
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:46:35.6
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:48:47.3
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:09:27.9
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:27:10.3
# 9 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:59:26.0
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:08:49.7
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:17:08.3


""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, applymode=applymode, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')


""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = LB_imsize,
              cellsize      = LB_cellsize,
              robust        = 1.0, # raise this from 0.5 to increase SNR; we can image with lower robust later
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 5
# 2022-12-10 00:56:54	WARN	tclean::::casa	Warning! Non-zero values at the edge of the .pb image can cause unexpected aliasing effects! (found value 0.6844037175178528 at index [1996, 983, 0, 0])

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'BB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'BB_concat': {'beammajor': 0.377480983734, 'beamminor': 0.269453197717668, 'beampa': -9.022125244141, 'disk_flux': 113.11104993796016, 'peak_intensity': 4.445771686732769, 'rms': 28.57707804339274, 'SNR': 155.57124769656667}}
os.system('mv *.last '+dir4lasts)



"""
################ PERFORM 2nd ROUND OF SELFCAL #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['BB_concat']['contp0'].replace('p0.ms', 'p1.ms')
caltable          = inputvis.replace('p1.ms', 'p2.cal')
outputvis         = inputvis.replace('p1.ms', 'p2.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['BB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = 'scan, spw' # too many failed solns with combine='' or 'scan'
solint      = '360s' # a little over half an entire EB
spwmap      = [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35] # since we did combine=spw, still we need a spwmap.
applymode   = 'calonly' # default is 'calflag', which applies solns and flags those with SNR<minSNR. Here we follow DSHARP and don't flag those with SNR<mnSNR; use with extreme caution
minsnr      = 2.5 # DSHARP uses 1.5 for combined, but that seems too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable) # had a small accident here (12Dec22), had to regenerate it
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['BB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)
# 9 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:16:32.6
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:25:49.2
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:48:01.0
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:02:04.0
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:24:16.4
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:11:37.2
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:12:38.1
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:17:03.8
# 18 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:20:44.7
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:24:23.0
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:31:59.4
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:39:09.1
# 4 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:42:34.7
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:46:43.1
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:53:17.1
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:56:41.0
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:00:17.2
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:07:52.6
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:14:45.6
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:18:25.9
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:20:16.4
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:13:29.3
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:14:29.8
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:18:57.2
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:22:28.8
# 1 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:26:00.8
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:33:39.7
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:40:39.7
# 13 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:44:13.7
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:48:39.2
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:54:45.9
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:58:19.3
# 2 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:11:48.8
# 1 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:16:33.8
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:20:08.7
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:23:44.1
# 2 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:31:19.0
# 3 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:38:19.6
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:41:53.3
# 3 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:46:10.1
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:52:27.3
# 11 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:55:59.9
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:59:37.3
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:07:13.5
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:14:20.5
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:17:49.5
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:19:41.1
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:40:41.4
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:45:31.7
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:49:04.5
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:52:40.9
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:10:56.4
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:15:24.9
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:21:34.3
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:25:07.5
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:43:24.5
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:46:56.8
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:48:47.3
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:09:27.9
# 7 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:17:47.6
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:21:23.5
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:28:56.9
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:39:29.7
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:50:02.0
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:53:33.8
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:57:06.8
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:04:43.0
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:11:43.4
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:15:16.1
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:17:05.3


""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, applymode=applymode, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')


""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = LB_imsize,
              cellsize      = LB_cellsize,
              robust        = 1.0, # raise this from 0.5 to increase SNR; we can image with lower robust later
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 5
# 2022-12-12 15:28:35	WARN	tclean::::casa	Warning! Non-zero values at the edge of the .pb image can cause unexpected aliasing effects! (found value 0.6844037175178528 at index [1996, 983, 0, 0])

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'BB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'BB_concat': {'beammajor': 0.37748101353659996, 'beamminor': 0.26945322751998, 'beampa': -9.022125244141, 'disk_flux': 109.00217241807141, 'peak_intensity': 4.454145673662424, 'rms': 25.28222451231294, 'SNR': 176.17696858491104}}
os.system('mv *.last '+dir4lasts)





"""
################ PERFORM 3rd ROUND OF SELFCAL #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['BB_concat']['contp0'].replace('p0.ms', 'p2.ms')
caltable          = inputvis.replace('p2.ms', 'p3.cal')
outputvis         = inputvis.replace('p2.ms', 'p3.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['BB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = 'scan, spw' # too many failed solns with combine='' or 'scan'
solint      = '180s' # a little over half an entire EB
spwmap      = [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35] # since we did combine=spw, still we need a spwmap.
applymode   = 'calonly' # default is 'calflag', which applies solns and flags those with SNR<minSNR. Here we follow DSHARP and don't flag those with SNR<mnSNR; use with extreme caution
minsnr      = 2.5 # DSHARP uses 1.5 for combined, but that seems too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['BB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:16:32.6
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:31:53.5
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:49:50.6
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:53:55.4
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:08:08.7
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:15:09.4
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:12:04.8
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:13:04.9
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:15:44.4
# 9 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:19:40.1
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:23:58.8
# 11 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:27:54.1
# 10 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:30:10.6
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:33:20.7
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:37:39.5
# 4 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:41:31.2
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:44:24.9
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:48:28.9
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:51:54.3
# 12 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:55:35.0
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:59:54.3
# 10 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:03:50.5
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:06:05.8
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:09:12.1
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:13:27.0
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:17:20.4
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:20:16.4
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:13:29.3
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:14:29.8
# 3 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:17:31.9
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:21:25.1
# 3 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:25:38.5
# 13 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:29:38.5
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:31:52.8
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:35:00.3
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:39:15.5
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:43:10.4
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:46:03.3
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:50:10.0
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:53:27.5
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:57:15.3
# 2 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:11:48.8
# 1 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:15:10.9
# 3 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:19:04.5
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:23:21.4
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:27:15.2
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:29:33.0
# 3 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:32:40.7
# 4 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:36:55.3
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:40:50.1
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:43:43.1
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:47:48.7
# 6 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:51:07.1
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:54:56.0
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:59:14.2
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:03:08.1
# 10 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:05:26.9
# 6 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:08:36.3
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:12:50.3
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:16:46.3
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:19:41.1
# 8 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:39:56.2
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:40:56.6
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:44:05.6
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:48:01.4
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:52:17.9
# 6 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:56:14.7
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:58:31.9
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:01:40.3
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:05:58.8
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:09:52.9
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:12:47.7
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:16:55.7
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:24:03.5
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:32:17.5
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:34:34.4
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:37:43.0
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:41:59.1
# 1 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:45:53.2
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:48:47.3
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:08:54.6
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:09:55.1
# 3 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:16:43.0
# 1 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:21:00.3
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:24:53.1
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:27:10.3
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:30:17.4
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:38:26.1
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:41:19.2
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:56:44.6
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:00:40.5
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:02:56.4
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:06:03.4
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:10:18.3
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:14:13.3
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:17:05.3



""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, applymode=applymode, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')


""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = LB_imsize,
              cellsize      = LB_cellsize,
              robust        = 1.0, # raise this from 0.5 to increase SNR; we can image with lower robust later
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 5
# 2022-12-12 16:35:46	WARN	tclean::::casa	Warning! Non-zero values at the edge of the .pb image can cause unexpected aliasing effects! (found value 0.6844037175178528 at index [1996, 983, 0, 0])

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'BB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'BB_concat': {'beammajor': 0.37748101353659996, 'beamminor': 0.26945322751998, 'beampa': -9.022125244141, 'disk_flux': 108.94340607637508, 'peak_intensity': 4.509070888161659, 'rms': 25.058145680905398, 'SNR': 179.94431613499734}}

os.system('mv *.last '+dir4lasts)
os.system('cp -r '+outputvis+' '+outputvis.replace('.ms', '.keepsafe.ms'))




"""
################ PERFORM 4th ROUND OF SELFCAL #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['BB_concat']['contp0'].replace('p0.ms', 'p3.ms')
caltable          = inputvis.replace('p3.ms', 'p4.cal')
outputvis         = inputvis.replace('p3.ms', 'p4.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['BB_concat']['refants_list']
calmode     = 'p' # phase-only
combine     = 'scan, spw' # too many failed solns with combine='' or 'scan'
solint      = '60s' # a little over half an entire EB
spwmap      = [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35] # since we did combine=spw, still we need a spwmap.
applymode   = 'calonly' # default is 'calflag', which applies solns and flags those with SNR<minSNR. Here we follow DSHARP and don't flag those with SNR<mnSNR; use with extreme caution
minsnr      = 2.5 # DSHARP uses 1.5 for combined, but that seems too low
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant)
for observation in data_dict['BB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine)
# 9 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:16:32.6
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:18:03.3
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:19:52.0
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:21:52.9
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:23:45.3
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:24:45.7
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:25:46.2
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:27:59.9
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:29:00.4
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:31:53.4
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:35:56.9
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:39:49.5
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:42:03.8
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:44:04.7
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:45:57.1
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:46:57.5
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:47:58.0
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:49:50.6
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:53:14.5
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:54:15.0
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/14:58:07.5
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:00:00.0
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:01:00.5
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:04:15.1
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:05:15.6
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:08:08.7
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:11:11.1
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:14:03.0
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:15:03.4
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:16:03.9
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:18:18.3
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:19:18.8
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:20:19.3
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:22:12.4
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:23:12.9
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=10 at 2022/07/17/15:24:13.4
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:11:34.2
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:12:16.7
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:13:04.9
# 9 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:14:45.9
# 14 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:15:46.3
# 10 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:16:46.8
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:18:40.7
# 15 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:19:41.2
# 16 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:20:41.7
# 4 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:22:57.6
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:23:58.1
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:24:58.6
# 10 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:26:36.8
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:29:25.4
# 10 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:30:25.9
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:32:19.4
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:33:19.9
# 4 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:34:20.4
# 12 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:36:36.7
# 16 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:37:37.2
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:38:37.6
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:40:30.8
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:41:31.2
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:42:31.7
# 5 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:44:24.9
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:47:50.2
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:48:50.7
# 12 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:50:44.0
# 4 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:51:44.5
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:52:44.9
# 17 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:54:37.0
# 15 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:55:37.5
# 10 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:56:38.0
# 14 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:58:53.1
# 9 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/12:59:53.6
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:00:54.1
# 13 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:02:31.7
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:05:20.2
# 11 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:06:20.7
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:08:12.6
# 9 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:09:13.0
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:10:13.5
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:12:27.9
# 8 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:13:28.6
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:14:28.9
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:16:21.9
# 6 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:17:22.3
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:18:22.8
# 7 of 36 solutions flagged due to SNR < 2.5 in spw=15 at 2022/07/19/13:20:16.4
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:13:29.3
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:14:20.7
# 14 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:15:00.0
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:16:31.1
# 11 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:17:31.6
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:18:32.0
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:20:24.8
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:21:25.3
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:22:25.8
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:24:40.0
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:25:40.4
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:26:40.9
# 15 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:28:19.2
# 11 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:31:07.4
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:32:07.9
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:34:00.0
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:35:00.4
# 11 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:36:00.9
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:38:15.7
# 11 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:39:16.2
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:40:16.6
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:42:09.8
# 11 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:43:10.2
# 14 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:44:10.7
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:46:03.3
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:49:28.6
# 11 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:50:29.1
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:52:21.1
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:53:21.5
# 11 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:54:22.0
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:56:15.3
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:57:15.7
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=20 at 2022/07/20/14:58:16.2
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:11:18.5
# 2 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:12:18.9
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:14:10.9
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:15:11.4
# 3 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:16:11.8
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:18:04.7
# 6 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:19:05.2
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:20:05.7
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:22:20.7
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:23:21.2
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:24:21.7
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:25:59.5
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:28:47.7
# 3 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:29:48.2
# 12 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:31:40.1
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:32:40.6
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:33:41.0
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:35:55.5
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:36:56.0
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:37:56.5
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:39:49.3
# 8 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:40:49.8
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:41:50.3
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:43:43.1
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:47:07.8
# 6 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:48:08.3
# 12 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:50:00.2
# 9 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:51:00.7
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:52:01.2
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:53:55.9
# 6 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:54:56.4
# 11 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:55:56.8
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:58:13.5
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/14:59:14.0
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:00:14.4
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:01:53.9
# 11 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:04:41.7
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:05:42.2
# 11 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:07:35.3
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:08:35.8
# 9 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:09:36.3
# 10 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:11:51.0
# 10 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:12:51.5
# 13 of 41 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:13:52.0
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:15:45.5
# 12 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:16:46.0
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=25 at 2022/07/21/15:17:46.5
# 8 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:39:56.2
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:40:41.5
# 8 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:41:26.9
# 6 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:43:05.4
# 6 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:44:05.9
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:45:06.3
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:47:00.5
# 8 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:48:01.0
# 6 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:49:01.4
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:51:17.3
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:52:17.8
# 7 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:53:18.3
# 10 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:54:57.6
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:57:46.5
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/11:58:47.0
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:00:40.3
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:01:40.8
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:02:41.2
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:04:57.4
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:05:57.9
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:06:58.3
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:08:52.5
# 7 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:09:52.9
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:10:53.4
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:12:47.7
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:16:14.5
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:17:15.0
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:19:08.2
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:20:08.7
# 7 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:21:09.2
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:23:03.5
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:24:04.0
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:25:04.5
# 3 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:27:20.4
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:28:20.9
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:29:21.3
# 8 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:31:00.2
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:33:49.0
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:34:49.5
# 7 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:36:42.2
# 8 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:37:42.7
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:38:43.1
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:40:58.5
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:41:59.0
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:42:59.4
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:44:52.8
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:45:53.3
# 5 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:46:53.8
# 4 of 43 solutions flagged due to SNR < 2.5 in spw=30 at 2022/07/22/12:48:47.3
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:08:54.6
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:09:52.0
# 16 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:10:25.3
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:11:49.9
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:12:50.4
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:13:50.9
# 2 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:15:43.6
# 4 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:16:44.0
# 4 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:17:44.5
# 3 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:19:58.9
# 4 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:20:59.3
# 1 of 40 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:21:59.8
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:23:37.1
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:26:24.9
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:27:25.4
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:29:17.1
# 9 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:30:17.6
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:31:18.1
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:33:32.3
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:34:32.7
# 6 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:35:33.2
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:37:25.8
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:38:26.2
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:39:26.7
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:41:19.2
# 1 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:44:43.9
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:45:44.4
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:47:36.2
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:48:36.7
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:49:37.2
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:51:29.8
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:52:30.3
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:53:30.7
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:55:44.7
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:56:45.2
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:57:45.7
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/13:59:23.0
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:02:11.0
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:03:11.4
# 2 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:05:03.1
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:06:03.6
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:07:04.1
# 6 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:09:18.4
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:10:18.9
# 6 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:11:19.4
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:13:12.2
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:14:12.6
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:15:13.1
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=35 at 2022/07/22/14:17:05.3

""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, applymode=applymode, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')


""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = LB_imsize,
              cellsize      = LB_cellsize,
              robust        = 1.0, # raise this from 0.5 to increase SNR; we can image with lower robust later
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 7
# 2022-12-12 17:58:43	WARN	tclean::::casa	Warning! Non-zero values at the edge of the .pb image can cause unexpected aliasing effects! (found value 0.6844037175178528 at index [1996, 983, 0, 0])

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'BB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'BB_concat': {'beammajor': 0.377480983734, 'beamminor': 0.269453197717668, 'beampa': -9.022125244141, 'disk_flux': 106.10958340099475, 'peak_intensity': 4.532103426754475, 'rms': 23.12901527796564, 'SNR': 195.94882757814952}}

os.system('mv *.last '+dir4lasts)
os.system('cp -r '+outputvis+' '+outputvis.replace('.ms', '.keepsafe.ms'))




"""
################ PERFORM 5th ROUND OF SELFCAL (AMPLITUDE+PHASE) #################
"""

inputvis          = data_dict['NRAO_path']+data_dict['BB_concat']['contp0'].replace('p0.ms', 'p4.ms')
caltable          = inputvis.replace('p4.ms', 'ap.cal')
outputvis         = inputvis.replace('p4.ms', 'ap.ms')
dir4lasts         = outputvis.replace('.ms', '/')
os.system('mkdir '+dir4lasts)

vis         = inputvis # input; these should be the same
parentvis   = inputvis # input; these should be the same
caltable    = caltable
gaintable   = [caltable]
gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
spw         = ''
refant      = data_dict['BB_concat']['refants_list']
calmode     = 'ap' # amplitude+phase
combine     = 'spw' # might need to be spw,scan and then 900s solint
solint      = 'inf' # to go up to scan length
spwmap      = [0,0,0,0,0, 5,5,5,5,5, 10,10,10,10,10, 15,15,15,15,15, 20,20,20,20,20, 25,25,25,25,25, 30,30,30,30,30, 35,35,35,35,35] # since we did combine=spw, still we need a spwmap.
applymode   = 'calonly' # default is 'calflag', which applies solns and flags those with SNR<minSNR. Here we follow DSHARP and don't flag those with SNR<mnSNR; use with extreme caution
minsnr      = 3. # raise this
minblperant = 4
interp      = 'linearPD' # PD is overkill but it's fine
calwt       = True
solnorm     = False # I don't know why we're doing this but we are because DSHARP does; default is False anyway

""" Generate the solutions """
os.system('rm -rf ' + caltable)
gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
        calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
        minblperant=minblperant, solnorm=solnorm)
for observation in data_dict['BB_concat']['observations']:
    for quantity in ['phase', 'amp', 'SNR']:
        plot_gaincal_solutions_per_antenna(caltable=caltable, parentvis=vis, quantity=quantity,
                               plot_average_soln=bool, solint=solint, minsnr=minsnr,
                               spw=spw, observation=observation, combine=combine, calmode=calmode)
        for bool in [False, True]:
            plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
                                   plot_average_soln=bool, solint=solint, minsnr=minsnr,
                                   spw=spw, observation=observation, combine=combine, calmode=calmode)
# 3 of 41 solutions flagged due to SNR < 3 in spw=10 at 2022/07/17/14:31:53.5
# 2 of 41 solutions flagged due to SNR < 3 in spw=10 at 2022/07/17/14:49:50.6
# 3 of 41 solutions flagged due to SNR < 3 in spw=10 at 2022/07/17/14:53:44.8
# 1 of 41 solutions flagged due to SNR < 3 in spw=10 at 2022/07/17/15:08:08.7
# 1 of 41 solutions flagged due to SNR < 3 in spw=10 at 2022/07/17/15:15:03.3
# 1 of 41 solutions flagged due to SNR < 3 in spw=10 at 2022/07/17/15:19:18.8
# 7 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:12:23.5
# 7 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:15:44.3
# 10 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:19:40.1
# 7 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:23:58.9
# 10 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:26:36.7
# 7 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:29:55.4
# 5 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:33:20.8
# 8 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:37:39.5
# 5 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:41:31.2
# 6 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:44:24.9
# 6 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:48:18.8
# 5 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:51:49.2
# 12 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:55:34.8
# 7 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/12:59:54.3
# 12 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/13:02:31.8
# 8 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/13:05:51.1
# 6 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/13:09:12.1
# 5 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/13:13:27.1
# 6 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/13:17:20.4
# 7 of 36 solutions flagged due to SNR < 3 in spw=15 at 2022/07/19/13:20:16.3
# 8 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:14:08.5
# 7 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:17:32.0
# 6 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:21:25.1
# 5 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:25:38.3
# 17 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:28:19.2
# 8 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:31:37.8
# 6 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:35:00.5
# 7 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:39:15.3
# 8 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:43:10.7
# 13 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:46:03.2
# 8 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:49:59.2
# 7 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:53:21.4
# 7 of 42 solutions flagged due to SNR < 3 in spw=20 at 2022/07/20/14:57:15.4
# 5 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:11:48.8
# 5 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:15:10.9
# 6 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:19:04.5
# 6 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:23:21.5
# 11 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:25:59.5
# 8 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:29:17.8
# 5 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:32:40.7
# 7 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:36:55.4
# 7 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:40:50.1
# 11 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:43:43.1
# 7 of 41 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:47:37.9
# 7 of 41 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:51:01.3
# 6 of 41 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:54:55.9
# 7 of 41 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/14:59:14.3
# 9 of 41 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/15:01:53.9
# 8 of 41 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/15:05:11.5
# 7 of 41 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/15:08:36.3
# 9 of 41 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/15:12:50.3
# 10 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/15:16:46.4
# 10 of 42 solutions flagged due to SNR < 3 in spw=25 at 2022/07/21/15:19:41.1
# 1 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/11:40:41.4
# 3 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/11:44:05.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/11:48:01.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/11:52:17.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/11:54:57.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/11:58:16.8
# 1 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:01:40.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:05:58.8
# 1 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:09:52.9
# 6 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:12:47.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:16:44.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:20:09.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:24:03.5
# 1 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:28:21.2
# 8 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:31:00.2
# 2 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:34:19.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:37:43.0
# 2 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:41:59.1
# 2 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:45:53.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=30 at 2022/07/22/12:48:47.3
# 1 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:09:27.9
# 3 of 40 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:16:43.0
# 3 of 40 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:21:00.3
# 8 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:23:37.0
# 3 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:26:55.4
# 7 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:30:17.4
# 1 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:38:26.1
# 6 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:41:19.2
# 1 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:45:14.1
# 1 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:48:37.0
# 1 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:52:30.0
# 3 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:56:44.6
# 8 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/13:59:23.0
# 6 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/14:02:41.4
# 3 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/14:06:03.4
# 4 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/14:10:18.2
# 5 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/14:14:13.4
# 6 of 41 solutions flagged due to SNR < 3 in spw=35 at 2022/07/22/14:17:05.3


""" Apply the solutions - Note: since we did combine=spw in gaincal, now we need a spwmap"""
applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp, applymode=applymode, calwt=calwt, spwmap=spwmap)
os.system('rm -rf ' + outputvis)
split(vis=vis, outputvis=outputvis, datacolumn='corrected')


""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
tclean_wrapper(vis          = outputvis,
              imagename     = outputvis.replace('.ms', ''), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = LB_imsize,
              cellsize      = LB_cellsize,
              robust        = 1.0, # raise this from 0.5 to increase SNR; we can image with lower robust later
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0mJy',
              interactive   = True,
              savemodel     = 'modelcolumn')
# number (interactive) iterations performed: 7
# 2022-12-12 20:02:30	WARN	tclean::::casa	Warning! Non-zero values at the edge of the .pb image can cause unexpected aliasing effects! (found value 0.6844037175178528 at index [1996, 983, 0, 0])

imagename         = outputvis.replace('.ms', '.fits')
image_metrics     = estimate_image_metrics(imagename      = imagename,
                                           disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                           noise_mask     = noise_annulus)
selfcal_dict      = update_selfcal_dict(save_dir      = dir4lasts,
                                        EB            = 'BB_concat',
                                        image_metrics = image_metrics)
print("selfcal_dict: ", selfcal_dict)
# selfcal_dict:  {'BB_concat': {'beammajor': 0.36318364739424003, 'beamminor': 0.258013367652888, 'beampa': -10.19715881348, 'disk_flux': 102.58159634227913, 'peak_intensity': 4.15557948872447, 'rms': 20.21376952048342, 'SNR': 205.58162021751832}}
os.system('mv *.last '+dir4lasts)


"""
######################################################################################
###### FINAL CONINTUUM OUTPUTS: SELFCAL'ED SHORT-BASELINE + LONG-BASELINE EB'S #######
######################################################################################
"""


""" Save the final MS """
chosenvis          = data_dict['NRAO_path']+data_dict['BB_concat']['contp0'].replace('p0.ms', 'ap.ms')
final_cont_ms      = 'ABAur_continuum.bin30s.ms'
split(vis=chosenvis, outputvis=final_cont_ms, spw='', timebin='30s', datacolumn='data')
listobs(vis=final_cont_ms, listfile=final_cont_ms+'.listobs.txt')
os.system('tar cvzf backups/' + final_cont_ms+'.tgz ' + final_cont_ms)


""" Image the final MS """
continuum       = data_dict['NRAO_path']+data_dict['continuum']
LB_cellsize     = '0.01arcsec' # samples beam ~10 times
LB_imsize       = 2000 # to image FOV of 20 arcsec (radius)

""" Define simple masks and clean scales for imaging """
mask_diam   = 2. 	# diameter of mask in arcsec; for LB data (SB data needs to be bigger)
mask_ra     = '04h55m45.8524s' # roughly; estimated in CARTA/DS9
mask_dec    = '+30.33.03.755'  # roughly; estimated in CARTA/DS9
SB_mask     = 'circle[[%s, %s], %.1farcsec]' %(mask_ra, mask_dec, mask_diam+0.8)

tclean_wrapper(vis          = continuum,
              imagename     = continuum.replace('.ms', '_robust05'), # will become .image, .residual, etc
              mask          = SB_mask,
              imsize        = LB_imsize,
              cellsize      = LB_cellsize,
              robust        = 0.5, # raise this from 0.5 to increase SNR; we can image with lower robust later
              gain          = 0.05, # nice and low; stops us from accidentally cleaning too deep between successive major cycles
              contspws      = '', # we want the model to include all spws
              threshold     = '0.08mJy',
              interactive   = False,
              savemodel     = 'modelcolumn')



"""
#########################################################
###### INSPECTION OF THE SELF-CALIBRATION RESULTS #######
#########################################################
"""

""" Export MS contents into Numpy save files """

SB_p0          = data_dict['NRAO_path']+data_dict['SB_concat']['contp0']
SB_p1          = SB_p0.replace('p0.ms', 'p1.ms')
SB_p2          = SB_p0.replace('p0.ms', 'p2.ms')
SB_p3          = SB_p0.replace('p0.ms', 'p3.ms')
SB_p4          = SB_p0.replace('p0.ms', 'p4.ms')
SB_p5          = SB_p0.replace('p0.ms', 'p5.ms')
SB_p6          = SB_p0.replace('p0.ms', 'p6.ms')
SB_ap          = SB_p0.replace('p0.ms', 'ap.ms')

for msfile in [SB_p0, SB_p1, SB_p2, SB_p3, SB_p4, SB_p5, SB_p6, SB_ap]:
    print("Exporting "+msfile+" as .npz...")
    export_MS(msfile)

BB_p0          = data_dict['NRAO_path']+data_dict['BB_concat']['contp0']
BB_p1          = BB_p0.replace('p0.ms', 'p1.ms')
BB_p2          = BB_p0.replace('p0.ms', 'p2.ms')
BB_p3          = BB_p0.replace('p0.ms', 'p3.ms')
BB_p4          = BB_p0.replace('p0.ms', 'p4.ms')
BB_ap          = BB_p0.replace('p0.ms', 'ap.ms')

for msfile in [BB_p0, BB_p1, BB_p2, BB_p3, BB_p4, BB_ap]:
    print("Exporting "+msfile+" as .npz...")
    export_MS(msfile)

""" Assign rough emission geometry parameters. """
PA, incl = 54, 23

""" Plot deprojected visibility profiles for all data together """
# plot_deprojected([SB_p0.replace('.ms', '.vis.npz'),
#                   SB_p1.replace('.ms', '.vis.npz'),
#                   SB_p2.replace('.ms', '.vis.npz'),
#                   SB_p3.replace('.ms', '.vis.npz'),
#                   SB_p4.replace('.ms', '.vis.npz'),
#                   SB_p5.replace('.ms', '.vis.npz'),
#                   SB_p6.replace('.ms', '.vis.npz'),
#                   SB_ap.replace('.ms', '.vis.npz')],
#                   fignametemplate='./workflow/step3/ABAur_SB_concat_shited_contselfcal',
#                   fluxscale=[1.,1.,1.,1.,1.,1.,1., 1.], PA=PA, incl=incl, show_err=False)
# plot_deprojected([BB_p0.replace('.ms', '.vis.npz'),
#                   BB_p1.replace('.ms', '.vis.npz'),
#                   BB_p2.replace('.ms', '.vis.npz'),
#                   BB_p3.replace('.ms', '.vis.npz'),
#                   BB_p4.replace('.ms', '.vis.npz'),
#                   BB_ap.replace('.ms', '.vis.npz')],
#                   fignametemplate='./workflow/step3/ABAur_BB_concat_shited_contselfcal',
#                   fluxscale=[1.,1.,1.,1.,1.,1.], PA=PA, incl=incl, show_err=False)
# plot_deprojected([SB_ap.replace('.ms', '.vis.npz'),
#                   BB_ap.replace('.ms', '.vis.npz')],
#                   fignametemplate='./workflow/step3/ABAur_BB+SB_concat_shited_contselfcal_ap',
#                   fluxscale=[1.,1.], PA=PA, incl=incl, show_err=False)
# plot_deprojected([SB_p6.replace('.ms', '.vis.npz'),
#                   BB_p4.replace('.ms', '.vis.npz')],
#                   fignametemplate='./workflow/step3/ABAur_BB+SB_concat_shited_contselfcal_pfin',
#                   fluxscale=[1.,1.], PA=PA, incl=incl, show_err=False)
#
# estimate_flux_scale(reference=SB_ap.replace('.ms', '.vis.npz'),
#                     comparison=BB_ap.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_BB+SB_concat_shited_contselfcal_ap',
#                     PA=PA, incl=incl)
# estimate_flux_scale(reference=SB_p6.replace('.ms', '.vis.npz'),
#                     comparison=BB_p4.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_BB+SB_concat_shited_contselfcal_pfin',
#                     PA=PA, incl=incl)
#
# estimate_flux_scale(reference=SB_p0.replace('.ms', '.vis.npz'),
#                     comparison=SB_p1.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_SB_concat_shited_contselfcal_p1',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=SB_p0.replace('.ms', '.vis.npz'),
#                     comparison=SB_p2.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_SB_concat_shited_contselfcal_p2',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=SB_p0.replace('.ms', '.vis.npz'),
#                     comparison=SB_p3.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_SB_concat_shited_contselfcal_p3',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=SB_p0.replace('.ms', '.vis.npz'),
#                     comparison=SB_p4.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_SB_concat_shited_contselfcal_p4',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=SB_p0.replace('.ms', '.vis.npz'),
#                     comparison=SB_p5.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_SB_concat_shited_contselfcal_p5',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=SB_p0.replace('.ms', '.vis.npz'),
#                     comparison=SB_p6.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_SB_concat_shited_contselfcal_p6',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=SB_p0.replace('.ms', '.vis.npz'),
#                     comparison=SB_ap.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_SB_concat_shited_contselfcal_ap',
#                     incl=incl, PA=PA)
#
# estimate_flux_scale(reference=BB_p0.replace('.ms', '.vis.npz'),
#                     comparison=BB_p1.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_BB_concat_shited_contselfcal_p1',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=BB_p0.replace('.ms', '.vis.npz'),
#                     comparison=BB_p2.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_BB_concat_shited_contselfcal_p2',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=BB_p0.replace('.ms', '.vis.npz'),
#                     comparison=BB_p3.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_BB_concat_shited_contselfcal_p3',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=BB_p0.replace('.ms', '.vis.npz'),
#                     comparison=BB_p4.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_BB_concat_shited_contselfcal_p4',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference=BB_p0.replace('.ms', '.vis.npz'),
#                     comparison=BB_ap.replace('.ms', '.vis.npz'),
#                     fignametemplate='./workflow/step3/ABAur_BB_concat_shited_contselfcal_ap',
#                     incl=incl, PA=PA)






















# sys.exit()
