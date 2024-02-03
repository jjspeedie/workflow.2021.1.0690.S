'''
This script was written for CASA 6.2.1.7
(on NRAO computers: casa -r 6.2.1-7-pipeline-2021.2.0.128)

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
execfile('dictionary_data.py') # loads data_dict
execfile('step1_utils.py') # loads multiple functions
execfile('selfcal_utils.py') # necessary for an initial round of selfcal at the end

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
######## Spectral averaging, pseudo-continuum ########
######################################################
"""
'''
We need to exlude the channels that contain line emission from being included
in the spectral averaging. To find those channels, we use the DSHARP function
get_flagchannels, by specifying a velocity range around the spectral line.
DSHARP says they flag +/- 25 km/s around the 12CO 2-1 line.
MAPS says they flag channels between -10 to 20 km/s, regardless of line (with all
sources having systemic velocity around 5 km/s.) This suggests +/- 15 km/s around
the lines.

From the pipeline-generated cubes, AB Aur seems to have a system velocity of ~5.8 km/s.
After trying what MAPS did (velocity_range = -10 km/s to 20 km/s), that looks to
remove too much data (the range is too wide...), to my eyes. We want as large a
bandwidth as possible, for high SNR for per-spw selfcal later.
As such, I've done custom velocity ranges.

Possible improvement: Make dirty image cubes for all EBs and all spws, and identify
        channels containing line emission manually. Reason for having not done this:
        (1) Time considerations.
'''

# """ Print to terminal and a file at the same time """
# old_print   = print
# log_file    = open('./workflow/step1/printlog.log', "w")
# print       = lambda *args, **kw: old_print(*args, **kw) or old_print(*args, file=log_file, **kw)
# print('Executing step1_prepare_continuum.py')
#
#
# print('\nPerforming spectral averaging of the continuum...') # takes about 38 mins per EB
# for EB in data_dict['EBs']:
#     inputvis            = data_dict['NRAO_path']+data_dict[EB]['.ms.split.cal.source']
#     outputvis           = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
#     print('inputvis = ', inputvis)
#     print('outputvis = ', outputvis)
#
#     flagchannels_string = get_flagchannels(ms_file          = inputvis,
#                                            ms_dict          = data_dict[EB],
#                                            velocity_range   = data_dict[EB]['velocity_ranges'])
#     print('\nFor '+EB+', flagchannels_string = ', flagchannels_string)
#
#     avg_cont(msfile         = inputvis,
#             outputvis       = outputvis,
#             flagchannels    = flagchannels_string,
#             contspws        = data_dict[EB]['cont_spws'],
#             width_array     = data_dict[EB]['width_array'],
#             datacolumn      = 'data',
#             field           = 'AB_Aur')
#
# log_file.close()

"""
######################################################
############# Initial imaging, continuum #############
######################################################
"""
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
SB_imsize       = 500 # to image FOV of 20 arcsec
LB_imsize       = 2000 # to image FOV of 20 arcsec

""" Define simple masks and clean scales for imaging """
mask_diam   = 2. 	# diameter of mask in arcsec; for LB data (SB data needs to be bigger)
mask_ra     = '04h55m45.8524s' # roughly; estimated in CARTA/DS9
mask_dec    = '+30.33.03.755'  # roughly; estimated in CARTA/DS9
SB_mask     = 'circle[[%s, %s], %.1farcsec]' %(mask_ra, mask_dec, mask_diam+0.8)
LB_mask     = 'circle[[%s, %s], %.1farcsec]' %(mask_ra, mask_dec, mask_diam)
SB_scales   = [0, 5, 10, 20, 40] # max scale ~ LAS ~ outer ring edge ~ 1.5 arcsec; min scale = delta function
LB_scales   = [0, 5, 30, 80, 150] # max scale ~ LAS ~ outer ring edge ~ 1.5 arcsec; min scale = delta function

""" Image each execution block individually, as well as per-spw images """ # gain=0.1, which is a bit high for images that matter
# for EB in data_dict['SB_EBs']:
#     image_each_obs(ms_dict      = data_dict[EB],
#                    inputvis     = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms'],
#                    mask         = SB_mask,
#                    scales       = SB_scales,
#                    imsize       = SB_imsize,
#                    cellsize     = SB_cellsize,
#                    contspws     = data_dict[EB]['cont_spws'],
#                    # threshold    = '0.0mJy', # this is taken care of inside image_each_obs
#                    interactive  = False)
# for EB in data_dict['LB_EBs']:
#     image_each_obs(ms_dict      = data_dict[EB],
#                    inputvis     = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms'],
#                    mask         = LB_mask,
#                    scales       = LB_scales,
#                    imsize       = LB_imsize,
#                    cellsize     = LB_cellsize,
#                    contspws     = data_dict[EB]['cont_spws'],
#                    # threshold    = '0.0mJy', # this is taken care of inside image_each_obs
#                    interactive  = False)


"""
######################################################
######## INITIAL SELF-CAL, ON EACH EXECUTION #########
######################################################

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

"""
################ MAKE STARTING MODEL #################

Before we start self calibration, we need to decide on our initial model. There
are many options, such as:
- A (not very deeply) cleaned image
- A model made from analytic fit to the visibilities
- An axisymmetric model, ie. azimuthal average of a (not very deeply) cleaned image
- Jess idea: The (not very deeply) cleaned image of the best EB, aligned? For
    example, for LB2-6, use LB1?
For now, we will go with a (not very deeply, but interactively) cleaned image.

There are also other considerations:
- What deconvolver do you use? For now, we will use hogbom.
- What beam (robust, uvtaper?) do you use? For now, we will use robust 0.5 and no
    uvtaper to keep it simple, and have something to improve upon later.
"""

""" Record the image metric results of each round of selfcal """
selfcal_dict = {}

# Create copies of _initcont.ms called _initcont_model.ms, whose model column
# will become the starting model for the self calibration
# for EB in data_dict['EBs']:
#     initcont        = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
#     initcont_model  = initcont.replace('.ms', '_model.ms')
#     print("Making a copy of "+initcont+" and saving it to "+initcont_model)
#     os.system('rm -rf '+initcont_model)
#     os.system('cp -r '+initcont+' '+initcont_model)


# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
# for EB in data_dict['SB_EBs']:
#     initcont        = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
#     initcont_model  = initcont.replace('.ms', '_model.ms')
#     tclean_wrapper(vis          = initcont_model,
#                   imagename     = initcont_model.replace('.ms', ''), # will become .image, .residual, etc
#                   mask          = SB_mask,
#                   imsize        = SB_imsize,
#                   cellsize      = SB_cellsize,
#                   gain          = 0.05, # nice and low this time; stops us from accidentally cleaning too deep between successive major cycles
#                   contspws      = '', # we want the model to include all spws
#                   threshold     = '0mJy',
#                   interactive   = True,
#                   savemodel     = 'modelcolumn')
# SB_EB1: number (interactive) iterations performed: 2
# SB_EB2: number (interactive) iterations performed: 2
# for EB in data_dict['LB_EBs']:
#     initcont        = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
#     initcont_model  = initcont.replace('.ms', '_model.ms')
#     tclean_wrapper(vis          = initcont_model,
#                   imagename     = initcont_model.replace('.ms', ''), # will become .image, .residual, etc
#                   mask          = LB_mask,
#                   imsize        = LB_imsize,
#                   cellsize      = LB_cellsize,
#                   gain          = 0.05, # nice and low this time; stops us from accidentally cleaning too deep between successive major cycles
#                   contspws      = '', # we want the model to include all spws
#                   threshold     = '0mJy',
#                   interactive   = True,
#                   savemodel     = 'modelcolumn')
# LB_EB1: number (interactive) iterations performed: 2
# LB_EB2: number (interactive) iterations performed: 1
# LB_EB3: number (interactive) iterations performed: 1
# LB_EB4: number (interactive) iterations performed: 2
# LB_EB5: number (interactive) iterations performed: 2
# LB_EB6: number (interactive) iterations performed: 2

""" Define a noise annulus, measure the peak SNR in map """
noise_annulus = "annulus[[%s, %s],['%.2farcsec', '9.00arcsec']]" % (mask_ra, mask_dec, 1.1*mask_diam) # outer annulus radius (9arcsec) almost = FOV

selfcal_dict['starting_model']    = {} # prepare to save this round's results
for EB in data_dict['EBs']:
    initcont        = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
    initcont_model  = initcont.replace('.ms', '_model.ms')
    imagename       = initcont_model.replace('.ms', '.fits')
    image_metrics   = estimate_image_metrics(imagename      = imagename,
                                             disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                             noise_mask     = noise_annulus)
    update_selfcal_dict(selfcal_dict  = selfcal_dict,
                        round         = 'starting_model',
                        EB            = EB,
                        image_metrics = image_metrics)
print("selfcal_dict['starting_model']: ", selfcal_dict['starting_model'])


"""
################ PERFORM INITIAL ROUND OF SELFCAL #################
"""

for EB in data_dict['EBs']:
    initcont          = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
    initcont_model    = initcont.replace('.ms', '_model.ms')
    initcont_selfcal  = initcont.replace('.ms', '_selfcal.ms')
    caltable          = initcont.replace('.ms', '_model.p1')

    vis         = initcont_model
    parentvis   = initcont_model
    caltable    = caltable
    gaintable   = [caltable]
    outputvis   = initcont_selfcal
    gaintype    = 'T' # could've used 'G', we only have 1 polarization mode
    spw         = ''
    refant      = data_dict[EB]['refants_list']
    calmode     = 'p' # phase-only
    combine     = 'scan' # ideally we want to combine nothing, but for a single EB, if you don't combine scans the solutions are poor quality
    solint      = 'inf' #
    minsnr      = 2.5 # don't want this to be too low
    minblperant = 4
    interp      = 'linearPD' # PD is overkill but it's fine
    calwt       = True

    """ Generate the solutions """
    # os.system('touch ' + caltable + '.printlog.log') # can paste printout of number of failed solutions in here later
    # os.system('rm -rf ' + caltable)
    # gaincal(vis=vis,caltable=caltable,gaintype=gaintype, spw=spw, refant=refant,
    #         calmode=calmode, combine=combine, solint=solint, minsnr=minsnr,
    #         minblperant=minblperant)
    # for quantity in ['phase', 'amp', 'SNR']:
    #     for bool in [False, True]:
    #         plot_gaincal_solutions(caltable=caltable, parentvis=vis, quantity=quantity,
    #                                plot_average_soln=bool, solint=solint, minsnr=minsnr,
    #                                spw=spw, observation='0', combine=combine)
    """ Apply the solutions """
    # applycal(vis=vis, spw=spw, gaintable=gaintable, interp=interp,  calwt=calwt)
    # os.system('rm -rf ' + outputvis)
    # split(vis=vis, outputvis=outputvis, datacolumn='corrected')

# SB_EB1: No failed solns
# SB_EB2:
# 1 of 38 solutions flagged due to SNR < 2.5 in spw=2 at 2022/05/15/18:30:36.2
# 3 of 38 solutions flagged due to SNR < 2.5 in spw=3 at 2022/05/15/18:30:40.9
# 2 of 38 solutions flagged due to SNR < 2.5 in spw=4 at 2022/05/15/18:30:41.6
# LB_EB1:
# 3 of 41 solutions flagged due to SNR < 2.5 in spw=1 at 2022/07/17/14:50:07.2
# 4 of 41 solutions flagged due to SNR < 2.5 in spw=2 at 2022/07/17/14:50:07.8
# 6 of 40 solutions flagged due to SNR < 2.5 in spw=3 at 2022/07/17/14:50:00.6
# 5 of 41 solutions flagged due to SNR < 2.5 in spw=4 at 2022/07/17/14:49:54.6
# LB_EB2:
# 9 of 41 solutions flagged due to SNR < 2.5 in spw=0 at 2022/07/19/12:46:01.7
# 19 of 41 solutions flagged due to SNR < 2.5 in spw=1 at 2022/07/19/12:45:56.1
# 16 of 41 solutions flagged due to SNR < 2.5 in spw=2 at 2022/07/19/12:46:04.4
# 12 of 41 solutions flagged due to SNR < 2.5 in spw=3 at 2022/07/19/12:45:40.3
# 18 of 41 solutions flagged due to SNR < 2.5 in spw=4 at 2022/07/19/12:45:58.8
# LB_EB3:
# 2 of 43 solutions flagged due to SNR < 2.5 in spw=0 at 2022/07/20/14:34:45.5
# 14 of 43 solutions flagged due to SNR < 2.5 in spw=1 at 2022/07/20/14:35:04.1
# 13 of 43 solutions flagged due to SNR < 2.5 in spw=2 at 2022/07/20/14:34:49.8
# 18 of 43 solutions flagged due to SNR < 2.5 in spw=3 at 2022/07/20/14:35:07.9
# 9 of 43 solutions flagged due to SNR < 2.5 in spw=4 at 2022/07/20/14:34:43.1
# LB_EB4:
# 9 of 42 solutions flagged due to SNR < 2.5 in spw=1 at 2022/07/21/14:42:23.1
# 7 of 42 solutions flagged due to SNR < 2.5 in spw=2 at 2022/07/21/14:42:23.5
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=3 at 2022/07/21/14:42:30.3
# 10 of 42 solutions flagged due to SNR < 2.5 in spw=4 at 2022/07/21/14:42:04.0
# LB_EB5:
# 7 of 43 solutions flagged due to SNR < 2.5 in spw=1 at 2022/07/22/12:15:04.8
# 7 of 43 solutions flagged due to SNR < 2.5 in spw=2 at 2022/07/22/12:15:13.9
# 5 of 42 solutions flagged due to SNR < 2.5 in spw=3 at 2022/07/22/12:14:38.7
# 7 of 43 solutions flagged due to SNR < 2.5 in spw=4 at 2022/07/22/12:15:11.7
# LB_EB6:
# 7 of 41 solutions flagged due to SNR < 2.5 in spw=1 at 2022/07/22/13:42:19.4
# 8 of 41 solutions flagged due to SNR < 2.5 in spw=2 at 2022/07/22/13:42:22.0
# 7 of 40 solutions flagged due to SNR < 2.5 in spw=3 at 2022/07/22/13:42:17.2
# 11 of 41 solutions flagged due to SNR < 2.5 in spw=4 at 2022/07/22/13:42:10.0


""" Image the results; check the resulting map """
# Interactively shallowly clean, and save the result to the _initcont_model.ms model column
# for EB in data_dict['SB_EBs']:
#     initcont          = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
#     initcont_selfcal  = initcont.replace('.ms', '_selfcal.ms')
#     tclean_wrapper(vis          = initcont_selfcal,
#                   imagename     = initcont_selfcal.replace('.ms', ''), # will become .image, .residual, etc
#                   mask          = SB_mask,
#                   imsize        = SB_imsize,
#                   cellsize      = SB_cellsize,
#                   gain          = 0.05, # nice and low this time; stops us from accidentally cleaning too deep between successive major cycles
#                   contspws      = '', # we want the model to include all spws
#                   threshold     = '0mJy',
#                   interactive   = True,
#                   savemodel     = 'modelcolumn')
# # SB_EB1: number (interactive) iterations performed: 3
# # SB_EB2: number (interactive) iterations performed: 3
# for EB in data_dict['LB_EBs']:
#     initcont          = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
#     initcont_selfcal  = initcont.replace('.ms', '_selfcal.ms')
#     tclean_wrapper(vis          = initcont_selfcal,
#                   imagename     = initcont_selfcal.replace('.ms', ''), # will become .image, .residual, etc
#                   mask          = LB_mask,
#                   imsize        = LB_imsize,
#                   cellsize      = LB_cellsize,
#                   gain          = 0.05, # nice and low this time; stops us from accidentally cleaning too deep between successive major cycles
#                   contspws      = '', # we want the model to include all spws
#                   threshold     = '0mJy',
#                   interactive   = True,
#                   savemodel     = 'modelcolumn')
# LB_EB1: number (interactive) iterations performed: 3
# LB_EB2: number (interactive) iterations performed: 1
# LB_EB3: number (interactive) iterations performed: 2
# LB_EB4: number (interactive) iterations performed: 2
# LB_EB5: number (interactive) iterations performed: 3
# LB_EB6: number (interactive) iterations performed: 3

selfcal_dict['initial_round']    = {} # prepare to save this round's results
for EB in data_dict['EBs']:
    initcont          = data_dict['NRAO_path']+data_dict[EB]['_initcont.ms']
    initcont_selfcal  = initcont.replace('.ms', '_selfcal.ms')
    imagename         = initcont_selfcal.replace('.ms', '.fits')
    image_metrics     = estimate_image_metrics(imagename      = imagename,
                                               disk_mask      = SB_mask, # a bit big for the LB EBs, but it's fine
                                               noise_mask     = noise_annulus)
    update_selfcal_dict(selfcal_dict  = selfcal_dict,
                        round         = 'initial_round',
                        EB            = EB,
                        image_metrics = image_metrics)
print("selfcal_dict['initial_round']: ", selfcal_dict['initial_round'])


# selfcal_dict['starting_model']:  {'SB_EB1': {'beammajor': 0.98555845022208, 'beamminor': 0.6502565145492001, 'beampa': -3.267037153244, 'disk_flux': 96.37119767635464, 'peak_intensity': 15.64478687942028, 'rms': 110.25500110126073, 'SNR': 141.89639221037916}, 'SB_EB2': {'beammajor': 0.84730178117736, 'beamminor': 0.55687248706824, 'beampa': -11.97756671906, 'disk_flux': 92.91614621722621, 'peak_intensity': 10.486830025911331, 'rms': 111.60027276620909, 'SNR': 93.96778131429969}, 'LB_EB1': {'beammajor': 0.25133076310158, 'beamminor': 0.158669948577888, 'beampa': -22.46494483948, 'disk_flux': 102.21922520103227, 'peak_intensity': 2.22735945135355, 'rms': 59.83876249396578, 'SNR': 37.22268573950138}, 'LB_EB2': {'beammajor': 0.2146227657795, 'beamminor': 0.143906086683288, 'beampa': 22.16491508484, 'disk_flux': 64.31607986480157, 'peak_intensity': 0.9940143208950758, 'rms': 74.2708706723907, 'SNR': 13.383636301770037}, 'LB_EB3': {'beammajor': 0.222091585397712, 'beamminor': 0.15012736618518, 'beampa': -20.48298072815, 'disk_flux': 98.33766323017612, 'peak_intensity': 1.7070583999156952, 'rms': 85.4101918477119, 'SNR': 19.986588988810784}, 'LB_EB4': {'beammajor': 0.222729876637452, 'beamminor': 0.152399435639364, 'beampa': -28.83377456665, 'disk_flux': 94.4923150756974, 'peak_intensity': 1.6992578748613596, 'rms': 66.02184294954915, 'SNR': 25.73781341062918}, 'LB_EB5': {'beammajor': 0.230228304862968, 'beamminor': 0.14739404618739602, 'beampa': 29.80331993103, 'disk_flux': 103.8761273112917, 'peak_intensity': 1.949697034433484, 'rms': 58.18971866991047, 'SNR': 33.505868029598496}, 'LB_EB6': {'beammajor': 0.204741835594164, 'beamminor': 0.154411122202884, 'beampa': -2.572570800781, 'disk_flux': 95.24123482288928, 'peak_intensity': 1.741013489663601, 'rms': 59.81558891546768, 'SNR': 29.106350388425138}}

# selfcal_dict['initial_round']:  {'SB_EB1': {'beammajor': 0.98555845022208, 'beamminor': 0.6502565145492001, 'beampa': -3.267035722733, 'disk_flux': 97.37744327442131, 'peak_intensity': 15.772012993693352, 'rms': 81.13996805896602, 'SNR': 194.38032046341846}, 'SB_EB2': {'beammajor': 0.8503259420393999, 'beamminor': 0.55800050497056, 'beampa': -12.06888866425, 'disk_flux': 93.93793275833693, 'peak_intensity': 10.558287613093853, 'rms': 89.668852844978, 'SNR': 117.747548653793}, 'LB_EB1': {'beammajor': 0.25463050603866, 'beamminor': 0.160491451621068, 'beampa': -22.56776428223, 'disk_flux': 101.04137510466822, 'peak_intensity': 2.245823619887233, 'rms': 52.304309448011786, 'SNR': 42.9376401980698}, 'LB_EB2': {'beammajor': 0.47473317384732, 'beamminor': 0.206704244017584, 'beampa': 25.47427558899, 'disk_flux': 56.16015007566346, 'peak_intensity': 2.806705655530095, 'rms': 139.53317041782125, 'SNR': 20.114970849767356}, 'LB_EB3': {'beammajor': 0.26145595312117204, 'beamminor': 0.158696442842472, 'beampa': -10.980427742, 'disk_flux': 96.63220026886405, 'peak_intensity': 1.9986922852694988, 'rms': 86.05921215271813, 'SNR': 23.22461750780008}, 'LB_EB4': {'beammajor': 0.22487176954746002, 'beamminor': 0.15411971509455602, 'beampa': -29.1448764801, 'disk_flux': 96.10158140277751, 'peak_intensity': 1.7663217149674892, 'rms': 66.62911868938559, 'SNR': 26.509756540557014}, 'LB_EB5': {'beammajor': 0.23194327950478802, 'beamminor': 0.14885918796063602, 'beampa': 29.62239074707, 'disk_flux': 104.86656000266917, 'peak_intensity': 1.9505757372826338, 'rms': 50.55783934790038, 'SNR': 38.581073923279504}, 'LB_EB6': {'beammajor': 0.207019567489608, 'beamminor': 0.15601308643818, 'beampa': -2.546737670898, 'disk_flux': 95.5095751180149, 'peak_intensity': 1.7425380647182465, 'rms': 54.00177322236003, 'SNR': 32.268163816456486}}

# import json
# with open('selfcal_dict.txt', 'w') as file:
#     file.write(json.dumps(selfcal_dict))
# with open('selfcal_dict.txt') as f:
#     data = f.read()
# js = json.loads(data)

# sys.exit()
