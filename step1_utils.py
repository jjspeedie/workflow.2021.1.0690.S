"""
Functions useful for Step 1: Prepare the continuum
Based off of reduction_utils.py by the DSHARP program
Adapted for CASA version 6.1.2.7 and the AB Aur program
Adapter: J. Speedie
"""
import os
import numpy as np
import sys

def LSRKvel_to_chan(msfile, field, spw, restfreq, LSRKvelocity):
    """
    Identifies the channel(s) corresponding to input LSRK velocities.
    Useful for choosing which channels to split out or flag if a line is expected to be present
    CASA tasks used:
        tb.open, tb.getcol, tb.close
        ms.open, ms.close, ms.cvelfreqs

    Args:
        msfile (string): Name of measurement set
        spw (int): Spectral window number
        field (string): Field name
        restfreq (float): Rest frequency of the spectral line to be flagged, in Hz
        LSRKvelocity (float or array of floats): input velocity in LSRK frame in km/s to be converted into channel numbers

    Returns:
        Channel number most closely corresponding to input LSRK velocity
    """
    cc = 299792458. # speed of light in m/s

    tb.open(msfile)
    spw_col = tb.getcol('DATA_DESC_ID')
    obs_col = tb.getcol('OBSERVATION_ID')
    tb.close()

    obsid = np.unique(obs_col[np.where(spw_col==spw)])

    tb.open(msfile+'/SPECTRAL_WINDOW')
    chanfreqs = tb.getcol('CHAN_FREQ', startrow = spw, nrow = 1)
    tb.close()

    tb.open(msfile+'/FIELD')
    fieldnames = tb.getcol('NAME')
    fieldid = np.where(fieldnames==field)[0][0]
    tb.close()

    tb.open(msfile+'/OBSERVATION')
    obstime = np.squeeze(tb.getcol('TIME_RANGE', startrow = obsid, nrow = 1))[0]
    tb.close()

    nchan = len(chanfreqs)

    ms.open(msfile)
    # Take the spectral grid of a given spectral window, tranform and regrid it as prescribed by the given grid parameters (same as in cvel and clean) and return the transformed values as a list.
    lsrkfreqs = ms.cvelfreqs(spwids=[spw], fieldids=[fieldid], mode='channel', nchan=nchan, obstime=str(obstime)+'s', start=0, width=1, outframe='LSRK')
    chanvelocities = (restfreq-lsrkfreqs)/restfreq*cc/1.e3 #converted to LSRK velocities in km/s
    ms.close()

    if type(LSRKvelocity)==np.ndarray:
        outchans = np.zeros_like(LSRKvelocity)
        for i in range(len(LSRKvelocity)):
            outchans[i] = np.argmin(np.abs(chanvelocities - LSRKvelocity[i]))
        return outchans
    else:
        return np.argmin(np.abs(chanvelocities - LSRKvelocity))



def get_flagchannels(ms_file, ms_dict, velocity_range = np.array([-20,20])):
    """
    Identify channels to flag based on provided velocity range of the line emission
    Calls:
        LSRKvel_to_chan

    Args:
        ms_file (string): The measurement set you wish to investigate
        ms_dict (dictionary): Dictionary of information about measurement set
        velocity_range np.array([min_velocity, max_velocity]): Velocity range (in km/s) over which line emission has been identified

    Returns:
        String of channels to be flagged, in a format that can be passed to the spw parameter in CASA's flagdata task.
    """
    flagchannels_string = ''
    for j,spw in enumerate(ms_dict['line_spws']):
        chans = LSRKvel_to_chan(ms_file, ms_dict['field'], spw, ms_dict['line_freqs'][j], velocity_range[j])
        if j==0:
            flagchannels_string+='%d:%d~%d' % (spw, np.min([chans[0], chans[1]]), np.max([chans[0], chans[1]]))
        else:
            flagchannels_string+=', %d:%d~%d' % (spw, np.min([chans[0], chans[1]]), np.max([chans[0], chans[1]]))
    print("Flagchannels input string for %s: \'%s\'" % (ms_dict['name'], flagchannels_string) )
    return flagchannels_string



def avg_cont(msfile, outputvis='avg_cont.ms', flagchannels='', datacolumn='data',
            contspws=None, width_array=None, field='AB_Aur'):
    '''
    Produces spectrally averaged (pseudo-) continuum measurement sets.
    CASA tasks used:
        flagmanager
        flagdata
        split

    Args:
        msfile (string): Name of measurement set
        outputvis (string): Name of the spectrally averaged continuum measurement set
        flagchannels (string): Argument to be passed for flagchannels parameter in flagdata task
        datacolumn (string): Column to pull from for continuum averaging (usually will be 'data', but may sometimes be 'corrected' if there was flux rescaling applied)
        contspws (string): Argument to be passed to CASA for the spw parameter in split. If not set, all SPWs will be selected by default.
        width_array (array): Argument to be passed to CASA for the width parameter in split. If not set, all SPWs will be selected by default.
    '''
    # Troubleshoot "Waiting for read-lock on file" https://casaguides.nrao.edu/index.php/Waiting_for_read-lock_on_file
    # clearstat
    # clearstat
    # clearstat
    # clearstat
    # clearstat

    # Before doing anything, plot the original amp vs. chan/freq, so can compare with post-flagging
    for i, _ in enumerate(width_array):
        plotms(vis=msfile, yaxis='amp', xaxis='channel', avgchannel='1', plotrange=[0,0,0,1.8], avgtime='1e8', avgscan=True, spw=str(i), plotfile=outputvis+'_before_lines_flagged_channel_spw'+str(i)+'.png', showgui=False, overwrite=True)
        plotms(vis=msfile, yaxis='amp', xaxis='freq', avgchannel='1', plotrange=[0,0,0,1.8], avgtime='1e8', avgscan=True, spw=str(i), plotfile=outputvis+'_before_lines_flagged_freq_spw'+str(i)+'.png', showgui=False, overwrite=True)

    # clear out old versions of the flags:
    if os.path.isdir(msfile+'.flagversions/flags.before_cont_flags'):
        flagmanager(vis=msfile, mode='delete', versionname='before_cont_flags')

    # Use the flagmanager task to save the state of the data before any flagging is applied:
    flagmanager(vis=msfile, mode='save', versionname='before_cont_flags', comment='Flag states before spectral lines are flagged')

    # Initialize the per-channel (or spectral) weights in the ms using initweights. This step ensures that when the flagged and
    # unflagged channels are combined, the appropriate weighting is given to the final set of averaged channels.
    initweights(vis=msfile, wtmode='weight', dowtsp=True)

    # Use the task flagdata to apply the spectral line flags:
    flagdata(vis=msfile, mode='manual', spw=flagchannels, flagbackup=False, field=field)

    # Check that the flags were applied correctly by using plotms to inspect the flagged ms
    for i, _ in enumerate(width_array):
        plotms(vis=msfile, yaxis='amp', xaxis='channel', avgchannel='1', plotrange=[0,0,0,1.8], avgtime='1e8', avgscan=True, spw=str(i), plotfile=outputvis+'_check_lines_flagged_correctly_channel_spw'+str(i)+'.png', showgui=False, overwrite=True)
        plotms(vis=msfile, yaxis='amp', xaxis='freq', avgchannel='1', plotrange=[0,0,0,1.8], avgtime='1e8', avgscan=True, spw=str(i), plotfile=outputvis+'_check_lines_flagged_correctly_freq_spw'+str(i)+'.png', showgui=False, overwrite=True)

    # Use the task split to average the channels together to produce the spectrally averaged continuum data set:
    os.system('rm -rf '+outputvis)
    split(vis         = msfile,
          field       = field,
          spw         = contspws,
          outputvis   = outputvis,
          width       = width_array, # number of channels to average together
          datacolumn  = datacolumn,
          keepflags   = False)
    print("Spetrally averaged continuum dataset saved to %s" % outputvis)
    listobs(vis=outputvis, listfile=outputvis+'listobs.txt')

    # "Now you should check the weights of the new continuum measurement set. The ratio of the weights value
    # between Time Domain Mode (TDM) and Frequency Domain Mode (FDM) windows should be approximately equal to
    # the ratio of the channel widths." (NRAO Imaging Guide)
    plotms(vis=msfile,yaxis='wtsp',xaxis='freq',avgchannel='1',avgtime='1e8',avgscan=True, plotfile=outputvis+'_check_weights.png', showgui=False)

    # Finally, we need to use the flagmanager tasks to restore the ms file to its original unflagged state, so that
    # later we can do continuum subtraction and line imaging:
    flagmanager(vis=msfile, mode='restore', versionname='before_cont_flags')



def image_each_obs(ms_dict, inputvis, scales, smallscalebias=0.6, mask='', contspws='',
                    threshold='0.5mJy', imsize=None, cellsize=None, interactive=False,
                    robust=0.5, gain=0.1, niter=50000, cycleniter=300):
    """
    Wrapper for tclean that will loop through all the observations in a measurement set and image them individually,
        as well as make images for each spectral window alone
    CASA tasks used:
        tclean
        exportfits

    Args:
        ms_dict (dictionary): Dictionary of information about measurement set
        inputvis (string): Measurement set to image
        See the CASA 6.2.1.7 documentation for tclean to get the definitions of all other parameters
    """
    tb.open(inputvis+'/OBSERVATION')
    num_observations = (tb.getcol('TIME_RANGE')).shape[1] #picked an arbitrary column to count the number of observations
    tb.close()

    tb.open(inputvis)
    spws = np.unique(tb.getcol('DATA_DESC_ID')) # get the spws [0, 1, 2, 3, 4] and not [0,0,0,0...4,4,4,4]
    tb.close()

    for i in range(num_observations):
        observation = '%d' % i
        print('We do not account for multiple EBs within this EB; check observation = 0, does it? observation = ', observation)
        imagename = inputvis.replace('.ms', '')
        for ext in ['.image', '.mask', '.model', '.pb', '.psf', '.residual', '.sumwt']:
            os.system('rm -rf '+ imagename + ext)
        tclean(vis              = inputvis, # ms to image
               imagename        = imagename, # file names preceding .image, .residual, etc.
               observation      = observation, # to be imaging each observation individually
               specmode         = 'mfs', # to make a continuum image
               deconvolver      = 'multiscale',
               scales           = scales,
               weighting        = 'briggs',
               robust           = robust,
               imsize           = imsize,
               cell             = cellsize,
               mask             = mask,
               spw              = contspws,
               niter            = niter, #we want to end on the threshold
               threshold        = '%.2fmJy'%(ms_dict['pipeline_cont_cleanthresh']),
               interactive      = interactive,
               cycleniter       = cycleniter,
               cyclefactor      = 1,
               smallscalebias   = smallscalebias, #set to CASA's default of 0.6 unless manually changed
               gain             = gain,
               nterms           = 1) # Number of Taylor coefficients in the spectral model; nterms=1 : Assume flat spectrum source
        os.system('rm -rf '+ imagename+'.image.fits')
        exportfits(imagename+'.image', imagename+'.image.fits')
        print('Done! Saved fits file: ', imagename+'.image.fits')

        for i in spws:
            spw = '%d' % i # need it to be a string
            print('Now imaging spectral window '+spw+' by itself')
            imagename = inputvis.replace('.ms', '')+'_spw'+spw
            for ext in ['.image', '.mask', '.model', '.pb', '.psf', '.residual', '.sumwt']:
                os.system('rm -rf '+ imagename + ext)
            tclean(vis              = inputvis, # msfile to image
                   imagename        = imagename, # file names preceding .image, .residual, etc.
                   observation      = observation, # to be imaging each observation individually
                   specmode         = 'mfs', # to make a continuum image
                   deconvolver      = 'multiscale',
                   scales           = scales,
                   weighting        = 'briggs',
                   robust           = robust,
                   imsize           = imsize,
                   cell             = cellsize,
                   mask             = mask,
                   spw              = spw, # one spectral window at a time
                   niter            = niter, #we want to end on the threshold
                   threshold        = '%.2fmJy'%(ms_dict['pipeline_cont_cleanthresh_perspw'][0][i]),
                   interactive      = interactive,
                   cycleniter       = cycleniter,
                   cyclefactor      = 1,
                   smallscalebias   = smallscalebias, #set to CASA's default of 0.6 unless manually changed
                   gain             = gain,
                   nterms           = 1) # Number of Taylor coefficients in the spectral model; nterms=1 : Assume flat spectrum source
            os.system('rm -rf '+ imagename+'.image.fits')
            exportfits(imagename+'.image', imagename+'.image.fits')
            print('Done! Saved fits file: ', imagename+'.image.fits')
