"""
Functions useful for self-calibration and investigating self-cal results
Written for CASA 6.2.1.7 and the AB Aur program
Author: J. Speedie
Version 0 November 21 2022

Also contains adaptations of functions from DSHARP LP's reduction_utils.py
"""
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
plt.rc('font', family='sans-serif', style='normal', weight='normal')
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
mpl.rcParams['axes.linewidth'] = 1.4

def update_selfcal_dict(save_dir='./', EB='EB', image_metrics=[]):
    """
    Creates and sets entries in a dictionary to keep an organized record of the
    results of each round of self calibration. Saves to a file.

    Args:
        save_dir (string): Directory to save the dict as a txt file.
        EB (string): First layer of the dictionary: e.g. "LB_EB1", or "LBs"
        image_metrics (list, len=7): The output of estimate_image_metrics().
    Returns:
        selfcal_dict (dictionary)
    """
    import json

    selfcal_dict = {}

    assert len(image_metrics)==7

    selfcal_dict[EB]                     = {}
    selfcal_dict[EB]['beammajor']        = image_metrics[0]
    selfcal_dict[EB]['beamminor']        = image_metrics[1]
    selfcal_dict[EB]['beampa']           = image_metrics[2]
    selfcal_dict[EB]['disk_flux']        = image_metrics[3]
    selfcal_dict[EB]['peak_intensity']   = image_metrics[4]
    selfcal_dict[EB]['rms']              = image_metrics[5]
    selfcal_dict[EB]['SNR']              = image_metrics[6]

    with open(save_dir+'selfcal_dict.txt', 'w') as file:
        file.write(json.dumps(selfcal_dict))

    print("selfcal_dict saved for "+EB+"!")

    return selfcal_dict


def tclean_wrapper(vis, imagename, smallscalebias=0.6, mask='', contspws='',
                    threshold='0.2mJy', imsize=None, cellsize=None, interactive=False,
                    robust=0.5, gain=0.05, niter=50000, uvtaper=[], cycleniter=300,
                    savemodel='none'):
    """
    *** Adapted from reduction_utils.py, by the DSHARP Large Program
    Wrapper for tclean with keywords set to values we desire for self calibration.
    CASA tasks used:
        tclean
        exportfits

    Args:
        See the CASA 6.2.1.7 documentation for tclean to get the definitions of all other parameters
    """

    for ext in ['.image', '.mask', '.model', '.pb', '.psf', '.residual', '.sumwt']:
        os.system('rm -rf '+ imagename + ext)

    tclean(vis              = vis, # msfile to image
           imagename        = imagename, # file names preceding .image, .residual, etc.
           specmode         = 'mfs', # to make a continuum image
           deconvolver      = 'hogbom', # better than multiscale, for continuum rings
           # scales           = scales,
           weighting        = 'briggs',
           robust           = robust,
           imsize           = imsize,
           cell             = cellsize,
           mask             = mask,
           spw              = contspws,
           niter            = niter, # we want to end on the threshold
           threshold        = threshold,
           interactive      = interactive,
           cycleniter       = cycleniter,
           cyclefactor      = 1,
           smallscalebias   = smallscalebias, # set to CASA's default of 0.6 unless manually changed
           gain             = gain,
           nterms           = 1, # Number of Taylor coefficients in the spectral model; nterms=1 : Assume flat spectrum source
           uvtaper          = uvtaper,
           savemodel        = savemodel) # VERY IMPORTANT ARG FOR SELF-CALIBRATION! MUST BE SET TO 'modelcolumn'

    os.system('rm -rf '+ imagename+'.fits')
    exportfits(imagename+'.image', imagename+'.fits')

    os.system('rm -rf '+ imagename+'.residual.fits')
    exportfits(imagename+'.residual', imagename+'.residual.fits')

    os.system('rm -rf '+ imagename+'.model.fits')
    exportfits(imagename+'.model', imagename+'.model.fits')

    os.system('rm -rf '+ imagename+'.psf.fits')
    exportfits(imagename+'.psf', imagename+'.psf.fits')

    os.system('rm -rf '+ imagename+'.pb.fits')
    exportfits(imagename+'.pb', imagename+'.pb.fits')



def estimate_image_metrics(imagename, disk_mask, noise_mask):
    """
    *** Adapted from estimate_SNR() from reduction_utils.py, by the DSHARP Large Program
    Return and print estimates for: beam dimensions, flux inside disk mask, peak intensity,
        rms noise, and peak SNR. Helpful for determining whether to continue with further rounds
        of self-calibration.
    CASA tasks used:
        imhead, imstat

    Args:
        imagename: Image name ending in '.image' (string)
        disk_mask: , in the CASa region format, e.g.
        noise_mask: Annulus to measure image rms, in the CASA region format, e.g. 'annulus[[500pix, 500pix],["1arcsec", "2arcsec"]]' (string)

    Returns:
        image_metrics (list): List of beammajor (arcsec), beamminor (arcsec),
                              beampa (deg), disk_flux (mJy), peak_intensity (mJy/beam),
                              rms (microJy/beam), SNR
    """
    headerlist  = imhead(imagename, mode = 'list')
    beammajor   = headerlist['beammajor']['value']
    beamminor   = headerlist['beamminor']['value']
    beampa      = headerlist['beampa']['value']
    print("# %s" % imagename)
    print("# Beam %.3f arcsec x %.3f arcsec (%.2f deg)" %(beammajor, beamminor, beampa))

    disk_stats  = imstat(imagename = imagename, region = disk_mask)
    disk_flux   = disk_stats['flux'][0]
    print("# Flux inside disk mask: %.2f mJy" %(disk_flux*1000))

    peak_intensity = disk_stats['max'][0]
    print("# Peak intensity of source: %.2f mJy/beam" %(peak_intensity*1000))

    rms = imstat(imagename = imagename, region = noise_mask)['rms'][0]
    print("# rms: %.2f microJy/beam" %(rms*1e6))

    SNR = peak_intensity/rms
    print("# Peak SNR: %.2f" %(SNR))

    image_metrics = [beammajor, beamminor, beampa, disk_flux*1e3, peak_intensity*1e3, rms*1e6, SNR]

    return image_metrics

def retrieve_from_caltable(caltable=None, xaxis='time', yaxis='phase', spw='',
                            observation='0', field='', timerange='', antenna='',
                            uvrange='', intent='', scan='', correlation=''):
    """
    Retrieves 2 columns from a caltable (those specified by xaxis and yaxis)
    using plotms(plotfile='something.txt') and subsequently sourcing the values
    from the something.txt file. It's a hack but it works. (CASA's 'calanalysis'
    is an alternative method, but the dictionaries are hard to interpret.)

    Args:
        caltable (string): The caltable from you wish to retrieve the data.
        xaxis (string): The first column you wish to retrieve.
        yaxis (string): The second column you wish to retrieve.
        All other arguments are like CASA's plotms.
    Returns:
        xaxis (array): 1D array of the caltable column specified by xaxis
            input arg. Possibilities are only: 'time' or 'ant1'
        yaxis (array): 1D array of the caltable column specified by yaxis
            input arg.
        xaxis_str (array): 1D array of the antenna1 column of the caltable.
            Ignore if xaxis='time'
    """
    if caltable is None:
        raise ValueError('You need to specify a caltable')

    plotms(vis=caltable, xaxis=xaxis, yaxis=yaxis, spw=spw, observation=observation,
            field=field, timerange=timerange, antenna=antenna, uvrange=uvrange,
            intent=intent, scan=scan, correlation=correlation,
            plotfile='temporary_file_for_retrieve_from_caltable.txt',
            showgui=False, overwrite=True)

    if xaxis=='time':
        caltable_txt = np.loadtxt('temporary_file_for_retrieve_from_caltable.txt',
                                comments='#', usecols=(1,9))
        xaxis = caltable_txt[:,1]
        yaxis = caltable_txt[:,0]

    elif xaxis=='antenna1':
        xaxis = np.loadtxt('temporary_file_for_retrieve_from_caltable.txt',
                                comments='#', usecols=(5))
        yaxis = np.loadtxt('temporary_file_for_retrieve_from_caltable.txt',
                                comments='#', usecols=(1))
    xaxis_str = np.loadtxt('temporary_file_for_retrieve_from_caltable.txt',
                            comments='#', usecols=(7), dtype=type(''))

    os.system('rm temporary_file_for_retrieve_from_caltable.txt')

    return xaxis, yaxis, xaxis_str


def get_scan_start_and_end_times(vis=None, observation='0'):
    """
    Retrieves the start and end times of the scans in an execution block.

    Args:
        vis (string): The measurement set whose scans you wish to get.
        observation (string): If the measurement set is a concatenation of
            execution blocks, then observation will identify their indices (like
            CASA's 'observation' parameter).
    Returns:
        scan_start_and_end_times (array): An array of the start and end times
            (in units of MJD seconds) of every scan in obs. The array has shape:
            (number of scans)x(2).
        num_scans (int): The number of scans, for convenience.
        scans (array): A list of the scans, with original index, for convenience.
    """
    if vis is None:
        raise ValueError('You need to specify a measurement set')

    tb.open(vis)
    scan_col_all    = tb.getcol('SCAN_NUMBER')
    obs_col_all     = tb.getcol('OBSERVATION_ID')
    time_col_all    = tb.getcol('TIME')
    tb.close()

    scan_col        = scan_col_all[np.where(obs_col_all==int(observation))]
    time_col        = time_col_all[np.where(obs_col_all==int(observation))]

    scans           = np.unique(scan_col)
    num_scans       = len(scans)

    scan_start_and_end_times = np.zeros((num_scans, 2))
    for i,scan in enumerate(scans):
        scan_start_and_end_times[i, 0] = time_col[np.where(scan_col==scan)][0]
        scan_start_and_end_times[i, 1] = time_col[np.where(scan_col==scan)][-1]

    return scan_start_and_end_times, num_scans, scans


def plot_gaincal_solutions(caltable=None, parentvis=None, quantity='phase',
                           plot_average_soln=False, solint='120', minsnr=2.5,
                           spw='', observation='0', combine='', field='',
                           timerange='', antenna='', uvrange='', intent='',
                           scan='', correlation='', calmode='p'):
    """
    Plots a quantity of the calibration table vs. time and saves the figure.

    Args:
        caltable (string): The caltable whose solutions you wish to plot.
        parentvis (string): The measurement set from which caltable was generated.
        quantity (string): The quantity you wish to plot on the y-axis, against
            time. Only available options: 'phase', 'amp', 'SNR'.
        plot_average_soln (bool): Whether or not to plot the average of quantity,
            ie., one point at each time step.
        solint (string): The solution interval in gaincal on which the phase solutions
            were generated.
        observation (string): If the measurement set is a concatenation of
            execution blocks, then observation will identify their indices (like
            CASA's 'observation' parameter). *NOTE* This must NOT be ''. It must
            be a single observation at a time, purely because otherwise the xaxis
            of the plot will be too stretched to understand.
        All other arguments are like CASA's plotms and gaincal.
    Returns:
        Figure (png): Saves the figure to a png file, named from the input args.
    """
    if ((quantity!='phase') & (quantity!='amp') & (quantity!='SNR')):
        raise ValueError("Sorry, you can't plot "+quantity+" with 'plot_gaincal_solutions'.")
    if caltable is None:
        raise ValueError('You need to specify a caltable')
    if parentvis is None:
        raise ValueError('You need to specify a measurement set')

# *** To be written
# check that scan, spw, are present in the caltable
# raise RuntimeError("The caltable does not contain spw ")

    # Times of start and end of each scan inside observation (EB) inside parentvis
    scan_start_and_end_times, num_scans, scans = get_scan_start_and_end_times(vis=parentvis, observation=observation)

    # To be able to plot the solutions with different colours for each spectral window:
    if spw=='':
        tb.open(parentvis)
        spw_col = tb.getcol('DATA_DESC_ID')
        obs_col = tb.getcol('OBSERVATION_ID')
        tb.close()
        spws = np.unique(spw_col[np.where(obs_col==int(observation))])
        print("Spectral windows in EB "+observation+" are:", spws)
    else:
        spws = np.array([int(spw)]) # for plotting purposes; even if there's just 1 spw, we still loop over a list
        print("Spectral windows in EB "+observation+" are:", spws)
    spws_reduced = spws - 5*int(observation) # for color indexing purposes
    print("Your parentvis may contain more than one EB. The 'reduced' form of your plotted spws is: ", spws_reduced)

    # Colors in which solutions of each spectral window will be plotted
    colors   = ['r', 'orange', 'mediumseagreen', 'steelblue', 'purple']

    fig = plt.figure(figsize=(9, 7))
    gs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
    ax = fig.add_subplot(gs[0,0])
    ax.tick_params(axis='both', direction='in', which='major', color='k', width=1.5, labelsize=16, length=8, top=True, bottom=True,right=True,left=True)
    ax.tick_params(axis='both', direction='in', which='minor', color='k', width=1.2, labelsize=16, length=4, top=True, bottom=True,right=True,left=True)
    ax.xaxis.set_minor_locator(MultipleLocator(1)) # 1 minute
    ax.xaxis.set_major_locator(MultipleLocator(5)) # 5 minutes
    ax.set_xlabel(r'time (MJD minutes)', fontsize=20)

    ax.text(0.97, 0.92, 'combine='+combine, color='k', transform=ax.transAxes, fontsize=15, horizontalalignment='right', verticalalignment='top')
    ax.text(0.97, 0.89, 'solint='+(solint), color='k', transform=ax.transAxes, fontsize=15, horizontalalignment='right', verticalalignment='top')
    ax.text(0.97, 0.86, 'calmode='+calmode, color='k', transform=ax.transAxes, fontsize=15, horizontalalignment='right', verticalalignment='top')
    EB = str(int(observation)+1)
    ax.text(0.97, 0.8, 'EB'+EB, color='k', transform=ax.transAxes, fontsize=17, horizontalalignment='right', verticalalignment='top')

    ax.set_title(caltable, fontsize=12)

    # Plot the parentvis's scan ranges in the background (ie. where gaincal solutions should be within)
    for s in range(num_scans):
        ax.axvspan(scan_start_and_end_times[s, 0]/60, scan_start_and_end_times[s, 1]/60, alpha=0.3, color='skyblue')

    for i,spw_i in enumerate(spws):
        time, quantity_vals, _ = retrieve_from_caltable(caltable=caltable, xaxis='time', yaxis=quantity, spw=str(spw_i),
                                    observation=observation, field=field, timerange=timerange, antenna=antenna,
                                    uvrange=uvrange, intent=intent, scan=scan, correlation=correlation)

        if plot_average_soln:
            time_avg = '_time-averaged'
            time_points = np.unique(time)
            quantity_vals_avg = np.zeros_like(time_points)
            for t, time_point in enumerate(time_points):
                quantity_vals_avg[t] = np.mean(quantity_vals[np.where(time==time_point)])
            ax.plot(time_points/60, quantity_vals_avg, markersize=40, lw=1.5, color=colors[spws_reduced[i]], zorder=1000000)
            ax.scatter(time_points/60, quantity_vals_avg, s=5, color=colors[spws_reduced[i]], zorder=1000000)
        else:
            time_avg = ''
            if solint!='inf':
                solint_str = solint.replace('s', '')
                solint_float = float(solint_str)
                ax.errorbar(time/60, quantity_vals, xerr=(solint_float/2)/60, linestyle='', lw=0.4, color=colors[spws_reduced[i]])
            ax.scatter(time/60, quantity_vals, s=5, color=colors[spws_reduced[i]], zorder=1000000)
        ax.text(0.03, 0.925-(i*0.035), 'spw '+str(spw_i)+' ('+str(spws_reduced[i])+')', color=colors[spws_reduced[i]], transform=ax.transAxes, fontsize=15, horizontalalignment='left', verticalalignment='top')

    if quantity=='SNR':
        ax.set_ylim(0, 100)
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.grid(alpha=0.5, which='both')
        ax.set_ylabel(r'SNR', fontsize=20)
        ax.axhline(minsnr, color='k', lw=0.8)
        for s in range(num_scans):
            ax.text(np.mean([scan_start_and_end_times[s, 0]/60, scan_start_and_end_times[s, 1]/60]), 0.95*100, str(scans[s]), horizontalalignment='center', color='k', fontsize=16)
    if quantity=='phase':
        ax.set_ylim(-180, 180)
        ax.yaxis.set_major_locator(MultipleLocator(30))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax.grid(alpha=0.5, which='both')
        ax.set_ylabel(r'gain phase (degrees)', fontsize=20)
        for s in range(num_scans):
            ax.text(np.mean([scan_start_and_end_times[s, 0]/60, scan_start_and_end_times[s, 1]/60]), 0.9*(180), str(scans[s]), horizontalalignment='center', color='k', fontsize=16)
    if quantity=='amp':
        ax.set_ylim(0, 3)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.grid(alpha=0.5, which='both')
        ax.set_ylabel(r'amplitude', fontsize=20)
        for s in range(num_scans):
            ax.text(np.mean([scan_start_and_end_times[s, 0]/60, scan_start_and_end_times[s, 1]/60]), 0.95*3, str(scans[s]), horizontalalignment='center', color='k', fontsize=16)

    gs.update(wspace=0.0, hspace=0.0) # set the spacing between axes.
    plt.tight_layout()
    if len(spws)!=1:
        filename = caltable+'_'+quantity+'-vs-time_EB'+EB+'_allspws'+time_avg+'.png'
    else:
        filename = caltable+'_'+quantity+'-vs-time_EB'+EB+'_spw'+spw+time_avg+'.png'
    plt.savefig(filename, dpi=500, transparent=True, bbox_inches='tight',pad_inches=0.015)
    plt.clf()
    print("Figure saved! To: "+filename)


def plot_gaincal_solutions_per_antenna(caltable=None, parentvis=None, quantity='phase',
                           plot_average_soln=False, solint=120, minsnr=2.5,
                           spw='', observation='0', combine='', field='',
                           timerange='', antenna='', uvrange='', intent='',
                           scan='', correlation='', calmode='p'):
    """
    Plots a quantity of the calibration table vs. time and saves the figure.

    Args:
        caltable (string): The caltable whose solutions you wish to plot.
        parentvis (string): The measurement set from which caltable was generated.
        quantity (string): The quantity you wish to plot on the y-axis, against
            time. Only available options: 'phase', 'amp', 'SNR'.
        plot_average_soln (bool): Whether or not to plot the average of quantity,
            ie., one point at each time step.
        solint (string): The solution interval in gaincal on which the phase solutions
            were generated.
        observation (string): If the measurement set is a concatenation of
            execution blocks, then observation will identify their indices (like
            CASA's 'observation' parameter). *NOTE* This must NOT be ''. It must
            be a single observation at a time, purely because otherwise the xaxis
            of the plot will be too stretched to understand.
        All other arguments are like CASA's plotms and gaincal.
    Returns:
        Figure (png): Saves the figure to a png file, named from the input args.
    """
    if ((quantity!='phase') & (quantity!='amp') & (quantity!='SNR')):
        raise ValueError("Sorry, you can't plot "+quantity+" with 'plot_gaincal_solutions'.")
    if caltable is None:
        raise ValueError('You need to specify a caltable')
    if parentvis is None:
        raise ValueError('You need to specify a measurement set')

# *** To be written
# check that scan, spw, are present in the caltable
# raise RuntimeError("The caltable does not contain spw ")

    # Times of start and end of each scan inside observation (EB) inside parentvis
    scan_start_and_end_times, num_scans, scans = get_scan_start_and_end_times(vis=parentvis, observation=observation)

    # To be able to plot the solutions with different colours for each spectral window:
    if spw=='':
        tb.open(parentvis)
        spw_col = tb.getcol('DATA_DESC_ID')
        obs_col = tb.getcol('OBSERVATION_ID')
        tb.close()
        spws = np.unique(spw_col[np.where(obs_col==int(observation))])
        print("Spectral windows in EB "+observation+" are:", spws)
    else:
        spws = np.array([int(spw)]) # for plotting purposes; even if there's just 1 spw, we still loop over a list
        print("Spectral windows in EB "+observation+" are:", spws)
    spws_reduced = spws - 5*int(observation) # for color indexing purposes
    print("Your parentvis may contain more than one EB. The 'reduced' form of your plotted spws is: ", spws_reduced)

    # Colors in which solutions of each spectral window will be plotted
    colors   = ['r', 'orange', 'mediumseagreen', 'steelblue', 'purple']

    fig = plt.figure(figsize=(9, 7))
    gs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
    ax = fig.add_subplot(gs[0,0])
    ax.tick_params(axis='both', direction='in', which='major', color='k', width=1.5, labelsize=16, length=8, top=True, bottom=True,right=True,left=True)
    ax.tick_params(axis='x', direction='in', which='major', color='k', width=1.5, labelsize=8, length=8, top=True, bottom=True,right=True,left=True)
    ax.tick_params(axis='both', direction='in', which='minor', color='k', width=1.2, labelsize=16, length=4, top=True, bottom=True,right=True,left=True)
    # ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_major_locator(MultipleLocator(1)) # 1 antenna
    ax.set_xlabel(r'antenna1', fontsize=20)

    ax.text(0.97, 0.92, 'combine='+combine, color='k', transform=ax.transAxes, fontsize=15, horizontalalignment='right', verticalalignment='top')
    ax.text(0.97, 0.89, 'solint='+(solint), color='k', transform=ax.transAxes, fontsize=15, horizontalalignment='right', verticalalignment='top')
    ax.text(0.97, 0.86, 'calmode='+calmode, color='k', transform=ax.transAxes, fontsize=15, horizontalalignment='right', verticalalignment='top')
    EB = str(int(observation)+1)
    ax.text(0.97, 0.8, 'EB'+EB, color='k', transform=ax.transAxes, fontsize=17, horizontalalignment='right', verticalalignment='top')

    ax.set_title(caltable, fontsize=12)

    # Plot the parentvis's scan ranges in the background (ie. where gaincal solutions should be within)
    # for s in range(num_scans):
    #     ax.axvspan(scan_start_and_end_times[s, 0]/60, scan_start_and_end_times[s, 1]/60, alpha=0.3, color='skyblue')

    for i,spw_i in enumerate(spws):
        antennas, quantity_vals, ant_names = retrieve_from_caltable(caltable=caltable, xaxis='antenna1', yaxis=quantity, spw=str(spw_i),
                                    observation=observation, field=field, timerange=timerange, antenna=antenna,
                                    uvrange=uvrange, intent=intent, scan=scan, correlation=correlation)



        # if plot_average_soln:
        #     time_avg = '_time-averaged'
        #     time_points = np.unique(time)
        #     quantity_vals_avg = np.zeros_like(time_points)
        #     for t, time_point in enumerate(time_points):
        #         quantity_vals_avg[t] = np.mean(quantity_vals[np.where(time==time_point)])
        #     ax.plot(time_points/60, quantity_vals_avg, markersize=40, lw=1.5, color=colors[spws_reduced[i]], zorder=1000000)
        #     ax.scatter(time_points/60, quantity_vals_avg, s=5, color=colors[spws_reduced[i]], zorder=1000000)
        # else:
            # time_avg = ''
        # if solint!='inf':
        #     solint_str = solint.replace('s', '')
        #     solint = float(solint_str)
        #     ax.errorbar(antennas, quantity_vals, xerr=(solint/2)/60, linestyle='', lw=0.4, color=colors[spws_reduced[i]])
        ax.scatter(antennas, quantity_vals, s=5, color=colors[spws_reduced[i]], zorder=1000000)
        ax.text(0.03, 0.925-(i*0.035), 'spw '+str(spw_i)+' ('+str(spws_reduced[i])+')', color=colors[spws_reduced[i]], transform=ax.transAxes, fontsize=15, horizontalalignment='left', verticalalignment='top')


    if quantity=='SNR':
        ax.set_ylim(0, 100)
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(2))
        ax.grid(alpha=0.5, which='both')
        ax.set_ylabel(r'SNR', fontsize=20)
        ax.axhline(minsnr, color='k', lw=0.8)
        for a,ant in enumerate(ant_names):
            if a>np.max(antennas):
                continue
            ax.text(a, 70+((a%2)*3), ant[0:4], fontsize=8, horizontalalignment='center', verticalalignment='bottom')
    if quantity=='phase':
        ax.set_ylim(-180, 180)
        ax.yaxis.set_major_locator(MultipleLocator(30))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax.grid(alpha=0.5, which='both')
        ax.set_ylabel(r'gain phase (degrees)', fontsize=20)
        for a,ant in enumerate(ant_names):
            if a>np.max(antennas):
                continue
            ax.text(a, -170+((a%2)*7), ant[0:4], fontsize=8, horizontalalignment='center', verticalalignment='bottom')
    if quantity=='amp':
        ax.set_ylim(0, 3)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.grid(alpha=0.5, which='both')
        ax.set_ylabel(r'amplitude', fontsize=20)
        for a,ant in enumerate(ant_names):
            if a>np.max(antennas):
                continue
            ax.text(a, 0.1+((a%2)*0.05), ant[0:4], fontsize=8, horizontalalignment='center', verticalalignment='bottom')

    gs.update(wspace=0.0, hspace=0.0) # set the spacing between axes.
    plt.tight_layout()
    if len(spws)!=1:
        filename = caltable+'_'+quantity+'-vs-antenna_EB'+EB+'_allspws'+'.png'
    else:
        filename = caltable+'_'+quantity+'-vs-antenna_EB'+EB+'_spw'+spw+'.png'
    plt.savefig(filename, dpi=500, transparent=True, bbox_inches='tight',pad_inches=0.015)
    plt.clf()
    print("Figure saved! To: "+filename)


def export_MS(msfile):
    """
    *** Adapted from reduction_utils.py, by the DSHARP Large Program (Jane Huang)
    Spectrally averages visibilities to a single channel per SPW and exports to .npz file

    msfile: Name of CASA measurement set, ending in '.ms' (string)
    """
    filename = msfile
    if filename[-3:]!='.ms':
        print("MS name must end in '.ms'")
        return
    # strip off the '.ms'
    MS_filename = filename.replace('.ms', '')

    # get information about spectral windows
    tb.open(MS_filename+'.ms/SPECTRAL_WINDOW')
    num_chan = tb.getcol('NUM_CHAN').tolist()
    tb.close()

    # spectral averaging (1 channel per SPW)
    os.system('rm -rf %s' % MS_filename+'_spavg.ms')
    split(vis=MS_filename+'.ms', width=num_chan, datacolumn='data',
    outputvis=MS_filename+'_spavg.ms')

    # get the data tables
    tb.open(MS_filename+'_spavg.ms')
    data   = np.squeeze(tb.getcol("DATA"))
    flag   = np.squeeze(tb.getcol("FLAG"))
    uvw    = tb.getcol("UVW")
    weight = tb.getcol("WEIGHT")
    spwid  = tb.getcol("DATA_DESC_ID")
    tb.close()

    # get frequency information
    tb.open(MS_filename+'_spavg.ms/SPECTRAL_WINDOW')
    freqlist = np.squeeze(tb.getcol("CHAN_FREQ"))
    tb.close()

    # get rid of any flagged columns
    good   = np.squeeze(np.any(flag, axis=0)==False)
    data   = data[:,good]
    weight = weight[:,good]
    uvw    = uvw[:,good]
    spwid = spwid[good]

    # compute spatial frequencies in lambda units
    get_freq = lambda ispw: freqlist[ispw]
    freqs = get_freq(spwid) #get spectral frequency corresponding to each datapoint
    u = uvw[0,:] * freqs / 2.9979e8
    v = uvw[1,:] * freqs / 2.9979e8

    #average the polarizations
    Re  = np.sum(data.real*weight, axis=0) / np.sum(weight, axis=0)
    Im  = np.sum(data.imag*weight, axis=0) / np.sum(weight, axis=0)
    Vis = Re + 1j*Im
    Wgt = np.sum(weight, axis=0)

    #output to npz file and delete intermediate measurement set
    os.system('rm -rf %s' % MS_filename+'_spavg.ms')
    os.system('rm -rf '+MS_filename+'.vis.npz')
    np.savez(MS_filename+'.vis', u=u, v=v, Vis=Vis, Wgt=Wgt)
    print("#Measurement set exported to %s" % (MS_filename+'.vis.npz',))


def deproject_vis(data, bins=np.array([0.]), incl=0., PA=0., offx=0., offy=0.,
                  errtype='mean'):
    """
    *** Adapted from reduction_utils.py, by the DSHARP Large Program (Jane Huang)
    Deprojects and azimuthally averages visibilities

    Parameters
    ==========
    data: Length-4 tuple of u,v, visibilities, and weight arrays
    bins: 1-D array of uv distance bins (kilolambda)
    incl: Inclination of disk (degrees)
    PA: Position angle of disk (degrees)
    offx: Horizontal offset of disk center from phase center (arcseconds)
    offy: Vertical offset of disk center from phase center (arcseconds)

    Returns
    =======
    uv distance bins (1D array), visibilities (1D array), errors on averaged visibilities (1D array)
    """

    # - read in, parse data
    u, v, vis, wgt = data
    # - convert keywords into relevant units
    inclr = np.radians(incl)
    PAr = 0.5*np.pi-np.radians(PA)
    offx *= -np.pi/(180.*3600.)
    offy *= -np.pi/(180.*3600.)

    # - change to a deprojected, rotated coordinate system
    uprime = (u*np.cos(PAr) + v*np.sin(PAr))
    vprime = (-u*np.sin(PAr) + v*np.cos(PAr)) * np.cos(inclr)
    rhop = np.sqrt(uprime**2 + vprime**2)

    # - phase shifts to account for offsets
    shifts = np.exp(-2.*np.pi*1.0j*(u*-offx + v*-offy))
    visp = vis*shifts
    realp = visp.real
    imagp = visp.imag

    # - if requested, return a binned (averaged) representation
    if (bins.size > 1.):
        avbins = 1e3*bins	# scale to lambda units (input in klambda)
        bwid = 0.5*(avbins[1]-avbins[0])
        bvis = np.zeros_like(avbins, dtype='complex')
        berr = np.zeros_like(avbins, dtype='complex')
        for ib in np.arange(len(avbins)):
            inb = np.where((rhop >= avbins[ib]-bwid) & (rhop < avbins[ib]+bwid))
            if (len(inb[0]) >= 5):
                bRe, eRemu = np.average(realp[inb], weights=wgt[inb],
                                        returned=True)
                eRese = np.std(realp[inb])
                bIm, eImmu = np.average(imagp[inb], weights=wgt[inb],
                                        returned=True)
                eImse = np.std(imagp[inb])
                bvis[ib] = bRe+1j*bIm
                if (errtype == 'scat'):
                    berr[ib] = eRese+1j*eImse
                else: berr[ib] = 1./np.sqrt(eRemu)+1j/np.sqrt(eImmu)
            else:
                bvis[ib] = 0+1j*0
                berr[ib] = 0+1j*0
        parser = np.where(berr.real != 0)
        output = avbins[parser], bvis[parser], berr[parser]
        return output

    # - if not, returned the unbinned representation
    output = rhop, realp+1j*imagp, 1./np.sqrt(wgt)

    return output

def plot_deprojected(filelist, fignametemplate='./out', incl=0, PA=0,
                     offx=0, offy=0, fluxscale=None, uvbins=None, show_err=True):
    """
    *** Adapted from reduction_utils.py, by the DSHARP Large Program (Jane Huang)
    Plots real and imaginary deprojected visibilities from a list of .npz files

    Args:
        filelist: List of names of .npz files storing visibility data
        fignametemplate (string): Figure will be save as 'fignametemplate_plot_deprojected.png'
        incl: Inclination of disk (degrees)
        PA: Position angle of disk (degrees)
        offx: Horizontal offset of disk center from phase center (arcseconds)
        offy: Vertical offset of disk center from phase center (arcseconds)
        fluxscale: List of scaling factors to multiply the visibility values by before plotting. Default value is set to all ones.
        uvbins: Array of bins at which to plot the visibility values, in lambda. By default, the range plotted will be from 10 to 1000 kilolambda
        show_err: If True, plot error bars.
    Returns:
        Figure (png): Saves the figure to a png file, named fignametemplate_plot_deprojected.png
    """
    if fluxscale is None:
        fluxscale = np.ones(len(filelist))
    assert len(filelist)==len(fluxscale)

    if uvbins is None:
        uvbins = 10.+10.*np.arange(250)
        # uvbins = 10.+5.*np.arange(250*2)

    minvis = np.zeros(len(filelist))
    maxvis = np.zeros(len(filelist))

    fig = plt.figure(figsize=(9, 14))
    gs = fig.add_gridspec(ncols=1, nrows=2, width_ratios=[1], height_ratios=[1,1])
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    for ax in [ax1, ax2]:
        ax.tick_params(axis='both', direction='in', which='major', color='k', width=1.5, labelsize=16, length=8, top=True, bottom=True,right=True,left=True)
        ax.tick_params(axis='both', direction='in', which='minor', color='k', width=1.2, labelsize=16, length=4, top=True, bottom=True,right=True,left=True)

    for i, filename in enumerate(filelist):

        # read in the data
        inpf = np.load(filename)
        u    = inpf['u']
        v    = inpf['v']
        vis  = fluxscale[i]*inpf['Vis']
        wgt  = inpf['Wgt']

        # deproject the visibilities and do the annular averaging
        vp   = deproject_vis([u, v, vis, wgt], bins=uvbins, incl=incl, PA=PA,
                         offx=offx, offy=offy)
        vp_rho, vp_vis, vp_sig = vp

        # calculate min, max of deprojected, averaged reals (for visualization)
        minvis[i] = np.min(vp_vis.real)
        maxvis[i] = np.max(vp_vis.real)

        # plot the profile
        directory,name = os.path.split(filename)
        if show_err:
            ax1.errorbar(1e-3*vp_rho, vp_vis.real, yerr = vp_sig.real, label = name, fmt = '.')
            ax2.errorbar(1e-3*vp_rho, vp_vis.imag, yerr = vp_sig.imag, label = name, fmt = '.')
        else:
            ax1.plot(1e-3*vp_rho, vp_vis.real, 'o', markersize=2.8, label = name)
            ax2.plot(1e-3*vp_rho, vp_vis.imag, 'o', markersize=2.8, label = name)

    allmaxvis = np.max(maxvis)
    allminvis = np.min(minvis)
    if ((allminvis < 0) or (allminvis-0.1*allmaxvis < 0)):
        ax1.axis([np.min(uvbins)*0.9 , np.max(uvbins), allminvis-0.1*allmaxvis, 1.1*allmaxvis])
        ax2.axis([np.min(uvbins)*0.9 , np.max(uvbins), 0.2*(allminvis-0.1*allmaxvis), 0.2*(1.1*allmaxvis)])
    else:
        ax1.axis([np.min(uvbins)*0.9 , np.max(uvbins), 0., 0.2*(1.1*allmaxvis)])
        ax2.axis([np.min(uvbins)*0.9 , np.max(uvbins), 0., 0.2*(1.1*allmaxvis)])

    ax1.axhline(0, color='k', linestyle='dashed')
    ax2.axhline(0, color='k', linestyle='dashed')
    ax2.set_xlabel(r'deprojected baseline length [kilo$\lambda$]', fontsize=20)
    ax1.set_ylabel('average real [Jy]', fontsize=20)
    ax2.set_ylabel('average imag [Jy]', fontsize=20)
    ax1.legend()
    ax1.grid(alpha=0.5, which='both')
    ax2.grid(alpha=0.5, which='both')
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    gs.update(wspace=0.0, hspace=0.1) # set the spacing between axes.
    plt.tight_layout()
    plt.savefig(fignametemplate+'_plot_deprojected.png', dpi=300, transparent=True, bbox_inches='tight',pad_inches=0.015)
    plt.clf()
    print("Figure saved! To: ", fignametemplate+'_plot_deprojected.png')


def estimate_flux_scale(reference, comparison, incl=0, PA=0, uvbins=None, offx=0,
                        offy=0, make_plot=True, fignametemplate='./out'):
    """
    *** Adapted from reduction_utils.py, by the DSHARP Large Program (Jane Huang)
    Calculates the weighted average of the flux ratio between two observations of a source
    The minimum baseline compared is the longer of the minimum baselines in the individual datasets
    The longest baseline compared is either the shorter of the longest baselines in the individual datasets, or 800 kilolambda

    Really, for our purposes, what it does it show you the phase decoherence on long baselines.
    Useful for diagnostics after phase alignment and during self calibration.

    Args:
        reference: Name of .npz file holding the reference dataset (with the "correct" flux")
        comparison: Name of .npz file holding the comparison dataset (with the flux ratio being checked)
        incl: Inclination of disk (degrees)
        PA: Position angle of disk (degrees)
        offx: Horizontal offset of disk center from phase center (arcseconds)
        offy: Vertical offset of disk center from phase center (arcseconds)
        uvbins: Array of bins at which to compare the visibility values, in lambda.
                By default, the minimum baseline compared is the longer of the minimum baselines in the individual datasets.
                The longest baseline compared is either the shorter of the longest baselines in the individual datasets, or 800 kilolambda, whichever comes first.
        make_plot (bool): Whether or not to make the figure.
        fignametemplate (string): Figure will be save as 'fignametemplate_estimate_flux_scale.png'
    Returns:
        ratio_avg (float): The ratio of the fluxes between the comparison and reference (comp/ref)
        np.sqrt(ratio_avg) (float): The scaling factor for gencal for your comparison measurement set
        Figure (png): Saves the figure to a png file, named fignametemplate_plot_deprojected.png
    """

    inpf = np.load(reference)
    u_ref    = inpf['u']
    v_ref    = inpf['v']
    vis_ref  = inpf['Vis']
    wgt_ref  = inpf['Wgt']

    inpf = np.load(comparison)
    u_comp    = inpf['u']
    v_comp    = inpf['v']
    vis_comp  = inpf['Vis']
    wgt_comp  = inpf['Wgt']

    uvdist_ref = np.sqrt(u_ref**2+v_ref**2)
    uvdist_comp = np.sqrt(u_comp**2+v_comp**2)

    mindist = np.max(np.array([np.min(uvdist_ref), np.min(uvdist_comp)]))
    maxdist = np.min(np.array([np.max(uvdist_ref), np.max(uvdist_ref), 8e5])) #the maximum baseline we want to compare is the longest shared baseline or 800 kilolambda, whichever comes first (we don't want to go out to a baseline that's too long because phase decorrelation becomes a bigger issue at longer baselines.

    if uvbins is None:
        uvbins = mindist/1.e3+10.*np.arange(np.floor((maxdist-mindist)/1.e4))

    # deproject the visibilities and do the annular averaging
    vp   = deproject_vis([u_ref, v_ref, vis_ref, wgt_ref], bins=uvbins, incl=incl, PA=PA,
                         offx=offx, offy=offy)
    ref_rho, ref_vis, ref_sig = vp

    # deproject the visibilities and do the annular averaging
    vp   = deproject_vis([u_comp, v_comp, vis_comp, wgt_comp], bins=uvbins, incl=incl, PA=PA,
                         offx=offx, offy=offy)
    comp_rho, comp_vis, comp_sig = vp

    # maxlen = np.min(np.array([len(comp_rho), len(ref_rho)])) # not sure what this is for; never gets used

    rho_intersection = np.intersect1d(ref_rho, comp_rho) #we only want to compare overlapping baseline intervals

    comp_sig_intersection = comp_sig[np.where(np.in1d(comp_rho, rho_intersection))].real #they're the same for the real and imaginary components
    comp_vis_intersection = comp_vis[np.where(np.in1d(comp_rho, rho_intersection))]
    ref_sig_intersection = ref_sig[np.where(np.in1d(ref_rho, rho_intersection))].real
    ref_vis_intersection = ref_vis[np.where(np.in1d(ref_rho, rho_intersection))]

    ratio = np.abs(comp_vis_intersection)/np.abs(ref_vis_intersection)
    err = ratio*np.sqrt((comp_sig_intersection/np.abs(comp_vis_intersection))**2+(ref_sig_intersection/np.abs(ref_vis_intersection))**2)

    w = 1/err**2
    ratio_avg =  np.sum(w*ratio)/np.sum(w)
    print("#The ratio of the fluxes of %s to %s is %.5f" % (comparison, reference, ratio_avg))
    print("#The scaling factor for gencal is %.3f for your comparison measurement" % (np.sqrt(ratio_avg)))
    print("#The error on the weighted mean ratio is %.3e, although it's likely that the weights in the measurement sets are off by some constant factor" % (1/np.sqrt(np.sum(w)),))

    if make_plot==True:
        fig = plt.figure(figsize=(9, 7))
        gs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
        ax = fig.add_subplot(gs[0,0])
        ax.tick_params(axis='both', direction='in', which='major', color='k', width=1.5, labelsize=16, length=8, top=True, bottom=True,right=True,left=True)
        ax.tick_params(axis='both', direction='in', which='minor', color='k', width=1.2, labelsize=16, length=4, top=True, bottom=True,right=True,left=True)
        ax.grid(alpha=0.5, which='both')
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.set_ylim(0, 4)
        ax.axhline(1, color='k', linestyle='dashed')

        ax.errorbar(1e-3*rho_intersection, ratio, yerr=err, fmt='.', label='binned ratios')
        ax.plot(1e-3*rho_intersection, np.ones_like(ratio)*ratio_avg, label='weighted average: %.3f'%(ratio_avg))
        ax.text(0.04, 0.95, 'Ratio = comp/ref', transform=ax.transAxes, fontsize=12, horizontalalignment='left', verticalalignment='top')
        directory,refname = os.path.split(reference)
        directory,compname = os.path.split(comparison)
        ax.text(0.04, 0.905, 'ref = '+refname, transform=ax.transAxes, fontsize=12, horizontalalignment='left', verticalalignment='top')
        ax.text(0.04, 0.86, 'comp = '+compname, transform=ax.transAxes, fontsize=12, horizontalalignment='left', verticalalignment='top')

        ax.set_ylabel('Visibility amplitude ratios', fontsize=20)
        ax.set_xlabel(r'UV distance [kilo$\lambda$]', fontsize=20)
        plt.legend(loc='lower left')
        gs.update(wspace=0.0, hspace=0.1) # set the spacing between axes.
        plt.tight_layout()
        plt.savefig(comparison+'_estimate_flux_scale.png', dpi=300, transparent=True, bbox_inches='tight',pad_inches=0.015)
        plt.clf()
        print("Figure saved! To: ", comparison+'_estimate_flux_scale.png')

    return ratio_avg, np.sqrt(ratio_avg)






















# estimate_flux_scale(reference='./workflow/step2_noselfcal/ABAur_LB_EB1_initcont_shift.vis.npz',
#                     comparison='./workflow/step2_noselfcal/ABAur_LB_EB2_initcont_shift.vis.npz',
#                     fignametemplate='./workflow/step2_noselfcal/ABAur_LB_EB2_initcont_shift',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference='./workflow/step2_noselfcal/ABAur_LB_EB1_initcont_shift.vis.npz',
#                     comparison='./workflow/step2_noselfcal/ABAur_LB_EB3_initcont_shift.vis.npz',
#                     fignametemplate='./workflow/step2_noselfcal/ABAur_LB_EB3_initcont_shift',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference='./workflow/step2_noselfcal/ABAur_LB_EB1_initcont_shift.vis.npz',
#                     comparison='./workflow/step2_noselfcal/ABAur_LB_EB4_initcont_shift.vis.npz',
#                     fignametemplate='./workflow/step2_noselfcal/ABAur_LB_EB4_initcont_shift',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference='./workflow/step2_noselfcal/ABAur_LB_EB1_initcont_shift.vis.npz',
#                     comparison='./workflow/step2_noselfcal/ABAur_LB_EB5_initcont_shift.vis.npz',
#                     fignametemplate='./workflow/step2_noselfcal/ABAur_LB_EB5_initcont_shift',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference='./workflow/step2_noselfcal/ABAur_LB_EB1_initcont_shift.vis.npz',
#                     comparison='./workflow/step2_noselfcal/ABAur_LB_EB6_initcont_shift.vis.npz',
#                     fignametemplate='./workflow/step2_noselfcal/ABAur_LB_EB6_initcont_shift',
#                     incl=incl, PA=PA)

# estimate_flux_scale(reference='./workflow/step1/ABAur_LB_EB1_initcont.vis.npz',
#                     comparison='./workflow/step1/ABAur_LB_EB2_initcont.vis.npz',
#                     fignametemplate='./workflow/step1/ABAur_LB_EB2_initcont',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference='./workflow/step1/ABAur_LB_EB1_initcont.vis.npz',
#                     comparison='./workflow/step1/ABAur_LB_EB3_initcont.vis.npz',
#                     fignametemplate='./workflow/step1/ABAur_LB_EB3_initcont',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference='./workflow/step1/ABAur_LB_EB1_initcont.vis.npz',
#                     comparison='./workflow/step1/ABAur_LB_EB4_initcont.vis.npz',
#                     fignametemplate='./workflow/step1/ABAur_LB_EB4_initcont',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference='./workflow/step1/ABAur_LB_EB1_initcont.vis.npz',
#                     comparison='./workflow/step1/ABAur_LB_EB5_initcont.vis.npz',
#                     fignametemplate='./workflow/step1/ABAur_LB_EB5_initcont',
#                     incl=incl, PA=PA)
# estimate_flux_scale(reference='./workflow/step1/ABAur_LB_EB1_initcont.vis.npz',
#                     comparison='./workflow/step1/ABAur_LB_EB6_initcont.vis.npz',
#                     fignametemplate='./workflow/step1/ABAur_LB_EB6_initcont',
#                     incl=incl, PA=PA)

#The ratio of the fluxes of ./reduced_data/ABAur_SB1_initcont.vis.npz to ./reduced_data/ABAur_LB1_initcont.vis.npz is 0.93143
#The scaling factor for gencal is 0.965 for your comparison measurement
#The error on the weighted mean ratio is 8.644e-04, although it's likely that the weights in the measurement sets are off by some constant factor


#
#
# plot_gaincal_solutions(caltable='./reduced_data/ABAur_SB.ap',
#                        parentvis='./reduced_data/ABAur_SB_contp4.ms',
#                        quantity='phase',
#                        plot_average_soln=False,
#                        solint=None,
#                        spw='', # up to you to make sure chosen spw (if individual) is inside observation
#                        observation='0')








#
# image the ms with tclean
# (amplitude vs time plot will have drop outs)
#
# plot the old/new image (plot the model image?)
#
#
# estimate the SNR of the image











# sys.exit()
