"""
Version0 by Ryan Loomis (?/?/22)
Version1 with updated by Rich Teague / Ryan Loomis (8/8/22)

All functions required to align several EBs.
"""
import casatasks
import casatools
import numpy as np
from numba import njit
from astropy import units as u
from scipy.optimize import minimize
from astropy.coordinates import SkyCoord, FK5
import os
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings

import faulthandler
faulthandler.enable()

# standard constants

arcsec = np.pi / 180.0 / 3600.0  # radians
cc = 2.99792458e10  # cm/s

skycoord_frames = {'ICRS':'icrs','J2000':FK5(equinox='J2000')}


def calculate_phase_shift(grid_vis, grid_nvis, grid_uu, grid_vv, mu_RA,
                          mu_DEC):
    """
    Apply a phase shifted to a single channel of gridded visibilities.

    Args:
        grid_vis (array): 2D array of gridded visibilities.
        grid_nvis (array): 2D array of number of gridded visibilities in cell.
        grid_uu (array): 2D array of grid u coordinates.
        grid_vv (array): 2D array of grid v coordinates.
        mu_RA (float): Offset in right ascension from phase center in [arcsec].
        mu_DEC (float): Offset in declination from phase center in [arcsec].

    Returns:
        shifted_grid_vis (array): 2D array of phase-shifted gridded
            visibilities.
    """
    phase_shifts = grid_uu * mu_RA * arcsec + grid_vv * mu_DEC * arcsec
    phase_shifts = np.exp(-2.0 * np.pi * 1.0j * phase_shifts)
    shifted_grid_vis = grid_vis * phase_shifts
    shifted_grid_vis[grid_nvis < 1] = 0.0 + 0.0j
    return shifted_grid_vis


def calculate_phase_difference(grid_vis1, grid_vis2, grid_wgts1, grid_wgts2,
                               grid_nvis):
    """
    Calculate the phase difference (no amplitude) between two sets of gridded
    visibilities, ``grid_vis1`` and ``grid_vis2``.

    Args:
        grid_vis1 (array): 2D array of gridded visibilities.
        grid_vis2 (array): 2D array of gridded visibilities.
        grid_wgts1 (array): 2D array of gridded weights for ``grid_vis1``.
        grid_wgts2 (array): 2D array of gridded weights for ``grid_vis2``.
        grid_nvis (array): 2D array of number of gridded visibilites in cell.

    Returns:
        phase_difference (array): 2D array of the phase differences between
            ``grid_vis1`` and ``grid_vis2``.
    """
    angle1, angle2 = np.angle(grid_vis1), np.angle(grid_vis2)
    phase_difference = np.minimum(2.0 * np.pi - np.abs(angle2 - angle1),
                                  np.abs(angle2 - angle1))
    phase_difference = 2.0 * np.sin(phase_difference / 2.0) * np.abs(grid_vis1)
    phase_difference *= np.mean([grid_wgts1, grid_wgts2], axis=0)
    return phase_difference


def calculate_full_phase_difference(grid_vis1, grid_vis2, grid_wgts1,
                                    grid_wgts2, grid_nvis):
    """
    Calculate the full phase difference (including amplitude) between two sets
    of gridded visibilities, ``grid_vis1`` and ``grid_vis2``.

    Args:
        grid_vis1 (array): 2D array of gridded visibilities.
        grid_vis2 (array): 2D array of gridded visibilities.
        grid_wgts1 (array): 2D array of gridded weights for ``grid_vis1``.
        grid_wgts2 (array): 2D array of gridded weights for ``grid_vis2``.
        grid_nvis (array): 2D array of number of gridded visibilites in cell.

    Returns:
        phase_difference (array): 2D array of the phase differences between
            ``grid_vis1`` and ``grid_vis2``.
    """
    phase_difference = np.abs(grid_vis2 - grid_vis1)
    phase_difference *= np.mean([grid_wgts1, grid_wgts2], axis=0)
    return phase_difference


@njit(fastmath=True) # This causes segmentation fault on CANFAR; without njit, the loop takes ~20 seconds
def grid(grid_vis, grid_nvis, grid_wgts, uu, vv, du, dv, npix, vis, wgts):
    """
    TBD

    Args:
        grid_vis (array): 2D array of gridded visibilities.
        grid_nvis (array): 2D array of number of gridded visibilities in cell.
        grid_wgts (array): 2D array of gridded weights.
        uu (array):
        vv (array):
        du (array):
        dv (array)
        npix (array):
        vis (array):
        wgts (array):

    Returns:
        grid_vis (array): 2D array of gridded visibilities.
        grid_nvis (array): 2D array of number of gridded visibilities in cell.
        grid_wgts (array): 2D array of gridded weights.
    """
    print('grid ---> Starting loop (uses up to ~1.3 GB of memory for our dataset)')
    for i in np.arange(uu.size):
        uidx_a = int(npix / 2.0 + uu[i] / du + 0.5)
        uidx_b = int(npix / 2.0 - uu[i] / du + 0.5)
        vidx_a = int(npix / 2.0 + vv[i] / dv + 0.5)
        vidx_b = int(npix / 2.0 - vv[i] / dv + 0.5)
        grid_vis[uidx_a, vidx_a] += vis[i]
        grid_vis[uidx_b, vidx_b] += np.conjugate(vis[i])
        grid_wgts[uidx_a, vidx_a] += wgts[i]
        grid_wgts[uidx_b, vidx_b] += wgts[i]
        grid_nvis[uidx_a, vidx_a] += 1
        grid_nvis[uidx_b, vidx_b] += 1
    print('grid ---> Finished loop!')
    return grid_vis, grid_nvis, grid_wgts


def ingest_ms(base_ms, npix, cell_size, grid_needs_to_cover_all_data, spwid=0):
    """
    Ingest a measurement set and grid onto a regular grid with ``npix`` cells,
    each with a size ``cell_size`` in [arcsec].

    Args:
        base_ms (str): Measurement set to ingest.
        npix (int): Number of pixels for the grid.
        cell_size (float): Cell size in [arcsec] for the grid.
        grid_needs_to_cover_all_data (bool): if True, make sure that grid cover all data
        spwid (int): The spectral window to ingest; defaults to 0.

    Returns:
        grid_vis (array): 2D array of gridded visibilities.
        grid_nvis (array): 2D array of number of gridded visibilities in cell.
        grid_uu (array): 2D array of grid u coordinates.
        grid_vv (array): 2D array of grid v coordinates.
        grid_wgts (array): 2D array of gridded weights.
    """

    # Use CASA table tools to get required columns.

    # this is an assumption that is valid for exoALMA data, but not in general
    data_desc_id = str(spwid)

    tb = casatools.table()
    tb.open(base_ms+"/SPECTRAL_WINDOW")
    chan_freqs_all = tb.getvarcol("CHAN_FREQ")
    tb.close()
    chan_freqs = chan_freqs_all["r"+str(spwid+1)]
    tb.open(base_ms)
    subt = tb.query("DATA_DESC_ID=="+data_desc_id)
    flag = subt.getcol("FLAG")
    uvw = subt.getcol("UVW")
    weight = subt.getcol("WEIGHT")
    data = subt.getcol("DATA")
    ant1 = subt.getcol("ANTENNA1")
    ant2 = subt.getcol("ANTENNA2")
    subt.close()
    tb.close()
    print('ingest_ms ---> Successfully grabbed what we need from the base_ms')

    # Define visibilities and weights.

    vis = (data[0, :] + data[1, :]) / 2.0
    wgts = weight[0, :] + weight[1, :]

    # Break out the u, v spatial frequencies, convert from m to lambda.

    uu = uvw[0, :][:, np.newaxis] * chan_freqs[:, 0] / (cc / 100.0)
    vv = uvw[1, :][:, np.newaxis] * chan_freqs[:, 0] / (cc / 100.0)

    # Toss out the autocorrelation placeholders.

    xc = np.where(ant1 != ant2)[0]
    uu, vv, vis, wgts = uu[xc], vv[xc], vis[:, xc], wgts[xc]

    # Remove flagged visibilities.

    flag = np.logical_not(np.prod(flag, axis=(0, 2)).T)
    uu, vv, vis = uu[:, flag], vv[:, flag], vis[flag]

    # Reshape the visibilities, weights and (u,v) points to a single list.
    print('ingest_ms ---> About to reshape the visibilities, weights and (u,v) points to a single list.')
    vis = np.ravel(np.broadcast_to(vis, (uu.shape[1], uu.shape[0])).T)
    wgts = np.ravel(np.broadcast_to(wgts, (uu.shape[1], uu.shape[0])).T)
    uu = np.ravel(uu)
    vv = np.ravel(vv)
    print('ingest_ms ---> Successfully reshaped the visibilities, weights and (u,v) points to a single list.')

    # Define grid in uv space.

    dl = cell_size * arcsec
    dm = cell_size * arcsec
    du = 1.0 / npix / dl
    dv = 1.0 / npix / dm

    # Empty arrays to hold gridded data.

    grid_vis = np.zeros((npix, npix)).astype('complex')
    grid_wgts = np.zeros((npix, npix))
    grid_nvis = np.zeros((npix, npix))
    print('About to do np.mgrid; apparently this could break the code')
    grid_uu, grid_vv = np.mgrid[(-npix/2.0+0.5)*du:(npix/2.0+0.5)*du:du,
                                (-npix/2.0+0.5)*dv:(npix/2.0+0.5)*dv:dv]
    #sometimes mgrid gives the wrong shape, i.e. one element too much
    #for example, cell_size=0.1 and npix=100 leads to grid_uu.shape=(101,101),
    #which then crashes the code
    print('Code not broken!')
    if not grid_uu.shape == (npix,npix) or not grid_vv.shape == (npix,npix):
        raise RuntimeError('please choose a slightly different npix')

    #toss away data that falls outside of the grid:
    min_uu,max_uu = np.min(grid_uu),np.max(grid_uu)
    min_vv,max_vv = np.min(grid_vv),np.max(grid_vv)
    inside_grid = (min_uu<uu) & (uu<max_uu) & (min_vv<vv) & (vv<max_vv)
    if not np.all(inside_grid):
        warnings.warn(f'some data of {base_ms} are outside your uv grid')
        if grid_needs_to_cover_all_data:
            raise ValueError('grid does not cover all data')
    vis = vis[inside_grid]
    wgts = wgts[inside_grid]
    uu = uu[inside_grid]
    vv = vv[inside_grid]

    # Grid the data and return.
    print('ingest_ms ---> About to call grid')
    grid_vis, grid_nvis, grid_wgts = grid(grid_vis,
                                          grid_nvis,
                                          grid_wgts,
                                          uu,
                                          vv,
                                          du,
                                          dv,
                                          npix,
                                          vis,
                                          wgts)
    print('ingest_ms ---> Woo leaving ingest_ms!')
    return grid_vis, grid_nvis, grid_uu, grid_vv, grid_wgts


def calculate_likelihood(x, data):
    """
    Calculate the likelihood using a full phase difference after phase shifting
    the second measurement set in ``data``.

    RICH -- Is this actually the likelihood, or are we just minimizing the
            aggregate phase difference between the two datasets?
         -- (Jess: Ryan told me it's not really)

    Args:
        x (list): A list containing the phase shift, ``(mu_RA, mu_DEC)``.
        data (list): A list containing a list of ``grid_vis``, ``grid_nvis``,
            ``grid_uu``, ``grid_vv`` and ``grid_wgts`` for the two datasets.

    Returns:
        likelihood (float): Likelihood value.
    """

    # Unpack the data.

    ms1_data, ms2_data = data

    # Apply a phase center shift to the second measurement set.

    shifted_data = calculate_phase_shift(grid_vis=ms2_data[0],
                                         grid_nvis=ms2_data[1],
                                         grid_uu=ms2_data[2],
                                         grid_vv=ms2_data[3],
                                         mu_RA=x[0],
                                         mu_DEC=x[1])

    # Calculate the likelihood and return.

    likelihood = calculate_full_phase_difference(grid_vis1=ms1_data[0],
                                                 grid_vis2=shifted_data,
                                                 grid_wgts1=ms1_data[4],
                                                 grid_wgts2=ms2_data[4],
                                                 grid_nvis=ms1_data[1])
    likelihood = np.sum(np.abs(likelihood))
    return likelihood


def find_offset(reference_ms, offset_ms, npix=1024, cell_size=0.01, spwid=0,
                fail_silently=False,verbose=False,plot_uv_grid=True, # Jess changed plot_uv_grid=True, previously False
                plot_filename=None):
    """
    Find the offset between ``offset_ms`` and ``reference_ms`` by minimizing
    the aggregate phase angle and amplitude.

    Args:
        reference_ms (str): The reference measurement set.
        offset_ms (str): The measurement set to use to derive an offset.
        npix (optional[int]): Number of pixels in the grid.
        cell_size (optional[float]): The grid cell size in [arcsec].
        spwid (optional[int]): The spectral window to evaluate; defaults to 0.
        fail_silently (optional[bool]): If ``True`` return a null offset if the
            minimization fails, otherwise raise a ``RuntimeError``.
        verbose (bool): whether to print out info
        plot_uv_grid (bool): whether to plot an overview of the uv grid
        plot_filename (str): filename out output uv grid plot
    Returns:
        offset (list): A list specifying the right ascension and declination
            offsets in [arcsec].
    """

    #to calculate the offset, we need ms in the same reference frame
    #thus, we convert to J2000 if necessary
    input_ms = {'ref':reference_ms,'offset':offset_ms}
    temporary_ms = []
    ms_for_offset_calculation = {}
    for ms_ID,ms in input_ms.items():
        print('find_offset ---> ms_ID, ms: ', ms_ID,ms)
        phase_center = get_phase_center(measurement_set=ms)
        print('find_offset ---> phase_center: ', phase_center)
        frame = get_coord_frame(phase_center)
        print('find_offset ---> frame: ', frame)
        if frame == 'J2000':
            if verbose:
                print(f'{ms_ID} ms {ms} is already in J2000, no need to transform')
            ms_for_offset_calculation[ms_ID] = ms
        elif frame == 'ICRS':
            if verbose:
                print(f'{ms_ID} ms {ms} is in ICRS, going to use copy in J2000 to calculate offset')
            assert ms[-3:] == '.ms'
            J2000_ms = ms[:-3] + '_J2000.ms'
            if os.path.isdir(J2000_ms):
                shutil.rmtree(J2000_ms)
            sky_coord_phase_center = get_skycoord(phase_center)
            sky_coord_phase_center = sky_coord_phase_center.transform_to(
                                                  skycoord_frames['J2000'])
            phase_center_J2000 = f'J2000 {sky_coord_phase_center.to_string("hmsdms")}'
            # casatasks.fixvis(vis=ms,outputvis=J2000_ms,phasecenter=phase_center_J2000) # Jess: since we have to use 6.4.4.31, we have to use phaseshift
            print('find_offset ---> About to do phaseshift task')
            os.system('rm -rf %s*' % J2000_ms) # Jess added this
            casatasks.phaseshift(vis=ms,outputvis=J2000_ms,phasecenter=phase_center_J2000)
            print('find_offset ---> Completed phaseshift task')
            temporary_ms.append(J2000_ms)
            ms_for_offset_calculation[ms_ID] = J2000_ms
        else:
            raise RuntimeError(f'unknown frame {frame}')


    # Ingest the required measurement sets. Will return a list of ``grid_vis``,
    # ``grid_nvis``, ``grid_uu``, ``grid_vv`` and ``grid_wgts``.
    print('About to do ingest_ms')
    ms1 = ingest_ms(base_ms=ms_for_offset_calculation['ref'], npix=npix,
                    cell_size=cell_size,grid_needs_to_cover_all_data=False,spwid=spwid)
    ms2 = ingest_ms(base_ms=ms_for_offset_calculation['offset'], npix=npix,
                    cell_size=cell_size, grid_needs_to_cover_all_data=True,spwid=spwid)

    # Define the overlap between the two measurement sets.
    overlap = np.logical_and(ms1[1].real >= 1, ms2[1].real >= 1).astype('int')

    # if plot_uv_grid:
        # vmin = 0
        # vmax = np.log10(np.max((ms1[1],ms2[1])))
        # images = []
        # for ID,ms_filename,ms,ax in zip(('ref','offset'),(reference_ms,offset_ms),
        #                                 (ms1,ms2),axes[0,:]):
        #     print('Plotting ID: ', ID)
        #     ax.set_title(f'{ms_filename} ({ID})',fontsize=8)
        #     #to avoid divide by zero warnings:
        #     log_vis_ngrid = np.log10(np.where(ms[1]>0,ms[1],1e-10))
        #     print('About to do pcolormesh...')
        #     # img = ax.pcolormesh(ms[2],ms[3],log_vis_ngrid,vmin=vmin,vmax=vmax, rasterized=True)
        #     # images.append(img)
        #     n_nonempty_grid_points = np.sum(ms[1]>0)
        #     print(f'{ID}: there are {n_nonempty_grid_points} uv grid points containing data')
        # # fig.colorbar(images[1],ax=axes,label='log(number of vis points)',
        # #              location='top',shrink=0.6)
        # axes[1,0].set_title('overlap')
        # axes[1,0].pcolormesh(ms1[2],ms1[3],overlap,cmap='Greys', rasterized=True)
        # for ax in axes.ravel():
        #     ax.set_xlabel('grid u [lambda]')
        #     ax.set_ylabel('grid v [lambda]')
        #     ax.set_aspect('equal')
        # axes[1,1].remove()
        # print('Ready to save figure to ', plot_filename)
        # if plot_filename is not None:
        #     print('plot_filename is not None')
        #     # plt.show()
        #     plt.savefig(plot_filename, dpi=50) # the memory usage gets up to 1.5 GB here
        # print('Saved figure! Or, did not, depending on whether the following is None or not: ', plot_filename)
        # print('Actually we skipped saving the figure, to be able to continue, and avoid another segmentation fault...')

    # Mask out all the cells where there is no overlap. Note that the np.clip()
    # is to avoid RuntimeWarnings when dividing by zero. These grid points will
    # be masked out by overlap anyway.
    print('find_offset ---> masking out the cells where there is no overlap')
    ms1 = [ms1[0] / np.clip(ms1[1], a_min=1.0, a_max=None) * overlap,
           ms1[1] * overlap, ms1[2], ms1[3], ms1[4]]
    ms2 = [ms2[0] / np.clip(ms2[1], a_min=1.0, a_max=None) * overlap,
           ms2[1] * overlap, ms2[2], ms2[3], ms2[4]]

    # Derive the offset. The starting point x0 is picked as a fraction (chosen as 1/6)
    #of the "resolution" expected from the longest baseline of the offset data set
    contains_data = ms2[1] > 0
    max_uu_vv = np.max(np.sqrt(ms2[2][contains_data]**2 + ms2[3][contains_data]**2))
    x0 = np.array((1/max_uu_vv/arcsec,1/max_uu_vv/arcsec)) / 6
    print('x0: ',x0)

    print('About to calculate the visibility residuals between the two EBs...')
    res = minimize(fun=calculate_likelihood,x0=x0,args=[ms1, ms2],method='L-BFGS-B')

    # Return the offset, otherwise if this fails, either raise an error or
    # returns a null offset.

    if not res.success:
        print('find_offset ---> res not a success')
        if fail_silently:
            print('find_offset ---> failing silently')
            return [0.0, 0.0]
        else:
            print(res)
            raise RuntimeError
    for temp_ms in temporary_ms:
        print('find_offset ---> temp_ms: ', temp_ms)
        if verbose:
            print(f'going to delete temporary ms {temp_ms}')
        print('find_offset ---> removing temp_ms')
        shutil.rmtree(temp_ms)
    print('Leaving find_offset!')
    return res.x


def assert_only_one_field_in_ms(ms):
    msmd = casatools.msmetadata()
    msmd.open(ms)
    field_names = msmd.fieldnames()
    msmd.close()
    if len(field_names) != 1:
        raise RuntimeError(f"Multiple fields found: {field_names}."
                           +" This script assumes that there is only one field")


def get_phase_center(measurement_set):
    """
    Read the phase center for the given measurement set.

    Args:
        measurement_set (str): Measurement to grab the phase center from.

    Returns:
        phase_center (str): Coordinates of the phase center.
    """

    assert_only_one_field_in_ms(ms=measurement_set)

    # Grab the measurement set metadata.

    msmd = casatools.msmetadata()
    msmd.open(measurement_set)

    # Check the returned epoch is known (either ICRS or J2000). If this can't
    # be determined, assumed that this is ICRS.

    ref_coord = [msmd.refdir()['m0']['value'], msmd.refdir()['m1']['value']]
    ref_epoch = msmd.refdir()['refer']
    msmd.close()
    if ref_epoch == "ICRS":
        frame = "icrs"
    elif ref_epoch == "J2000":
        frame = FK5(equinox="J2000")
    else:
        print("WARNING: Could not determine reference frame. Assuming ICRS")
        frame = "icrs"

    # Convert the reference coordinate into a parsable string.

    c = SkyCoord(ref_coord[0], ref_coord[1], frame=frame, unit=u.rad)
    phase_center = ref_epoch + " " + c.to_string('hmsdms')
    return phase_center


def get_coord_frame(coord):
    if coord[:4] == 'ICRS':
        frame = 'ICRS'
    elif coord[:5] == 'J2000':
        frame = 'J2000'
    else:
        raise RuntimeError(f'unable to determine reference frame of {coord}')
    return frame


def get_skycoord(coord):
    frame = get_coord_frame(coord)
    skycoord_input = coord.replace(f'{frame} ', '')
    splitted_skycoord_input = skycoord_input.split(" ")
    c = SkyCoord(splitted_skycoord_input[0],splitted_skycoord_input[1],
                 frame=skycoord_frames[frame])
    return c


def generate_shifted_J2000_coords(original_coord, offset):
    """
    For a given reference coordinate, ``original_coord``, and a list of RA and
    Dec offsets, ``offset``, calculate the new, shifted coordinate.

    NOTE: Due to quirks of CASA, some functions only work with J2000
          coordinates and so the input must first be convered from ICRS to
          J2000 coordinates, and then the offset applied.

    Args:
        original_coord (str): Initial coordinate in ICRS or J2000 'frame hmsdms' format, for
            example: ``"ICRS 12h43m12.159252s +85d52m12.952837s"``.
        offset (tuple): Tuple of RA and Dec offset in [arcsec], for example:
            ``[0.01, -0.02]``.

    Returns:
        shifted_J2000 (str): Shifted coordinate in J2000 'hmsdms' format.
    """

    # Read in the position and convert to J2000 coordinate if necessary.

    original_frame = get_coord_frame(original_coord)

    c = get_skycoord(original_coord)

    # Apply the RA and Dec offsets.
    ra_offset, dec_offset = offset
    c.data.lon[()] = c.ra + ra_offset / 3600.0 / np.cos(c.dec.rad) * u.degree
    c.data.lat[()] = c.dec + dec_offset / 3600.0 * u.degree
    c.cache.clear()

    # Extract the string and return.
    if original_frame == 'ICRS':
        shifted_J2000 = c.transform_to(skycoord_frames['J2000']).to_string('hmsdms')
    elif original_frame == 'J2000':
        shifted_J2000 = c.to_string('hmsdms')
    shifted_J2000 = 'J2000 ' + shifted_J2000
    return shifted_J2000


def update_phase_center(vis, new_phase_center, ref_phase_center,
                        suffix='_shift', move_to_folder=None):
    """
    Apply the updated phase center to the provided measurement set. This will
    update both the phase center using 'fixvis' (Jess: now 'phaseshift') with ``new_phase_center`` and
    the phase center coordinates using 'fixplanets' with ``ref_phase_center``.

    Args:
        vis (str): Visibility set to update the phase center of.
        new_phase_center (str): New phase center to apply in J2000.
        ref_phase_center (str): Reference phase center to update in J2000.
        suffix (optional[str]): Suffix to add prior to '.ms' for the new MS.
    """
    assert_only_one_field_in_ms(ms=vis)
    for pc in (new_phase_center,ref_phase_center):
        assert get_coord_frame(pc) == 'J2000'
    shifted_vis = vis.replace('.ms', suffix+'.ms')
    # casatasks.fixvis(vis=vis, outputvis=shifted_vis, field='0',
    #                  phasecenter=new_phase_center)        # Jess: since we have to use 6.4.4.31, we have to use phaseshift
    os.system('rm -rf %s*' % shifted_vis) # Jess added this
    casatasks.phaseshift(vis=vis, outputvis=shifted_vis, field='0',
                     phasecenter=new_phase_center)
    casatasks.fixplanets(vis=shifted_vis, field='0', fixuvw=False,
                         direction=ref_phase_center)
    if move_to_folder is not None: # Jess added
        os.system('mv %s %s' %(shifted_vis, move_to_folder)) # Jess added


def align_measurement_sets(reference_ms, align_ms, align_offsets=None,npix=1024,
                           cell_size=0.01,spwid=0,plot_uv_grid=False,
                           plot_file_template=None, ouput_folder=None):
    """
    Using ``reference_ms`` as the truth, align all meausrement sets in
    ``align_ms``. This includes calculating the RA and Dec offset between the
    two measurement sets, calcuating the updated phase center coordinate, and
    then appling this phase center shift to the data.

    Args:
        reference_ms (str): The MS to use as a the fixed reference point.
        align_ms (str or list): The MS to align to the reference MS.
        align_offsets (optional[list]): list of offsets to be used for the alignment.
            Each element corresponds to an element of align_ms.
            If None, offsets will be calculated.
        npix (optional[int]): Number of pixels in the grid.
        cell_size (optional[float]): Cell size in [arcsec] for the grid.
        spwid (optional[int]): The spectral window to align based on; defaults to 0.
        plot_uv_grid (bool): whether to plot an overview of the uv grid
        plot_file_template (str): template to produce output file of uv grid plot
    """

    # Use the reference MS as the phase center for all the shifted MSs. Note
    # the call to generate_shifted_J2000_coords() is to convert from ICRS to J2000 for
    # the fixplanets() call later which only uses J2000.

    source_phase_center = get_phase_center(measurement_set=reference_ms)
    print('---> Phase center of reference EB: ', source_phase_center)
    source_phase_center = generate_shifted_J2000_coords(
                              original_coord=source_phase_center,offset=[0.0,0.0])
    print('---> Phase center of reference EB: ', source_phase_center)

    # Cycle through each measurement set and find the offset and then update
    # the phase center and replace the coordinates to match that of the
    # reference MS.
    align_ms = np.atleast_1d(align_ms)
    if align_offsets is not None:
        assert len(align_offsets) == len(align_ms),\
                'number of provided offsets does not correspond to number of ms'
    zero_offset = np.zeros(2)
    for i,ms in enumerate(align_ms):
        print('---> Cycling; i, ms: ', i,ms)
        is_ref_ms = (ms == reference_ms)
        if is_ref_ms:
            if align_offsets is not None:
                assert np.all(align_offsets[i] == zero_offset),\
                                   'offset of ref ms has to be [0,0]'
            offset = zero_offset
        else:
            if align_offsets is not None:
                offset = align_offsets[i]
            else:
                if plot_file_template is None:
                    plot_filename = None
                else:
                    directory,file_template = os.path.split(plot_file_template)
                    plot_filename = os.path.join(directory,f'{ms}_{file_template}')
                    print('---> directory, file_template, plot_filename: ', directory,file_template,plot_filename)
                offset = find_offset(reference_ms=reference_ms,offset_ms=ms,npix=npix,
                                     cell_size=cell_size,spwid=spwid,
                                     plot_uv_grid=plot_uv_grid,
                                     plot_filename=plot_filename)
        ms_phase_center = get_phase_center(measurement_set=ms)
        print('---> ms_phase_center: ', ms_phase_center)
        shifted = generate_shifted_J2000_coords(original_coord=ms_phase_center,
                                                offset=offset)
        print('---> shifted: ', shifted)
        if align_offsets is None:
            print('#New coordinates for {}'.format(ms))
            if is_ref_ms:
                print('#no shift, reference MS.\n')
            else:
                print('#requires a shift of [{:.5g},{:.5g}]\n'.format(*offset))
        else:
            print(f'applying shift {offset} to {ms}')
        update_phase_center(vis=ms,new_phase_center=shifted,
                            ref_phase_center=source_phase_center, move_to_folder=ouput_folder)

# Jess commented this out:
if __name__ == '__main__':

    reference_ms = 'data/SY_Cha_LB_EB3_spw1.ms'
    align_ms = ['data/SY_Cha_{}_EB{:d}_spw1.ms'.format(b, i)
                for b, i in zip(['LB', 'LB', 'LB', 'SB'], [0, 1, 2, 0])]
    align_measurement_sets(reference_ms, align_ms, spwid=1)
