import copy
import os
import glob

from astropy.io import fits
import copy
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.signal import medfilt

""" Contains functions for identifying and reparing CR hits in spatially
    scanned data. Optimized for WFC3/UVIS scans that are nearly vertical
    or horizontal. Adapted from an IDL routine by S. Casertano for STIS.

    Author
    ------
    Clare Shanahan, May 2018

    Notes
    -----
    Currently, these routines are optimized for data that are nearly vertically
    or horizontally scanned. Functionality for arbitrary scan angle will be
    added later. If your scans are at a significant angle on the detector, you
    can rotate the image. Also keep in mind that have also been optimized and
    tested on a specific WFC3/UVIS data set, so they may not be universally
    useful and output should be inspected carefully.

    The main routine is in `mask_and_repair_flagged_pixels`. This function
    takes in a 2D data array and returns the corrected data and the mask
    indicating CR hits. This function calls make_cr_mask as a first pass
    to identify CRs in data, and then does a second step to identify CRs in
    the tail of scans.

    The function `make_crcorr_file_scan_wfc3` takes an FLC or FLT as input,
    uses the 'scan_ang' keyword to determine if scan was vertical or horizonal,
    runs mask_and_repair_flagged_pixels on the data in the specified fits
    extension, and finally writes out a single-extension fits file with the
    corrected data.

"""


def make_cr_mask(data, nlen=4, ncomp=10, mult=4.):
    """Create a boolean mask of image where cosmic ray hits have been detected.

    Identifies suspected CRs in data as pixels that are more than `mult`
    sigma above the median of ncomp pixels +`nlen` ablve and -`nlen` below
    them. It is expected that valid pixels may be above the median of pixels at
    higher Y *or* lower Y because a scan is starting or ending, but not
    both. The filter may mark pixels incorrectly in the presence of
    anomalous low pixels. The few tests so far suggest it may be a little
    too agressive. Returns a mask (2D array) of identified cosmic rays in
    data.

    This function tends to miss CRs in scan tails. A second pass is done along
    with this function in `mask_and_repair_flagged_pixels` to fully correct
    spatial scans.

    Parameters
    ----------
    data : `np.array`
        2D numpy array of floats.
    nlen : int
        Integer number of pixels above and below a given pixel that will
        be used to compute the median.
    ncomp : int
        Number of pixels used to compute surrounding median.
    mult : int
        Sigma, above the median, used as a threshold for a CR detection.

    Returns
    -------
    mask : `np.array`
        2D numpy boolean mask for data marking locations of detected CRs.
        Pixels with value 1 signify clean data, and 0 indicate a CR.

    """

    data = copy.copy(data)
    nrange = nlen + ncomp

    upmedian = np.zeros(data.shape)  # median of NCOMP pixels above Y+NLEN
    downmedian = np.zeros(data.shape)  # median of NCOMP pixels below Y-NLEN
    dispimage = np.zeros(data.shape)  # rms of these 2*NCOMP pixels
    updisp = np.zeros(data.shape)
    downdisp = np.zeros(data.shape)

    ilow = 0  # column 0
    ihigh = data.shape[1] - 1  # last column
    jlow = nrange
    jhigh = data.shape[0] - nrange - 1

    for i in range(ilow, ihigh+1):  # rows
        for j in range(jlow, jhigh+1):  # columns
            upvec = data[j+nlen:j+nlen+ncomp-1+1, i]
            downvec = data[j-nlen-ncomp+1:j-nlen+1, i]
            upmedian[j, i] = np.median(upvec)
            downmedian[j, i] = np.median(downvec)
            dispimage[j, i] = np.std(np.concatenate((upvec, downvec)))
            updisp[j, i] = np.std(upvec)
            downdisp[j, i] = np.std(downvec)
            updown = np.concatenate((upvec, downvec))

    mask = (data > ((upmedian + mult * dispimage))) & \
           (data > ((downmedian + mult * dispimage)))

    return mask


def _fill_nan_interp(im_arr):
    """ Helper function to interpolate data with NaNs, used to replace CR hits
       with surrounding 'good' pixels."""
    trans_im_arr = im_arr.T

    for c in range(0, im_arr.shape[1]):
        row = trans_im_arr[c]  # 10 rows on the edge are garbage, leave as NaNs
        if row[0] == np.nan:
            row[0] = row[~np.isnan(row)][0]
        if row[-1] == np.nan:
            row[-1] = row[~np.isnan(row)][-1]
        trans_im_arr[c] = pd.Series(row).interpolate()
    return trans_im_arr.T


def _determine_scan_orientation_wfc3(hdr):
    """Uses 'SCAN_ANG' header keyword to determine if scan is vertical or
       horizontal."""

    scan_ang = hdr['SCAN_ANG']
    if np.abs(scan_ang - 138.5) < 5:
        return 'H'
    else:
        return 'V'


def mask_and_repair_flagged_pixels(data, scan_orient):
    """Given a 2D array of data, as well as the orientation of the scan
        (vertical or horizontal), this function will identify and correct
        cosmic rays and return both the corrected data and the mask
        identifying CRs.

        This process consists of two steps. First, make_cr_mask is run on the
        data to identify CR hits in the image. Next, a second pass is done on
        the data to identify CR hits in the tails of the scan, which
        make_cr_mask tends to miss due to the strong gradient in the tails of
        scans.

        Pixels identified as CRs are replaced with a linear interpolation of
        the 2 nearest surrounding 'good' pixels above or below, in the
        direction of the scan.

        Parameters
        ----------
        data : `np.array`
            2D array of data
        scan_orient : str
            'V' for vertical scans or 'H' for horizontal.
        Returns
        -------
        (corrected_data, mask) : Tuple of 'np.array'
            Tuple of 2D arrays containing the repaired data, and mask,
            respectivley.
    """

    data = copy.copy(data)
    if scan_orient == 'H':
        data = data.T

    # pass #1
    mask = 1.*make_cr_mask(data)

    # replace masked pix with the linear interp. of the nearest 2 unmasked pix
    k = np.where(mask > 0)

    if len(k) > 0:
        data[k] = np.nan

    data = _fill_nan_interp(data)

    # pass #2
    # process the ends of the trails differently than the rest

    gr = data - np.roll(data, 1, 0)

    gr_min = np.nanmin(gr.flatten())
    k_min = np.where(gr == gr_min)  # y,x where min occurs
    if len(k_min[0]) > 1:
        # select only one of these, sometimes it will find two nearby
        k_min = [np.array(item[0]) for item in k_min]

    gr_max = np.nanmax(gr.flatten())
    k_max = np.where(gr == gr_max)  # y,x where max occurs
    if len(k_max[0]) > 1:
        # select only one of these, sometimes it will find two nearby
        k_max = [np.array(item[0]) for item in k_max]

    mult = 5
    dy = 5
    for j in range(min([k_min[1], k_max[1]])[0] - 3,
                   max([k_min[1], k_max[1]])[0] + 3 + 1):
        v = data[:, j]
        mv = medfilt(v, 5)
        v_mv = v - mv
        sig = np.nanstd(v_mv)

        k = (np.where(v_mv > (mult * sig)))[0]

        if len(k) > 0:
            for i, val in enumerate(k):
                if ((np.abs(k[i] - k_min[0]) < dy) or
                   (np.abs(k[i] - k_max[0]) < dy)):
                    mask[k, j] = 1

    k = np.where(mask > 0)

    if len(k) > 0:
        data[k] = np.nan

    corrected_data = _fill_nan_interp(data)

    if scan_orient == 'H':
        data = corrected_data.T
        mask = mask.T

    return (data, mask)


def _write_fcr(input_file, output_dir, masked_im, ext, file_type):
    """Writes out CR corrected data (masked_im) as single extension fits file
        in output_dir. The 0th and 'ext' headers are concatenated and written
        to the output file. The name of this file will be the same as the input
        but with file_type 'flt' or 'flc' replaced with 'fcr'. """

    # concatenate 0th and ['SCI', ext] headers
    hdr_out = fits.open(input_file)[0].header + \
        fits.open(input_file)['SCI', ext].header

    hdu_new = fits.PrimaryHDU(masked_im, header=hdr_out)
    output_path = output_dir+os.path.basename(input_file).\
        replace('{}.fits'.format(file_type), 'fcr.fits')

    if os.path.isfile(output_path):
        os.remove(output_path)
    print('Writing', output_path)
    hdu_new.writeto(output_path)


def _write_mask(input_file, output_dir, mask, ext, file_type):
    """Writes out CR mask as single_extension fits file in output_dir. See
        _write_fcr()."""

    # concatenate 0th and ['SCI', ext] headers
    hdr_out = fits.open(input_file)[0].header + \
        fits.open(input_file)['SCI', ext].header

    hdu_new = fits.PrimaryHDU(mask, header=hdr_out)

    output_path = output_dir + os.path.basename(input_file).\
        replace('{}.fits'.format(file_type), 'mask.fits')

    if os.path.isfile(output_path):
        os.remove(output_path)
    print('Writing', output_path)
    hdu_new.writeto(output_path)


def make_crcorr_file_scan_wfc3(input_file, output_dir=None, ext=1,
                               write_mask=True):
    """ Wrapper function that calls routine to identify and correct
     cosmic rays in spatially scanned HST flt.fits or flc.fits images and
     write out a corrected single-extension 'fcr.fits' file. Only the data in
     the input file at the specified extension 'ext' is corrected and written
     out, keep this in mind when working with multi-extension fits files.

    This function calls the main CR correction routine
    `mask_and_repair_flagged_pixels`, and writes the output of this to file.
    If 'write_mask' is set to True, the mask indentifying locations of CR hits
    in the data will be written out as well, to a 'mask.fits' file.

    Currently, this function is optimized for files that are nearly vertically
    or horizontally scanned. Functionality for arbitrary scan angle will be
    added later. If your scans are at a significant angle on the detector, you
    can rotate the image, padding it with a fill value, save, and pass that as
    input to this function.

    Parameters
    ----------
    input_file : str
        Full path to input fits file.
    output_dir : str
        Directory where corrected files should be output. Defaults to
        input_file directory.
    ext : int
        FITS extension of data.
    write_mask : bool
        If True, in addition to the CR corrected image the mask indicating the
        location of CR hits will be saved as well.

    Outputs
    -------
        A single extension CR-corrected fits file. The file name is the same
        as the input, but with 'flt' or 'flc' changed to 'fcr'. By default,
        corrected files are output in the same directory as the input unless a
        different 'output_dir' is specified.) If write_mask, a 'mask.fits' file
        will be written out as well.
    """

    print('Running CR rejection on {}'.format(input_file))
    # define output directory, same directory as input
    if output_dir is None:
        output_dir = input_file.replace(os.path.basename(input_file), '')
    output_dir = os.path.join(output_dir, '')  # ensure trailing slash.
    # check that input file is either an flt or flc
    if 'flt' in os.path.basename(input_file):
        file_type = 'flt'
    elif 'flc' in os.path.basename(input_file):
        file_type = 'flc'
    else:
        raise ValueError('Input file must be flt.fits or flc.fits')

    # open file and get data, 0th header
    hdu = fits.open(input_file)
    hdr0 = hdu[0].header
    data = hdu[ext].data
    scan_orient = _determine_scan_orientation_wfc3(hdr0)

    # call CR rejection routine
    masked_im, mask = mask_and_repair_flagged_pixels(data, scan_orient)

    _write_fcr(input_file, output_dir, masked_im, ext, file_type)

    if write_mask:
        _write_mask(input_file, output_dir, mask, ext, file_type)
