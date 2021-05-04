"""

This module contains functions to perform source detection and aperture
photometry on spatial scan data.

Author
------
Clare Shanahan, Dec 2019

"""

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clipped_stats
import copy
import matplotlib.pyplot as plt
import numpy as np
from photutils import (detect_sources, detect_threshold, source_properties,
                       RectangularAperture, aperture_photometry)


def detect_sources_scan(data, snr_threshold=3.0, sigma_kernel=3,
                        size_kernel=(3, 3), n_pixels=500, show=False):

    """Uses image segmentation to detect sources in spatially scanned images.

        A pixel-wise threshold image used to detect sources is generated based
        on the data and the snr_threshold provided. Data is then convolved with
        a 2D Gaussian kernel, of width sigma_kernel (default 3.0) and x, y size
        given by size_kernel (default 3 pixels x 3 pixels) to smooth out some
        of the background noise.

        A segmentation map of the convolved image is generated using the
        threshold image and npixels, the lower limit on the number of connected
        pixels that represent a true source (default is 1000., since scans
        cover a larger area of the detector).

        Optionally, a plot showing the detected source(s) can be shown.

        Parameters
        ----------
        data : `~np.array`
            2D array of data (floats)
        snr_threshold : int or float
            For creation of the threshold image, the signal-to-noise ratio per
            pixel above the background used for which to consider a pixel as
            possibly being part of a source. The background is calculated for
            the entire image using sigma-clipped statistics.
        sigma_kernel : float or int
            Width of 2D gaussian kernel, in pixels.
        size_kernel : tuple
            x, y size in pixels of kernel.
        n_pixels : int
            The (positive) integer number of connected pixels, each greater
            than the threshold, that an object must have to be detected.
        show : bool
            If True, the image will be displayed with all of the identified
            sources marked.
        Returns
        -------
        properties_tbl : '~astropy.table.QTable'
            Table containing properties of detected sources()

        """
    # make threshold image
    threshold = detect_threshold(data, nsigma=snr_threshold)

    # construct gaussian kernel to smooth image
    sigma = sigma_kernel * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=size_kernel[0],
                              y_size=size_kernel[1])
    kernel.normalize()

    # pass in data, convolution kernel to make segmentation map

    segm = detect_sources(data, threshold, npixels=n_pixels,
                          filter_kernel=kernel)

    props = source_properties(data, segm)
    properties_tbl = props.to_table()

    if show:
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 7))
        ax[0].imshow(segm.data, origin='lower')
        ax[1].scatter(properties_tbl['xcentroid'],
                      properties_tbl['ycentroid'], marker='x', c='r')
        z1, z2 = (-12.10630989074707, 32.53888328838081)
        ax[1].imshow(data, origin='lower', cmap='Greys_r', vmin=z1, vmax=z2)
        title_str = '{} detected sources'.format(len(properties_tbl))
        if len(properties_tbl) == 1:
            title_str = '{} detected source'.format(len(properties_tbl))
        ax[0].set_title(title_str)
        plt.tight_layout()
        plt.show()
        plt.close()

    return properties_tbl


def calc_sky(data, x_pos, y_pos, source_mask_len, source_mask_width, n_pix,
             method='median'):
    """
    Calculates sky level in a rectangular annulus around source. The source is
    first masked with a rectangle of dimensions source_mask_width x
    source_mask_length. Then, the background is computed in a rectangular
    rind around the source mask of span n_pix. The background method defaults
    to a sigmia-clipped median, but the mean can also be returned.

    Parameters
    ----------
    data : `np.array`
        2D array of floats
    x_positions : float
        X position of source.
    y_positions : float
        Y position of source.
    source_mask_len : int
        Length of rectangle (along y axis) used to mask source.
    source_mask_width : int
        Width of rectangle (along x axis) used to mask source.
    n_pix : int
        Number of pixels around source masking rectangle that define the rind
        regionused to measure the background.
    method : str
        'Median' or 'Mean', sigma clipped.
     """
    temp_data = copy.deepcopy(data)

    # mask rect. aperture around source, to exclude these pix in sky calc.
    temp_data[int(y_pos-(source_mask_len/2.)):int(y_pos+(source_mask_len/2.)),
              int(x_pos-(source_mask_width/2.)):
              int(x_pos+(source_mask_width/2.))] = np.nan

    flat_dat = temp_data.flatten()

    flat_masked_dat = flat_dat[~np.isnan(flat_dat)]
    mean, median, backrms = sigma_clipped_stats(flat_masked_dat)
    if method == 'median':
        back = median
    if method == 'mean':
        back = mean

    return back, backrms


def aperture_photometry_scan(data, x_pos, y_pos, ap_width, ap_length,
                             theta=0.0, show=False, plt_title=None):
    """Aperture photometry on source located on x_pos, y_pos with
    rectangular aperture of dimensions specified by ap_length, ap_width
    is used. Aperture sums are NOT sky subtracted.

    Parameters
    ----------
    data : `np.array`
        2D array of floats
    x_pos : float
        X position of source.
    y_pos : float
        Y position of source
    ap_width : int
        Width (along x axis) of photometric aperture.
    ap_length : int
        Length (along y axis) of photometric aperture.
    theta : float
        Angle of orientation (from x-axis) for aperture, in radians.
        Increases counter-clockwise.
    show : bool, optional
        If true, plot showing aperture(s) on source will pop up. Defaults to F.
    plt_title : str or None, optional
        Only used if `show` is True. Title for plot. Defaults to None.
    Returns
    -------
    phot_tab : `astropy.table`
        Table containing
    """

    copy_data = copy.copy(data)

    rect_ap = RectangularAperture((x_pos, y_pos), w=ap_width,
                                  h=ap_length, theta=theta)

    phot_table = aperture_photometry(copy_data, rect_ap,
                                     method='exact')
    if show:
        mask = rect_ap.to_mask(method='center')
        data_cutout = mask.cutout(data)
        plt.title(plt_title)
        z1, z2 = (-12.10630989074707, 32.53888328838081)
        plt.imshow(data, origin='lower', vmin=z1, vmax=z2)
        rect_ap.plot(color='white', lw=2)
        plt.show()
        plt.close()

    return phot_table
