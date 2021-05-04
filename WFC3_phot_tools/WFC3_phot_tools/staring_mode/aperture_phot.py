from astropy.stats import mad_std
from ginga.util import zscale
from photutils import DAOStarFinder
from photutils import aperture_photometry, CircularAperture, CircularAnnulus
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import (DAOStarFinder, aperture_photometry,  deblend_sources,
                       detect_sources, detect_threshold, source_properties)

import warnings
warnings.filterwarnings("ignore")


def show_source_detection_plot(data, coo_tab):
    """
    Displays a plot of `data` (z-scaled), with sources in `coo_tab` overplotted.

    Parameters
    ----------
    data : array
        Image array.
    coo_tab : `astropy.table.Table`
        Table with source positions as columns 'xcentroid' and 'ycentroid'.

    Returns
    -------
    Window with plot.

    """

    z1, z2 = zscale.zscale(data)

    plt.imshow(data, origin='lower', cmap='Greys_r', vmin=z1, vmax=z2)
    if coo_tab:
        plt.scatter(coo_tab['xcentroid'], coo_tab['ycentroid'], c='r')
        plt.title('{} Sources Detected'.format(len(coo_tab)))
    else:
        plt.title('No Sources Detected')
    plt.show()


def detect_sources_segmap(data, threshold, npixels, kernel_fwhm=1.8,
                          show=False):
    """
    Runs image segmentation to detect sources in `data`.

    Parameters
    ----------
    data : array
        Image array.
    threshold : float or array
        Detection threshold value, or pixel-wise threshold image (must be same
        shape as `data`.)
    npixels : int
        Positive integer number of connected pixels, each greater that
        `threshold` that an object must have to be detected.
    kernel_fwhm : float
        FWHM of gaussian kernel used to smooth image before segmentation.
    show : bool
        Show a plot of detected source(s).

    Returns
    -------
    coo_tab : `astropy.table.Table` or int
        Table with detected source(s). Returns '0' if no sources are detected.

    """

    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(data, threshold=threshold, npixels=npixels,
                          filter_kernel=kernel)
    if segm:
        coo_tab = source_properties(data, segm).to_table()
        if show:
            show_source_detection_plot(data, coo_tab)
        return coo_tab

    if not segm:
        if show:
            show_source_detection_plot(data, None)
        return 0


def calc_sky_annulus(data, x, y, r_in, r_out, sky_method='median',
                     sigma_clip=True):
    """
        Calculates the background level in a circular annulus centered
        at `x, y` with an inner radius of `r_in` pixels and an outer radius of
        `r_out` pixels.

        Sky level and standard deviation in annulus are returned.
        Options for `sky_method` are `mean`, `median`, or `mode`.
        Uses a sigma clip with sigma=`n_sigma_clip` to compute stats.
        If n_sigma_clip is set to 0, sigma clipping won't be done.

        Parameters
        ----------
            data : array
                Science data array.
            x, y : float
                Pixel position for annulus center.
            r_in, r_out : int
                Inner, outer circular annulus radii
            sky_method : str
                'mean', 'median', 'mode', or 'gaussian'
            n_sigma_clip : int
                Used for sigma clipping. 0 to turn off sigma clipping

        Returns
        -------
            (back, backstd) : tuple
                Sky level in annulus, std of sky.

    """
    from .background import aperture_stats_tbl, calc_1d_gauss_background

    sky_apertures = CircularAnnulus((x, y), r_in, r_out)

    if sky_method in ['mean', 'median', 'mode']:
        sky_tbl = aperture_stats_tbl(data, sky_apertures, sigma_clip=sigma_clip)
        return(sky_tbl['aperture_'+sky_method].item(),
               sky_tbl['aperture_std'].item())

    elif sky_method == 'gaussian':
        print('fitting gaussian to data for background.')
        try:
            A, mu, sigma = calc_1d_gauss_background(data, bins=100)
        except RuntimeError:
            try:
                A, mu, sigma = calc_1d_gauss_background(data, bins=100,
                                                        hist_range=(-20, 10))
            except RuntimeError:
                try:
                    A, mu, sigma = calc_1d_gauss_background(data, bins=100,
                                                            hist_range=(-10, 20))
                except RuntimeError:
                    try:
                        A, mu, sigma = calc_1d_gauss_background(data, bins=100,
                                                                hist_range=(-5, 5))
                    except RuntimeError:
                        A, mu, sigma = -999, -999, -999
                        print("Couldn't fit gaussian")

        return(mu, sigma)


def circular_aperture_photometry(data, xc, yc, aperture_radii):
    """
    Aperture photometry with circular apertures on single source in `data`
    located at (xc, yc).

    Returns a table of aperture sums for every aperture size in
    `aperture_radii`.

    Parameters
    ----------
    data : array
        Image array.
    xc, yc : float
        x, y location of source in `data`.
    aperture_radii : float or list of floats
        Desired aperture radii, in pixels.

    Returns
    -------
    phot_table : `astropy.table.Table`
        Table with columns for x&y source position, and sum in every circular
        aperture in `aperture_radii`.

    """
    if type(aperture_radii) == float:
        aperture_radii = [aperture_radii]
    aps = [CircularAperture((xc, yc), r=rad) for rad in aperture_radii]
    phot_table = Table(aperture_photometry(data, aps))
    phot_table.remove_column('id')

    table_order = []
    for i, rad in enumerate(aperture_radii):
        rad = str(rad)
        phot_table.rename_column('aperture_sum_{}'.format(str(i)),
                                 'countrate_{}'.format(rad))
        table_order.append('countrate_{}'.format(rad))

    table_order = ['xcenter', 'ycenter'] + table_order

    phot_table = phot_table[table_order]

    return(phot_table)
