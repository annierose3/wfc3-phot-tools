"""
Tools for computing backgroun statistics within photutils apertures.

Authors
-------
    - Varun Bajaj, December 2017
    - Clare Shanahan, December 2019

"""

from astropy.table import Table
import numpy as np
from scipy.stats import sigmaclip
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def _gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def calc_1d_gauss_background(data, bins=100, hist_range=(-10, 10)):
    """ 
    Fits a 1D gaussian distribution to pixels in `data`, returns mean of fit.

    Parameters
    -----------
    data : array
        Science array.
    bins : int
        Number of bins for histogram fit.
    hist_range : tuple of ints
        Range for fit (min, max).
    Returns
    -------
    coeff : tuple
        A, mu, sigma of gaussian fit. 

    """
    data = data.flatten()
    h, b = np.histogram(data, range=hist_range, bins=bins)
    bd = b[1]-b[0]
    bdist = int(1.1/bd)
    locs = find_peaks(h, distance=bdist)[0]
    p0 = [1., 0., 1.]
    centers = .5*b[:-1] + .5*b[1:]
    coeff, var_matrix = curve_fit(_gauss, centers[locs], h[locs], p0=p0)
    # Get the fitted curve
    hist_fit = _gauss(centers, *coeff)
    print(data[0:5])
    print(coeff)
    return coeff

def aperture_stats_tbl(data, apertures,
                       method='exact', sigma_clip=True):
    """Computes mean/median/mode/std in Photutils apertures.
    Compute statistics for custom local background methods.
    This is primarily intended for estimating backgrounds
    via annulus apertures.  The intent is that this falls easily
    into other code to provide background measurements.
    Parameters
    ----------
    data : array
        The data for the image to be measured.
    apertures : photutils PixelAperture object (or subclass)
        The phoutils aperture object to measure the stats in.
        i.e. the object returned via CirularAperture,
        CircularAnnulus, or RectangularAperture etc.
    method: str
        The method by which to handle the pixel overlap.
        Defaults to computing the exact area.
        NOTE: Currently, this will actually fully include a
        pixel where the aperture has ANY overlap, as a median
        is also being performed.  If the method is set to 'center'
        the pixels will only be included if the pixel's center
        falls within the aperture.
    sigma_clip: bool
        Flag to activate sigma clipping of background pixels
    Returns
    -------
    stats_tbl : astropy.table.Table
        An astropy Table with the colums X, Y, aperture_mean,
        aperture_median, aperture_mode, aperture_std, aperture_area
        and a row for each of the positions of the apertures.
    """

    # Get the masks that will be used to identify our desired pixels.
    masks = apertures.to_mask(method=method)

    # Compute the stats of pixels within the masks
    aperture_stats = [calc_aperture_mmm(data, mask, sigma_clip)
                      for mask in masks]

    aperture_stats = np.array(aperture_stats)

    # Place the array of the x y positions alongside the stats
    stacked = np.hstack([apertures.positions, aperture_stats])
    # Name the columns
    names = ['X', 'Y', 'aperture_mean', 'aperture_median', 'aperture_mode',
            'aperture_std', 'aperture_area']
    # Make the table
    stats_tbl = Table(data=stacked, names=names)

    return stats_tbl


def calc_aperture_mmm(data, mask, sigma_clip):
    """Helper function to actually calculate the stats for pixels
        falling within some Photutils aperture mask on some array
        of data.
    """
    cutout = mask.cutout(data, fill_value=np.nan)
    if cutout is None:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    else:
        values = cutout * mask.data / mask.data
        values = values[~np.isnan(values)]
        if sigma_clip:
            values, clow, chigh = sigmaclip(values, low=3, high=3)

        mean = np.mean(values)
        median = np.median(values)
        std = np.std(values)

        mode = 3 * median - 2 * mean
        actual_area = (~np.isnan(values)).sum()
        return (mean, median, mode, std, actual_area)