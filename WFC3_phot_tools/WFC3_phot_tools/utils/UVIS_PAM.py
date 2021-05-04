from astropy.io import fits
import copy
import numpy as np


def make_PAMcorr_image_UVIS(data, prihdr, scihdr, pamdir):
    """Creates the Pixel Area Map (PAM) image.
    Parameters:
        data : array
            Name of FITS file.
        pri : header
            Primary header of file for data.
        scihdr : header
            Header from science extension of data.
        pamdir : str
            Path to where pixel area maps for UVIS1 and/or UVIS2 are located.
    Returns:
        pamcorr_data : array
            PAM corrected data
    """

    data = copy.copy(data)
    x0 = int(np.abs(scihdr['LTV1']))
    y0 = int(np.abs(scihdr['LTV2']))
    x1 = int(x0 + scihdr['NAXIS1'])
    y1 = int(y0 + scihdr['NAXIS2'])

    if scihdr['CCDCHIP'] == 1:
        pam = fits.getdata(pamdir + 'UVIS1wfc3_map.fits')
        pamcorr_data = data * pam[y0:y1, x0:x1]

    elif scihdr['CCDCHIP'] == 2:
        pam = fits.getdata(pamdir + 'UVIS2wfc3_map.fits')
        pamcorr_data = data * pam[y0:y1, x0:x1]
    else:
        raise Exception('Chip case not handled.')

    return pamcorr_data