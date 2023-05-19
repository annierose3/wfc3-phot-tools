"""
This module contains utilities for performing aperture
photometry for staring mode observations.

Authors
-------
    Clare Shanahan, December 2019
    Mariarosa Marinelli, April 2023
"""

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astroquery.simbad import Simbad

def query_simbad_by_targname(targname):
    """Construct `SkyCoord` using SIMBAD query of target.

    Queries SIMBAD for `targname`, and constructs a
    `astropy.coordinates.SkyCoord` object containing the
    RA, Dec, proper motion, radial velocity, and distance
    to the object.

    Paramters
    ---------
    targname : str
        Target name with which to query SIMBAD.

    Returns
    -------
    c : `astropy.coordinates.SkyCoord`
        SkyCoord object with query results.
    """

    customSimbad = Simbad()
    customSimbad.add_votable_fields('pmra', 'pmdec', 'distance', 'rv_value')
    result_table = customSimbad.query_object(targname)

    simbad_ra = str(result_table['RA'].item())
    ra = f'{simbad_ra[:2]}h{simbad_ra[3:5]}m{simbad_ra[6:]}s'

    simbad_dec = str(result_table['DEC'].item())
    dec = f'{simbad_dec[:3]}h{simbad_dec[4:6]}m{simbad_dec[7:]}s'

    pm_ra = result_table['PMRA'].item()*u.mas/u.yr
    pm_dec = result_table['PMDEC'].item()*u.mas/u.yr

    c = SkyCoord(ra, dec, frame='icrs', obstime=Time('J2000'),
                 distance=result_table['Distance_distance'].item()*u.pc,
                 pm_ra_cosdec=pm_ra,
                 pm_dec=pm_dec)

    return c


def apply_proper_motion_targ(targname, mjd):
    """Apply PM to RA/Dec of target.

    Queries Simbad for ra, dec, and proper motions for
    'targname'. Returns (RA, and Dec) at date `mjd`
    considering the proper motion.

    Parameters
    ---------
    targname : str
        Target name. Simbad is queried by target name.
    mjd : float
        MJD for which to calculate RA, and Dec.

    Returns
    -------
    (ra_new, dec_new) : tuple of floats
        RA and Dec for `targname` after applying proper
        motion. In degrees.
    """
    c = query_simbad_by_targname(targname)

    delta_t_yr = (mjd - 51544) / 365. * u.yr
    ra_new = c.ra.deg * u.deg + (c.pm_ra_cosdec * delta_t_yr)
    dec_new = c.dec.deg * u.deg + (c.pm_dec * delta_t_yr)

    return ra_new.value, dec_new.value
