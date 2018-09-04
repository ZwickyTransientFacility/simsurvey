#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for drawing random redshifts and sky coordinates"""

import numpy as np
import random

_d2r = np.pi / 180 

__all__ = ["radec", "redshift",
           "simulate_lb", "simulate_z"]

# ============================== #
# = High Level Function        = #
# ============================== #
def radec(npoints=1,
          ra_range=(-180,180),dec_range=(-90,90),
          mw_exclusion=10,**kwargs):
    """
    """
    return np.asarray(simulate_lb(npoints,MW_exclusion=mw_exclusion,
                       ra_range=ra_range,dec_range=dec_range,
                       output_frame="j2000",**kwargs))

def redshift(npoints, zrange,
             pdfkind="flat",
             **kwargs):
    """
    """
    # Note for Uly: This redshift function could do much more
    # to parse easily the z_pdf and z_pdf_bins.
    # This through "pdfkind"
    # We can imagine having a string parser that is None
    
    if pdfkind.lower() in ["flat","None"]:
        pdfkind = None
        
    if pdfkind is None:
        # - The default Stuff
        z_pdf = kwargs.pop("z_pdf",None)
        z_pdf_bins = kwargs.pop("z_pdf_bins",None)
    elif type(pdfkind) is str:
        raise NotImplementedError("the 'pdfkind' could only be 'flat' or None")
    else:
        raise NotImplementedError("the 'pdfkind' could only be 'flat' or None")

    
    return np.asarray(simulate_z(npoints,zrange,z_pdf=z_pdf,z_pdf_bins=z_pdf_bins))


# ============================== #
# = Low level Functions        = #
# ============================== #
def simulate_lb(Npoints,MW_exclusion=10,ra_range=(-180,180),dec_range=(-90,90),
                output_frame='galactic',radius=None):
    """
    Draw a set of coordinates for particular RA and Dec range with MW exclusion 

    Arguments:
    Npoints -- number of coordinates to draw
    
    Keyword arguments:
    MW_exclusion -- redraw coordinates with b < this valu in degrees (default: 10)
    ra_range     -- range of RA distribution
    dec_range    -- range of DEC distribution
    output_frame -- output coordinate system ('galactic' or 'j2000')
    radius -- (r, l, b) force coordinates to be within r degrees of (l, b)
              Only works in galactic coordinates so far
    """
    # ----------------------- #
    # --                   -- #
    # ----------------------- #
    def _draw_radec_(Npoints_,ra_range_,dec_sin_range_):
        """
        """
        ra = np.random.random(Npoints_)*(ra_range_[1] - ra_range_[0]) + ra_range_[0]
        dec = np.arcsin(np.random.random(Npoints_)*(dec_sin_range_[1] - dec_sin_range_[0]) + dec_sin_range_[0]) / _d2r

        return ra,dec

    def _draw_without_MW_(Npoints_,ra_range_,dec_sin_range_,MW_exclusion_,radius_):
        """
        """
        
        l,b = np.array([]),np.array([])
        while( len(l) < Npoints_ ):
            ra,dec = _draw_radec_(Npoints_ - len(l),ra_range_,dec_sin_range_)
            l_,b_ = radec2gcs(ra,dec)

            if radius is not None:
                as_mask = ang_sep(radius[1], radius[2], l_, b_) < radius[0]
            else:
                as_mask = np.ones(len(l_), dtype=bool)

            mask = as_mask & (np.abs(b_)>MW_exclusion_)
            if output_frame == 'galactic':
                l = np.concatenate((l,l_[mask]))
                b = np.concatenate((b,b_[mask]))
            else:
                l = np.concatenate((l,ra[mask]))
                b = np.concatenate((b,dec[mask]))                

        return l,b

    # ----------------------- #
    # --                   -- #
    # ----------------------- #

    if output_frame not in ['galactic','j2000']:
        raise ValueError('output_frame must "galactic" or "j2000"')

    if ra_range[0] < -180 or ra_range[1] > 360 or ra_range[0] > ra_range[1]:
        raise ValueError('ra_range must be contained in [-180,360]')

    if dec_range[0] < -90 or dec_range[1] > 90 or dec_range[0] > dec_range[1]:
        raise ValueError('dec_range must be contained in [-90,90]')

    dec_sin_range = (np.sin(dec_range[0]*_d2r),np.sin(dec_range[1]*_d2r)) 

    if MW_exclusion > 0. or radius is not None:
        return _draw_without_MW_(Npoints, ra_range, dec_sin_range,
                                 MW_exclusion, radius)
    else:
        ra,dec = _draw_radec_(Npoints, ra_range, dec_sin_range)
        if output_frame == 'galactic':
            return radec2gcs(ra,dec)
        else:
            return ra,dec

def simulate_z(NPoints,z_range,z_pdf=None,z_pdf_bins=None):
    """
    Draw redshifts from distribution based on histogram

    Arguments:
    NPoints -- number of redshifts to draw
    z_range -- redshift range (tuple of length 2)

    Keyword arguments:
    z_pdf      -- redshift histogramm values (need not be normalized)
    z_pdf_bins -- redshift bins for z_pdf (must contain one more element 
                  than z_pdf)
    """
    if (len(z_range) != 2 or z_range[0] > z_range[1]):
        raise ValueError('Invalid z_range')
        
    if z_pdf is None:
        if z_pdf_bins is None:
            z_pdf = np.ones(1)
            z_pdf_bins = np.array(z_range)
            widths = np.array([z_range[1]-z_range[0]])
        else:
            z_pdf_bins = np.array(z_pdf_bins)
            z_pdf = np.ones(len(z_pdf_bins)-1)/(len(z_pdf_bins)-1)
    else:
        if z_pdf_bins is None:
            z_pdf_bins = np.linspace(z_range[0],z_range[1],len(z_pdf)+1)
        elif (np.abs(z_pdf_bins[0] - z_range[0]) / z_range[0] > 1e-9 
              or np.abs(z_pdf_bins[-1] - z_range[1]) / z_range[1] > 1e-9 
              or True in [a>b for a,b in zip(z_pdf_bins[:-1],z_pdf_bins[1:])]):
            print(np.abs(z_pdf_bins[0] - z_range[0]) / z_range[0] > 1e-9)
            print(np.abs(z_pdf_bins[-1] - z_range[1]) / z_range[1] > 1e-9) 
            print([a>b for a,b in zip(z_pdf_bins[:-1],z_pdf_bins[1:])])
            print(True in [a>b for a,b in zip(z_pdf_bins[:-1],z_pdf_bins[1:])])
            raise ValueError('Invalid z_pdf_bins')
        else:
            z_pdf_bins = np.array(z_pdf_bins)

    widths = z_pdf_bins[1:]-z_pdf_bins[:-1]
    z_pdf = np.array(z_pdf,dtype=float)/np.sum(np.array(z_pdf*widths))

    if len(z_pdf) > 1:
        z_cdf = np.cumsum(z_pdf*widths)
        val_uni = np.random.random(NPoints)
        val_bins = np.array([np.where(z_cdf > val)[0][0] for val in val_uni])
        val_rem = ((val_uni - z_cdf[val_bins-1])%1)/((z_cdf[val_bins]-z_cdf[val_bins-1])%1)

        z = z_pdf_bins[val_bins] + (z_pdf_bins[val_bins+1] - z_pdf_bins[val_bins]) * val_rem
    else:
        z = np.random.random(NPoints) * (z_range[1]-z_range[0]) + z_range[0]

    return z

# ----------------------------------------------------- #
# -- Required functions that might go somewhere else -- #
# ----------------------------------------------------- #

def ang_sep(l1,b1,l2,b2):
    """
    Angular separation between two positions on the sky 
    (l1,b1) and (l2,b2) in degrees.
    """
    sin_theta = np.sqrt((np.cos(b2 * _d2r) * np.sin((l1 - l2) * _d2r)) ** 2 +
                        (np.cos(b1 * _d2r) * np.sin(b2 * _d2r) - 
                         np.sin(b1 * _d2r) * np.cos(b2 * _d2r) * np.cos((l1 - l2) * _d2r)) ** 2)
    cos_theta = (np.cos(b1 * _d2r) * np.cos(b2 * _d2r) *
                 np.cos((l1 - l2) * _d2r) +
                 np.sin(b1 * _d2r) * np.sin(b2 * _d2r))
    return np.arctan2(sin_theta,cos_theta) / _d2r

# -------------------------------- #
# ----  FROM THE SNf ToolBox ----- #
# -------------------------------- #

def radec2gcs(ra, dec, deg=True):
    """
    Authors: Yannick Copin (ycopin@ipnl.in2p3.fr)
    
    Convert *(ra,dec)* equatorial coordinates (J2000, in degrees if
    *deg*) to Galactic Coordinate System coordinates *(lII,bII)* (in
    degrees if *deg*).

    Sources:

    - http://www.dur.ac.uk/physics.astrolab/py_source/conv.py_source
    - Rotation matrix from
      http://www.astro.rug.nl/software/kapteyn/celestialbackground.html

    .. Note:: This routine is only roughly accurate, probably at the
              arcsec level, and therefore not to be used for
              astrometric purposes. For most accurate conversion, use
              dedicated `kapteyn.celestial.sky2sky` routine.

    >>> radec2gal(123.456, 12.3456)
    (210.82842704243518, 23.787110745502183)
    """

    if deg:
        ra  =  ra * _d2r
        dec = dec * _d2r

    rmat = np.array([[-0.054875539396, -0.873437104728, -0.48383499177 ],
                    [ 0.494109453628, -0.444829594298,  0.7469822487  ],
                    [-0.867666135683, -0.198076389613,  0.455983794521]])
    cosd = np.cos(dec)
    v1 = np.array([np.cos(ra)*cosd,
                  np.sin(ra)*cosd,
                  np.sin(dec)])
    v2 = np.dot(rmat, v1)
    x,y,z = v2

    c,l = rec2pol(x,y)
    r,b = rec2pol(c,z)

    assert np.allclose(r,1), "Precision error"

    if deg:
        l /= _d2r
        b /= _d2r

    return l, b

def rec2pol(x,y, deg=False):
    """
    Authors: Yannick Copin (ycopin@ipnl.in2p3.fr)
    
    Conversion of rectangular *(x,y)* to polar *(r,theta)*
    coordinates
    """

    r = np.hypot(x,y)
    t = np.arctan2(y,x)
    if deg:
        t /= _d2r

    return r,t
