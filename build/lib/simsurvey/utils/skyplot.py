#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Assorted functions for setting up the skyplots, e.g. prepating an axis with
projection, transforming coordinates from degrees to radians and reverse
longitude direction to get a celestial plot.
"""

import numpy as np
import matplotlib.pyplot as mpl
import warnings
from copy import deepcopy

_d2r = np.pi / 180 

__all__ = ["ax_skyplot"]

# ============================== #
# = Axis setup                 = #
# ============================== #
def ax_skyplot(fig=None, figsize=(12, 6), rect=[0.1, 0.1, 0.8, 0.8], 
               projection='mollweide', xlabelpad=None, xlabelmode='hist'): 
    """
    Initialize axis for skyplot and make grid and labels nicer.
    [Fill in kwargs]
    """
    # Work around to get correct xlabel position in older matplotlib
    import matplotlib as m
    mv = m.__version__.split('.')
    if int(mv[0]) < 2 and int(mv[1]) < 5:
        xlabeladjust = {'show': 180, 
                        'hist': 165}
        
        if xlabelmode not in xlabeladjust.keys():
            raise ValueError('Unknown xlabelmode. Know values: %2'%
                             (', '.join(xlabeladjust.keys())))
            
        warnings.warn('You are using matplotlib version < 1.5.0. '+
                      'The padding of the xlabel has been adjusted. '+
                      'You can use the option "xlabelpad" to adjust.')

        if xlabelpad is None:
            xlabelpad = xlabeladjust[xlabelmode]
        else:
            xlabelpad += xlabeladjust[xlabelmode]

    allowed_proj = ['mollweide', 'hammer']

    if fig is None:
        fig = mpl.figure(figsize=figsize)

    if projection not in allowed_proj:
        raise ValueError("Projection not supported; allowed values: %s"
                         %','.join(allowed_proj))

    ax = fig.add_axes(rect, projection=projection)
    ax.grid(True)
    xlabels = [u'%i\xb0'%ra for ra in range(150,-1,-30) + range(330,209,-30)]
    ax.set_xticklabels(xlabels)

    ax.set_xlabel(r"$\mathrm{RA\ [deg]}$", fontsize="x-large", 
                  labelpad=xlabelpad)
    ax.set_ylabel(r"$\mathrm{Dec\ [deg]}$", fontsize="x-large")

    return fig, ax

# ============================== #
# = Conversion                 = #
# ============================== #
def convert_radec_azel(ra, dec, edge=0):
    """
    Convert ra, dec to azimuth and elavation in radians as used in matplotlib 
    projections and switch sign of ra to get celestial plot.
    
    edge -- can be used to set move points at exactly ra = -180 or 180
            slightly off that
    
    [This could be extended to also convert between coordinate systems.]
    """
    #Make sure RA is between -180 and 180, then invert axis
    if edge > 0:
        if type(ra) == float:
            if ra < -180 + edge:
                ra = -180 + edge
            elif ra > 180 - edge:
                ra = 180 - edge
        else:
            ra[ra < -180 + edge] = -180 + edge
            ra[ra > 180 - edge] = 180 - edge

    ra = ((ra + 180) % 360) - 180
    ra *= -1

    az = _d2r * ra
    el = _d2r * dec

    return az, el

def cart2sph(vec, cov=None):
    """
    Convert vector in Cartesian coordinates to spherical coordinates 
    (angles in degrees). Convariance matrix can be converted as well
    if it is stated.
    """
    x = vec[0]
    y = vec[1] 
    z = vec[2] 

    v = np.sqrt(x**2 + y**2 + z**2)
    v_sph = np.array([v, (np.arctan2(y,x) / _d2r + 180) % 360 - 180, 
                          np.arcsin(z/v) / _d2r])
    
    if cov is None:
        return v_sph
    else:
        jacobian = np.zeros((3,3))
        jacobian[0,0] = x / v
        jacobian[1,0] = - y / (x**2 + y**2)
        jacobian[2,0] = - x * z / (v**2 * np.sqrt(x**2 + y**2))
        jacobian[0,1] = y / v
        jacobian[1,1] = x / (x**2 + y**2)
        jacobian[2,1] = - y * z / (v**2 * np.sqrt(x**2 + y**2))
        jacobian[0,2] = z / v
        jacobian[1,2] = 0
        jacobian[2,2] = np.sqrt(x**2 + y**2) / (v**2)

        cov_sph = (jacobian.dot(cov)).dot(jacobian.T)
        cov_sph[1,1] /= _d2r**2
        cov_sph[2,2] /= _d2r**2
        cov_sph[2,1] /= _d2r**2
        cov_sph[1,2] /= _d2r**2
    
        cov_sph[0,1] /= _d2r
        cov_sph[0,2] /= _d2r
        cov_sph[1,0] /= _d2r
        cov_sph[2,0] /= _d2r    

        return v_sph, cov_sph

def sph2cart(vec, cov=None):
    """
    Convert vector in spherical coordinates (angles in degrees)
    to Cartesian coordinates. Convariance matrix can be converted as well
    if it is stated.
    """
    v = vec[0]
    l = vec[1]*_d2r
    b = vec[2]*_d2r

    v_cart = np.array([v*np.cos(b)*np.cos(l), v*np.cos(b)*np.sin(l), 
                       v*np.sin(b)])    

    if cov is None:
        return v_cart
    else:
        cov_out = deepcopy(cov)
        cov_out[1,1] *= _d2r**2
        cov_out[2,2] *= _d2r**2
        cov_out[2,1] *= _d2r**2
        cov_out[1,2] *= _d2r**2
        cov_out[0,1] *= _d2r
        cov_out[0,2] *= _d2r
        cov_out[1,0] *= _d2r
        cov_out[2,0] *= _d2r

        jacobian = np.zeros((3,3))
        jacobian[0,0] = np.cos(b) * np.cos(l)
        jacobian[1,0] = np.cos(b) * np.sin(l)
        jacobian[2,0] = np.sin(b)
        jacobian[0,1] = - v * np.cos(b) * np.sin(l)
        jacobian[1,1] = v * np.cos(b) * np.cos(l)
        jacobian[2,1] = 0
        jacobian[0,2] = - v * np.sin(b) * np.cos(l)
        jacobian[1,2] = - v * np.sin(b) * np.sin(l)
        jacobian[2,2] = v * np.cos(b)

        cov_cart = (jacobian.dot(cov_out)).dot(jacobian.T)

        return v_cart, cov_cart

def rot_xz(v, theta):
    """
    Rotate Cartesian vector v by angle theta around axis (0,1,0)
    """
    return np.array([v[0]*np.cos(theta*_d2r) - v[2]*np.sin(theta*_d2r),
                     v[1],
                     v[2]*np.cos(theta*_d2r) + v[0]*np.sin(theta*_d2r)])

def rot_xz_sph(l, b, theta):
    """
    Rotate Spherical coordinate (l,b) by angle theta around axis (0,1,0)
    """
    v_cart = sph2cart([1,l,b])
    v_rot = rot_xz(v_cart, theta)
    return cart2sph(v_rot)[1:]
