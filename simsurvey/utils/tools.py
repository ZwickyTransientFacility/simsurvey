#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module gather the basic tools used in a lot of methods"""

import warnings
import numpy as np

__all__ = ["kwargs_update", "kwargs_extract",
           "range_args", "range_length",
           "load_pkl", "dump_pkl"]


def kwargs_update(default,**kwargs):
    """
    """
    k = default.copy()
    for key,val in kwargs.items():
        k[key] = val
        
    return k

def kwargs_extract(default,**kwargs):
    """
    like kwargs_update but extracts keys of default from kwargs

    Returns:
    k -- dictionary based on default update for kwargs
    l -- kwargs without keys defined in default
    """
    k = default.copy()
    l = {}
    for key,val in kwargs.iteritems():
        if key in k.keys():
            k[key] = val
        else:
            l[key] = val

    return k, l

def range_args(n_max, *args):
    """Process args similar to those range(), i.e. (start, end, step)
    but provide n_max to limit end
    """
    if len(args) == 0:
        start, end, step = 0, n_max, 1
    elif len(args) == 1:
        start, end, step = 0, args[0], 1
    elif len(args) == 2:
        start, end, step = args[0], args[1], 1
    elif len(args) == 3:
        start, end, step = args
    else:
        raise TypeError('range_args requires 1-4 int arguments')
        
    if end > n_max:
        end = n_max
        warnings.warn('only %i items were available'%n_max)

    return start, end, step

def range_length(start, end, step):
    if step > 0:
        lo, hi = start, end
    else:
        hi, lo = start, end
        step = -step

    if lo >= hi:
        return 0
    else:
        return (hi - lo - 1) // step + 1
    
    
# --------------------------- #
# - I/O Tools               - #
# --------------------------- #
def ipython_info():
    import sys
    return 'notebook' if 'ipykernel' in sys.modules \
      else 'terminal' if 'Ipython' in sys.modules \
      else None

def load_pkl(filename):
    """
    """
    import cPickle as pkl
    try:
        pkl_file = open(filename,'rb')
    except:
        raise IOError("The given file does not exist %s"%filename)
    
    return pkl.load(pkl_file)


def dump_pkl(data,filename,**kwargs):
    """
    """
    from cPickle import dump
    if len(filename.split("."))>1 and filename.split(".")[-1]=="pkl":
        outfile =  open(filename,"wb")
    else:
        outfile =  open(filename+".pkl","wb")
    
    dump(data, outfile,**kwargs)
    outfile.close()

def fitsrec_to_dict(data):
    fields = data.dtype.fields.keys()
    dico = {}
    for f in fields:
        dico[f] = data[f]
    return dico

# --------------------------- #
# - Conversion Tools        - #
# --------------------------- #
def flux_to_mag(flux,dflux,wavelength):
    """
    """
    F_Lxlambda2  = flux * wavelength**2
    if dflux is None:
        return -2.5*np.log10(F_Lxlambda2) - 2.406, None
    
    err = -2.5/np.log(10) * dflux / flux
    
    return -2.5*np.log10(F_Lxlambda2) - 2.406, np.abs(err)

def mag_to_flux(mag,magerr,wavelength):
    """
    mag must be ABmag
    wavelength in Angstrom
    return Flux in erg/s/cm2/A
    """
    F_Lxlambda2 = 10**(-(mag+2.406)/2.5)
    flux = F_Lxlambda2/wavelength**2
    if magerr is None:
        return flux
    
    dflux = np.abs(flux*(-magerr/2.5*np.log(10))) # df/f = dcount/count
    return flux,dflux

def hourangle_2_degree(ra_hours,dec_hours,obstime="J2000"):
    """given the coordinate of the target in hours units
    (e.g. 10:58:59.072 +46:40:25.23) this will retour them
    in degree"""
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    c = SkyCoord("%s %s"%(ra_hours,dec_hours),
                 unit=(u.hourangle, u.deg), obstime=obstime)
    return c.ra.value,c.dec.value


# --------------------------- #
# - Array Tools             - #
# --------------------------- #
def shape_ajustment(X,Y,model_X,k=4,s=0,
                    verbose=False):
    """ DOC TO BE DONE
    return  Y /w the same binning as model_X  in their commun wavelength
    """
    def commun_wavelength(x1,x2):
        """
        """
        flagx1 = (x1>=x2.min()) & (x1<=x2.max())
        flagx2 = (x2>=x1.min()) & (x2<=x1.max())
        return flagx1, flagx2

    from scipy.interpolate import UnivariateSpline
    flagX,flagmodel = commun_wavelength(X,model_X)
    Yrebin = UnivariateSpline(X[flagX], Y[flagX],k=k,s=s)(model_X)

    if len(Yrebin)==len(model_X):
        return Yrebin
    else:
        if verbose:
            print('WARNING [shape_adjustment] non-mached shape ... I am fixing that')

        YrebinOK = np.empty((len(Yrebin)+1),)
        YrebinOK[1:] = Yrebin
        YrebinOK[0]  = Yrebin[0]

        return YrebinOK

# --------------------------- #
# - Array Tools             - #
# --------------------------- #
def get_progressbar(gen, notebook=False):
    """
    """
    from astropy.utils.console import ProgressBar

    if not notebook:
        gen = ProgressBar(gen)
    else:
        try:
            gen = ProgressBar(gen, ipython_widget=True)
        except ImportError as e:
            warnings.warn('ProgressBar in notebook not working. Is ipywidgets installed?')
            raise e

    return gen
