#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" Based on provided simumation this module enables to retreive the input values """

"""
Simulated values input for the fitting techniques:

* pixel coordinates: x, y
* observed magnitudes: g, r, g_err, r_err (i, i_err)
* ref magnitudes (gaia): gref, rref, gref_err, rref_err (iref,iref_err)
* observational condition: airmass
* Bandpasses: ZTF, Gaia
"""

import numpy as np
import matplotlib.pyplot as mpl
import pandas as pd

from astrobject.astrobject.collections.photospatial import PhotoMap



def get_simulphotomap(filecsv):
    """ to be defined """


###########################################
#                                         #
#  Spatial Variations                     #
#                                         #
###########################################
def get_polycorr(poly,params,x,y):
    """ get a 2d polynomial for the given parameters
    (the degree to amount of parameters is the following
    {1=>3 ; 2=>6 ; 3=>10 etc.})
    """
    
    return None


###########################################
#                                         #
#  Atmosphere                             #
#                                         #
###########################################
def get_atmtransmission(extmodel, ozoneint, aerosol, aerosolexp,
                        pressure,airmass):
    """ returns the atmosphere transmisson in the given condition
    for the given extinction model (Buton et al.)
    (Transmission = 10**(-0.4* _extinction_ * airmass ))

    Parameters:
    ----------
    ozoneint, aerosol, aerosolexp: floats
        parameters of the Buton et al. model:
        - ozoneint   (oi): ozone intensity [Dobson units]
                           typically 260 [220,300] at Mauna Kea (season dependent)
        - aerosol    (ai): aerosol optical depth at reference wavelength
                           typically 0.008 [0.005,0.02] at Mauna Kea
        - aerosolexp (ap): aerosol angstrom exponent
                           typically 1.8 [-0.5,4] at Mauna Kea
                           
    pressure: float
        last parameters of the Buton et al model. This might be given
        - pressure   (pr): surface pressure [mbar]
        
    airmass: float
        The airmass of the observation

    Return
    -------
    Transmission of the atmosphere (vector similar to wavelength)
    """
    em.setParams(pressure, ozoneint, aerosol, aerosolexp)
    return 10**(-0.4*em.extinction()*airmass)

    
def atmophere_transmission( ozoneint, aerosol, aerosolexp,
                     pressure,airmass, wavelength):
    """ The transmission of the atmosphere in the given condition
    at the given wavelength.
    (Transmission = 10**(-0.4* _extinction_ * airmass );
    see get_atmtransmission)

    _extinction_ model is based on Buton et al model.
    This requieres pyExtinction
    http://snfactory.in2p3.fr/soft/atmosphericExtinction/

    Parameters:
    ----------
    ozoneint, aerosol, aerosolexp: floats
        parameters of the Buton et al. model:
        - ozoneint   (oi): ozone intensity [Dobson units]
                           typically 260 [220,300] at Mauna Kea (season dependent)
        - aerosol    (ai): aerosol optical depth at reference wavelength
                           typically 0.008 [0.005,0.02] at Mauna Kea
        - aerosolexp (ap): aerosol angstrom exponent
                           typically 1.8 [-0.5,4] at Mauna Kea
                           
    pressure: float
        last parameters of the Buton et al model. This might be given
        - pressure   (pr): surface pressure [mbar]
        
    airmass: float
        The airmass of the observation

    wavelength: float/vector
        wavelength of definition of the atm model

    Return
    -------
    Transmission of the atmosphere (vector similar to wavelength)
    """
    try:
        import pyExtinction
    except ImportError:
        raise ImportError("please install pyExtinction"+\
                " (http://snfactory.in2p3.fr/soft/atmosphericExtinction/)")

    em = pyExtinction.AtmosphericExtinction.ExtinctionModel(wavelength)
    return get_atmtransmission(em,ozoneint, aerosol, aerosolexp,
                               pressure,airmass)


