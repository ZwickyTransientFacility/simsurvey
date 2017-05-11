#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Useful sncosmo model/source definitions"""

import warnings
import numpy as np
import cPickle
from copy import deepcopy
from collections import OrderedDict as odict
from itertools import izip

from scipy.constants import speed_of_light, Planck, Boltzmann
from scipy.interpolate import (InterpolatedUnivariateSpline as Spline1d,
                               RectBivariateSpline as Spline2d)

import sncosmo

__all__ = ["BlackBodySource", "MultiSource"]

#######################################
#                                     #
# New source classes                  #
#                                     #
#######################################
class BlackBodySource(sncosmo.Source):
    """A simple spectral source for a transient with a 
    black-body spectrum of constant temperature
    The spectral flux density of this model is given by
    .. math::
       F(t, \lambda) = A \\times LC(t / s) \\times B(\lambda, T)
    where _A_ is the amplitude, _s_ is the "stretch" and 
    _T_ is the black-body temperature; LC is a dimensionless 
    lightcurve that defines the time evolution of the flux 
    while B is the black-body spectrum.
    Parameters
    ----------
    phase : `~numpy.ndarray`
        Phases in days.
    flux : `~numpy.ndarray`
        Model spectral flux density in arbitrary units.
        Must have shape `(num_phases)`.
    """

    _param_names = ['amplitude', 's', 'T']
    param_names_latex = ['A', 's', 'T']

    def __init__(self, phase, flux, name=None, version=None):
        self.name = name
        self.version = version
        self._phase = phase
        self._parameters = np.array([1., 1., 2e4])
        self._model_lc = Spline1d(phase, flux, k=2)

    def minwave(self):  
        return 1e-100

    def maxwave(self):  
        return 1e100
    
    def minphase(self):
        return self._parameters[1] * self._phase[0]

    def maxphase(self):
        return self._parameters[1] * self._phase[-1]

    def _flux(self, phase, wave):
        return np.outer(self._parameters[0] *
                        self._model_lc(phase / self._parameters[1]),
                        blackbody(wave, self._parameters[2]))

class MultiSource(sncosmo.Source):
    """A simple spectral source for a transient with a 
    black-body spectrum of constant temperature
    The spectral flux density of this model is given by
    .. math::
       F(t, \lambda) = A \\times LC(t / s) \\times B(\lambda, T)
    where _A_ is the amplitude, _s_ is the "stretch" and 
    _T_ is the black-body temperature; LC is a dimensionless 
    lightcurve that defines the time evolution of the flux 
    while B is the black-body spectrum.
    Parameters
    ----------
    phase : `~numpy.ndarray`
        Phases in days.
    flux : `~numpy.ndarray`
        Model spectral flux density in arbitrary units.
        Must have shape `(num_phases)`.
    """   
    _param_names = ['amplitude', 's', 'template_index']
    param_names_latex = ['A', 's', 'k']

    def __init__(self, phase, wave, flux, zero_before=False, name=None,
                 version=None, kx=3, ky=3):
        self.name = name
        self.version = version
        self._loaded_index = None

        self._parameters = np.array([1., 1., 0.])
        self._zero_before = zero_before
        self._spline_k = (kx, ky)

        self._phase = phase
        self._wave = wave

        self._model_flux = [Spline2d(p_, w_, f_, kx=kx, ky=ky)
                            for p_, w_, f_ in zip(phase, wave, flux)]

    def minwave(self):
        return self._parameters[1] * self._wave[self._k][0]

    def maxwave(self):
        return self._parameters[1] * self._wave[self._k][-1]

    def minphase(self):
        return self._parameters[1] * self._phase[self._k][0]

    def maxphase(self):
        return self._parameters[1] * self._phase[self._k][-1]
        
    def _flux(self, phase, wave):
        #f = self._parameters[0] * self._model_flux[k](phase, wave)
        
        #if self._zero_before:
        #    mask = np.atleast_1d(phase) < self.minphase()
        #    f[mask, :] = 0.

        #return f

        return (self._parameters[0] *
                self._model_flux[self._k](phase / self._parameters[1], wave))

    @property
    def _k(self):
        """template index"""
        return int(self._parameters[2])

#######################################
#                                     #
# Auxilary functions                  #
#                                     #
#######################################
def blackbody(wl, T=2e4):
    # wl in Ångströms
    # T in Kelvin
    # output is not normalized because sncosmo amplitude will be scaled to some
    # arbitrary value
    return 1.19104295262e+34 * wl**-5 / (np.exp(1.43877735383e+8 / (wl*T)) - 1)
