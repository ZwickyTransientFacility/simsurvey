#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Useful sncosmo model/source definitions"""

import warnings
import numpy as np
import pickle
from copy import deepcopy
from collections import OrderedDict as odict

try:
    from itertools import izip
except ImportError: #Python3.x
    izip = zip

from scipy.constants import speed_of_light, Planck, Boltzmann
from scipy.interpolate import (InterpolatedUnivariateSpline as Spline1d,
                               RectBivariateSpline as Spline2d)

import sncosmo

__all__ = ["BlackBodySource", "ExpandingBlackBodySource",
           "MultiSource", "SpectralIndexSource"]

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
    """
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

class ExpandingBlackBodySource(sncosmo.Source):
    """
    """
    def __init__(self, name=None, version=None, minphase=0.1, maxphase=30.,
                 tempfunc=(lambda p, x: p[0] * np.exp(p[1]*x)),
                 radiusfunc=(lambda p, x: p[0] + p[1] * x),
                 tempparam=[6e4, 1.], radiusparam=[5000, 1000]):
        self.name = name
        self.version = version
        self._minphase = minphase
        self._maxphase = maxphase

        self._tempfunc = tempfunc
        self._radiusfunc = radiusfunc

        self._n_tempparam = len(tempparam)
        self._n_radiusparam = len(radiusparam)

        self._parameters = np.array([1e-5])
        self._parameters = np.append(self._parameters, tempparam)
        self._parameters = np.append(self._parameters, radiusparam)

        self._param_names = ['d']
        self.param_names_latex = ['d']

        for k in range(self._n_tempparam):
            self._param_names.append('T%i'%k)
            self.param_names_latex.append('T_%i'%k)

        for k in range(self._n_radiusparam):
            self._param_names.append('R%i'%k)
            self.param_names_latex.append('R_%i'%k)

    def minwave(self):  
        return 1e2

    def maxwave(self):  
        return 1e6
    
    def minphase(self):
        return self._minphase

    def maxphase(self):
        return self._maxphase

    def temperature(self, phase):
        param = self._parameters[1:self._n_tempparam+1]
        T = self._tempfunc(param, phase)
        try:
            T[(T < 1.)] = 1.
        except TypeError:
            if T < 1.:
                return 1.
        return T
        
    def radius(self, phase):
        param = self._parameters[self._n_tempparam+1:
                                 self._n_tempparam+self._n_radiusparam+1]
        R = self._radiusfunc(param, phase)
        try:
            R[(phase < self.minphase())
              | (phase > self.maxphase()) | (R < 0.)] = 0.
        except TypeError:
            if phase < self.minphase() or phase > self.maxphase() or R < 0.:
                return 0.
        return R
    
    def luminosity(self, phase):
        return (0.5670367 * self.temperature(phase)**4
                * 4 * np.pi * (self.radius(phase)*6.957e8) **2)
    
    def _flux(self, phase, wave):
        wave = np.array(wave)
        out = [np.pi*blackbody(wave, T=self.temperature(p_))
               * (self.radius(p_)/self._parameters[0]*2.25459e-14)**2 
               for p_ in phase]
        
        return np.array(out)
    
class SpectralIndexSource(sncosmo.Source):
    """
    """
    def __init__(self, minphase=0.1, maxphase=30,
                 fluxfunc=(lambda p, t: t ** (-0.4 * p[0])),
                 specfunc=(lambda p, t: p[0]),
                 fluxparam=[1.], specparam=[1.],
                 name=None, version=None):
        self.name = name
        self.version = version

        self._minphase = minphase
        self._maxphase = maxphase

        self._fluxfunc = fluxfunc
        self._specfunc = specfunc

        self._n_fluxparam = len(fluxparam)
        self._n_specparam = len(specparam)

        self._parameters = np.array([1.])
        self._parameters = np.append(self._parameters, fluxparam)
        self._parameters = np.append(self._parameters, specparam)

        self._param_names = ['amplitude']
        self.param_names_latex = ['A']

        for k in range(self._n_fluxparam):
            self._param_names.append('f%i'%k)
            self.param_names_latex.append('f_%i'%k)

        for k in range(self._n_specparam):
            self._param_names.append('a%i'%k)
            self.param_names_latex.append('\alpha_%i'%k)


    def minwave(self):  
        return 1e-100

    def maxwave(self):  
        return 1e100
    
    def minphase(self):
        return self._minphase

    def maxphase(self):
        return self._maxphase

    def _flux(self, phase, wave):
        wave = np.array(wave)
        fluxparam = self._parameters[1:self._n_fluxparam+1]
        specparam = self._parameters[self._n_fluxparam+1:
                                     self._n_fluxparam+self._n_specparam+1]
        return np.array([((self._parameters[0]
                          * self._fluxfunc(fluxparam, phase)
                          * wave ** self._specfunc(specparam, phase))
                          if p_ >= self.minphase and p_ <= self.maxphase
                          else np.zeros(wave.shape))
                         for p_ in phase])

#######################################
#                                     #
# Auxilary functions                  #
#                                     #
#######################################
def blackbody(wl, T=6e3):
    # wl in Ångströms
    # T in Kelvin
    # R in solar radii
    # d in Mpc
    # output in erg s^-1 cm^-2 Angstrom^-1
    return 1.19104295262e+27 * wl**-5 / (np.exp(1.43877735383e+8 / (wl*T)) - 1)
    
