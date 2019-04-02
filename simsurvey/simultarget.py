#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains generators for simulated transients"""

import os
import warnings
import numpy as np
from numpy.random import uniform, normal
import sncosmo
import pickle

from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from scipy.stats import truncnorm
from astropy.cosmology import Planck15

from propobject import BaseObject

from .models import ExpandingBlackBodySource#, SpectralIndexSource, MultiSource

from .utils       import random
from .utils.tools import kwargs_extract, kwargs_update, range_args

_d2r = np.pi / 180

__all__ = ["get_transient_generator", "load_transient_generator",
           "generate_transients", "generate_lightcurves"]

def get_transient_generator(zrange,ratekind="basic",ratefunc=None,
                            ra_range=[0,360], dec_range=[-90,90],
                            ntransient=None,
                            **kwargs):
    """
    This model returns the object that enables to create and change
    the kind of transient you wish to set in the sky.

    # - HERE COPY PASTE THE TransientGenerator INIT - #
    # - TO BE DONE

    """
    return TransientGenerator(ratekind=ratekind,ratefunc=ratefunc,
                              ntransient=ntransient,zrange=zrange,
                              ra_range=ra_range,dec_range=dec_range,
                              **kwargs)

def load_transient_generator(filename):
    """
    """
    return TransientGenerator(load=filename)

def generate_transients(zrange,**kwargs):
    """
    This module calls get_transient_generator to create the
    TransientGenerator object and then returns the associated
    TransientGenerator.transients

    # - HERE COPY PASTE the transient_generator docstring
    """
    return get_transient_generator(zrange,**kwargs).transients

def generate_lightcurves(zrange, obs, trim_observations=True, **kwargs):
    """
    This module calls get_transient_generator to create the
    TransientGenerator object and then generates lightcurves based
    on the random transients 

    # - HERE COPY PASTE the transient_generator docstring
    """
    tr =  get_transient_generator(zrange,**kwargs)

    return tr.get_lightcurves(obs, trim_observations=trim_observations)


#######################################
#                                     #
# Generator: Any Transient            #
#                                     #
#######################################
class TransientGenerator( BaseObject ):
    """
    """
    __nature__ = "TransientGenerator"

    PROPERTIES         = ["transient_coverage",
                          "event_coverage"]
    SIDE_PROPERTIES    = ["sfd98_dir", "ratefunc", "model", "err_mwebv", "apply_mwebv"]
    DERIVED_PROPERTIES = ["simul_parameters", "mwebv", "mwebv_sfd98", 
                          "has_mwebv_sfd98", "lightcurve_parameters"]

    def __init__(self, empty=False, **kwargs):
        """
        """
        self.__build__()
        if empty:
            return

        self.create(**kwargs)

    def create(self, zrange=(0.0, 0.2), ratekind="basic", ratefunc=None,
               ntransient=None, transient=None, template=None, load=False,
               mjd_range=(57754.0,58849.0),
               ra_range=(0,360), dec_range=(-90,90), apply_mwebv=True,
               mw_exclusion=0, sfd98_dir=None, transientprop=None, err_mwebv=0.01):
        """
        """
        # == Add the Input Test == #
        #   TO BE DONE

        # *************** #
        # * Create      * #
        # *************** #
        # -- This will be directly used as random.radec inputs
        if transientprop is None:
            transientprop = {}

        self.set_sfd98_dir(sfd98_dir)
        self.set_apply_mwebv(apply_mwebv)
        
        if not load:
            self.set_event_parameters(update=False,
                                      **{"ra_range":ra_range, "dec_range":dec_range,
                                         "zcmb_range":zrange, "mjd_range":mjd_range,
                                         "mw_exclusion":mw_exclusion})

            self.set_transient_parameters(ratekind=ratekind, ratefunc=ratefunc,
                                          transient=transient, template=template,
                                          ntransient=ntransient,
                                          update=False, **transientprop)

            self.set_err_mwebv(err_mwebv)
            self._update_()
        else:
            self.load(load)
            self.set_transient_parameters(ratekind=None, ratefunc=ratefunc,
                                          type_=None, template=template,
                                          ntransient=None,
                                          update=False, **transientprop)
    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def load(self, filename):
        """
        """
        loaded = pickle.load(open(filename, 'rb'))

        for k, v in loaded['properties'].items():
            self._properties[k] = v
        for k, v in loaded['side_properties'].items():
            self._side_properties[k] = v
        for k, v in loaded['derived_properties'].items():
            self._derived_properties[k] = v

    def save(self, filename):
        """
        """
        prop_save = ["transient_coverage", "event_coverage"]
        side_save = ["err_mwebv", "apply_mwebv"]
        deri_save = ["simul_parameters", "mwebv", "mwebv_sfd98", 
                     "lightcurve_parameters", "has_mwebv_sfd98"]

        out = {
            "properties": {k: self._properties[k] for k in prop_save},
            "side_properties": {k: self._side_properties[k] for k in side_save},
            "derived_properties": {k: self._derived_properties[k] for k in deri_save},
        }
        out['properties']['transient_coverage']['lightcurve_prop'] = None

        pickle.dump(out, open(filename, 'wb'))

    # --------------------------- #
    # - Set Methods             - #
    # --------------------------- #
    def set_event_parameters(self,update=True,**kwargs):
        """
        Change the properties associated to the transient events.
        Known properties: "ra_range","dec_range","zcmb_range","mjd_range",
                          "mw_exclusion"

        Set update to True to update the derived properties
        """
        known_event_prop = ["ra_range","dec_range","zcmb_range",
                            "mw_exclusion","mjd_range"]

        for k in kwargs.keys():
            if k not in known_event_prop:
                raise ValueError("'%s' is not an 'event' parameter."%k +\
                                 " These are: "+", ".join(known_event_prop))

        self._properties["event_coverage"] = kwargs_update(self.event_coverage,
                                                           **kwargs)
        if update:
            self._update_()

    def set_transient_parameters(self,ratekind="basic", ratefunc=None,
                                 ntransient=None, update=True,
                                 transient=None,
                                 template=None, **kwargs):
        """
        This method will define the transient properties.
        """
        if self._properties["transient_coverage"] is None:
            self._properties["transient_coverage"] = {}

        # - you are good to fill it
        if transient is not None:
            self._properties["transient_coverage"]["transienttype"] = transient
        if template is not None:
            self._properties["transient_coverage"]["template"] = template
        if ratekind is not None:
            self._properties["transient_coverage"]["ratekind"] = ratekind
        if ntransient is not None:
            self._properties["transient_coverage"]["ntransient"] = ntransient

        # -- if this works, you are good to go
        f = RateGenerator().get_ratefunc(transient=self.transienttype,
                                         ratekind=self.ratekind,
                                         ratefunc=ratefunc)

        self._side_properties["ratefunction"] = f

        self.set_lightcurve_prop(**kwargs)

        # if "lcmodel" in kwargs.keys():
        #     self.set_model(kwargs["lcmodel"])
        #     self.set_lightcurve_prop(model=self.model,**kwargs)
        # elif "lcsource" in kwargs.keys():
        #     self.set_lightcurve_prop(source=kwargs["lcsource"],**kwargs)
        # else:
        #     self.set_lightcurve_prop(source=kwargs["lcsource"],**kwargs)

        if update:
            self._update_()

    def set_lightcurve_prop(self, lcmodel=None, lcmodel_prop=None,
                            lcsimul_func='basic', lcsimul_prop=None, **kwargs):
        """
        lcsimul_func must be function with redshift and sncosmo.model as arguments
        lcsimul_prop can be used for options of lcsimul_func
        """
        if lcmodel_prop is None:
            lcmodel_prop = {}
        if lcsimul_prop is None:
            lcsimul_prop = {}

        # TODO: Try to obtain lcsimul_func first
        if type(lcsimul_func) is str:
            lcsimul_func = LightCurveGenerator().get_lightcurve_func(
                transient=self.transienttype,
                template=self.template,
                simulation=lcsimul_func
            )

        props = {
            "param_func": lcsimul_func,
            "param_func_kwargs": lcsimul_prop
        }

        if lcmodel is None:
            self.set_model(
                LightCurveGenerator().get_model(
                    transient=self.transienttype,
                    template=self.template,
                    **lcmodel_prop
                )
            )
        else:
            self.set_model(lcmodel)

        self._properties["transient_coverage"]["lightcurve_prop"] = props

    # --------------------------- #
    # - Get Methods             - #
    # --------------------------- #
    def get_bandmag(self, band='bessellb', magsys='vega', t=0):
        """
        Returns the magnitudes of transient according to lightcurve parameters
        """
        # Save old params, so you can restore them
        param0 = {name: value for name, value
                  in zip(self.model.param_names,
                         self.model.parameters)}
        out = []
        for param in self.get_lightcurve_full_param():
            p = {k: param[k] for k in self.model.param_names if k != 'mwr_v'}
            self.model.set(**p)
            out.append(self.model.bandmag(band, magsys, p['t0'] + t))
        self.model.set(**param0)

        return np.array(out)

    def get_lightcurve_full_param(self, *args, **kwargs):
        """Transient lightcurve parameters"""
        full_out = kwargs.get("full_out", True)
        for i in range(*range_args(self.ntransient, *args)):
            out = dict(z=self.zcmb[i], t0=self.mjd[i],
                       **{p: v[i] for p, v in self.lightcurve.items()})

            if self.apply_mwebv:
                if self.has_mwebv_sfd98:
                    out["mwebv"] = self.mwebv[i]
                else:
                    out["mwebv"] = 0
                
            if full_out:
                out["ra"] = self.ra[i]
                out["dec"] = self.dec[i]
                if self.has_mwebv_sfd98:
                    out["mwebv_sfd98"] = self.mwebv_sfd98[i]
                else:
                    out["mwebv_sfd98"] = 0

            yield out

    def get_lightcurves(self, obs, trim_observations=True, **kwargs):
        """Realize lightcurves based on the randomized lightcurve parameters
        and a single set of observations"""
        params = self.get_lightcurve_full_param(full_out=False)
        return sncosmo.realize_lcs(obs, self.model, params,
                                   trim_observations=trim_observations, **kwargs)

    # --------------------------- #
    # - Plots Methods           - #
    # --------------------------- #
    def show_skycoverage(self, ax=None, savefile=None, show=True, cscale=None, 
                         cblabel=None, cmargin=5, mask=None, **kwargs):
        """This function enable to draw on the sky the position of the
        transients

        Parameters:
        -----------
        ax [matplotlib axis]       This axis must have Mollweide or Hammer projection

        - color options -

        cscale: [None/array]       array used to set a color scale to the markers

        cblabel: [string]          label of the colorbar (if any)

        cmargin: [float<100]       to avoid outlier issue in the colorbars, this
                                   enable to set vmin and vmax value, being the
                                   `cmargin`% and 100-`cmargin`% limits of the array
                                   (set 0 to effectively remove this option)

        vmin/vmax: [float]         set vmin/vmax of the colorbar manually;
                                   overrides results cmargin but if e.g. only
                                   vmin is given vmax from cmargin is still used

        mask: [None/bool. array]   mask for the scatter plot

        - output option -

        savefile, show [string, bool] Output options

        -- kwargs goes to skyscatter -> mpl's scatter --

        Returns:
        --------
        dict of the plot parameters

        """
        import matplotlib.pyplot as mpl
        from utils.mpladdon import figout, skyplot
        from utils.skyplot import ax_skyplot
        self._plot = {}

        # ------------------
        # - Color Scale
        if cscale == 'zcmb':
            c = np.asarray([t['zcmb'] for t in self.transientsources])
            if cblabel is None:
                cblabel = r"$\mathrm{Redshift}$"
        elif type(cscale) != str and hasattr(cscale, '__iter__'):
            c = cscale
        elif cscale is not None:
            raise ValueError('cscale must be array or predefined string, '+\
                             ' e.g. "redshift"')

        # ------------------
        # - Mask
        if mask is None:
            mask = np.ones(self.ntransient, dtype=bool)

        # ------------------
        # - Axis definition
        if ax is None:
            ax_default = dict(fig=None, figsize=(12, 6),
                              rect=[0.1, 0.1, 0.8, 0.8],
                              projection='mollweide',
                              xlabelpad=None,
                              xlabelmode='show')
            if cscale is not None:
                ax_default['figsize'] = (12,8)

            ax_kw, kwargs = kwargs_extract(ax_default, **kwargs)
            fig, ax = ax_skyplot(**ax_kw)
        elif ("MollweideTransform" not in dir(ax) and
              "HammerTransform" not in dir(ax)):
            raise TypeError("The given 'ax' most likely is not a matplotlib axis "+\
                            "with Mollweide or Hammer projection. Transform "+\
                            "function not found.")
        else:
            fig = ax.fig

        # ------------------
        # - Actual plotting
        if cscale is None:
            pl = ax.skyplot(self.ra[mask], self.dec[mask], **kwargs)
            cb = None
        else:
            # --------------------
            # - To avoid outliers
            vmin = kwargs.pop("vmin", np.percentile(c[mask], cmargin))
            vmax = kwargs.pop("vmax", np.percentile(c[mask], 100-cmargin))
            # - Da plo
            pl = ax.skyscatter(self.ra[mask], self.dec[mask], c=c[mask], 
                               vmin=vmin, vmax=vmax, **kwargs)
            cb = fig.colorbar(pl, orientation='horizontal', shrink=0.85, pad=0.08)
            if cblabel is not None:
                cb.set_label(cblabel, fontsize="x-large") 

        # ------------------- #
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["plot"]   = pl
        if cb is not None:
            self._plot["cbar"] = cb

        fig.figout(savefile=savefile,show=show)
        
        return self._plot

    def hist_skycoverage(self, ax=None, savefile=None, show=True, 
                         cblabel=r"$N_{SNe}$", **kwargs):
        """This function draws a sky histogram of the transient coverage"""
        import matplotlib.pyplot as mpl
        from utils.mpladdon import figout, skyhist
        from utils.skyplot import ax_skyplot
        self._plot = {}

        if ax is None:
            ax_default = dict(fig=None, figsize=(12, 6),
                              rect=[0.1, 0.1, 0.8, 0.8],
                              projection='mollweide',
                              xlabelpad=None,
                              xlabelmode='hist')
            ax_kw, kwargs = kwargs_extract(ax_default, **kwargs)
            fig, ax = ax_skyplot(**ax_kw)
        elif ("MollweideTransform" not in dir(ax) and
              "HammerTransform" not in dir(ax)):
            raise TypeError("The given 'ax' most likely is not a matplotlib axis "+\
                        "with Mollweide or Hammer projection. Transform "+\
                        "function not found.")
        else:
            fig = ax.fig

        collec, cb = ax.skyhist(self.ra, self.dec, cblabel=cblabel, **kwargs)
        cb.set_label(cblabel, fontsize="x-large") 

        # ------------------- #
        # -- Save the data -- #
        self._plot["figure"] = fig
        self._plot["ax"]     = ax
        self._plot["collection"] = collec
        self._plot["cbar"] = cb

        fig.figout(savefile=savefile,show=show)
        return self._plot 

    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _update_simulation_(self):
        """
        """
        # -----------------------
        # - Redshift from Rate
        if "ntransient" not in self.transient_coverage.keys():
            self.simul_parameters["zcmb"] = \
                list(sncosmo.zdist(self.zcmb_range[0], self.zcmb_range[1],
                                   time=self.timescale, area=self.coveredarea,
                                   ratefunc=self.ratefunc))
        else:
            self.simul_parameters["zcmb"] = \
                list(zdist_fixed_nsim(self.transient_coverage["ntransient"],
                                      self.zcmb_range[0], self.zcmb_range[1],
                                      ratefunc=self.ratefunc))
            
        self.simul_parameters["mjd"] = self._simulate_mjd_()
        self.simul_parameters["ra"], self.simul_parameters["dec"] = \
          random.radec(self.ntransient,
                       ra_range=self.ra_range,
                       dec_range=self.dec_range,
                       mw_exclusion=self._get_event_property_("mw_exclusion"))
        self._derived_properties['mwebv'] = None

        if "lightcurve_prop" in self.transient_coverage.keys():
            lc = self.transient_coverage["lightcurve_prop"]
            param = lc["param_func"](self.zcmb,self.model,
                                     **lc["param_func_kwargs"])
            self._derived_properties["simul_parameters"]["lightcurve"] = param 

    def _update_mwebv_sfd98_(self):
        try:
            import sfdmap
            self._derived_properties["mwebv_sfd98"] = sfdmap.ebv(
                self.ra, self.dec,
                mapdir=self._sfd98_dir
            )
        except ImportError:
            warnings.warn("sfdmap is not installed. "
                          "MW E(B-V) will be set to zero.")
        except IOError:
            warnings.warn("SFD98 dust map files not found. "
                          "MW E(B-V) will be set to zero.")

        if self._derived_properties["mwebv_sfd98"] is not None:
            self._derived_properties["has_mwebv_sfd98"] = True
        else:
            self._derived_properties["has_mwebv_sfd98"] = False

    def _update_mwebv_(self):
        if self.has_mwebv_sfd98:
            self._derived_properties["mwebv"] = self.mwebv_sfd98.copy()
            off = self.err_mwebv * np.random.randn(self.ntransient)
            self._derived_properties["mwebv"] += off
        
    def _update_(self):
        """This module create the derived values based on the
        fundamental ones"""
        # --------------
        # - update the actual simulation
        self._update_simulation_()
        self._update_mwebv_()

    def _reset_mwebv_(self):
        self._derived_properties['mwebv_sfd98'] = None
        self._derived_properties['mwebv'] = None
        
    def _simulate_mjd_(self):
        """
        Be default, this is a random flat time distribution returning a float
        per transient.
        Simple overwrite this function in a child-class to have more advanced
        properties.
        The Simulated mjd will be stored in transientsources.
        """
        return np.random.rand(self.ntransient)*self.timescale + self.mjd_range[0]

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    def _get_event_property_(self,key):
        """
        """
        if self.transient_coverage is None:
            raise AttributeError("'transient_coverage' has not been defined")
        
        return self.event_coverage[key]
    # -----------------
    # - Ranges
    @property
    def zcmb_range(self):
        """zcmb range used to draw transient"""
        return self._get_event_property_("zcmb_range")
    
    @property
    def ra_range(self):
        return self._get_event_property_("ra_range")
    
    @property
    def dec_range(self):
        return self._get_event_property_("dec_range")
    
    @property
    def mjd_range(self):
        """zcmb range used to draw transient"""
        return self._get_event_property_("mjd_range")
    # -----------------
    # - Rates
    @property
    def ratefunc(self):
        return self._side_properties["ratefunction"]

    # -------------------------------
    # - Derived Transient Properties
    @property
    def timescale(self):
        return self.mjd_range[1] - self.mjd_range[0]
       
    @property
    def ntransient(self):
        """number of transient requested"""
        return len(self.zcmb)

    @property
    def coveredarea(self):
        """Covered area in degree squared"""
        mw_exclusion = self._get_event_property_("mw_exclusion")
        
        # Area in steradians without accounting for MW exclusion
        area_sr = ((np.sin(self.dec_range[1] * _d2r) 
                    - np.sin(self.dec_range[0] * _d2r)) 
                   * (self.ra_range[1] - self.ra_range[0]) * _d2r)

        if mw_exclusion > 0:
            if self.ra_range == [-180, 180] and self.dec_range == [-90, 90]:
                area_sr -= 4 * np.pi * np.sin(mw_exclusion * _d2r) 
            else:
                # Make sure the warning is issued every time
                warnings.simplefilter('always', UserWarning)
                warnings.warn("MW exclusion was ignored when calculating covered area.")

        return area_sr / _d2r ** 2

    # --------------------
    # - Target Coverage
    @property
    def transient_coverage(self):
        if self._properties["transient_coverage"] is None:
            self._properties["transient_coverage"] = {}
        return self._properties["transient_coverage"]
        
    @property
    def event_coverage(self):
        """where and when the transient could be drawn"""
        if self._properties["event_coverage"] is None:
            self._properties["event_coverage"] = {}
            
        return self._properties["event_coverage"]

    # --------------------
    # - Derived Properties
    @property
    def simul_parameters(self):
        if self._derived_properties["simul_parameters"] is None:
            self._derived_properties["simul_parameters"] = {}
        return self._derived_properties["simul_parameters"]

    # -----------------
    # - transient info
    @property
    def zcmb(self):
        """Simulated zcmb based on the given rate function"""
        return self.simul_parameters["zcmb"]

    @property
    def mjd(self):
        """Loop over the transient sources to get the mjd"""
        return np.asarray(self.simul_parameters['mjd'])
    
    @property
    def ra(self):
        """Loop over the transient sources to get the ra"""
        return np.asarray(self.simul_parameters['ra'])
    
    @property
    def dec(self):
        """Loop over the transient sources to get the ra"""        
        return np.asarray(self.simul_parameters['dec'])

    # ------------------------
    # - derived transient info

    @property
    def mwebv(self):
        """Return MW E(B-V) 
        (if None or not up to date fetch from SFD98 map, and perturb)"""
        if self._derived_properties['mwebv'] is None:
            self._update_mwebv_()
        
        # if it is still None after update, some thing went wrong
        # likely map files were missing or in wrong directory
        if self._derived_properties['mwebv'] is None:
            return None
        else:
            return np.asarray(self._derived_properties['mwebv'])

    @property
    def has_mwebv_sfd98(self):
        if self._derived_properties['has_mwebv_sfd98'] is None:
            self._update_mwebv_sfd98_()
            
        return self._derived_properties['has_mwebv_sfd98']
    
    @property
    def mwebv_sfd98(self):
        """Return 'true' MW E(B-V) 
        (if None or not up to date fetch from SFD98 map)"""
        if self._derived_properties['mwebv_sfd98'] is None:
            self._update_mwebv_sfd98_()
        
        # if it is still None after update, some thing went wrong
        # likely map files were missing or in wrong directory
        if self._derived_properties['mwebv_sfd98'] is None:
            return None
        else:
            return np.asarray(self._derived_properties['mwebv_sfd98'])
    
    # ------------------
    # - Side properties

    @property
    def _sfd98_dir(self):
        """Directory where the dust maps are located"""
        return self._side_properties["sfd98_dir"]

    def set_sfd98_dir(self, value):
        """Clears fetched values for MW E(B-V)"""
        self._side_properties['sfd98_dir'] = value
        self._reset_mwebv_()

    @property
    def model(self):
        """Light curve model (derived from source if not set)"""
        
        return self._side_properties["model"]

    def set_model(self, model):
        """
        Set the transient model.
        If it does not have MW dust effect, the effect is added.
        """
        # if model.__class__ is not sncosmo.models.Model:
        #     raise TypeError("model must be sncosmo.model.Model")

        if "mwebv" not in model.param_names and self.apply_mwebv:
            model.add_effect(sncosmo.CCM89Dust(), 'mw', 'obs')

        self._side_properties["model"] = model

    @property
    def apply_mwebv(self):
        """Option to switch of MW dust completely (need at wavelengths > 3.3 micron)"""
        return self._side_properties["apply_mwebv"]

    def set_apply_mwebv(self, apply_mwebv):
        """
        """
        self._side_properties["apply_mwebv"] = apply_mwebv
        
    @property
    def err_mwebv(self):
        """Assumed error of dustmap; will be applied to lightcurve creation"""
        return self._properties["err_mwebv"]

    def set_err_mwebv(self, err):
        """
        """
        self._properties['err_mwebv'] = err

    @property
    def transienttype(self):
        """
        """
        return (self._properties["transient_coverage"]["transienttype"]
                if "transienttype"
                in self._properties["transient_coverage"].keys()
                else None)

    @property
    def ratekind(self):
        """
        """
        return (self._properties["transient_coverage"]["ratekind"]
                if "ratekind"
                in self._properties["transient_coverage"].keys()
                else None)

    @property
    def template(self):
        """
        """
        return (self._properties["transient_coverage"]["template"]
                if "template"
                in self._properties["transient_coverage"].keys()
                else None)

    # -----------------------
    # - LightCuve Properties
    @property
    def lightcurve(self):
        """Check if there is a lightcurve simul_parameter defined"""
        return None if "lightcurve" not in self.simul_parameters.keys()\
          else self.simul_parameters['lightcurve']
          
    @property
    def lightcurve_param_names(self):
        """ """
        return self.lightcurve.keys()

    def has_lightcurves(self):
        return False if self.lightcurve is None \
          else True

    @property
    def lightcurve_properties(self):
        """
        """
        return self.transient_coverage["lightcurve_prop"]
        
#######################################
#                                     #
# Generator: Transient Rate           #
#                                     #
#######################################
class _PropertyGenerator_(BaseObject):
    """
    """
    __nature__ = "PropertyGenerator"

    # ========================== #
    # = Internal               = #
    # ========================== #
    def _parse_rate_(self,key="rate"):
        return [m.split(key+"_")[-1] for m in dir(self)
                if m.startswith(key)]
    
class RateGenerator( _PropertyGenerator_ ):
    """
    This follows the SN cosmo ratefunc
    
    ratefunc : callable
        A callable that accepts a single float (redshift) and returns the
        comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.
        The default is a function that returns ``1.e-4``.
    
    """
    __nature__ = "RateGenerator"
    #
    # User: Add any rate_TTT_BLA(self,z) to define a new rate for the TTT transient.
    #       where this return the comoving volumetric rate at each redshift in
    #        units of yr^-1 Mpc^-3.
    # 
    #       This rate will be accessible from
    #       get_ratefunc(transient="TTT",ratekind="BLA")
    #

    def get_ratefunc(self, transient="Ia", ratekind="basic", ratefunc=None):
        """
        Parameters
        ----------

        Return
        ------
        function(z)
        """
        # ---------------
        # - Rate Kinds
        if ratefunc is not None:
            return ratefunc

        name = '%s_%s'%(transient, ratekind)
        if name in self.known_rates:
            return eval("self.rate_%s"%(name))
        else:
            raise ValueError("no '%s' rate found"%(name)+\
                             "Available rates: "+",".join(self.known_rates))

    # ========================== #
    # = Rates                  = #
    # ========================== #
    def rate_basic(self,z):
        """
        Basic default rate function in sncosmo: returns ``1.e-4``.
        (comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.)
        """
        return 1.e-4

    # ----------------- #
    # - Ia rate       - #
    # ----------------- #
    def rate_Ia_basic(self,z):
        """
        More realistic value for low-z than sncosmo default
        (comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.)
        """
        return 3.e-5 * (1 + z)

    
    def rate_Ia_basiclow(self,z):
        """
        """
        return self.rate_basic(z) / 10

    # ----------------- #
    # - CC rates       - #
    # ----------------- #
    def rate_Ibc_basic(self,z):
        """
        [TODO: Add source]
        (comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.)
        """
        return 2.25e-5 * (1 + z)

    def rate_IIn_basic(self,z):
        """
        [TODO: Add source]
        (comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.)
        """
        return 7.5e-6 * (1 + z)

    def rate_IIP_basic(self,z):
        """
        [TODO: Add source]
        (comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.)
        """
        return 1.2e-4 * (1 + z)

    # ========================== #
    # = *** Rates              = #
    # ========================== #

    
    # ========================== #
    # = Properties             = #
    # ========================== #


    @property
    def known_rates(self):
        return self._parse_rate_("rate")
    
    @property
    def known_Ia_rates(self):
        return self._parse_rate_("rate_Ia")
    
#######################################
#                                     #
# Generator: Light Curve              #
#                                     #
#######################################
class LightCurveGenerator( _PropertyGenerator_ ):
    """
    """

    # ========================== #
    # = Method                 = #
    # ========================== #
    def get_lightcurve_func(self, transient="Ia", template="hisao", simulation="basic"):
        """
        Parameters
        ----------
        Return
        ------
        a lightcurve generator function
        """
        # ---------------
        # - Transients
        # if transient is None or transient == "":
        #     available_lc = self.known_lightcurve_simulations
        #     if len(avialable_lc) == 0:
        #         raise NotImplementedError("no lightcurve simulation implemented")
            
        #     transient = None            
        # # elif transient == "Ia":
        # #     avialable_lc = self.known_Ia_lightcurve_simulation
        # else:
        #     raise ValueError("'%s' is not a known transient"%transient)

        name = "%s_%s_%s"%(transient, template, simulation)
        # ---------------
        # - Rate Kinds   
        if name in self.known_lightcurve_simulations:
            return eval("self.lightcurve_%s"%name)
        else:
            raise ValueError("No lightcurve randomizer available for '%s'"%name)


    def get_model(self, transient="Ia",  template="hisao", **kwargs):
        """
        """
        name = "%s_%s"%(transient, template)
        if name in self.known_models:
            return eval("self.model_%s(**kwargs)"%name)
        else:
            raise ValueError("No lightcurve model available for '%s'"%name)

    def set_model(self, model):
        """
        """
        self._side_properties["model"] = model

    # ============================ #
    # = Model definitions        = #
    # ============================ #

    def model_Ia_salt2(self):
        """
        """
        return sncosmo.Model(source='salt2')

    def model_Ia_hsiao(self):
        """
        """
        sncosmo.get_source('hsiao', version='3.0')
        p, w, f = sncosmo.read_griddata_fits(
            os.path.join(sncosmo.builtins.get_cache_dir(),
                         'sncosmo/models/hsiao/Hsiao_SED_V3.fits')
        )

        return sncosmo.Model(
            source=sncosmo.StretchSource(p, w, f, name='hsiao-stretch'),
            effects=[sncosmo.CCM89Dust()],
            effect_names=['host'],
            effect_frames=['rest']
        )

    def model_Ibc_nugent(self):
        """
        """
        return sncosmo.Model(source='nugent-sn1bc',
                             effects=[sncosmo.CCM89Dust()],
                             effect_names=['host'],
                             effect_frames=['rest'])

    def model_Ibc_snana(self):
        """
        """
        filenames = get_snana_filenames('Ib')
        filenames.append(get_snana_filenames('Ic'))

        return self.model_generic_MultiSource(files=filenames)

    def model_IIn_nugent(self):
        """
        """
        sncosmo.get_source('nugent-sn2n', version='2.1')
        p, w, f = sncosmo.read_griddata_ascii(
            os.path.join(sncosmo.builtins.get_cache_dir(),
                         'sncosmo/models/nugent/sn2n_flux.v2.1.dat')
        )
        mask = (p < 150)

        return sncosmo.Model(source=sncosmo.TimeSeriesSource(p[mask], w, f[mask]),
                             effects=[sncosmo.CCM89Dust()],
                             effect_names=['host'],
                             effect_frames=['rest'])

    def model_IIn_snana(self):
        """
        """
        filenames = get_snana_filenames('IIn')

        return self.model_generic_MultiSource(files=filenames)

    def model_IIP_nugent(self):
        """
        """
        sncosmo.get_source('nugent-sn2p', version='1.2')
        p, w, f = sncosmo.read_griddata_ascii(
            os.path.join(sncosmo.builtins.get_cache_dir(),
                         'sncosmo/models/nugent/sn2p_flux.v1.2.dat')
        )
        mask = (p < 130)

        return sncosmo.Model(source=sncosmo.TimeSeriesSource(p[mask], w, f[mask]),
                             effects=[sncosmo.CCM89Dust()],
                             effect_names=['host'],
                             effect_frames=['rest'])

    def model_IIP_snana(self):
        """
        """
        filenames = get_snana_filenames('IIP')

        return self.model_generic_MultiSource(files=filenames)

    def model_generic_ExpandingBlackBody(self, **kwargs):
        """
        """
        return sncosmo.Model(source=ExpandingBlackBodySource(**kwargs),
                             effects=[sncosmo.CCM89Dust()],
                             effect_names=['host'],
                             effect_frames=['rest'])

    def model_generic_SpectralIndex(self, **kwargs):
        """
        """
        return sncosmo.Model(source=SpectralIndexSource(**kwargs),
                             effects=[sncosmo.CCM89Dust()],
                             effect_names=['host'],
                             effect_frames=['rest'])

    def model_generic_MultiSource(self, files=None, **kwargs):
        """
        """
        data_sed = {'p': [], 'w': [], 'f': []}
        for sed_file in sed_files[0]:
            p_, w_, f_ = sncosmo.read_griddata_ascii(sed_file)
            data_sed['p'].append(p_)
            data_sed['w'].append(w_)
            data_sed['f'].append(f_)

        source = simsurvey.MultiSource(data_sed['p'], data_sed['w'], data_sed['f'])
        return sncosmo.Model(source=source,
                             effects=[sncosmo.CCM89Dust()],
                             effect_names=['host'],
                             effect_frames=['rest'])

    # ============================ #
    # = LC Parameter randomizers = #
    # ============================ #

    # ----------------- #
    # - Ia LC         - #
    # ----------------- #

    # Functions to be implemented
    # Lightcurve parameter randomizers:
    # - Ibc snana basic
    # - IIn snana basic (including mag distro trunction)
    # - IIP snana basic
    # - SpectralIndexSource

    def lightcurve_Ia_salt2_basic(self, redshifts, model,
                                  color_mean=0, color_sigma=0.1,
                                  stretch_mean=0, stretch_sigma=1,
                                  alpha=0.13, beta=3.,
                                  cosmo=Planck15):
        """
        """
        ntransient = len(redshifts)
            
        x1 = normal(stretch_mean, stretch_sigma, ntransient)
        c = normal(color_mean, color_sigma, ntransient)
            
        x0 = []
        for z, x1_, c_ in zip(redshifts, x1, c):
            model.set(z=z, x1=x1_, c=c_)
            mabs = normal(-19.3, 0.1)
            mabs -= alpha*x1_ - beta*c_            
            model.set_source_peakabsmag(mabs, 'bessellb', 'ab', cosmo=cosmo)
            x0.append(model.get('x0'))
        
        return {'x0': np.array(x0),
                'x1': x1,
                'c': c}
    
    def lightcurve_Ia_salt2_realistic(self, redshifts, model,
                                      mabs_mean=-19.3, mabs_sigma=0.12,
                                      color_mean=0, color_sigma=0.1,
                                      stretch_mean=(0.5, -1), stretch_sigma=(1, 1),
                                      stretch_thr=0.75, alpha=0, beta=3,
                                      cosmo=Planck15):
        """
        stretch parameters assume bimodal distribution
        stretch_thr is a threshold for uniformly drawn number used to determine 
        in which mode the SNe is  
        """
        ntransient = len(redshifts)

        c = normal(color_mean, color_sigma, ntransient)
        x1 = normal(0,1,ntransient)
        
        # Shift x1 values to bimodal distribution
        x1_thr = uniform(size=ntransient)
        x1[x1_thr <= 0.75] *= stretch_sigma[0]
        x1[x1_thr <= 0.75] += stretch_mean[0]
        x1[x1_thr > 0.75] *= stretch_sigma[1]
        x1[x1_thr > 0.75] += stretch_mean[1]

        x0 = []
        for k,z in enumerate(redshifts):
            model.set(z=z, x1=x1[k], c=c[k])
            mabs = normal(mabs_mean, mabs_sigma)
            mabs -= alpha*x1[k] - beta*c[k]
            model.set_source_peakabsmag(mabs, 'bessellb', 'ab', cosmo=cosmo)
            x0.append(model.get('x0'))
            
        ntransient = len(redshifts)
        
        return {'x0': np.array(x0),
                'x1': np.array(x1),
                'c': np.array(c)}

    def lightcurve_Ia_hsiao_basic(self, redshifts, model,
                                  mag=(-19.3, 0.1),
                                  r_v=2., ebv_rate=0.11,
                                  alpha=1.3,
                                  **kwargs):
        """
        """
        out = lightcurve_scaled_to_mag(redshifts, model, mag=mag,
                                       r_v=r_v, ebv_rate=ebv_rate, **kwargs)

        out['s'] = np.random.normal(1., 0.1, len(redshifts))
        out['amplitude'] *= 10 ** (0.4 * alpha * (out['s'] - 1))

        return out

    # ----------------- #
    # - CC LC         - #
    # ----------------- #
    def lightcurve_Ibc_nugent_basic(self, redshifts, model,
                                    mag=(-17.5, 1.2),
                                    r_v=2., ebv_rate=0.11,
                                    **kwargs):
        """
        """
        return lightcurve_scaled_to_mag(redshifts, model, mag=mag,
                                        r_v=r_v, ebv_rate=ebv_rate, **kwargs)

    def lightcurve_Ibc_snana_basic(self, redshifts, model,
                                   mag=(-17.5, 1.2),
                                   r_v=2., ebv_rate=0.11,
                                   **kwargs):
        """
        """
        return self.lightcurve_generic_MultiSource_basic(redshifts, model,
                                                         mag=mag, r_v=r_v,
                                                         ebv_rate=ebv_rate,
                                                         **kwargs)

    def lightcurve_IIn_nugent_basic(self, redshifts, model,
                                    mag=(-18.5, 1.4),
                                    mag_dist_trunc=(-1, 1e6),
                                    r_v=2., ebv_rate=0.11,
                                    **kwargs):
        """
        """
        return lightcurve_scaled_to_mag(redshifts, model, mag=mag,
                                        mag_dist_trunc=mag_dist_trunc,
                                        r_v=r_v, ebv_rate=ebv_rate, **kwargs)

    def lightcurve_IIn_snana_basic(self, redshifts, model,
                                   mag=(-18.5, 1.4),
                                   mag_dist_trunc=(-1, 1e6),
                                   r_v=2., ebv_rate=0.11,
                                   **kwargs):
        """
        """
        return self.lightcurve_generic_MultiSource_basic(redshifts, model,
                                                         mag=mag, r_v=r_v,
                                                         ebv_rate=ebv_rate,
                                                         **kwargs)

    def lightcurve_IIP_nugent_basic(self, redshifts, model,
                                    mag=(-16.75, 1.),
                                    r_v=2., ebv_rate=0.11,
                                    **kwargs):
        """
        """
        return lightcurve_scaled_to_mag(redshifts, model, mag=mag,
                                        r_v=r_v, ebv_rate=ebv_rate, **kwargs)
    def lightcurve_IIP_snana_basic(self, redshifts, model,
                                   mag=(-16.75, 1.),
                                   r_v=2., ebv_rate=0.11,
                                   **kwargs):
        """
        """
        return self.lightcurve_generic_MultiSource_basic(redshifts, model,
                                                         mag=mag, r_v=r_v,
                                                         ebv_rate=ebv_rate,
                                                         **kwargs)

    def lightcurve_generic_ExpandingBlackBody_basic(self, redshifts, model,
                                                    sig_mag=0.1,
                                                    r_v=2., ebv_rate=0.11,
                                                    cosmo=Planck15,
                                                    **kwargs):
        """
        """
        return {
            'd': (cosmo.luminosity_distance(redshifts).value
                  * 10**(-0.4*np.random.normal(0, sig_mag))),
            'hostr_v': r_v * np.ones(len(redshifts)),
            'hostebv': np.random.exponential(ebv_rate, len(redshifts))
        }

    def lightcurve_generic_SpectralIndex_basic(self, redshifts, model, **kwargs):
        """
        """
        return lightcurve_scaled_to_mag(redshifts, model, **kwargs)

    def lightcurve_generic_MultiSource_basic(self, redshifts, model, **kwargs):
        """
        """
        return lightcurve_scaled_to_mag(redshifts, model,
                                        n_templates=len(model._source._model_flux),
                                        **kwargs)

    def lightcurve_Ia_salt2_hostdependent():
        raise NotImplementedError("To be done")

    # ========================== #
    # = Properties             = #
    # ========================== #
    @property
    def model(self):
        return self._side_properties['model']
    
    @property
    def known_lightcurve_simulations(self):
        return self._parse_rate_("lightcurve")
    
    @property
    def known_models(self):
        return self._parse_rate_("model")
    

#######################################
#                                     #
# Utilities                           #
#                                     #
#######################################
def lightcurve_scaled_to_mag(redshifts, model,
                             mag=(-18, 0.1),
                             mag_dist_trunc=None,
                             r_v=2., ebv_rate=0.11,
                             t_scale=None, cosmo=Planck15,
                             n_templates=1,
                             **kwargs):
    """
    """
    out = {}
    if n_templates > 1:
        out['template_index'] = np.random.randint(0, n_temp, len(redshifts))
        model_kw = [{'template_index': k} for k in out['template_index']]
    elif n_templates == 1:
        model_kw = [{} for z in redshifts]

    # Amplitude
    amp = []
    for z, model_kw_ in zip(redshifts, model_kw):
        if mag_dist_trunc is None:
            mabs = np.random.normal(mag[0], mag[1])
        else:
            mabs = truncnorm.rvs(mag_dist_trunc[0],
                                 mag_dist_trunc[1],
                                 mag[0],
                                 mag[1])

        if t_scale is None:
            model.set(z=z, **model_kw_)
            model.set_source_peakabsmag(mabs, 'bessellb', 'vega', cosmo=cosmo)
            amp.append(model.get('amplitude'))
        else:
            model.set(z=0, amplitude=1, **model_kw_)
            mag_abs = np.random.normal(mag[0], mag[1])
            mag_current = model.bandmag('sdssr', 'ab', 1)
            dm = mag_current - mag_abs
            amp.append(10**(0.4*(dm-cosmo.distmod(z).value)))

    out['amplitude'] = np.array(amp)
    out['hostr_v'] = r_v * np.ones(len(redshifts))
    out['hostebv'] =  np.random.exponential(ebv_rate, len(redshifts))

    return out


def zdist_fixed_nsim(nsim, zmin, zmax, 
                     ratefunc=lambda z: 1.,
                     cosmo=Planck15):
    """Generate a distribution of redshifts.

    Generates redshifts for a given number of tranisents with the correct
    redshift distribution, given the input volumetric SN rate function and
    cosmology.

    (Adapted from sncosmo.zdist)
    
    Parameters
    ----------
    nsim : int
        Number of transient redshifts to be simulated. 
    zmin, zmax : float
        Minimum and maximum redshift.
    ratefunc : callable
        A callable that accepts a single float (redshift) and returns the
        comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.
        The default is a function that returns ``1.``.
    cosmo : `~astropy.cosmology.Cosmology`, optional
        Cosmology used to determine volume. The default is Planck15.

    Examples
    --------

    TBA
    """

    # Get comoving volume in each redshift shell.
    z_bins = 100  # Good enough for now.
    z_binedges = np.linspace(zmin, zmax, z_bins + 1)
    z_binctrs = 0.5 * (z_binedges[1:] + z_binedges[:-1])
    sphere_vols = cosmo.comoving_volume(z_binedges).value
    shell_vols = sphere_vols[1:] - sphere_vols[:-1]

    # SN / (observer year) in shell
    shell_snrate = np.array([shell_vols[i] *
                             ratefunc(z_binctrs[i]) / (1.+z_binctrs[i])
                             for i in range(z_bins)])

    # SN / (observer year) within z_binedges
    vol_snrate = np.zeros_like(z_binedges)
    vol_snrate[1:] = np.add.accumulate(shell_snrate)

    # Create a ppf (inverse cdf). We'll use this later to get
    # a random SN redshift from the distribution.
    snrate_cdf = vol_snrate / vol_snrate[-1]
    snrate_ppf = Spline1d(snrate_cdf, z_binedges, k=1)

    for i in range(nsim):
        yield float(snrate_ppf(uniform()))

def get_snana_filenames(sntype):
    """
    """
    if not sntype.startswith('SN '):
        sntype = 'SN %s'%sntype

    reg = sncosmo.registry._get_registry(sncosmo.Source)
    source_tuples = [(v['name'], v['version'])
                     for v in reg.get_loaders_metadata()
                     if v['name'].startswith('snana')
                     and v['type'] == sntype]

    filenames = []
    for name, version in source_tuples:
        filepath = os.path.join(sncosmo.builtins.get_cache_dir(),
                                'sncosmo',
                                reg._loaders[(name, version)][1])
        if not os.exists(filepath):
            sncosmo.get_source(name, version=version)
        filenames.append(filepath)

    return filenames
