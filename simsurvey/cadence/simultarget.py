#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This pupose of this method is the have a generator of fake astrotarget, like SN"""

import warnings
import numpy as np
from numpy.random import uniform, normal
import sncosmo

from astrobject                      import BaseObject
from astrobject                      import get_target
from astrobject.utils.tools          import kwargs_extract,kwargs_update
from astrobject.utils                import random

_d2r = np.pi / 180

__all__ = ["get_sn_generator","get_transient_generator","generate_transients"]


def get_sn_generator(zrange,ratekind="basic",**kwargs):
    """
    This function enables to load the generator of Type Ia supernovae
    """
    return SNIaGenerator(ratekind=ratekind,zrange=zrange,
                         **kwargs)

def get_transient_generator(zrange,ratekind="basic",ratefunc=None,
                        ra_range=[-180,180], dec_range=[-90,90],
                        ntransients=None,
                        **kwargs):
    """
    This model returns the object that enables to create and change
    the kind of transient you wish to set in the sky.

    # - HERE COPY PASTE THE TransientGenerator INIT - #
    # - TO BE DONE
    
    """
    return TransientGenerator(ratekind=ratekind,ratefunc=ratefunc,
                              ntransients=ntransients,zrange=zrange,
                              ra_range=ra_range,dec_range=dec_range,
                              **kwargs)
    
def generate_transients(zrange,**kwargs):
    """
    This module calls transient_generator to create the
    TransientGenerator object and then returns the associated
    TransientGenerator.transients

    # - HERE COPY PASTE the transient_generator docstring
    """
    return transient_generator(zrange,**kwargs).transients


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
    SIDE_PROPERTIES    = ["sfd98_dir", "ratefunc", "model", "err_mwebv"]
    DERIVED_PROPERTIES = ["simul_parameters", "mwebmv", "mwebmv_sfd98", 
                          "lightcurve_parameters"]

    def __init__(self,zrange=[0.0,0.2], ratekind="basic", ratefunc=None,# How deep
                 mjd_range=[57754.0,58849.0],
                 ra_range=(-180,180),dec_range=(-90,90), # Where, see also kwargs
                 ntransients=None,empty=False,sfd98_dir=None,**kwargs):
        """
        """
        self.__build__()
        if empty:
            return
        
        self.create(zrange,
                    ratekind=ratekind, ratefunc=ratefunc, 
                    ntransients=ntransients,
                    ra_range=ra_range, dec_range=dec_range,
                    mjd_range=mjd_range,sfd98_dir=sfd98_dir,
                    **kwargs)

    def create(self,zrange,ratekind="basic",ratefunc=None,
               ntransients=None,type_=None,
               mjd_range=[57754.0,58849.0],
               ra_range=(-180,180),dec_range=(-90,90),
               mw_exclusion=0,sfd98_dir=None,transientprop={},err_mwebv=0.01):
        """
        """
        # == Add the Input Test == #
        #   TO BE DONE
        
        # *************** #
        # * Create      * #
        # *************** #
        # -- This will be directly used as random.radec inputs
        self.set_event_parameters(update=False,
                                  **{"ra_range":ra_range,"dec_range":dec_range,
                                   "zcmb_range":zrange,"mjd_range":mjd_range,
                                   "mw_exclusion":mw_exclusion})
        
        self.set_transient_parameters(ratekind=ratekind,ratefunc=ratefunc,
                                      type_=type_,
                                      update=False,**transientprop)
        
        self.set_err_mwebv(err_mwebv)

        self._update_()
        
    # =========================== #
    # = Main Methods            = #
    # =========================== #
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

    def set_transient_parameters(self,ratekind="basic",ratefunc=None,
                                 update=True,type_=None,**kwargs):
        """
        This method will define the transient properties.
        """
        if self._properties["transient_coverage"] is None:
            self._properties["transient_coverage"] = {}
            
        # -- if this works, you are good to go
        f = RateGenerator().get_ratefunc(transient=type_,ratekind=ratekind,
                                         ratefunc=ratefunc)
        
        # - you are good to fill it
        self._properties["transient_coverage"]["transienttype"] = type_
        self._properties["transient_coverage"]["ratekind"] = ratekind
        self._side_properties["ratefunction"] = f
        
        if "lcmodel" in kwargs.keys():
            self.set_model(kwargs["lcmodel"])
            self.set_lightcurve_prop(model=self.model,**kwargs)
        elif "lcsource" in kwargs.keys():
            self.set_lightcurve_prop(source=kwargs["lcsource"],**kwargs)

        if update:
            self._update_()
            
    def set_lightcurve_prop(self,source=None, model=None, lcsimul_func=None,
                            lcsimul_prop={}, **kwargs):
        """
        lcsimul_func must be function with redshift and sncosmo.model as arguments
        lcsimul_prop can be used for options of lcsimul_func
        """
        if lcsimul_func is None:
            raise ValueError("please set lcsimul_func")

        props = { 
            "param_func": lcsimul_func,
            "param_func_kwargs": lcsimul_prop
        }

        if model is not None:
            props["source"] = None
            props["model"] = model
        elif source is not None:
            accepted = [str, sncosmo.models.SALT2Source, 
                        sncosmo.models.TimeSeriesSource]
            if source.__class__ not in accepted:
                raise TypeError("source must be string, " + 
                                "sncosmo.models.TimeSeriesSource or" +
                                "sncosmo.models.SALT2Source")    
            props["source"] = source
            props["model"] = None
        else:
            raise ValueError("Please provide a model or a source")

        self._properties["transient_coverage"]["lightcurve_prop"] = props

    # --------------------------- #
    # - Get Methods             - #
    # --------------------------- #
    def get_transients(self,index=None,pass_mwebmv=True):
        """loops over the transientsources to load the transients objects.
        This method could be a bit slow..."""
        return [get_target(**s) for s in self.get_transientsource(index, pass_mwebmv)]
    
    def get_transientsource(self,index=None,pass_mwebmv=True):
        """dictionary containing the fundamental parameters that enable to
        load the transient objects"""
        if index is not None and "__iter__" not in dir(index):
            index = [index]

        # If self.mwebmv remains None, something went wrong and
        # pass_mwebmb is set to None. This should already have
        # resulted in a warning, no need for another one. 
        if pass_mwebmv and self.mwebmv is None:
            pass_mwebmv = False

        return [dict(name="simul%d"%i,ra=self.ra[i],dec=self.dec[i], zcmb=self.zcmb[i],
                     mjd=self.mjd[i],type_=self.transient_coverage["transienttype"],
                     lightcurve=None, 
                     forced_mwebmv=(self.mwebmv[i] if pass_mwebmv else None))
                for i in xrange(self.ntransient) if index is None or i in index]

    def get_bandmag(self, band='bessellb', magsys='vega', t=0):
        """
        Returns the magnitudes of transient according to lightcurve parameters
        """
        # Save old params, so you can restore them 
        param0 = {name: value for name, value 
                  in zip(self.model.param_names,
                         self.model.parameters)}
        out = []
        for param in self.lightcurve_full_param:
            self.model.set(**param)
            out.append(self.model.bandmag(band, magsys, param['t0'] + t))
        self.model.set(**param0)

        return np.array(out)

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
        from astrobject.utils.mpladdon import figout, skyplot
        from astrobject.utils.plot.skyplot import ax_skyplot
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
        from astrobject.utils.mpladdon import figout, skyhist
        from astrobject.utils.plot.skyplot import ax_skyplot
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
        self.simul_parameters["zcmb"] = \
          list(sncosmo.zdist(self.zcmb_range[0], self.zcmb_range[1],
                             time=self.timescale, area=self.coveredarea,
                             ratefunc=self.ratefunc))

        self.simul_parameters["mjd"] = self._simulate_mjd_()
        self.simul_parameters["ra"], self.simul_parameters["dec"] = \
          random.radec(self.ntransient,
                       ra_range=self.ra_range,
                       dec_range=self.dec_range,
                       mw_exclusion=self._get_event_property_("mw_exclusion"))
        self._derived_properties['mwebmv'] = None

        if "lightcurve_prop" in self.transient_coverage.keys():
            lc = self.transient_coverage["lightcurve_prop"]
            param = lc["param_func"](self.zcmb,self.model,
                                     **lc["param_func_kwargs"])
            self._derived_properties["simul_parameters"]["lightcurve"] = param 

    def _update_mwebmv_sfd98_(self):
        try:
            m = sncosmo.SFD98Map(mapdir=self._sfd98_dir)
            self._derived_properties["mwebmv_sfd98"] = m.get_ebv((self.ra, self.dec))
        except IOError:
            warnings.warn("MW E(B-V) could not be fetched. Please set sfd98_dir to the map directory.")

    def _update_mwebmv_(self):
        self._derived_properties["mwebmv"] = self.mwebmv_sfd98.copy()
        off = self.err_mwebmv * np.random.randn(self.ntransient)
        self._derived_properties["mwebmv"] += off
                 
    def _update_(self):
        """This module create the derived values based on the
        fundamental ones"""
        # --------------
        # - update the actual simulation
        self._update_simulation_()        

        
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
    def mwebmv(self):
        """Return MW E(B-V) 
        (if None or not up to date fetch from SFD98 map, and perturb)"""
        if self._derived_properties['mwebmv'] is None:
            self._update_mwebmv_()
        
        # if it is still None after update, some thing went wrong
        # likely map files were missing or in wrong directory
        if self._derived_properties['mwebmv'] is None:
            return None
        else:
            return np.asarray(self._derived_properties['mwebmv'])

    @property
    def mwebmv_sfd98(self):
        """Return 'true' MW E(B-V) 
        (if None or not up to date fetch from SFD98 map)"""
        if self._derived_properties['mwebmv_sfd98'] is None:
            self._update_mwebmv_sfd98_()
        
        # if it is still None after update, some thing went wrong
        # likely map files were missing or in wrong directory
        if self._derived_properties['mwebmv_sfd98'] is None:
            return None
        else:
            return np.asarray(self._derived_properties['mwebmv_sfd98'])
    
    # ------------------
    # - Side properties

    @property
    def _sfd98_dir(self):
        """Director where the maps are. Default option set"""
        if self._side_properties["sfd98_dir"] is None:
            from astrobject.utils.io import get_default_sfd98_dir
            self._side_properties["sfd98_dir"] = get_default_sfd98_dir(download_if_needed=True)
    
        return self._side_properties["sfd98_dir"]

    def set_sfd98_dir(self, value):
        self._side_properties['sfd98_dir'] = value
        self._update_mwebmv_()

    @property
    def model(self):
        """Light curve model (derived from source if not set)"""
        if self._side_properties["model"] is None:
            self.set_model(sncosmo.Model(source=self.lightcurve_source))
        
        return self._side_properties["model"]

    def set_model(self, model):
        """
        Set the transient model.
        If it does not have MW dust effect, the effect is added.
        """
        if model.__class__ is not sncosmo.models.Model:
            raise TypeError("model must be sncosmo.model.Model")

        if "mwebv" not in model.param_names:
            model.add_effect(sncosmo.CCM89Dust(), 'mw', 'obs')

        self._side_properties["model"] = model

    def reset_model(self):
        """
        Resets model to None. 
        Next time self.model is called it will be rederived from 
        self.lightcurve_source
        """
        self._side_properties["model"] = None

    @property
    def err_mwebv(self):
        """Assumed error of dustmap; will be applied to lightcurve creation"""
        return self._properties["err_mwebv"]

    def set_err_mwebv(self, err):
        """
        """
        self._properties['err_mwebv'] = err

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
        
    @property
    def lightcurve_source(self):
        if "lightcurve_prop" not in self.transient_coverage:
            raise AttributeError("no 'lightcurve_prop' defined")
        return self.lightcurve_properties["source"]

    @property
    def lightcurve_full_param(self):
        """Transient lightcurve parameters"""

        return [dict(z=self.zcmb[i], t0=self.mjd[i], 
                     ra=self.ra[i], dec=self.dec[i],
                     mwebv_sfd98=self.mwebv_sfd98[i], mwebv=self.mwebmv[i], 
                     **{p: v[i] for p, v in self.lightcurve.items()})
                for i in range(self.ntransient)]

    
#######################################
#                                     #
# Generator: SN Ia                    #
#                                     #
#######################################
class SNIaGenerator( TransientGenerator ):
    """
    This child Class enable to add the SN properties to the
    transient generator.
    """
    
    # -------------------- #
    # - Hacked Methods   - #
    # -------------------- #
    def set_transient_parameters(self,ratekind="basic",ratefunc=None,
                                 lcsimulation="basic",
                                 lcsource="salt2",type_=None,
                                 update=True,lcsimul_prop={}):
        """
        Add to the TransientGenerator the SN Ia properties
        """
        # - If this works, you are good to go
        type_= "Ia" if type_ is None or "Ia" not in type_ else type_
        f = LightCurveGenerator().get_lightcurve_func(transient=type_,
                                                      simulation=lcsimulation)

        super(SNIaGenerator,self).set_transient_parameters(type_=type_,
                                                           ratekind=ratekind,
                                                           ratefunc=ratefunc,
                                                           lcsource=lcsource,
                                                           lcsimul_func=f,
                                                           lcsimul_prop=lcsimul_prop,
                                                           update=update)
        
    # ----------------- #
    # - Set Methods   - #
    # ----------------- #
    
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def color(self):
        if not self.has_lightcurves() or "c" not in self.lightcurve:
            raise AttributeError("no 'lightcurve' defined or no color ('c') on it")
        return np.asarray(self.lightcurve["c"])

    @property
    def x1(self):
        if not self.has_lightcurves() or "x1" not in self.lightcurve:
            raise AttributeError("no 'lightcurve' defined or no x1 on it")
        return np.asarray(self.lightcurve["x1"])

    @property
    def host(self):
        if not self.has_host_param():
            return np.asarray([None]*self.ntransient)
        return np.asarray(self.lightcurve["host"])

    def has_host_param(self):
        if not self.has_lightcurves():
            raise AttributeError("no 'lightcurve' defined")
        return True if "host" in self.lightcurve_param_names \
          else False
    

    
    
#######################################
#                                     #
# Generator: SN Rate                  #
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

    def get_ratefunc(self,transient="Ia",ratekind="basic",ratefunc=None):
        """
        Parameters
        ----------

        Return
        ------
        function(z)
        """
        # ---------------
        # - Transients
        if transient is None or transient == "":
            avialable_rates = self.known_rates
            transient = None
        elif transient == "Ia":
            avialable_rates = self.known_Ia_rates
        elif ratekind != "custom":
            raise ValueError("'%s' is not a known transient"%transient)
    
        # ---------------
        # - Rate Kinds
        if ratekind == "custom":
            if ratefunc is None:
                raise ValueError("Ratekind 'custom' requires ratefunc")
            return ratefunc
        elif ratekind not in avialable_rates:
            raise ValueError("not '%s' rate kind for '%s'"%(ratekind,transient)+\
                             "These are: "+",".join(avialable_rates))
        # -- the output
        if transient is None:
            return eval("self.rate_%s"%(ratekind))
        return eval("self.rate_%s_%s"%(transient,ratekind))

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
        return 3.e-5

    
    def rate_Ia_basiclow(self,z):
        """
        """
        return self.rate_basic(z) / 10
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
    def get_lightcurve_func(self,transient="Ia",simulation="basic"):
        """
        Parameters
        ----------

        Return
        ------
        a lightcurve generator function
        """
        # ---------------
        # - Transients
        if transient is None or transient == "":
            avialable_lc = self.known_lightcurve_simulation
            if len(avialable_lc) == 0:
                raise NotImplementedError("no lightcurve simulation implemented")
            
            transient = None
            
        elif transient == "Ia":
            avialable_lc = self.known_Ia_lightcurve_simulation

        else:
            raise ValueError("'%s' is not a known transient"%transient)
    
        # ---------------
        # - Rate Kinds
        if simulation not in avialable_lc:
            raise ValueError("not '%s' rate kind for '%s'"%(simulation,transient)+\
                             "These are: "+",".join(avialable_lc))
        # -- the output
        if transient is None:
            return eval("self.lightcurve_%s"%(simulation))
        return eval("self.lightcurve_%s_%s"%(transient,simulation))
    
    def set_model(self,model):
        """
        """
        self._side_properties["model"] = model
        
    # ========================== #
    # = LC Kind                = #
    # ========================== #
    
    # ----------------- #
    # - Ia LC         - #
    # ----------------- #
    def lightcurve_Ia_basic(self,redshifts,model=None,
                            color_mean=0,color_sigma=0.1,
                            stretch_mean=0,stretch_sigma=1,
                            #source="salt2",
                            ):
        """
        """
        # ----------------
        # - Models        
        #self.set_model(sncosmo.Model(source=source))
        if model is None:
            self.set_model(sncosmo.Model(source='salt2'))
        else:
            self.set_model(model)
        x0 = []
        for z in redshifts:
            self.model.set(z=z)
            mabs = normal(-19.3, 0.3)
            self.model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
            x0.append(self.model.get('x0'))
            
        ntransient = len(redshifts)
        
        return {'x0':np.array(x0),
                'x1':normal(stretch_mean, stretch_sigma, ntransient),
                'c':normal(color_mean, color_sigma, ntransient)}

    def lightcurve_Ia_random(self,redshifts,model=None,
                            color_mean=0,color_sigma=0.1,
                            stretch_mean=0,stretch_sigma=1,
                            x0_mean=1e-5,x0_sigma=1.e-5,
                            #source="salt2",
                            ):
        """
        """
        # ----------------
        # - Models        
        ntransient = len(redshifts)
        
        return {'x0':normal(x0_mean, x0_sigma, ntransient),
                'x1':normal(stretch_mean, stretch_sigma, ntransient),
                'c':normal(color_mean, color_sigma, ntransient)}

    
    def lightcurve_Ia_realistic(self,redshifts,model=None,
                            mabs_mean = -19.3, mabs_sigma=0.12,
                            color_mean=0,color_sigma=0.1,
                            stretch_mean=(0.5,-1),stretch_sigma=(1,1),
                            stretch_thr=0.75, alpha=0.13, beta=3
                            ):
        """
        stretch parameters assume bimodal distribution
        stretch_thr is a threshold for uniformly drawn number used to determine 
        in which mode the SNe is  
        """
        # ----------------
        # - Models        
        #self.set_model(sncosmo.Model(source=source))
        if model is None:
            self.set_model(sncosmo.Model(source='salt2'))
        else:
            self.set_model(model)
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
            self.model.set(z=z)
            mabs = normal(mabs_mean, mabs_sigma)
            mabs -= alpha*x1[k] - beta*c[k]
            self.model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
            x0.append(self.model.get('x0'))
            
        ntransient = len(redshifts)
        
        return {'x0':np.array(x0),
                'x1':np.array(x1),
                'c':np.array(c)}

    
    def lightcurve_Ia_hostdependent():
        raise NotImplementedError("To be done")
    # ========================== #
    # = Properties             = #
    # ========================== #
    @property
    def model(self):
        return self._side_properties['model']
    
    @property
    def known_lightcurve_simulation(self):
        return self._parse_rate_("lightcurve")
    
    @property
    def known_Ia_lightcurve_simulation(self):
        return self._parse_rate_("lightcurve_Ia")
    
    
