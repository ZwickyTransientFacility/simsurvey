#! /usr/bin/env python
# -*- coding: utf-8 -*-


import warnings
import numpy as np
from copy import deepcopy
from collections import OrderedDict as odict

import sncosmo
from astropy.table import Table, vstack

from astrobject.astrobject.baseobject import BaseObject
from astrobject.utils.tools           import kwargs_update
from astrobject.utils.plot.skybins    import SurveyField, SurveyFieldBins 


__all__ = ["SimulSurvey", "SurveyPlan"] # to be changed

#######################################
#                                     #
# Survey: Simulation Base             #
#                                     #
#######################################
class SimulSurvey( BaseObject ):
    """
    Basic survey object
    (far from finished)
    """
    _properties_keys         = ["generator","instruments","plan"]
    _side_properties_keys    = ["cadence"]
    _derived_properties_keys = ["observations"]
    
    def __init__(self,generator=None, plan=None,
                 instprop=None,
                 empty=False):
        """
        Parameters:
        ----------
        generator: [simultarget.transient_generator or derived child like sn_generator]

        
        """
        self.__build__()
        if empty:
            return

        self.create(generator, plan, instprop)

    def create(self, generator, plan, instprop):
        """
        """
        if generator is not None:
            self.set_target_generator(generator)

        if plan is not None:
            self.set_plan(plan)

        if instprop is not None:
            self.set_instruments(instprop)

    # =========================== #
    # = Main Methods            = #
    # =========================== #

    # ---------------------- #
    # - Get Methods        - #
    # ---------------------- #
    def get_lightcurves(self):
        """
        """
        if not self.is_set():
            raise AttributeError("plan, generator or instrument not set")

        return [(sncosmo.realize_lcs(obs, self.generator.model, [p])[0] 
                 if obs is not None else None)
                for p, obs in zip(self.generator.lightcurve_full_param, 
                                  self.observations)]

    # ---------------------- #
    # - Setter Methods     - #
    # ---------------------- #
    
    # -------------
    # - Targets
    def set_target_generator(self, generator):
        """
        """
        if "__nature__" not in dir(generator) or\
          generator.__nature__ != "TransientGenerator":
            raise TypeError("generator must be an astrobject TransientGenerator")

        if not generator.has_lightcurves():
            warnings.warn("No lightcurves set in the given transient generator")

        self._properties["generator"] = generator

    # -------------
    # - SurveyPlan
    def set_plan(self,plan):
        """
        """
        # ----------------------
        # - Load cadence here
        if "__nature__" not in dir(plan) or \
          plan.__nature__ != "SurveyPlan":
            raise TypeError("the input 'plan' must be an astrobject SurveyPlan")
        self._properties["plan"] = plan
        
        # ----------------------------
        # - Set back the observations
        self._reset_observations_()

    # -------------
    # - Instruments
    def set_instruments(self,properties):
        """
        properties must be a dictionary containing the
        instruments' information (bandname,gain,zp,zpsys) related
        to each bands
        
        example..
        ---------
        properties = {"desg":{"gain":1,"zp":30,"zpsys":'ab'},
                      "desr":{"gain":1,"zp":30,"zpsys":'ab'}}
        """
        prop = deepcopy(properties)
        for band,d in prop.items():
            gain,zp,zpsys=d.pop("gain"),d.pop("zp"),d.pop("zpsys","ab")
            if gain is None or zp is None:
                raise ValueError('gain or zp is None or not defined for %s'%band)
            self.add_instrument(band,gain,zp,zpsys,
                                update=False,**d)
        
        self._reset_observations_()

    # ---------------------- #
    # - Add Stuffs         - #
    # ---------------------- #
    def add_instrument(self,bandname,gain,zp,zpsys="ab",
                       force_it=True,update=True,**kwargs):
        """
        kwargs could be any properties you wish to save with the instrument
        """
        if self.instruments is None:
            self._properties["instruments"] = {}
            
        if bandname in self.instruments.keys() and not force_it:
            raise AttributeError("%s is already defined."+\
                                 " Set force_it to True to overwrite it. ")
                                 
        instprop = {"gain":gain,"zp":zp,"zpsys":zpsys}
        self.instruments[bandname] = kwargs_update(instprop,**kwargs)
        
        if update:
            self._reset_observations_()
            
    # ---------------------- #
    # - Recover Methods    - #
    # ---------------------- #
    #def recover_targets(self):
    #    """
    #    bunch threshold...
    #    """
    #
    #def recover_lightcurves(self):
    #    """
    #    """

    # =========================== #
    # = Internal Methods        = #
    # =========================== #
    def _update_lc_(self):
        """
        """
        # -----------------------------
        # -- Do you have all you need ?
        if not self.is_set():
            return

    def _reset_observations_(self):
        """
        """
        self._derived_properties["observations"] = None
        
    def _load_observations_(self):
        """
        """
        # -------------
        # - Input test
        if self.plan is None or self.instruments is None:
            raise AttributeError("Plan or Instruments is not set.")
        
        # -----------------------
        # - Check if instruments exists
        all_instruments = np.unique(self.cadence["band"])
        if not np.all([i in self.instruments.keys() for i in all_instruments]):
            raise ValueError("Some of the instrument in cadence have not been defined."+"\n"+
                             "given instruments :"+", ".join(all_instruments.tolist())+"\n"+
                             "known instruments :"+", ".join(self.instruments.keys()))
            
        # -----------------------
        # - Based on the model get a reasonable time scale for each transient
        mjd = self.generator.mjd
        z = np.array(self.generator.zcmb)
        mjd_range = np.array([mjd + self.generator.model.mintime() * (1 + z), 
                              mjd + self.generator.model.maxtime() * (1 + z)])
        
        # -----------------------
        # - Lets build the tables
        self.plan.observe(self.generator.ra, self.generator.dec,
                          mjd_range=mjd_range)

        self._derived_properties["observations"] = [(Table(
            {"time": obs["time"],
             "band": obs["band"],
             "skynoise": obs["skynoise"],
             "gain":[self.instruments[b]["gain"] for b in obs["band"]],
             "zp":[self.instruments[b]["zp"] for b in obs["band"]],
             "zpsys":[self.instruments[b]["zpsys"] for b in obs["band"]]
            }) if len(obs) > 0 else None) for obs in self.plan.observed]
    
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def instruments(self):
        """The basic information relative to the instrument used for the survey"""
        return self._properties["instruments"]
    
    @property
    def generator(self):
        """The instance that enable to create fake targets"""
        return self._properties["generator"]

    @property
    def plan(self):
        """This is the survey plan including field definitions and telescope pointings"""
        return self._properties["plan"]

    def is_set(self):
        """This parameter is True if this has cadence, instruments and genetor set"""
        return not (self.instruments is None or \
                    self.generator is None or \
                    self.plan is None)

    # ------------------
    # - Side properties
    @property
    def cadence(self):
        """This is a table containing where the telescope is pointed with which band."""
        if self._properties["plan"] is not None:
            return self._properties["plan"].cadence
        else:
            raise ValueError("Property 'plan' not set yet")

                    
    # ------------------
    # - Derived values
    @property
    def observations(self):
        """Observations derived from cadence and instrument properties.
        Note that the first time this is called, observations will be recorded"""
        
        if self._derived_properties["observations"] is None:
            self._load_observations_()
            
        return self._derived_properties["observations"]
                                          

#######################################
#                                     #
# Survey: Plan object                 #
#                                     #
#######################################
class SurveyPlan( BaseObject ):
    """
    Survey Plan
    contains the list of observation times, bands and pointings and
    can return that times and bands, which a transient is observed at/with.
    A list of fields can be given to simplify adding observations and avoid 
    lookups whether an object is in a certain field.
    Currently assumes a single instrument, especially for FoV width and height.
    [This may be useful for the cadence property of SimulSurvey]
    """
    __nature__ = "SurveyPlan"

    _properties_keys         = ["cadence", "width", "height"]
    _side_properties_keys    = ["fields"]
    _derived_properties_keys = ["observed"]
    
    def __init__(self, time=None, ra=None, dec=None, band=None, skynoise=None, 
                 obs_field=None, width=7., height=7., fields=None, empty=False,
                 load_opsim=None):
        """
        Parameters:
        ----------
        TBA
        
        """
        self.__build__()
        if empty:
            return
    
        self.create(time=time,ra=ra,dec=dec,band=band,skynoise=skynoise,
                    obs_field=obs_field,fields=fields, load_opsim=load_opsim)

    def create(self, time=None, ra=None, dec=None, band=None, skynoise=None, 
               obs_field=None, width=7., height=7., fields=None, 
               load_opsim=None):
        """
        """
        self._properties["width"] = float(width)
        self._properties["height"] = float(height)
        
        if fields is not None:
            self.set_fields(**fields)

        if load_opsim is None:
            self.add_observation(time,band,skynoise,ra=ra,dec=dec,field=obs_field)
        else:
            self.load_opsim(load_opsim)

    # =========================== #
    # = Main Methods            = #
    # =========================== #

    # ---------------------- #
    # - Get Methods        - #
    # ---------------------- #
    
    # ---------------------- #
    # - Setter Methods     - #
    # ---------------------- #
    def set_fields(self, ra=None, dec=None, **kwargs):
        """
        """
        kwargs["width"] = kwargs.get("width", self.width)
        kwargs["height"] = kwargs.get("height", self.height)

        self._side_properties["fields"] = SurveyFieldBins(ra, dec, **kwargs)

        if self.cadence is not None and np.any(np.isnan(self.cadence['field'])):
            warnings.warning("cadence was already set, field pointing will be updated")
            self._update_field_radec()

    def add_observation(self, time, band, skynoise, ra=None, dec=None, field=None):
        """
        """
        if ra is None and dec is None and field is None:
            raise ValueError("Either field or ra and dec must to specified.")
        elif ra is None and dec is None:
            if self.fields is None:
                raise ValueError("Survey fields not defined.")
            else:
                ra = self.fields.ra[field]
                dec = self.fields.dec[field]
        elif field is None:
            field = np.array([np.nan for r in ra])

        new_obs = Table({"time": time,
                         "band": band,
                         "skynoise": skynoise,
                         "RA": ra,
                         "Dec": dec,
                         "field": field})

        if self._properties["cadence"] is None:
            self._properties["cadence"] = new_obs
        else:
            self._properties["cadence"] = vstack((self._properties["cadence"], 
                                                  new_obs))

    # ---------------------- #
    # - Load Method        - #
    # ---------------------- #            
    def load_opsim(self, filename, table_name="ptf", band_prefix="ptf"):
        """
        see https://confluence.lsstcorp.org/display/SIM/Summary+Table+Column+Descriptions
        for format description
        
        Currently only the used columns are loaded

        table_name -- name of table in SQLite DB (deafult "ptf" because of 
                      Eric's example)
        band_prefix -- string prefix for sncosmo band name
        """
        import sqlite3
        connection = sqlite3.connect(filename)

        # define columns name and keys to be fetched
        to_fetch = odict()
        to_fetch['time'] = 'expMJD'
        to_fetch['band'] = 'filter' # Currently is a float not a string in Eric's example
        to_fetch['skynoise'] = 'filtSkyBrightness' # in mag/arcsec^2 [What do we need to convert this?]
        to_fetch['ra'] = 'fieldRA'
        to_fetch['dec'] = 'fieldDec'
        to_fetch['field'] = 'fieldID'
        
        loaded = odict()
        for key, value in to_fetch.items():
            # This is not safe against injection (but should be OK)
            # TODO: Add function to sanitize input
            cmd = 'SELECT %s from %s;'%(value, table_name)
            tmp = connection.execute(cmd)
            loaded[key] = [a[0] for a in tmp]

        connection.close()

        self.add_observation(loaded['time'],loaded['band'],loaded['skynoise'],
                             ra=loaded['ra'],dec=loaded['dec'],
                             field=loaded['field'])

    # ================================== #
    # = Observation time determination = #
    # ================================== #
    def observe(self, ra, dec, mjd_range=None):
        """
        """
        self._derived_properties["observed"] = self.observed_on(ra, dec,
                                                                mjd_range)

        return self._derived_properties["observed"]

    def observed_on(self, ra, dec, mjd_range=None):
        """
        mjd_range must be (2,N)-array 
        where N is the length of ra and dec
        """
        single_coord = None

        # first get the observation times and bands for pointings without a
        # field number use this to determine whether ra and dec were arrays or
        # floats (since this is done in SurveyField.coord_in_field there is no
        # need to redo this)
        for k, obs in enumerate(self.cadence[np.isnan(self.cadence["field"])]):
            tmp_f = SurveyField(obs["RA"], obs["Dec"], 
                                self.width, self.height)
            b = tmp_f.coord_in_field(ra, dec)
            
            # Setup output as dictionaries that can be converted to Tables and
            # sorted later
            if k == 0:
                if type(b) is np.bool_:
                    single_coord = True
                    out = {'time': [], 'band': [], 'skynoise': []}
                else:
                    single_coord = False
                    out = [{'time': [], 'band': [], 'skynoise': []} for r in ra]

            if single_coord:
                if b:
                    out['time'].extend(obs['time'].quantity.value)
                    out['band'].extend(obs['band'])
                    out['skynoise'].extend(obs['skynoise'].quantity.value)
            else:
                for l in np.where(b)[0]:
                    out[l]['time'].extend(obs['time'].quantity.value)
                    out[l]['band'].extend(obs['band'])
                    out[l]['skynoise'].extend(obs['skynoise'].quantity.value)

        # Now get the other observations (those with a field number)
        if (self.fields is not None and 
            not np.all(np.isnan(self.cadence["field"]))):
            b = self.fields.coord2field(ra, dec)
            
            # if all pointings were in fields create new dicts, otherwise append
            if single_coord is None:
                if type(b) is not list:
                    single_coord = True
                    out = {'time': [], 'band': [], 'skynoise': []}
                else:
                    single_coord = False
                    out = [{'time': [], 'band': [], 'skynoise': []} for r in ra]
            
            if single_coord:
                for l in b:
                    mask = (self.cadence['field'] == l)
                    out['time'].extend(self.cadence['time'][mask].quantity.value)
                    out['band'].extend(self.cadence['band'][mask])
                    out['skynoise'].extend(self.cadence['skynoise']
                                           [mask].quantity.value)
            else:
                for k, idx in enumerate(b):
                    for l in idx:
                        mask = (self.cadence['field'] == l)
                        out[k]['time'].extend(self.cadence['time'][mask].quantity.value)
                        out[k]['band'].extend(self.cadence['band'][mask])
                        out[k]['skynoise'].extend(self.cadence['skynoise']
                                                  [mask].quantity.value)

        # Make Tables and sort by time
        if single_coord:
            table = Table(out, meta={'RA': ra, 'Dec': dec})
            idx = np.argsort(table['time'])
            if mjd_range is None:
                return table[idx]
            else:
                t = table[idx]
                return t[(t['time'] >= mjd_range[0]) &
                         (t['time'] <= mjd_range[1])]
        else:
            tables = [Table(a, meta={'RA': r, 'Dec': d}) for a, r, d 
                      in zip(out, ra, dec)]
            idx = [np.argsort(t['time']) for t in tables]
            if mjd_range is None:
                return [t[i] for t, i in zip(tables, idx)]
            else:
                ts = [t[i] for t, i in zip(tables, idx)]
                return [t[(t['time'] >= mjd_range[0][k]) &
                          (t['time'] <= mjd_range[1][k])] 
                        for k, t in enumerate(ts)]

    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def cadence(self):
        """Table of observations"""
        return self._properties["cadence"]

    @property
    def width(self):
        """field width"""
        return self._properties["width"]

    @property
    def height(self):
        """field height"""
        return self._properties["height"]

    # ------------------
    # - Side properties                    
    @property
    def fields(self):
        """Observation fields"""
        return self._side_properties["fields"]

    # ------------------
    # - Derived values
    @property
    def observed(self):
        """Saved observation times per object"""
        return self._derived_properties["observed"]
