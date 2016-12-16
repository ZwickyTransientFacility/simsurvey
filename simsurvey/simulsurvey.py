#! /usr/bin/env python
# -*- coding: utf-8 -*-


import warnings
import numpy as np
import cPickle
from copy import deepcopy
from collections import OrderedDict as odict
from itertools import izip

import sncosmo
from sncosmo.photdata      import dict_to_array
from astropy.table         import Table, vstack, hstack

from propobject import BaseObject

from utils.tools   import kwargs_update, range_args, range_length, get_progressbar
from utils.skybins import SurveyField, SurveyFieldBins 

_d2r = np.pi/180

__all__ = ["SimulSurvey", "SurveyPlan", "LightcurveCollection"]

#######################################
#                                     #
# Survey: Simulation Base             #
#                                     #
#######################################
class SimulSurvey( BaseObject ):
    """
    Basic survey object
    """
    __nature__ = "SimulSurvey"
    
    PROPERTIES         = ["generator", "instruments","plan"]
    SIDE_PROPERTIES    = ["cadence", "blinded_bias"]
    DERIVED_PROPERTIES = ["obs_fields", "obs_ccd",
                          "non_field_obs", "non_field_obs_ccd",
                          "non_field_obs_exist"]

    def __init__(self,generator=None, plan=None,
                 instprop=None, blinded_bias=None,
                 empty=False):
        """
        Parameters:
        ----------
        generator: [simultarget.transient_generator or derived child like sn_generator]

        """
        self.__build__()
        if empty:
            return

        self.create(generator, plan, instprop, blinded_bias)

    def create(self, generator, plan, instprop, blinded_bias):
        """
        """
        if generator is not None:
            self.set_target_generator(generator)

        if plan is not None:
            self.set_plan(plan)

        if instprop is not None:
            self.set_instruments(instprop)

        if blinded_bias is not None:
            self.set_blinded_bias(blinded_bias)

    # =========================== #
    # = Main Methods            = #
    # =========================== #

    # ---------------------- #
    # - Get Methods        - #
    # ---------------------- #
    def get_lightcurves(self, *args, **kwargs):
        """
        """
        args = range_args(self.generator.ntransient, *args)
        progress_bar = kwargs.pop('progress_bar', False)

        if not self.is_set():
            raise AttributeError("plan, generator or instrument not set")

        lcs = LightcurveCollection(empty=True)
        gen = izip(xrange(*args),
                   self.generator.get_lightcurve_full_param(*args),
                   self._get_observations_(*args))

        progress_bar_success = False
        
        if progress_bar:
            self._assign_obs_fields_(progress_bar=True)
            self._assign_non_field_obs_(progress_bar=True)
            
            try:
                print 'Generating lightcurves'
                ntransient = range_length(*args)
                bar = get_progressbar(ntransient)

                with get_progressbar(ntransient) as bar:
                    for k, p, obs in gen:
                        if obs is not None:
                            lcs.add(self._get_lightcurve_(p, obs, k))
                            bar.update()
                            
                progress_bar_success = True
            except ImportError:
                progress_bar_success = False
            except IOError:
                progress_bar_success = False

        if not progress_bar_success:
            for k, p, obs in gen:
                if obs is not None:
                    lcs.add(self._get_lightcurve_(p, obs, k))

        return lcs
            
    def _get_lightcurve_(self, p, obs, idx_orig=None):
        """
        """        
        if obs is not None:
            ra, dec, mwebv_sfd98 = p.pop('ra'), p.pop('dec'), p.pop('mwebv_sfd98')

            # Get unperturbed lc from sncosmo
            lc = sncosmo.realize_lcs(obs, self.generator.model, [p],
                                     scatter=False)[0]

            if 'field' in obs.colnames:
                lc = hstack((lc, obs[('field',)]))
            if 'ccd' in obs.colnames:
                lc = hstack((lc, obs[('ccd',)]))
                
            # Replace fluxerrors with covariance matrix that contains
            # correlated terms for the calibration uncertainty
            fluxerr = np.sqrt(obs['skynoise']**2 +
                              np.abs(lc['flux']) / obs['gain'])
            
            fluxcov = np.diag(fluxerr**2)
            save_cov = False
            for band in set(obs['band']):
                if self.instruments[band]['err_calib'] is not None:
                    save_cov = True
                    idx = np.where(obs['band'] == band)[0]
                    err = self.instruments[band]['err_calib']
                    for k0 in idx:
                        for k1 in idx:
                            fluxcov[k0,k1] += (lc['flux'][k0] * 
                                               lc['flux'][k1] *
                                               err**2)

            # Add random (but correlated) noise to the fluxes
            fluxchol = np.linalg.cholesky(fluxcov)
            flux = lc['flux'] + fluxchol.dot(np.random.randn(len(lc)))

            # Apply blinded bias if given
            if self.blinded_bias is not None:
                bias_array = np.array([self.blinded_bias[band]
                                       if band in self.blinded_bias.keys() else 0
                                       for band in obs['band']])
                flux *= 10 ** (-0.4*bias_array)

            lc['flux'] = flux
            lc['fluxerr'] = np.sqrt(np.diag(fluxcov))

            # Additional metadata for the lc fitter
            lc.meta['ra'] = ra
            lc.meta['dec'] = dec
            if save_cov:
                lc.meta['fluxcov'] = fluxcov
            lc.meta['mwebv_sfd98'] = mwebv_sfd98
            if idx_orig is not None:
                lc.meta['idx_orig'] = idx_orig
        else:
            lc = None

        return lc

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
        # self._reset_observations_()
        self._reset_obs_fields_()
        self._reset_non_field_obs_()

    # -------------
    # - Instruments
    def set_instruments(self,properties):
        """
        properties must be a dictionary containing the
        instruments' information (bandname,gain,zp,zpsys,err_calib) related
        to each bands


        example..
        ---------
        properties = {"desg":{"gain":1,"zp":30,"zpsys":'ab',"err_calib":0.005},
                      "desr":{"gain":1,"zp":30,"zpsys":'ab',"err_calib":0.005}}
        """
        prop = deepcopy(properties)
        for band,d in prop.items():
            gain,zp,zpsys = d.pop("gain"),d.pop("zp"),d.pop("zpsys","ab")
            err_calib = d.pop("err_calib", None)
            if gain is None or zp is None:
                raise ValueError('gain or zp is None or not defined for %s'%band)
            self.add_instrument(band,gain,zp,zpsys,err_calib,
                                update=False,**d)

        #self._reset_observations_()

    # -----------------------
    # - Blinded bias in bands
    def set_blinded_bias(self, bias):
        """Expect input dict of band and bounds maximum bias
        Bias will be drawn from uniform distribution
        """
        self._side_properties['blinded_bias'] = {k: np.random.uniform(-v, v) 
                                            for k, v in bias.items()}

    # ---------------------- #
    # - Add Stuffs         - #
    # ---------------------- #
    def add_instrument(self,bandname,gain,zp,zpsys="ab",err_calib=None,
                       force_it=True,update=True,**kwargs):
        """
        kwargs could be any properties you wish to save with the instrument
        """
        if self.instruments is None:
            self._properties["instruments"] = {}

        if bandname in self.instruments.keys() and not force_it:
            raise AttributeError("%s is already defined."+\
                                 " Set force_it to True to overwrite it. ")

        instprop = {"gain":gain,"zp":zp,"zpsys":zpsys,"err_calib":err_calib}
        self.instruments[bandname] = kwargs_update(instprop,**kwargs)

        if update:
            # self._reset_observations_()
            pass

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
            
    def _get_observations_(self, *args):
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
        mjd_range = [mjd + self.generator.model.mintime() * (1 + z), 
                     mjd + self.generator.model.maxtime() * (1 + z)]

        # -----------------------
        # - Let's build the tables
        for k in xrange(*range_args(self.generator.ntransient, *args)):
            obs = self.plan.observed_on(self.obs_fields[k],
                                        (self.obs_ccds[k]
                                         if self.obs_ccds is not None
                                         else None),
                                        self.non_field_obs[k],
                                        (self.non_field_obs_ccds[k]
                                         if self.non_field_obs_ccds is not None
                                         else None),
                                        (mjd_range[0][k], mjd_range[1][k]))

            
            data = [[self.instruments[b]["gain"] for b in obs["band"]],
                    [self.instruments[b]["zp"] for b in obs["band"]],
                    [self.instruments[b]["zpsys"] for b in obs["band"]]]
            names = ["gain", "zp", "zpsys"]
            
            if len(obs) > 0:
                yield hstack((obs,Table(data=data, names=names)))
            else:
                yield None

    def _assign_obs_fields_(self, progress_bar=False):
        """
        """
        f, c = self.plan.get_obs_fields(
            self.generator.ra,
            self.generator.dec,
            field_id=np.unique(self.cadence['field']),
            progress_bar=progress_bar
        )
        self._derived_properties["obs_fields"] = f
        self._derived_properties["obs_ccds"] = c

    def _reset_obs_fields_(self):
        """
        """
        self._derived_properties["obs_fields"] = None
        self._derived_properties["obs_ccds"] = None

    def _assign_non_field_obs_(self, progress_bar=False):
        """
        """
        f, c = self.plan.get_non_field_obs(
            self.generator.ra,
            self.generator.dec,
            progress_bar=progress_bar
        )
        self._derived_properties["non_field_obs"] = f
        self._derived_properties["non_field_obs_ccds"] = c

    def _reset_non_field_obs_(self):
        """
        """
        self._derived_properties["non_field_obs"] = None
        self._derived_properties["non_field_obs_ccd"] = None
        self._derived_properties["non_field_obs_exist"] = None
    
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

    @property
    def blinded_bias(self):
        """Blinded bias applied to specific bands for all observations"""
        return self._side_properties["blinded_bias"]

    # ------------------
    # - Derived values
    @property
    def obs_fields(self):
        """Transients are assigned fields that they are found"""
        if self._derived_properties["obs_fields"] is None:
            self._assign_obs_fields_()

        return self._derived_properties["obs_fields"]

    @property
    def obs_ccds(self):
        """Transients are assigned fields that they are found"""
        if (self._derived_properties["obs_fields"] is None and
            self._derived_properties["obs_ccds"] is None and
            self.plan.ccds is not None):
            self._assign_obs_fields_()

        return self._derived_properties["obs_ccds"]

    @property
    def non_field_obs(self):
        """If the plan contains pointings with field id, prepare a list of those."""
        if (self._derived_properties["non_field_obs"] is None
            and self.non_field_obs_exist is False):
            self._assign_non_field_obs_()
            
        if self._derived_properties["non_field_obs"] is None:
            self._derived_properties["non_field_obs_exist"] = False
        else:
            self._derived_properties["non_field_obs_exist"] = True

        if self.non_field_obs_exist is False:
            return [None for k in xrange(self.generator.ntransient)]
        return self._derived_properties["non_field_obs"]

    @property
    def non_field_obs_ccds(self):
        """If the plan contains pointings with field id, prepare a list of those."""
        if self.non_field_obs_exist is False:
            return [None for k in xrange(self.generator.ntransient)]
        return self._derived_properties["non_field_obs_ccds"]
        
    @property
    def non_field_obs_exist(self):
        """Avoid checking for non-field pointings more than once."""
        return self._derived_properties["non_field_obs_exist"]

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
    """
    __nature__ = "SurveyPlan"

    PROPERTIES         = ["cadence", "width", "height"]
    SIDE_PROPERTIES    = ["fields", "ccds"]
    DERIVED_PROPERTIES = []

    def __init__(self, time=None, ra=None, dec=None, band=None, skynoise=None, 
                 obs_field=None, width=7.295, height=7.465, fields=None, empty=False,
                 load_opsim=None, **kwargs):
        """
        Parameters:
        ----------
        TBA

        """
        self.__build__()
        if empty:
            return

        self.create(time=time,ra=ra,dec=dec,band=band,skynoise=skynoise,
                    obs_field=obs_field,fields=fields, load_opsim=load_opsim, **kwargs)

    def create(self, time=None, ra=None, dec=None, band=None, skynoise=None, 
               obs_field=None, width=7.295, height=7.465, fields=None, 
               load_opsim=None, **kwargs):
        """
        """        
        self._properties["width"] = float(width)
        self._properties["height"] = float(height)
        self._side_properties["ccds"] = kwargs.pop('ccds', None)
        
        if fields is not None:
            self.set_fields(**fields)

        if load_opsim is None:
            self.add_observation(time, band, skynoise, ra=ra, dec=dec,
                                 field=obs_field)
        else:
            self.load_opsim(load_opsim, **kwargs)

    # =========================== #
    # = Main Methods            = #
    # =========================== #

    # ---------------------- #
    # - Get Methods        - #
    # ---------------------- #

    # ---------------------- #
    # - Setter Methods     - #
    # ---------------------- #
    def set_fields(self, ra=None, dec=None, ccds=None, **kwargs):
        """
        """
        kwargs["width"] = kwargs.get("width", self.width)
        kwargs["height"] = kwargs.get("height", self.height)

        self._side_properties["fields"] = SurveyFieldBins(ra, dec, ccds=self.ccds,
                                                          **kwargs)

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
                idx = self.fields.field_id_index[field]
                ra = self.fields.ra[idx]
                dec = self.fields.dec[idx]
        elif field is None:
            field = np.array([np.nan for r in ra])

        new_obs = Table(data=[time, band, skynoise, ra, dec, field],
                        names=['time', 'band', 'skynoise', 'RA', 'Dec', 'field'])

        if self._properties['cadence'] is None:
            self._properties['cadence'] = new_obs
        else:
            self._properties['cadence'] = vstack((self._properties['cadence'], 
                                                  new_obs))

    # ---------------------- #
    # - Load Method        - #
    # ---------------------- #
    def load_opsim(self, filename, survey_table='Summary', field_table='Field',
                   band_dict=None, default_depth=21, zp=30):
        """
        see https://confluence.lsstcorp.org/display/SIM/Summary+Table+Column+Descriptions
        for format description

        Currently only the used columns are loaded

        table_name -- name of table in SQLite DB (deafult "ptf" because of 
                      Eric's example)
        band_dict -- dictionary for converting filter names 
        zp -- zero point for converting sky brightness from mag to flux units
              (should match the zp used in instprop for SimulSurvey)
        """        
        import sqlite3
        connection = sqlite3.connect(filename)

        def _fetch(keys, table):
            loaded = odict()
            for key in keys:
                # This is not safe against injection (but should be OK)
                # TODO: Add function to sanitize input
                cmd = 'SELECT %s from %s;'%(key, table)
                tmp = connection.execute(cmd)
                loaded[key] = np.array([a[0] for a in tmp])

            return loaded
        
        loaded = _fetch(['expMJD', 'filter', 'fieldRA', 'fieldDec',
                         'fieldID', 'fiveSigmaDepth'],
                        survey_table)
        fields = _fetch(['fieldID', 'fieldRA', 'fieldDec'],
                       field_table)
        
        connection.close()

        loaded['fieldRA'] /= _d2r
        loaded['fieldDec'] /= _d2r

        loaded['fiveSigmaDepth'] = np.array([(d if d is not None else default_depth)
                                    for d in loaded['fiveSigmaDepth']])
        
        loaded['skynoise'] = 10 ** (-0.4 * (loaded['fiveSigmaDepth']-zp)) / 5
        
        if band_dict is not None:
            loaded['filter'] = [band_dict[band] for band in loaded['filter']]
        else:
            loaded['filter'] = loaded['filter']
 
        self.add_observation(loaded['expMJD'],loaded['filter'],loaded['skynoise'],
                             ra=loaded['fieldRA'],dec=loaded['fieldDec'],
                             field=loaded['fieldID'])

        self.set_fields(ra=fields['fieldRA'], dec=fields['fieldDec'],
                        field_id=fields['fieldID'], ccds=self.ccds)
        
    # ================================== #
    # = Observation time determination = #
    # ================================== #
    def get_obs_fields(self, ra, dec, field_id=None,
                       progress_bar=False):
        """
        """
        if (self.fields is not None and 
            not np.all(np.isnan(self.cadence["field"]))):
            return self.fields.coord2field(ra, dec, field_id=field_id,
                                           progress_bar=progress_bar)
        else:
            return None, None
        
    def get_non_field_obs(self, ra, dec, progress_bar=False):
        """
        """
        observed = False
        gen = self.cadence[np.isnan(self.cadence["field"])]
        
        if progress_bar and len(gen) > 0:
            try:
                print "Finding transients observed in custom pointings"
                gen = get_progressbar(gen)
            except ImportError:
                pass
            except IOError:
                pass
                
        for k, obs in enumerate(gen):
            tmp_f = SurveyField(obs["RA"], obs["Dec"], 
                                self.width, self.height)
            b, c = tmp_f.coord_in_field(ra, dec, ccds=self.fields.ccds)

            # Setup output as dictionaries that can be converted to Tables and
            # sorted later
            if k == 0:
                if type(b) is np.bool_:
                    single_coord = True
                    out = np.array([], dtype=int)
                    ccd = np.array([], dtype=int)
                else:
                    out = [np.array([], dtype=int) for r in ra]
                    ccd = [np.array([], dtype=int) for r in ra]

            if single_coord:
                if b:
                    observed = True
                    out = np.append(out, [k])
                    ccd = np.append(ccd, [c])
            else:
                for l in np.where(b)[0]:
                    observed = True
                    out[l] = np.append(out[l], [k])
                    ccd[l] = np.append(ccd[l], [c])

        if observed:
            return out, ccd
        else:
            return None, None

    def observed_on(self, fields=None, ccds=None,
                    non_field=None, non_field_ccds=None,
                    mjd_range=None):
        """
        mjd_range must be 2-tuple
        fields and non_field np.arrays
        """
        if fields is None and non_field is None:
            raise ValueError("Provide arrays of fields and/or other pointings") 

        out = {'time': [], 'band': [], 'skynoise': [], 'field': []}

        if ccds is not None:
            out['ccd'] = []
            
        if fields is not None:
            for k, l in enumerate(fields):
                mask = (self.cadence['field'] == l)
                out['time'].extend(self.cadence['time'][mask].quantity.value)
                out['band'].extend(self.cadence['band'][mask])
                out['skynoise'].extend(self.cadence['skynoise']
                                       [mask].quantity.value)
                out['field'].extend(l*np.ones(np.sum(mask), dtype=int))
                if 'ccd' in out.keys():
                    out['ccd'].extend(ccds[k]*np.ones(np.sum(mask), dtype=int))
                
        if non_field is not None:
            mask = np.isnan(self.cadence["field"])
            out['time'].extend(self.cadence['time'][mask][non_field].quantity.value)
            out['band'].extend(self.cadence['band'][mask][non_field])
            out['skynoise'].extend(self.cadence['skynoise']
                                   [mask][non_field].quantity.value)
            out['field'].extend(np.nan*np.ones(np.sum(mask), dtype=int))
            if 'ccd' in out.keys():
                out['ccd'].extend(non_field_ccds[k][mask])
            
        table = Table(out, meta={})
        idx = np.argsort(table['time'])
        if mjd_range is None:
            return table[idx]
        else:
            t = table[idx]
            return t[(t['time'] >= mjd_range[0]) &
                     (t['time'] <= mjd_range[1])]

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

    @property
    def ccds(self):
        """Camera CCDs"""
        return self._side_properties["ccds"]

#######################################
#                                     #
# LigthcurveCollecion object          #
#                                     #
#######################################
class LightcurveCollection( BaseObject ):
    """
    LightcurveCollection
    Collects and organizes lightcurves (e.g. simulated by a Survey object)
    for easy access and serialization while try to avoid excessive memory
    use by Astropy Tables. Superficially acts like a list of tables but
    creates them on the fly from structured numpy arrays
    """
    __nature__ = "LightcurveCollection"

    PROPERTIES         = ['lcs','meta']
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    def __init__(self, lcs=None, empty=False, load=None):
        """
        Parameters:
        ----------
        TBA

        """
        self.__build__()
        if empty:
            return

        self.create(lcs=lcs, load=load)

    def create(self, lcs=None, load=None):
        """
        """
        if load is None:
            self.add(lcs)
        else:
            self.load(load)

    # =========================== #
    # = Main Methods            = #
    # =========================== #
    def add(self, lcs):
        """
        """
        if type(lcs) is list:
            meta = [lc.meta for lc in lcs]
        else:
            meta = lcs.meta

        self._add_lcs_(lcs)
        self._add_meta_(meta)    

    def load(self, filename):
        """
        """
        loaded = cPickle.load(open(filename))
        self._properties['lcs'] = loaded['lcs']
        self._properties['meta'] = loaded['meta']

    def save(self, filename):
        """
        """
        cPickle.dump({'lcs': self._properties["lcs"],
                      'meta': self._properties["meta"]},
                     open(filename, 'w'))

    # ---------------------- #
    # - Get Methods        - #
    # ---------------------- #
    def __getitem__(self, given):
        """
        """
        if isinstance(given, slice):
            return [Table(data=data,
                          meta={k: v for k, v in zip(meta.dtype.names, meta)})
                    for data, meta in
                    zip(self.lcs[given], self.meta[given])]
        else:
            meta = self.meta[given]
            return Table(data=self.lcs[given],
                         meta={k: v for k, v in zip(meta.dtype.names, meta)}) 
            
    # ---------------------- #
    # - Add Methods        - #
    # ---------------------- #

    def _add_lcs_(self, lcs):
        """
        """
        if self.lcs is None:
            self._properties['lcs'] = []

        if type(lcs) is list:
            for lc in lcs:
                self._add_lc_(lc)
        else:
            self._add_lc_(lcs)

    def _add_lc_(self, lc):
        """
        """
        self._properties['lcs'].append(lc.as_array())
            
    def _add_meta_(self, meta):
        """
        """
        if type(meta) is list:
            if self.meta is None:
                keys = [k for k in meta[0].keys()]
                dtypes = [type(v) for v in meta[0].values()]
                self._create_meta_(keys, dtypes)
                
            for meta_ in meta:
                for k in self.meta.dtype.names:
                    self._properties['meta'][k] = np.append(
                        self._properties['meta'][k],
                        meta_[k]
                    )
        else:
            if self.meta is None:
                keys = [k for k in meta.keys()]
                dtypes = [type(v) for v in meta.values()]
                self._create_meta_(keys, dtypes)
                
            for k in self.meta.dtype.names:
                self._properties['meta'][k] = np.append(
                    self._properties['meta'][k],
                    meta[k]
                )            

    def _create_meta_(self, keys, dtypes):
        """
        Create the ordered ditcionary of meta parameters based of first item
        """
        self._properties['meta'] = odict()
        for k, t in zip(keys, dtypes):
            self._properties['meta'][k] = np.array([], dtype=t)
                
    # =========================== #
    # = Properties and Settings = #
    # =========================== #
    @property
    def lcs(self):
        """List of lcs as numpy structured arrays without meta parameters"""
        return self._properties["lcs"]

    @property
    def meta(self):
        """numpy structured array with of meta parameters"""
        if self._properties["meta"] is None:
            return None
        return dict_to_array(self._properties["meta"])

    