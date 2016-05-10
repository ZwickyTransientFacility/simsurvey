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

from astrobject.astrobject.baseobject import BaseObject
from .iomock import read_mockdata



class RetrieveGaia( BaseObject ):
    """
    """
    __nature__ = "Retriever"
    
    _properties_keys = ["data"]
    _side_properties_keys = []
    _derived_properties_keys = []


    freeparameters = ["zeropoint", "a_gaiaGmB","a_gaiaBmR",
                       "a_airmass","a_airmass_gaiaGmB","a_airmass_gaiaBmR"]
        
    def __init__(self,filename=None):
        """
        """
        self.__build__()
        if filename is not None:
            self.read_data(filename)

    
    # ================= #
    # = Init          = #
    # ================= #
    def read_data(self,filename, select_ztfband=None):
        """ read the mock data in the given filename
        and create the data base"""
        # --------------------
        # -- Create the data
        dataframe = read_mockdata(filename)
        if select_ztfband is not None:
            dataframe = dataframe.mask( ~(dataframe["ZTF_filter"]==select_ztfband)).dropna(axis=0,how="all")
        
        self._properties["data"] = dataframe
        # --------------------
        # -- and the Derived Colors
        self.data["GPmBP_gaia"] = self.data["G_gaia"] - self.data["BP_gaia"]
        self.data["GPmBP_gaia_err"] = np.sqrt(self.data["G_gaia_err"]**2 + self.data["BP_gaia_err"]**2) # CAUTION NO COV TERM HERE
        self.data["BPmRP_gaia"] = self.data["BP_gaia"] - self.data["RP_gaia"]
        self.data["BPmRP_gaia_err"] = np.sqrt(self.data["BP_gaia_err"]**2 + self.data["RP_gaia_err"]**2) # CAUTION NO COV TERM HERE

    # ================= #
    # = Method        = #
    # ================= #
    def get_model(self,parameter):
        """ return the expected magnitude correction based on the input parameters """

        zeropoint, a_gaiaGmB,a_gaiaBmR,a_airmass,a_airmass_gaiaGmB,a_airmass_gaiaBmR = parameter
        
        return zeropoint +  \
          a_gaiaGmB*self.data["GPmBP_gaia"] + a_gaiaBmR*self.data["BPmRP_gaia"] + \
          a_airmass*self.data['airmass']+ \
          a_airmass_gaiaGmB*self.data['airmass']*self.data["GPmBP_gaia"] + a_airmass_gaiaBmR*self.data['airmass']*self.data["BPmRP_gaia"]    
        
    
    def get_residual(self,parameter):
        """ """
        model = self.get_model(parameter)+self.data["BP_gaia"]
        return model - self.data["ZTF_mag"]
    
    def get_chi2(self,parameter):

        res = self.get_residual(parameter) / self.data["ZTF_mag_err"]
        
        return np.sum(res**2)


    # ========================== #
    # = Bayesian stuff         = #
    # ========================== #
    def lnprob(self,parameter):
        """ This is the Bayesian posterior function (in log).
        it returns  lnprior - 0.5*Chi2
        (Assuming Chi2 = -2logLikelihood)
        """
        priors = self.lnprior(parameter)
        if not np.isfinite(priors):
            return -np.inf
        # not sure it works with the _minuit_chi2_/_scipy_chi2_  tricks
        return priors - 0.5*self.get_chi2(parameter)

        
    def lnprior(self,parameter):
        """ flat prior so far """
        return 0
    
    # ================= #
    # = Properties    = #
    # ================= #
    @property
    def data(self):
        """ pandas DataFrame containing the data """
        return self._properties["data"]

    def has_data(self):
        """ test if the current instance has data, True means yes """
        return self.data is not None

    # = Model
    @property
    def nparam(self):
        return len(self.freeparameters)
    
    # ===================== #
    # = MCMC Stuffs       = #
    # ===================== #
    def run_mcmc(self,init_guess,init_guesserr,
                 nrun=2000, walkers_per_dof=3):
        """
        """
        import emcee

        # -- set up the mcmc
        self.mcmc["ndim"], self.mcmc["nwalkers"] = \
          self.nparam, self.nparam*walkers_per_dof
        self.mcmc["nrun"] = nrun
        
        # -- init the walkers
            
        self.mcmc["pos_init"] = init_guess
        self.mcmc["pos"] = [self.mcmc["pos_init"] + np.random.randn(self.mcmc["ndim"])*init_guesserr for i in range(self.mcmc["nwalkers"])]
        # -- run the mcmc        
        self.mcmc["sampler"] = emcee.EnsembleSampler(self.mcmc["nwalkers"], self.mcmc["ndim"], self.lnprob)
        _ = self.mcmc["sampler"].run_mcmc(self.mcmc["pos"], self.mcmc["nrun"])
    
    def show_mcmc_corner(self, savefile=None, show=True,
                         truths=None,**kwargs):
        """
        **kwargs goes to corner.corner
        """
        import corner
        from astrobject.utils.mpladdon import figout
        fig = corner.corner(self.mcmc_samples, labels=self.freeparameters, 
                        truths=self.mcmc["pos_init"] if truths is None else truths,
                        show_titles=True,label_kwargs={"fontsize":"xx-large"})

        fig.figout(savefile=savefile, show=show)
        
    def show_mcmcwalkers(self, savefile=None, show=True,
                        cwalker=None, cline=None, truths=None, **kwargs):
        """ Show the walker values for the mcmc run.

        Parameters
        ----------

        savefile: [string]
            where to save the figure. if None, the figure won't be saved

        show: [bool]
            If no figure saved, the function will show it except if this is set
            to False

        cwalker, cline: [matplotlib color]
            Colors or the walkers and input values.
        """
        # -- This show the 
        import matplotlib.pyplot as mpl
        from astrobject.utils.mpladdon import figout
        if not self.has_mcmc_ran():
            raise AttributeError("you must run mcmc first")
        
        fig = mpl.figure(figsize=[7,3*self.mcmc["ndim"]])
        # -- inputs
        if cline is None:
            cline = mpl.cm.Blues(0.4,0.8)
        if cwalker is None:
            cwalker = mpl.cm.binary(0.7,0.2)
        
        # -- ploting            
        for i, name, fitted in zip(range(self.mcmc["ndim"]), self.freeparameters, self.mcmc["pos_init"] if truths is None else truths):
            ax = fig.add_subplot(self.mcmc["ndim"],1,i+1, ylabel=name)
            _ = ax.plot(np.arange(self.mcmc["nrun"]), self.mcmc["sampler"].chain.T[i],
                        color=cwalker,**kwargs)
            
            ax.axhline(fitted, color=cline, lw=2)

        fig.figout(savefile=savefile, show=show)

    # ================ #
    # = Properties   = #
    # ================ #
    @property
    def mcmc(self):
        """ dictionary containing the mcmc parameters """
        if self._derived_properties["mcmc"] is None:
            self._derived_properties["mcmc"] = {}
        return self._derived_properties["mcmc"]

    @property
    def mcmc_samples(self):
        """ the flatten samplers after burned in removal, see set_mcmc_samples """
        if not self.has_mcmc_ran():
            raise AttributeError("run mcmc first.")
        
        if "burnin" not in self.mcmc.keys():
            raise AttributeError("You did not specified the burnin value. see 'set_mcmc_burnin")
        
        return self.mcmc["sampler"].chain[:, self.mcmc["burnin"]:, :].reshape((-1, self.mcmc["ndim"]))
    
    def set_mcmc_burnin(self, burnin):
        """ """
        if burnin<0 or burnin>self.mcmc["nrun"]:
            raise ValueError("the mcmc burnin must be greater than 0 and lower than the amount of run.")
        
        self.mcmc["burnin"] = burnin
        
    def _set_mcmc_(self,mcmcdict):
        """ Advanced methods to avoid rerunning an existing mcmc """
        self._derived_properties["mcmc"] = mcmcdict
        
    def has_mcmc_ran(self):
        """ return True if you ran 'run_mcmc' """
        return "sampler" in self.mcmc.keys()


    


############################
#                          #
# Models                   #
#                          #
############################
def PhotoCalibrationModel( object ):
    """
    """
    def __init__(self):
        print "tobedone"


    
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


