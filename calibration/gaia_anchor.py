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

# - modefit dependency
from modefit.fitter.baseobjects import BaseFitter,BaseModel
# - local dependency
from .iomock import read_mockdata




__all__ = ["get_gaiafit"]

def get_gaiafit(gaiamock, modelname="Basic", **kwargs):
    """ load a structure to fit the given Model on the mockdata.

    Parameters
    ----------
    gaiamock: [string]
        Datafile of the mock data. Some example provided in ./data

    modelname: [string] - default 'Basic' -
        name of the model you wish to use to fit the data.

    **kwargs goes to GaiaFitter

    Return
    ------
    GaiaFitter object.
    """
    return GaiaFitter(gaiamock, modelname=modelname, **kwargs)




# ========================== #
#                            #
#     Gaia Model             #
#                            #
# ========================== #
"""
# Create a new Model.

## How to call it?
You can implement any model named e.g. ModelNewOne. This model will then
be accessible from the GaiaFitter using the keywork 'NewOne'. (`modelname='NewOne'`)

## What must be its structure?
There is only a few methods and parameter you need to implement:

* `FREEPARAMETERS` the list of name of *all* the parameters

* `setup` the methods that read the parameter for the model.
   Its argument must be `parameters`

* `get_loglikelihood` the must return the loglikelihood (could be -0.5*chi2) of the
   model given the data (the parameters has been set to the model before).
   Its argument must be `data`: the `GaiaFitter.data`
   
* `_minuit_chi2_` a stupid function that explicitly takes in argument *all* the
  freeparameters. the function feed them to parameters and returns 
  self.get_chi2(parameter).
  For example, if FREEPARAMETERS = ["a","b","toto"]
  then add this method:
  ```python
  def _minuit_chi2_(self,a,b,toto):
      parameter = a,b,toto
      return self.get_chi2(parameter)
   ```
[optional]
* `lnprior` (optional but better): method that measure the log of the prior
   value for the parameters. This methods could be as trivial as a flat
   'no-informatio' prior:
   ```python
   def lnprior(self,parameters):
       return 0
    ```
[optional]
* `get_set_param_input`: method with no argument that returns a dictionary
  containing part or the totality of the parameter fit properties:
   * x`_guess`
   * x`_boundaries`
   * x`_fixed`
  where x is any of the freeparameters.
  For example, if FREEPARAMETERS = ["a","b","toto"]
  and you know `b` have to be positive around, say 2, you could do somethings like:
  ```python
  def get_set_param_input(self):
      return {"b_guess":2,"b_boundaries":[0,None]}
  ```
  You can overwrite anything when calling the `fit()` method. Its kwargs does just so
  and overwrite anything returned by this default value setter.
  
## Done
"""
class ModelBasic( BaseModel ):
    """
    """
    # -- unrequested if empty
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    FREEPARAMETERS  = ["zeropoint", "a_gaiaGmB","a_gaiaBmR",
                       "a_airmass","a_airmass_gaiaGmB","a_airmass_gaiaBmR"]

    # ================= #
    # = Method        = #
    # ================= #
    def setup(self,parameters):
        """
        """
        self._parameters = parameters
                
    def get_model(self,data):
        """ return the expected magnitude correction based on the input parameters """

        zeropoint, a_gaiaGmB,a_gaiaBmR,a_airmass,a_airmass_gaiaGmB,a_airmass_gaiaBmR = self._parameters
        
        return zeropoint +  \
          a_gaiaGmB*data["GPmBP_gaia"] + a_gaiaBmR*data["BPmRP_gaia"] + \
          a_airmass*data['airmass']+ \
          a_airmass_gaiaGmB*data['airmass']*data["GPmBP_gaia"] +\
           a_airmass_gaiaBmR*data['airmass']*data["BPmRP_gaia"]    

          
    def get_loglikelihood(self,data):
        """ """
        res = (self.get_model(data) + data["BP_gaia"] - data["ZTF_mag"]) / data["ZTF_mag_err"]
        
        return -0.5*np.sum(res**2)

    def lnprior(self,parameters):
        """ Overwritting the priors, flat here"""
        return 0
    
    # ----------------------- #
    # - Model Particularity - #
    # ----------------------- #
    def _minuit_chi2_(self,zeropoint, a_gaiaGmB,a_gaiaBmR,
                       a_airmass,a_airmass_gaiaGmB,a_airmass_gaiaBmR):
        """
        """
        parameter = zeropoint, a_gaiaGmB,a_gaiaBmR,a_airmass,a_airmass_gaiaGmB,a_airmass_gaiaBmR
        return self.get_chi2(parameter)



# ================================= #
#                                   #
#   Gaia Data Structure Base        #
#                                   #
# ================================= #
class GaiaFitter( BaseFitter ):
    """
    """
    __nature__ = "GaiaFitter"

    PROPERTIES         = ["data"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    def __init__(self,filename, modelname="Basic", select_ztfband=None, select_maxwidth=None ):
        """
        """
        self.__build__()

        self.set_data(filename, select_ztfband=None, select_maxwidth=None)
        self.set_model(eval("Model%s()"%modelname))
    
    # ================= #
    # = Main          = #
    # ================= #
    def set_data(self,filename, select_ztfband=None, select_maxwidth=None):
        """ read the mock data in the given filename
        and create the data base"""
        # --------------------
        # -- Create the data
        dataframe = read_mockdata(filename)
        if select_ztfband is not None:
            dataframe = dataframe.mask( ~(dataframe["ZTF_filter"]==select_ztfband)).dropna(axis=0,how="all")
        if select_maxwidth is not None:
            dataframe = dataframe.mask( ~(dataframe["x"]<=select_maxwidth) | ~(dataframe["y"]<=select_maxwidth)).dropna(axis=0,how="all")
        
        self._properties["data"] = dataframe
        # --------------------
        # -- and the Derived Colors
        self.data["GPmBP_gaia"]     = self.data["G_gaia"] - self.data["BP_gaia"]
        self.data["GPmBP_gaia_err"] = np.sqrt(self.data["G_gaia_err"]**2 + self.data["BP_gaia_err"]**2) # CAUTION NO COV TERM HERE
        self.data["BPmRP_gaia"]     = self.data["BP_gaia"] - self.data["RP_gaia"]
        self.data["BPmRP_gaia_err"] = np.sqrt(self.data["BP_gaia_err"]**2 + self.data["RP_gaia_err"]**2) # CAUTION NO COV TERM HERE


    def get_modelchi2(self,parameters):
        """ Parses the parameters and return the associated -2 log Likelihood
        Both the parser and the log Likelohood functions belongs to the
        model.
        This should usually be passed to the model with loading it. 
        
        parameters: [array]
            a list of parameter as they could be understood
            by self.model.setup to setup the current model.
                                   
        Returns
        -------
        float (chi2 defines as -2*log_likelihood)
        """
        self.model.setup(parameters)
        return -2*self.model.get_loglikelihood(self.data)

    def get_fitresidual(self):
        """ Return the residuals """
        model = self.model.get_model(self.data)
        return model - self.data["ZTF_mag"]

        

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




# ================================= #
#                                   #
#   Internal Functions              #
#                                   #
# ================================= #

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


