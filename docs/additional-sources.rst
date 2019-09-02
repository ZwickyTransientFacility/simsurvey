**************************
Additional sncosmo Sources
**************************

Several additional ``sncosmo.Source`` classes are included in
simsurvey. These classes simplify defining a transient model e.g. if a
general model is not available and classes like ``TimeSeriesSource``
are insufficient.

.. _multisource:

MultiSource
===========

This source allows loading multiple SEDs for ``TimeSeriesSource``
objects into a single source. The index of the source is added as a
parameter. An example of its use are the built-in sources using the
SNANA SEDs. These sources can be loaded as described in
:doc:`builtin-transients` but the examples below shows the individual
step to reproduce this.

::

   import numpy as np
   import sncosmo
   import simsurvey
   from astropy.cosmology import Planck15
   
   # Assume that sed_files contains the paths to the SEDs used for the sources
   data_sed = {'p': [], 'w': [], 'f': []}
   for sed_file in files:
       p_, w_, f_ = sncosmo.read_griddata_ascii(sed_file)
       data_sed['p'].append(p_)
       data_sed['w'].append(w_)
       data_sed['f'].append(f_)
   
   source = simsurvey.MultiSource(data_sed['p'], data_sed['w'], data_sed['f'])
   
   # Define model using MultiSource and CCM89 extinction in rest-frame
   model = sncosmo.Model(
       source=source,
       effects=[sncosmo.CCM89Dust()],
       effect_names=['host'],
       effect_frames=['rest']
   )
   
   # Define the function that generates the lightcurve parameters
   # (Note how the value for mags is not in the typical range for SNe Ia.
   #  This will be fixed by passing new value
   def random_parameters(redshifts, model,
                         mag=(-18., 1.),
                         r_v=2., ebv_rate=0.11,
                         cosmo=Planck15,
                         n_templates=len(model._source._model_flux),
                         **kwargs):
       """
       """
       out = {}
       if n_templates > 1:
           out['template_index'] = np.random.randint(0, n_templates, len(redshifts))
           model_kw = [{'template_index': k} for k in out['template_index']]
       elif n_templates == 1:
           model_kw = [{} for z in redshifts]
   
       # Amplitude
       amp = []
       for z, model_kw_ in zip(redshifts, model_kw):
           mabs = np.random.normal(mag[0], mag[1])
           model.set(z=z, **model_kw_)
           model.set_source_peakabsmag(mabs, 'bessellb', 'vega', cosmo=cosmo)
           amp.append(model.get('amplitude'))
   
       out['amplitude'] = np.array(amp)
       out['hostr_v'] = r_v * np.ones(len(redshifts))
       out['hostebv'] =  np.random.exponential(ebv_rate, len(redshifts))
   
       return out
   
   
   transientprop = {
       'lcmodel': model,
       'lcsimul_func': random_parameters,
       'lcsimul_prop': {'mag': (-16.75, 1.), 'n_templates': 2}
   }
   
   tr = simsurvey.get_transient_generator((0.0, 0.05),
                                          ratefunc=lambda z: 1.2e-4,
                                          ra_range=(0,360),
                                          dec_range=(-30,90),
                                          mjd_range=(58178, 58543),
                                          transientprop=transientprop)
   
CompoundSource
==============

The ``CompoundSource`` class also combines multiple source. However,
it is not keep them separate but instead returns the sum of the fluxes
of each individual source. This can be use to e.g. model the combined
flux of the images of a strongly lensed supernova or add an outburst
of the supernova progenitor at a random time before explosion.

When creating the ``CompoundSource`` a list
of ``sncsomo.Source`` objects needs to be provided. These sources need
not be of the same class and can have different sets of
parameters. The parameters of the ``CompoundSource`` will be the ones
of the individual sources with an index added to the name, e.g. for
two SALT2 sources it would be ``x0_0``, ``x1_0``, ``c_0``, ``x0_1``,
``x1_1``, and ``c_1``.  Furthermore the source classe permits a time
delay between the sources, i.e. the first source will be defined
around phase 0 (or time ``t0`` in an ``sncosmo.Model`` based on it)
and all additional sources get an additional parameter ``dt_1``,
``dt_2``, and so on.

For the example of a strongly lensed type Ia supernova with two
images, the source would then be defined as follows:

::

   import os
   import sncosmo
   from simsurvey import CompoundSource
   
   p, w, f = sncosmo.read_griddata_fits(
       os.path.join(sncosmo.builtins.get_cache_dir(),
       'sncosmo/models/hsiao/Hsiao_SED_V3.fits')
   )
   
   source = CompoundSource((sncosmo.StretchSource(p, w, f, name='hsiao-stretch1'),
                            sncosmo.StretchSource(p, w, f, name='hsiao-stretch2')))
      

ExpandingBlackBodySource
========================

This source is based on Planck's law for the radiation emitted by a
black body. Unlike the other sources typically used in sncosmo, it
does not use a series of SED that are interpolated. Instead the user
must provide the radius and temperature of a black body (e.g. the
photosphere of a transient) as functions of time. These functions must
accept a list of parameter and an array of time values as input. For
instance the model used in `Kasliwal et al. (2017)
<https://arxiv.org/abs/1710.05436>`_ to describe the optical
counterpart of GW170817 can be implemented like this:

::

   import numpy as np
   from simsurvey import ExpandingBlackBodySource

   source = ExpandingBlackBodySource(
       minphase=0.5, maxphase=15.,
       tempfunc=(lambda p, t: p[0] * t**p[1]),
       tempparam=(6050, -0.62),
       radiusfunc=(lambda p, t: p[0] * (1-np.exp(-p[1]*t) + p[2]*t)),
       radiusparam=(24000, 0.42, 2500)
   )

The parameters ``minphase`` and ``maxphase`` simply set the phases for
which the source is defined, where ``tempparam`` and ``radiusparam``
defined the default values for the parameters of the source and also
set the number of parameters required for ``tempfunc`` and
``radiusfunc``, which will be numbered starting with 0. In this case
for example the source will have the following parameters for
temperature and radius: ``T0``, ``T1``, ``R0``, ``R1``, and
``R2``. The radius needs to be given in units of solar radii adn
temperature is in Kelvin.  Additional the distance ``d`` to the
transient is a model parameter and is given in units of Mpc. It should
be set based on the redshift according to a cosmological model.
   
SpectralIndexSource
===================

The ``SpectralIndexSource`` class is for cases of tranients for which
only very limited information is available. it behaves very similarly
two the ``ExpandingBlackBodySource`` class (see above), in that it
accepts two functions of time that define the evolution of its
spectral form. The keyword arguments for the functions are
``fluxfunc`` and ``specfunc`` whereas their default parameters are
``fluxparam`` and ``specparam``, respectively. The former function
defines the monochromatic evolution of the source's flux, e.g. based
on a single-band lightcurve, and the latter defines the evolution of
the spectral index, i.e. the exponent :math:`\alpha` of the spectral
shape :math:`F(\lambda)\propto\lambda^\alpha`. Again the lists of
parameters set the numbers of required paremeters, which are named
``f0``, ``f1``, ``a0``, ``a1`` and so on. The units of these functions
can be arbitrary. For the overall normalization the parameter
``amplitude`` can be set (similar to setting the amplitude of the
``MultiSource`` in an example above).

