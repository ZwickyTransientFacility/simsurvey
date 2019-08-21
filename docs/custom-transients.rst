*****************
Custom Transients
*****************

The use of simsurvey is not limited to the built-in SN types but can
be used to simulated the lightcurves of most types of non-periodic
extra-galactic optical transients. In the sections of this page the
components required for a ``TransientGenerator`` are described,
followed by an example that shows how the components fit
together. While all the listed components are required for creating a
generator for a completely custom model, they can be replaced on an
individual basis for the built-in models. Additional useful examples
can be found in the :doc:`additional-sources`.

sncosmo Model
=============

At the center of a ``TransientGenerator`` is an ``sncosmo.Model``
object, which will be uses to synthesize the photometry of the
transient. While these models are often based on time series of SEDs
that are interpolated that can practically be any function that
returns the model flux as a function of time and wavelength. For more
information on making a model, see the sncosmo documentation on
`Supernova Models
<https://sncosmo.readthedocs.io/en/latest/models.html>`_ and `Creating
a new Source class
<https://sncosmo.readthedocs.io/en/latest/examples/plot_custom_source.html>`_.

Note that by default the ``TransientGenerator`` will add a Milky Way
extinction effect to the model, hence it is not necessary to include
it yourself. This may however lead to problems if the simulated survey
uses filters outside the definiton range of the sncosmo implementation
of the CCM89 extinction law, e.g. for an IR survey (for which MW
extinction is negligible anyway). If this is the case, set the option
``apply_mwebv=False`` to disable this effect.

Volumetric Rate
===============

The number transients that will be simulated as well as their redshift
distribution is determined based on the volumetric rate as a function
of redshift. This function will accept only the redshift as input and
should return the rate in units of
:math:`\textrm{Mpc}^{-3}\textrm{yr}^{-1}`. The total number of
transients is determined by integrating the rate over the comoving
volume within the selected redshift range and multiplying by the
length of the survey. It is then adjusted for the fraction of the sky
that the survey covers. Lastly the number is assumed to be slightly
random and a random number drawn from a Gaussian with a standard
deviation equal to the square root of the "exact" number is added
before rounding the result.

In some simulations it may be useful to fix the number of transients
to a specific value. This can be done using the ``ntransient``
argument of the ``TransientGenerator``.

Luminosity Function and Lightcurve Parameter Distributions
==========================================================

The final component is funtion takes the list of redhsifts and the
``sncosmo.Model`` as input and needs to return a dictionary of arrays
that contain the lightcurve parameters for each transient. In most
cases this will include the ``amplitude`` parameter, which scales the
overall brightness of the model and thus needs to be based on an
absolute magnitude drawn from the luminosity function and adjust for
the luminosity distance, which can be calculated from the
redshift. Additional lightcurve parameters, such as extinction in the
host galaxy or the "stretch" parameter of SN Ia ligthcurve model, also
need to be determined for every simulated transient unless their
default values are acceptable for the simulation.

The function can further have any number of keyword argument, values
for which can be provided as a dictionary as shown in the example
below.

Example
=======

This example corresponds to the built-in SN Ia model based on the
Hsiao template, showing all steps required for it.

::

   import numpy as np
   import sncosmo
   import simsurvey
   from astropy.cosmology import Planck15

   # Make sure that the Hsiao template is downloaded
   # and then load it from its default location
   sncosmo.get_source('hsiao', version='3.0')
   p, w, f = sncosmo.read_griddata_fits(
       os.path.join(sncosmo.builtins.get_cache_dir(),
       'sncosmo/models/hsiao/Hsiao_SED_V3.fits')
   )

   # Define model using StretchSource and CCM89 extinction in rest-frame
   model = sncosmo.Model(
       source=sncosmo.StretchSource(p, w, f, name='hsiao-stretch'),
       effects=[sncosmo.CCM89Dust()],
       effect_names=['host'],
       effect_frames=['rest']
   )

   # Define the function that generates the lightcurve parameters
   # (Note how the value for mags is not in the typical range for SNe Ia.
   #  This will be fixed by passing new value 
   def random_parameters(self, redshifts, model,
                         mag=(-18., 1.),
                         r_v=2., ebv_rate=0.11,
                         alpha=1.3,
			 cosmo=Planck15
                         **kwargs):
        """
        """
        out = {}

	amp = []
	for z in redshifts:
	    model.set(z=z)
            model.set_source_peakabsmag(mabs, 'bessellb', 'vega', cosmo=cosmo)
            amp.append(model.get('amplitude'))

	out['amplitude'] = np.array(amp)
	out['hostr_v'] = r_v * np.ones(len(redshifts))
	out['hostebv'] =  np.random.exponential(ebv_rate, len(redshifts))
	    
        out['s'] = np.random.normal(1., 0.1, len(redshifts))
        out['amplitude'] *= 10 ** (0.4 * alpha * (out['s'] - 1))

        return out

   transientprop = {
       'lcmodel': model,
       'lcsimul_func': random_parameters,
       'lcsimul_prop': {'mag': (-19.3, 0.1)}
   }

   tr = simsurvey.get_transient_generator((0.0, 0.05),
                                          ratefunc=lambda z: 3e-5,
					  ra_range=(0,360),
                                          dec_range=(-30,90),
                                          mjd_range=(58178, 58543),
					  transientprop=transientprop)
