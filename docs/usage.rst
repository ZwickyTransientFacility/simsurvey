*****
Usage
*****

Two components are needed for simsurvey to simulate transient lightcurves:

* A ``SurveyPlan`` that contains the pointing schedule of the survey,
  including the times of observation, filter used, and a measure of
  the sky brightness limiting the observation. Furthermore the outline
  of the individual CCDs of the camera can be provided to simulate the
  losses due to chip gaps.
* A ``TransientGenerator`` that generates and stores the distributions
  of transient lightcurve parameters, along with their times of
  explosion (or peak depending on the model used) and their
  coordinates. If sfdmap is installed correctly, it will also retrieve
  the E(B-V) values from the Schlegel, Finkbeiner & Davis (1998) maps.

SurveyPlan
==========

Survey plans store all information about the survey strategy and the
telescope. The minimum required information is a set of lists for the
pointings, which needs to contain the times, filters, coordinates (RA,
Dec), and skynoise values. The definition of skynoise is the same as
for ``sncosmo.realize_lcs``, which is used to generate the
lightcurves.

Instead of giving individual coordinates for each pointing, it is
better to define a field grid (or assign repeated pointings field
numbers if there is no grid). In this case the coordinates need not be
passed to the ``SurveyPlan`` as a long list but can instead be
provided through a dictionary containing lists of the coordinates of
the field centers and (optionally) field numbers (in case numbering
does not start at 0).

The following example assumes that ``obs`` is an ``astropy.Table``
with the columns ``'time'``, ``'band'``, ``'field'``, and
``'skynoise'`` and that ``fields`` is dictionary containing list for
the keys ``'ra'`` and ``'dec'``.

::

   import simsurvey
   
   plan = simsurvey.SurveyPlan(time=obs['time'],
                               band=obs['band'],
			       obs_field=obs['field'],
                               skynoise=obs['skynoise'], 
                               fields=fields)

For more detailed information, see :doc:`surveyplan`.

TransientGenerator
==================

Transient generators for several supernova types have been built into
simsurvey. Based on the type, the volumetric rate and luminosity
functions (as well as distributions of color and stretch for type Ia
SNe) are set to realistic values. Each individual setting can be
adjusted to the user's liking. For more information on the built-in
transient type, see :doc:`builtin-transients`.

The minimum information required to create a ``TransientGenerator``
for a built-in transient type is the redshift range, a range of RA and
Dec values that will cover the full survey footprint, and a range of
MJD values (or any other times in units of days) that covers the whole
duration of the survey. The range of MJDs is the range for the ``t0``
value of the ``sncosmo.Model`` used to simulate the lightcurve,
i.e. typically the time of explosion or the time of peak, thus the
range should be set such that transients that would only be observed
while they are fading when the survey starts or while they are rising
and the survey is ending.

To create the ``TransientGenerator`` use the following command:

::
   
    tr = simsurvey.get_transient_generator((0.0, 0.05),
                                           transient='Ia',
                                           template='salt2',
					   ra_range=(0,360),
					   dec_range=(-30,90),
                                           mjd_range=(58178, 58543))

The resulting ``TransientGenerator`` will simulate type Ia supernovae
at redshifts between 0 and 0.05 using the SALT2 template. The survey
area covers the whole sky down to a declination of 30 degrees and the
time scale is one year.

Furthermore transient generators can be created for any type of real
or hypothtical extra-galactic transient if the temporal evolution of
its spectral energy density (SED) can be defined in an
``sncosmo.Source`` object. simsurvey includes several custom source
classes that make this easier. For more information, see
:doc:`custom-transients`.

Generating the Lightcurves
==========================

Once the ``SurveyPlan`` and the ``TransientGenerator`` have been
created, they can be combined in a ``SimulSurvey``. Then the
lightcurves can be realized:

::
   
   survey = simsurvey.SimulSurvey(generator=tr, plan=plan)

   lcs = survey.get_lightcurves()

For a simulations up to redhsift 0.2, this process will typically take
ten minutes to one hour depending on the supernova type. Progress bars
can be shown by setting the keyword argument ``progress_bar=True``
(also set ``notebook=True`` if running a jupyter notebook).

Since many simulated transients will be too faint to be detected by a
trypical survey, only lightcurves that pass a detection criterion will
be saved. By default two points with S/N > 5 are required, for more
information on modifying this criterion see :doc:`simulsurvey`.


Reading the Output
==================

The final output of the simulation is a ``LightcurveCollection``
object. Essentially it behaves like a list of ``astropy.Table``
objects, i.e. ``lcs[0]`` is the first lightcurve that passed the
detection criterion and ``lcs[0].meta`` contains the parameters for
the simulation as well as basic statistics of the transient, e.g. its
phase of detection. For easy selections of subsamples, the object
contains dictionaries of ``numpy.arrays`` that list the simulation
parameters and lightcurve statistics, e.g. ``lcs.meta['z']`` for the
redshifts or ``lcs.stats['p_det']`` for the phases of detection. For
more information, see :doc:`lightcurve-collections`.
