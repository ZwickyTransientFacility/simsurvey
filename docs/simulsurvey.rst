**************************
Generating the Lightcurves
**************************

Once the ``SurveyPlan`` and the ``TransientGenerator`` have been
created, they can be combined in a ``SimulSurvey``. Then the
lightcurves can be realized:

::
   
   survey = simsurvey.SimulSurvey(generator=tr, plan=plan)

   lcs = survey.get_lightcurves()

   # The lightcurves can further be saved as a pickle file
   lcs.save('lcs.pkl')

For a simulations up to redhsift 0.2, this process will typically take
ten minutes to one hour depending on the supernova type. Progress bars
can be shown by setting the keyword argument ``progress_bar=True``
(also set ``notebook=True`` if running a jupyter notebook).

Detection Criteria and Lightcurve Statistics
============================================

Since many simulated transients will be too faint to be detected by a
trypical survey, only lightcurves that pass a detection criterion will
be saved. By default two points with S/N > 5 are required. These
settings can be changed using the ``n_det`` and ``threshold`` keyword
argument, with the former specifying the required number of detections
and the latter the required S/N level.

Additionally, several statistics will be generated for each
lightcurve, e.g. the phase of detection, see
:doc:`lightcurve-collections`. Among those statistics are arrays with
the number of observations per filter in specific phase bins relative
to the model's ``t0`` parameter. By default the bins are set to
cover -30 to 70 days with a 5 day width but this can be changed by
passing the bin edges using the option ``p_bins``.
