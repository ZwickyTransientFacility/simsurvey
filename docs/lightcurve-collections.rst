**********************
Lightcurve Collections
**********************

The final output of the simulation is a ``LightcurveCollection``
object. Essentially it behaves like a list of ``astropy.Table``
objects, i.e. ``lcs[0]`` is the first lightcurve that passed the
detection criterion and ``lcs[0].meta`` contains the parameters for
the simulation as well as basic statistics of the transient, e.g. its
phase of detection. Here we will assume that lightcurves have
previously been generated and save as ``lcs.pkl``.

::

   In [1]: import simsurvey

   In [2]: lcs = simsurvey.LightcurveCollection(load='lcs.pkl')

   In [3]: lcs[0]
   Out[3]: 
   <Table length=18>
          time        band         flux             fluxerr          zp   zpsys field  ccd  comment
        float64       str4       float64            float64       float64  str2 int64 int64   str1 
   ------------------ ---- ------------------- ------------------ ------- ----- ----- ----- -------
   2458560.6563542006 ztfg  -170.8542867234545  5844.653315573051    30.0    ab   698    30        
      2458572.6273264 ztfr  13762.603125898202  8445.064350412049    30.0    ab   698    30        
      2458572.6323264 ztfr  10540.091531246484 10088.720419990608    30.0    ab   698    30        
      2458572.6373032 ztfr  29905.462158596303  19422.11566646984    30.0    ab   698    30        
   2458572.6422917005 ztfr   37152.11891220342  21244.33148986577    30.0    ab   698    30        
      2458575.6350926 ztfr -7988.6485422021615 31512.163787343936    30.0    ab   698    30        
      2458575.6395023 ztfr  27495.696490342685  37132.35617399112    30.0    ab   698    30        
      2458575.6439236 ztfr  53363.018490809874  45694.64685905065    30.0    ab   698    30        
      2458636.9612269 ztfr   38498.82429972523  2641.457143995491    30.0    ab   697    19        
      2458636.9616782 ztfr    37560.9124104196  2544.949781696994    30.0    ab   697    19        
   2458636.9635417005 ztfr  35801.486813393385  3258.595391324462    30.0    ab   697    19        
      2458639.9581944 ztfr  35013.732048571605  2356.152931435486    30.0    ab   697    19        
      2458639.9610069 ztfr   32697.02664079166 2505.4176771406887    30.0    ab   697    19        
      2458639.9684954 ztfr  30639.904113769848   2282.97758546629    30.0    ab   697    19        
      2458642.9571875 ztfr   30693.38412720665  2374.852856081374    30.0    ab   697    19        
      2458642.9609838 ztfr  28717.113554489457  2434.926783974113    30.0    ab   697    19        
      2458642.9643519 ztfr   31660.29677394924 2465.9182826338188    30.0    ab   697    19        
      2458642.9677199 ztfr  31685.250054557826 2412.9192082315303    30.0    ab   697    19

   In [4]: lcs[0].meta
   Out[4]: 
   {'z': 0.040268437228642796,
    't0': 2458593.8428128595,
    'x0': 0.0032943469122044553,
    'x1': 1.5399498609092175,
    'c': -0.053627722692004026,
    'mwebv': 0.05216634585651255,
    'ra': 33.126855403133355,
    'dec': 38.75016731904456,
    'mwebv_sfd98': 0.0475812059382565,
    'idx_orig': 0,
    'stats': {'p_det': 43.118865340482444,
     'p_last': 49.124907040502876,
     'dt_det': 61.31730329990387,
     'p_binned': {'all': array([0, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 7, 0, 0, 0, 0]),
      'ztfg': array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
             0., 0., 0.]),
      'ztfr': array([0., 4., 3., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3., 7., 0.,
             0., 0., 0.]),
      'ztfi': array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
             0., 0., 0.])},
     'mag_max': {'ztfg': 99.0, 'ztfr': 18.18189903108589, 'ztfi': 99.0}}}


The first items in the meta dictionary are the parameters of the
``sncosmo.Model`` followed by the coordinates. ``mwebv_sfd98`` is the
E(B-V) value at the transient's coordinates according to the SFD98
dust maps, which was automatically read from them. This is not the
value used in the simulation, however, which is ``mwebv``. In order to
simulate the uncertainty of the dust maps, Gaussian noise is added to
the E(B-V) values (:math:`\sigma=1\textrm{mag}`). ``idx_orig`` is the
original index of the transient in the ``TransientGenerator``.

The dictionary ``stats`` contains additional information about the
lightcurve:

- ``p_det``: time of detection relative to ``t0`` (in observer frame),
  specifically the time of the last observation required for the
  discovery criterion

- ``p_last``: time of the last observation relative to ``t0`` (in
  observer frame, no S/N requirement)
  
- ``dt_det``: time between the last observation without a detection
  before discovery and the first detection (both defined by the S/N
  threshold)
  
- ``p_binned``: arrays of the number of observations per filter in
  specific phase bins relative to the model's ``t0`` parameter (bin
  defintions in ``lcs.p_bins``)
  
- ``mag_max``: brightest observed magnitude for each filter (may be
  based on very low S/N point, 99.0 means observations in that filter)

Arrays containing information for all lightcurve can also be found in
the dictionaries ``lcs.meta`` and ``lcs.stats``. These contain only
information about the transients that passed the detection
criterion. The meta information of rejected lightcurves is stored in
``lcs.meta_rejected`` and the complete list of parameters can be
accessed via ``lcs.meta_full``.  This full list of parameters can the
e.g. be used to determine the survey efficiency. Note that these lists
only contain the parameters of transients, for which a lightcurves
were simulated. All other simulation parameters, i.e. those of
transients that did not lie in an observed region of the sky or that
were not bright when their fields were observed, can be found in
`lcs.meta_notobserved`.
  
Filtering the Lightcurves
=========================

In some cases it may be interesting to compare survey performance
based on small changes without rerunning the whole simulation,
e.g. assessing how many transients pass the detection criterion with
just one of the filters. To do this, a function that takes a
lightcurve is input and returns it can be passed to the function
``lcs.filter()``:

::

   def filter_function(lc):
       return lc[lc['band'] == 'ztfr']

   lcs_r = lcs.filter(filter_function)

Furthermore this function can also be used to just adjust the
detection criterion using the arguments ``n_det`` and ``threshold``
(required number of detections and S/N requirement for a detection,
respectively), as well a changing the phase bins ``p_bins`` used when
determining the lightcurve statistics.
