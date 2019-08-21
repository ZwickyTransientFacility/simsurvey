************
Survey Plans
************

Survey plans store and organize any inofrmation regarding the survey
strategy and the telescope.  The minimum required information is a set
of lists for the pointings, which needs to contain the times, filters,
coordinates (RA and Dec, both in degress), and skynoise values. The
definition of skynoise is the same as for ``sncosmo.realize_lcs``,
which is used to generate the lightcurves.

::

   import simsurvey

   obs = {'time': [56176.19, 56188.254, 56207.172],
          'band': ['desg', 'desr', 'desi'],
          'skynoise': [1261.0, 1261.0, 1261.0],
	  'ra': [269.8, 269.8, 269.8],
	  'dec': [40.5, 40.5, 40.5]}
   
   plan = simsurvey.SurveyPlan(time=obs['time'],
                               band=obs['band'],
			       ra=obs['ra'],
			       dec=obs['dec'],
                               skynoise=obs['skynoise'])
			       
Note that each pair of RA and Dec will be treated separately for
determining which transients are within the field of view for this
pointing. Thus this is only a good option for quick simulations of a
short observation run. Recurring pointings should best be assign field
numbers (see below), even if the simulated survey does not have a
pre-defined field grid.

Using a Field Grid
==================

Instead of passing the individual coordinates for each observation, a
set of fields can be defined. Using such a field grid, will allow
simsurvey to check only once, which simulated transient will be
observed when pointing the telescope at a given set of coordinates. In
addition to RA and Dec, a field ID number can be listed. Otherwise the
fields will be numbered started with 0.

::

   obs = {'time': [56176.19, 56188.254, 56207.172],
          'band': ['desg', 'desr', 'desi'],
          'skynoise': [1261.0, 1261.0, 1261.0],
	  'field': [725, 725, 725]}

   fields = {'ra': [269.8, 278.4],
             'dec': [40.5, 40.5].
	     'field_id': [725, 726]}

   plan = simsurvey.SurveyPlan(time=obs['time'],
                               band=obs['band'],
			       obs_field=obs['field'],
                               skynoise=obs['skynoise'], 
                               fields=fields)


Changing the Field of View and Adding Comments
==============================================

By default field of view is set such that it encompasses the whole
field of view of ZTF. When simulation another survey the size of the
field of view can be changed using the ``width`` and ``height``
arguments. For matching which simulated transients are observed, a
rectangle with the given width and height (centered on the
coordinates) will be projected onto the sky.

For a more realistic simulation, the layout of the camera,
e.g. individual CCD chips or their read-out channels, can also be
passed to the ``SurveyPlan``. For this a list of focal plane
coordinates relative to the field center needs to be provided, with
each entry containing four sets of coordinates in degrees. The order
of these coordinates does not matter; they will be sorted in order to
form a concave quadrilateral. The field will be numbered starting
with 0. Unlike the fields, this numbering can currently not be
changed.

To reflect that each CCD will produce a separate image, which may have
a different depth, the list of observations can also be provided for
individual CCDs. In this case, not all CCDs need to have information
for each observation (e.g. because reference images were only
completed for some of the CCDs on a specific field).

Lastly each observation can be given a comment, which allows the user
to pass additional information, e.g. which subsurvey of a larger
survey the observation was taken for.

::

   obs = {'time': [56176.19, 56176.19, 56188.254, 56188.254, 56207.172, 56207.172],
          'band': ['desg', 'desg', 'desr', 'desr', 'desi', 'desi],
          'skynoise': [1261.0, 1444.0, 1261.0, 1444.0, 1261.0, 1444.0],
	  'field': [725, 725, 725, 725, 725, 725],
	  'ccd': [0, 1, 0, 1, 0, 1],
	  'comments': ['all-sky', 'all-sky', 'all-sky', 'all-sky', 'i-band', 'i-band']}

   fields = {'ra': [269.8, 278.4],
             'dec': [40.5, 40.5].
	     'field_id': [725, 726]}

   ccds = [[[-1., -1.],   [-0.1, -1.], [-0.1, -0.1], [-1, -0.1]],
           [[0.1, -1.],   [1., -1.],   [1., -0.1],   [0.1, -0.1]],
	   [[-1., 0.1],   [-0.1, 0.1], [-0.1, 1.],   [-1, 1.]],
	   [[0.1., 0.1.], [1., 0.1.],  [1., 1.],     [0.1, 1.]]]
   

   plan = simsurvey.SurveyPlan(time=obs['time'],
                               band=obs['band'],
			       obs_field=obs['field'],,
			       obs_ccd=obs['ccd']
                               skynoise=obs['skynoise'], 
                               comments=obs['comments'],
			       fields=fields,
			       ccds=ccds)
