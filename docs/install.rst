************
Installation
************

simsurvey works on Python 3.6+ and depends on the following Python packages:

- `numpy <http://www.numpy.org/>`_
- `scipy <http://www.scipy.org/>`_
- `astropy <http://www.astropy.org>`_
- `sncomso <https://sncosmo.readthedocs.io>`_
- `extinction <http://extinction.readthedocs.io>`_
- `sfdmap <https://github.com/kbarbary/sfdmap>`_
- `propobject <https://github.com/MickaelRigault/propobject>`_

Install using pip
=================

It is advisable to first insall sncosmo and its dependencies according
to its instructions. Then::

  pip install simsurvey

In order to simulate lightcurves with Milky Way extinction based on
the transients coordinates and the SFD98 dust maps, you will further
need to download the `sfddata <https://github.com/kbarbary/sfddata>`_
files and either set ``$SFD_DIR`` to their path or pass the path in
your script running simsurvey.

Install latest development version
==================================

simsurvey is being developed `on github
<https://github.com/ufeindt/simsurvey>`_. To get the latest
development version using ``git``::

    git clone git://github.com/ufeindt/simsurvey.git
    cd simsurvey

then::

    python setup.py install
