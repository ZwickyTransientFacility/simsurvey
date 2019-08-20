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
SNANA SEDs.

CompoundSource
==============

ExpandingBlackBody
==================

SpectralIndex
=============
