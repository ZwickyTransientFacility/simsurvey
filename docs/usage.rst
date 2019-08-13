*****
Usage
*****

Two components are needed for simsurvey to simulate transient lightcurves:

* A ``SurveyPlan`` that contains the pointing schedule of the survey,
  including the times of observation, filter used, and a measure of
  the sky brightness limiting the observation. Furthermore the outline
  of the individual CCDs of the camera can be provided to simulate the
  losses due to chip gap.
* A ``TransientGenerator`` that generates and stores the distributions
  of transient lightcurve parameters, along with their times of
  explosion (or peak depending on the model used) and their
  coordinates. If sfdmap is installed correctly, it will also retrieve
  the E(B-V) values from the Schlegel, Finkbeiner & Davis (1998) maps.

SurveyPlan
==========

TransientGenerator
==================

Generating the Lightcurves
==========================

Reading the Output
==================
