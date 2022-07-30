*******************
Built-in Transients
*******************

Common supernova types can be simulated easily using built-in
settings.  The supported types are Ia, Ib/c, IIn, and IIP. Based on
the type, the volumetric rate and luminosity functions (as well as
distributions of color and stretch for type Ia SNe) are set to
realistic values. Each individual setting can be adjusted to the
user's liking, see :doc:`custom-transients` for details.

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
time scale is one year. The full list of options for ``transient`` and
``template``, as well as their associated distributions and parameters
are discussed below.

sncosmo Sources
===============

For built-in transients the models will be created using sources that
are built into sncosmo and will be downloaded at the first attempt of
loading them.

+---------------+--------------+---------------------+--------------------+--------------------------------+
| ``transient`` | ``template`` | sncosmo Source      | Source Class       | Notes                          |
+===============+==============+=====================+====================+================================+
| Ia            | salt2        | salt2               | SALT2Source        |                                |
+---------------+--------------+---------------------+--------------------+--------------------------------+
| Ia            | hsiao        | hsiao               | StretchSource      |                                |
+---------------+--------------+---------------------+--------------------+--------------------------------+
| Ibc           | nugent       | nugent-sn1bc        | TimeSeriesSource   |                                |
+---------------+--------------+---------------------+--------------------+--------------------------------+
| Ibc           | snana        | All 'snana' sources | :ref:`multisource` |                                |
|               |              | of types Ib and Ic  |                    |                                |
+---------------+--------------+---------------------+--------------------+--------------------------------+
| IIn           | nugent       | nugent-sn2n         | TimeSeriesSource   | Limited to first 150 days      |
|               |              |                     |                    | to avoid interpolations errors |
+---------------+--------------+---------------------+--------------------+--------------------------------+
| IIn           | snana        | All 'snana' sources | :ref:`multisource` |                                |
|               |              | of type IIn         |                    |                                |
+---------------+--------------+---------------------+--------------------+--------------------------------+
| IIP           | nugent       | nugent-sn2p         | TimeSeriesSource   | Limited to first 130 days      |
|               |              |                     |                    | to avoid interpolations errors |
+---------------+--------------+---------------------+--------------------+--------------------------------+
| IIP           | snana        | All 'snana' sources | :ref:`multisource` |                                |
|               |              | of type IIP         |                    |                                |
+---------------+--------------+---------------------+--------------------+--------------------------------+

For references to the SEDs used, see `the sncosmo documentation
<https://sncosmo.readthedocs.io/en/latest/source-list.html>`_.


Rates and Luminosity Functions
==============================

Based on the supernova type the volumetric rate is set and is then is
used to determine the number of SNe to be simulated as well as their
redshift distribution. Furthermore the absolute B-band magnitude at
peak of each individual supernova is set to a value drawn from a
Gaussian distribution centered on :math:`M_B` with a standard
deviation :math:`\sigma_M`. Lastly all models except SALT2 have had an
extinction effect at the transient redshift added, which simulates
extinction by dust in the host galaxy (SALT2 has its on color law for
that). For this effect :math:`R_V` is fixed to 2and :math:`E(B-V)` is
drawn from an exponential distribution with a scale of 0.11.

+---------+--------------------------------------------------+--------------------+----------------------+
| SN type | Rate [:math:`\textrm{Mpc}^{-3}\textrm{yr}^{-1}`] | :math:`M_B` (peak) | :math:`\sigma_M`     |
+=========+==================================================+====================+======================+
| Ia      | :math:`3 \times 10^{-5} (1+z)`                   | :math:`-19.3`      | :math:`0.1^\dagger`  |
+---------+--------------------------------------------------+--------------------+----------------------+
| Ib/c    | :math:`2.25 \times 10^{-5} (1+z)`                | :math:`-17.5`      | :math:`1.2`          |
+---------+--------------------------------------------------+--------------------+----------------------+
| IIn     | :math:`7.5 \times 10^{-6} (1+z)`                 | :math:`-18.5`      | :math:`1.4^\ddagger` |
+---------+--------------------------------------------------+--------------------+----------------------+
| IIP     | :math:`1.2 \times 10^{-4} (1+z)`                 | :math:`-16.75`     | :math:`1`            |
+---------+--------------------------------------------------+--------------------+----------------------+

:math:`\dagger` In addition to the intrinsic scatter of SN Ia peak
magnitudes, the Tripp relations [25] were used to simulate a realistic
population. In case of SALT2, the parameters :math:`x_1` and :math:`c`
are drawn from Gaussians centered around 0 with standard deviations 1
and 0.1, respectively. From the peak magnitude the code then subtracts
:math:`\alpha x_1 - \beta c` with :math:`\alpha=0.13` and
:math:`\beta=3`. In case of the Hsiao template, the "stretch" parameter
:math:`s` is drawn from a Gaussian centered around 1 with a standard
deviation of 0.1. From the peak magnitude the code then subtracts
:math:`\alpha (s-1)` with :math:`\alpha=1.3`.

:math:`\ddagger` To avoid simulating a large number of unrealistically
bright SNe IIn, the Gaussian distribution of peak magnitudes was
truncated at :math:`1\sigma` on the brighter-than-average side.
