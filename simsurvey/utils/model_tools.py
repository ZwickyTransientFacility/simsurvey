
import numpy as np

def read_possis_file(filename):
    """Read in a spectral model created by POSSIS (1906.04205), as appropriate
       for injestion as a simsurvey.models.AngularTimeSeriesSource. Model
       grids can be found here: https://github.com/mbulla/kilonova_models.
    Parameters
    ----------
    filename : str
        Path to the POSSIS file
    Returns
    -------
    phase : `~numpy.ndarray`
        Phases in days.
    wave : `~numpy.ndarray`
        Wavelengths in Angstroms.
    flux : `~numpy.ndarray`
        Model spectral flux density in arbitrary units.
        Must have shape `(num_phases)`.
    cos_theta : `~numpy.ndarray`
        cosine of viewing angle
    """

    f = open(filename)
    lines = f.readlines()

    nobs = int(lines[0])
    nwave = float(lines[1])
    line3 = (lines[2]).split(' ')
    ntime = int(line3[0])
    t_i = float(line3[1])
    t_f = float(line3[2])

    cos_theta = np.linspace(0, 1, nobs)  # 11 viewing angles
    phase = np.linspace(t_i, t_f, ntime)  # epochs

    file_ = np.genfromtxt(filename, skip_header=3)

    wave = file_[0:int(nwave),0]
    flux = []
    for i in range(int(nobs)):
        flux.append(file_[i*int(nwave):i*int(nwave)+int(nwave),1:])
    flux = np.array(flux).T

    phase = np.linspace(t_i, t_f, len(flux.T[0][0]))  # epochs

    return phase, wave, cos_theta, flux
