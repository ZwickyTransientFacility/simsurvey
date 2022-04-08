#! /usr/bin/env python
#
# Copyright (C) 2015-22 simsurvey Developers
import simsurvey.version as sv

DESCRIPTION = "simsurvey: Basic tools for simulating astronomical surveys"
LONG_DESCRIPTION = """\
Simulation tools for transient surveys (cadence, strategy).
"""

DISTNAME = 'simsurvey'
AUTHOR = 'Simsurvey Developers'
MAINTAINER = 'Ulrich Feindt, Ana Sagués Carracedo'
MAINTAINER_EMAIL = 'ufeindt@gmail.com, ana.sagues-carracedo@fysik.su.se'
URL = 'https://github.com/ZwickyTransientFacility/simsurvey/'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/ufeindt/ZwickyTransientFacility/tarball/'+sv.__VERSION__
VERSION = sv.__VERSION__

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup
    _has_setuptools = False

if __name__ == "__main__":
    if _has_setuptools:
        packages = find_packages()
    else:
        # This should be updated if new submodules are added
        packages = ['simsurvey', 'simsurvey.utils']

    requirements = [
        'propobject>=0.1.0',
        'astropy>=1.2.1',
        'sncosmo>=1.3.1',
        'sfdmap>=0.1.0',
    ]

    setup(name=DISTNAME,
          author=AUTHOR,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          packages=packages,
          classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 2.7',
              'Programming Language :: Python :: 3',
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'],
          install_requires=requirements,
    )
