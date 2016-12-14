#! /usr/bin/env python
#
# Copyright (C) 2015 SNSurvey Developers

DESCRIPTION = "simsurvey: Basic tools for simulating astronomical surveys"
LONG_DESCRIPTION = """\
Simulation tools for transient surveys (cadence, strategy).
"""

DISTNAME = 'simsurvey'
AUTHOR = 'Simsurvey Developers'
MAINTAINER = 'Ulrich Feindt' 
MAINTAINER_EMAIL = 'ulrich.feindt@fysik.su.se'
URL = 'https://github.com/MickaelRigault/simsurvey/'
LICENSE = 'BSD (3-clause)'
DOWNLOAD_URL = 'https://github.com/MickaelRigault/simsurvey/tarball/0.2.1'
VERSION = '0.2.1'

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup

def check_dependencies():
    install_requires = []

    # Just make sure dependencies exist, I haven't rigorously
    # tested what the minimal versions that will work are
    # (help on that would be awesome)
    try:
        import astropy
    except ImportError:
        install_requires.append('astropy')
    try:
        import sncosmo
    except ImportError:
        install_requires.append('sncosmo')
    try:
        import propobject
    except ImportError:
        install_requires.append('propobject')
    try:
        import sfdmap
    except ImportError:
        install_requires.append('sfdmap')
        
    return install_requires

if __name__ == "__main__":

    install_requires = check_dependencies()

    if _has_setuptools:
        packages = find_packages()
    else:
        # This should be updated if new submodules are added
        packages = ['simsurvey', 'simsurvey.utils']

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
          install_requires=install_requires,
          packages=packages,
          classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 2.7',
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'],
      )
