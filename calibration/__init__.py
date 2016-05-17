""" Simulate and retreive data assuming """


try:
    import gaia_anchor
except ImportError:
    import warnings
    warnings.warn("modefit not imported. the gaia_anchor tools are not available")

from iomock import *
