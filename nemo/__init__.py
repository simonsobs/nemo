"""

nemo - SZ cluster and source detection package

"""

from . import catalogs 
from . import photometry
from . import plotSettings
from . import startUp
from . import filters
from . import signals
from . import completeness
from . import maps
from . import pipelines

__all__ = ['startUp', 'catalogs', 'maps', 'filters', 'photometry', 
           'signals', 'completeness', 'pipelines', 'plotSettings']

__version__ = "git"
