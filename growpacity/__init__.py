from importlib import metadata as _md

__name__ = "growpacity"
__version__ = _md.version("growpacity")

# importing the c-extension
import growpacity_c
from .growpacity import *
