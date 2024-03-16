from .scConGraph import *

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata

package_name = "scConGraph"
__version__ = importlib_metadata.version(package_name)