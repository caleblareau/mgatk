from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("mgatk")
except PackageNotFoundError:  # running from a source tree without an install
    __version__ = "0.8.0"
