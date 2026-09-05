import inspect
import importlib
import os


def get_fuse202_dir():
	"""Returns the project root (the directory containing pyproject.toml), i.e.
	two levels above this file: fuse202/ lives under src/, which lives under
	the project root."""
	fuse202_module = importlib.import_module("fuse202")
	fuse202_dir = os.path.dirname(inspect.getabsfile(fuse202_module))
	return os.path.abspath(os.path.join(fuse202_dir, "..", ".."))



