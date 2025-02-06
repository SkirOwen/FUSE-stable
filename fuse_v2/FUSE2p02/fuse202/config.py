import inspect
import importlib
import os


def get_fuse202_dir():
	fuse202_module = importlib.import_module("fuse202")
	fuse202_dir = os.path.dirname(inspect.getabsfile(fuse202_module))
	return os.path.abspath(os.path.join(fuse202_dir, ".."))



