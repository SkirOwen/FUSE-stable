from __future__ import annotations

import os
import shutil

from fuse202.config import get_fuse202_dir


def cache_solutions_dir():
	return guarantee_existence(os.path.join(cache_dir(), "solutions"))


def cache_dir():
	"""./cache"""
	return guarantee_existence(os.path.join(get_fuse202_dir(), "cache"))


def guarantee_existence(path: str) -> str:
	"""Function to guarantee the existence of a path, and returns its absolute path.

	Parameters
	----------
	path : str
		Path (in str) to guarantee the existence.

	Returns
	-------
	str
		The absolute path.
	"""
	os.makedirs(path, exist_ok=True)
	return os.path.abspath(path)


def guarantee_empty_existence(path: str) -> str:
	if os.path.exists(path):
		shutil.rmtree(str(path))
	os.makedirs(path, exist_ok=True)
	return os.path.realpath(path, strict=True)
