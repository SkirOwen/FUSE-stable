from __future__ import annotations

import glob
import os
import shutil


def clean_files(dir_path: str, file_ext: str) -> None:
	for file in glob.glob(os.path.join(dir_path, f"*.{file_ext}")):
		os.remove(file)


def backup_restart_files() -> None:
	"""Copy the restart state into backup/ before relaxing another structure.

	FUSE keeps two copies of its resumable state: restart/ is rewritten as the
	run progresses, and backup/ holds the previous good copy to fall back on if
	a crash lands mid-write. This is called before each relaxation so the
	fallback is never more than one structure behind.

	Operates on relative paths from the run directory, matching the original
	inline code - the caller is expected to be in the run directory, and the
	restart/, backup/ and structures/ directories are expected to exist.

	run_fuse() inlined this twice and the copies had drifted: the initial
	population backed up the *.p state pickles, the search loop did not, so a
	crash during the search left an incomplete backup. Both paths now back up
	everything.
	"""
	os.chdir("restart")
	for name in glob.glob("*"):
		shutil.copy(name, "../backup/" + name)
	os.chdir("../")

	os.chdir("backup")
	if not os.path.isdir("structures"):
		os.mkdir("structures")
	os.chdir("../")

	os.chdir("structures")
	for name in glob.glob("*.cif"):
		shutil.copy(name, "../backup/structures/" + name)
	os.chdir("../")

	# the unit cell shape / search state pickles
	for name in glob.glob("*.p"):
		shutil.copy(name, "backup/.")
