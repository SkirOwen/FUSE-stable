"""Tests for the restart-file backup step.

FUSE keeps two copies of its resumable state: restart/ is rewritten as the run
progresses, and backup/ holds the previous good copy to fall back on if a crash
lands mid-write. run_fuse() inlined this dance twice and the copies had
drifted - the initial-population one backed up the *.p state pickles, the
search one did not, so a crash during the search left an incomplete backup.
"""
import os

import pytest

from fuse202.utils.file_ops import backup_restart_files


@pytest.fixture
def run_directory(tmp_path, monkeypatch):
	"""A run directory shaped the way run_fuse() sets one up."""
	(tmp_path / "restart").mkdir()
	(tmp_path / "backup").mkdir()
	(tmp_path / "structures").mkdir()
	(tmp_path / "restart" / "energies.p").write_bytes(b"restart-state")
	(tmp_path / "structures" / "S-0000001.cif").write_text("cif contents")
	(tmp_path / "moves.p").write_bytes(b"search-state")
	(tmp_path / "cubes.p").write_bytes(b"cell-shapes")
	monkeypatch.chdir(tmp_path)
	return tmp_path


def test_restart_state_and_structures_are_backed_up(run_directory):
	backup_restart_files()
	assert (run_directory / "backup" / "energies.p").read_bytes() == b"restart-state"
	assert (run_directory / "backup" / "structures" / "S-0000001.cif").read_text() == "cif contents"


def test_state_pickles_are_backed_up(run_directory):
	"""The drift this extraction fixed: the search loop skipped these, so a
	crash mid-search left backup/ without the search state."""
	backup_restart_files()
	assert (run_directory / "backup" / "moves.p").read_bytes() == b"search-state"
	assert (run_directory / "backup" / "cubes.p").read_bytes() == b"cell-shapes"


def test_returns_to_the_run_directory(run_directory):
	"""It navigates with os.chdir, so it must leave the caller where it found
	it - the calculators that run next all use relative paths."""
	before = os.getcwd()
	backup_restart_files()
	assert os.getcwd() == before


def test_is_idempotent_and_refreshes_the_backup(run_directory):
	"""Called before every relaxation, so it runs many times per job and must
	overwrite cleanly rather than accumulate or fail on the second pass."""
	backup_restart_files()
	(run_directory / "restart" / "energies.p").write_bytes(b"newer-state")
	backup_restart_files()
	assert (run_directory / "backup" / "energies.p").read_bytes() == b"newer-state"
