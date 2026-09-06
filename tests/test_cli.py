"""Tests for the `fuse202 check` command.

The command exists to answer "is my installation set up?" without starting a
calculation and waiting to see what breaks. What matters is that it
distinguishes a *broken* install (missing something FUSE cannot work without)
from a merely *incomplete* one (no external chemistry code, which is fine if
you only want CHGNet) - the first is an error, the second is not.
"""
import pytest

from fuse202.cli import check, main


def test_check_passes_on_a_working_install(capsys):
	assert check() == 0
	output = capsys.readouterr().out
	assert "FUSE installation check" in output
	assert "atomic simulation environment" in output


def test_missing_core_dependency_is_an_error(capsys, monkeypatch):
	"""spglib is not optional - without it FUSE cannot standardise cells.

	Patches importlib.import_module rather than builtins.__import__: by the time
	the suite runs, spglib is already in sys.modules, and import_module returns
	the cached module without going through __import__ at all.
	"""
	import importlib as importlib_module
	real_import_module = importlib_module.import_module

	def blocked(name, *args, **kwargs):
		if name == "spglib":
			raise ImportError("No module named spglib")
		return real_import_module(name, *args, **kwargs)

	monkeypatch.setattr("fuse202.cli.importlib.import_module", blocked)
	assert check() == 1
	output = capsys.readouterr().out
	assert "spglib is not installed" in output
	assert "NOT usable" in output


def test_no_external_calculator_is_not_an_error(capsys, monkeypatch):
	"""Having neither GULP nor VASP is a perfectly normal state - CHGNet alone
	is enough to run - so it must not be reported as a broken install."""
	for variable in ("ASE_GULP_COMMAND", "VASP_PP_PATH", "VASP_SCRIPT"):
		monkeypatch.delenv(variable, raising=False)
	monkeypatch.setattr("shutil.which", lambda name: None)

	assert check() == 0
	output = capsys.readouterr().out
	assert "NOT usable" not in output


def test_configured_gulp_is_reported_as_runnable(capsys, monkeypatch):
	monkeypatch.setenv("ASE_GULP_COMMAND", "gulp < gulp.gin > gulp.got")
	assert check() == 0
	output = capsys.readouterr().out
	assert "gulp" in output.split("Ready to run with ctype:")[1]


def test_bare_command_shows_help(capsys):
	"""A FUSE run is a Python script, not a subcommand - the help should point
	there rather than implying the CLI runs calculations."""
	assert main([]) == 0
	output = capsys.readouterr().out
	assert "run_fuse()" in output
