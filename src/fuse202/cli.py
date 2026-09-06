"""Command line entry point for FUSE.

Deliberately small. A FUSE calculation is an ordinary Python script that calls
run_fuse() (see examples/), which is the convention in this field and works
well, so this does not try to wrap that. What it does provide is the one thing
a script cannot: a way to answer "is my installation actually set up?" without
starting a calculation and waiting to see what breaks.

A fuller CLI is planned for v3.
"""
from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import os
import shutil
import sys

# (import name, what it is, why you need it)
# (import name, distribution name, label, why you need it)
CORE_DEPENDENCIES = [
	("ase", "ase", "atomic simulation environment", "structure handling, all calculators"),
	("numpy", "numpy", "numpy", "everywhere"),
	("pandas", "pandas", "pandas", "reading ranked structure lists"),
	("spglib", "spglib", "spglib", "symmetry standardisation of cells"),
	("func_timeout", "func-timeout", "func-timeout", "timing out stuck GULP runs"),
	("matplotlib", "matplotlib", "matplotlib", "the summary plot at the end of a run"),
]

ML_DEPENDENCIES = [
	("chgnet", "chgnet", "CHGNet", "the ML potential (ctype='chgnet')"),
	("torch", "torch", "PyTorch", "required by CHGNet"),
	("pymatgen", "pymatgen", "pymatgen", "required by CHGNet"),
]

# (variable, what it is for, whether FUSE itself reads it)
ENVIRONMENT_VARIABLES = [
	("ASE_GULP_COMMAND", "how ase invokes GULP", True),
	("GULP_LIB", "GULP potential library directory (may be empty)", False),
	("VASP_PP_PATH", "VASP pseudopotential directory", False),
	("VASP_SCRIPT", "script ase uses to run VASP", False),
	("SPP_PATH", "unpacked statistical proxy potential library", True),
	("GNBOSS_TEMP", "unpacked gnboss template", True),
]

TICK, CROSS, DASH = "ok  ", "MISSING", "--  "


def _check_import(module_name: str, distribution_name: str = "") -> tuple[bool, str]:
	"""Import a module and report its version.

	The version comes from the installed distribution metadata rather than a
	__version__ attribute, because not every package exposes one (pymatgen does
	not) and the distribution name is not always the import name
	(func-timeout / func_timeout).
	"""
	try:
		importlib.import_module(module_name)
	except Exception as error:  # ImportError, but a broken install can raise others
		return False, str(error).split("\n")[0][:60]
	try:
		return True, importlib.metadata.version(distribution_name or module_name)
	except importlib.metadata.PackageNotFoundError:
		return True, ""


def check(argv=None) -> int:
	"""Report what is installed and what FUSE can therefore run.

	Exit status is 0 when the core installation is usable, 1 when something the
	package cannot work without is missing. A missing external chemistry code is
	*not* an error - plenty of people only want CHGNet - so those are reported
	but do not affect the status.
	"""
	print("FUSE installation check")
	print("=" * 60)

	problems = []

	# the package itself
	ok, detail = _check_import("fuse202", "fuse202")
	if ok:
		print(f"[{TICK}] fuse202 {detail}")
	else:
		print(f"[{CROSS}] fuse202 could not be imported: {detail}")
		problems.append("the fuse202 package itself could not be imported")

	print(f"\nCore dependencies (needed for any run)")
	for module_name, distribution, label, purpose in CORE_DEPENDENCIES:
		ok, detail = _check_import(module_name, distribution)
		if ok:
			print(f"  [{TICK}] {label:<32} {detail}")
		else:
			print(f"  [{CROSS}] {label:<32} needed for: {purpose}")
			problems.append(f"{label} is not installed")

	print(f"\nML potential (optional - install with: uv sync --extra ml)")
	ml_ok = True
	for module_name, distribution, label, purpose in ML_DEPENDENCIES:
		ok, detail = _check_import(module_name, distribution)
		print(f"  [{TICK if ok else DASH}] {label:<32} {detail if ok else 'not installed'}")
		ml_ok &= ok

	print(f"\nExternal chemistry codes (optional - each enables one calculator)")
	gulp_command = os.environ.get("ASE_GULP_COMMAND", "")
	gulp_binary = shutil.which("gulp")
	gulp_ok = bool(gulp_command or gulp_binary)
	print(f"  [{TICK if gulp_ok else DASH}] GULP                             "
	      f"{'found' if gulp_ok else 'not found (set ASE_GULP_COMMAND)'}")
	vasp_ok = bool(os.environ.get("VASP_PP_PATH") and os.environ.get("VASP_SCRIPT"))
	print(f"  [{TICK if vasp_ok else DASH}] VASP                             "
	      f"{'configured' if vasp_ok else 'not configured (set VASP_PP_PATH, VASP_SCRIPT)'}")

	print(f"\nEnvironment variables")
	for name, purpose, read_by_fuse in ENVIRONMENT_VARIABLES:
		value = os.environ.get(name)
		marker = TICK if value else DASH
		shown = (value[:34] + "...") if value and len(value) > 37 else (value or f"unset - {purpose}")
		print(f"  [{marker}] {name:<20} {shown}")

	# what does all of that add up to
	print("\n" + "=" * 60)
	runnable = []
	if ml_ok:
		runnable.append("chgnet")
	if gulp_ok:
		runnable.append("gulp")
	if vasp_ok:
		runnable.append("vasp")
	if gulp_ok:
		runnable.append("mixed (gulp-based)")

	if problems:
		print("Installation is NOT usable:")
		for problem in problems:
			print(f"  - {problem}")
		print("\nTry: uv sync --extra ml")
		return 1

	if runnable:
		print(f"Ready to run with ctype: {', '.join(runnable)}")
	else:
		print("Core install is fine, but no calculator is available yet.")
		print("The quickest way to get a working run is the ML potential:")
		print("  uv sync --extra ml")

	if ml_ok:
		print("\nTry a real run now (takes a few seconds, needs no external code):")
		print("  python examples/quickstart_chgnet.py")
	return 0


def main(argv=None) -> int:
	parser = argparse.ArgumentParser(
		prog="fuse202",
		description="Flexible Unit Structure Engine. To run a calculation, write a "
		            "Python script that calls run_fuse() - see examples/.",
	)
	subparsers = parser.add_subparsers(dest="command")
	subparsers.add_parser("check", help="report what is installed and what FUSE can run")

	args = parser.parse_args(argv)
	if args.command == "check":
		return check()
	parser.print_help()
	return 0


if __name__ == "__main__":
	sys.exit(main())
