"""Single point at which FUSE hands a structure to an energy calculator.

Every backend shares one contract: given an Atoms object, return
`(atoms, energy, converged)`. Choosing between them is therefore a lookup
rather than anything structural, and this module is the one place that lookup
happens.
"""
from __future__ import annotations

from fuse202.calculators.run_chgnet import run_chgnet
from fuse202.calculators.run_gulp import run_gulp
from fuse202.calculators.run_multiple_calculators import run_calculators
from fuse202.calculators.run_qe import run_qe
from fuse202.calculators.run_vasp import run_vasp

# Energy assigned to a structure whose calculation failed. Large enough that
# the search will never prefer it.
FAILED_ENERGY = 1.e20

SUPPORTED_CALCULATORS = ('gulp', 'vasp', 'qe', 'chgnet', 'mixed')


def run_calculator(
		ctype,
		atoms,
		*,
		# GULP
		shel=None,
		kwds='',
		gulp_opts='',
		lib='',
		gulp_command='gulp < gulp.gin > gulp.got',
		gulp_timeout='',
		# VASP
		vasp_opts='',
		kcut='',
		dist_cutoff=1.0,
		# Quantum Espresso
		qe_opts='',
		# CHGNet
		n_opts=2,
		rel=None,
		relaxer_opts=None,
		opt_class=None,
		opt_device='cpu',
		mode='relax',
		use_spglib=True,
		# "mixed": the ordered list of backends to run
		calcs='',
) -> tuple:
	"""Relax `atoms` with the calculator named by `ctype`.

	A structure search is expected to propose structures that make calculators
	fail, so any backend exception rejects that one structure rather than
	aborting the run. An unrecognised `ctype` is treated the same way.

	Parameters
	----------
	ctype : str
		Which calculator to use. One of SUPPORTED_CALCULATORS.
	atoms : ase.Atoms
		The structure to relax.

	Returns
	-------
	tuple of (ase.Atoms, float, bool)
		The relaxed structure, its energy, and whether the calculation
		converged. On any failure the atoms are returned unchanged with
		FAILED_ENERGY and converged False.
	"""
	try:
		if ctype == 'gulp':
			return run_gulp(
				atoms=atoms, shel=shel, kwds=kwds, opts=gulp_opts, lib=lib,
				produce_steps=False, gulp_command=gulp_command, gulp_timeout=gulp_timeout,
			)

		if ctype == 'vasp':
			return run_vasp(
				atoms=atoms, vasp_opts=vasp_opts, kcut=kcut, produce_steps=False,
				dist_cutoff=dist_cutoff,
			)

		if ctype == 'qe':
			return run_qe(atoms=atoms, qe_opts=qe_opts, kcut=kcut, produce_steps=False)

		if ctype == 'chgnet':
			return run_chgnet(
				atoms, n_opts=n_opts, rel=rel, relaxer_opts=relaxer_opts,
				opt_class=opt_class, opt_device=opt_device, use_spglib=use_spglib,
				mode=mode,
			)

		if ctype == 'mixed':
			return run_calculators(
				atoms=atoms, vasp_opts=vasp_opts, kcut=kcut, produce_steps=None,
				shel=shel, kwds=kwds, gulp_opts=gulp_opts, lib=lib, calcs=calcs,
				dist_cutoff=dist_cutoff, qe_opts=qe_opts, gulp_command=gulp_command,
				gulp_timeout=gulp_timeout, n_opts=n_opts, rel=rel,
				relaxer_opts=relaxer_opts, opt_class=opt_class, opt_device=opt_device,
				mode=mode,
			)

	except Exception:
		# Deliberately broad: any backend failure rejects this one structure and
		# the search moves on. Exception rather than a bare except, so that
		# KeyboardInterrupt and SystemExit still stop the run.
		return atoms, FAILED_ENERGY, False

	return atoms, FAILED_ENERGY, False
