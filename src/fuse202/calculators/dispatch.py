"""Single place where FUSE hands a structure to an energy calculator.

Every backend shares one contract - given an Atoms object, return
`(atoms, energy, converged)` - so the choice between them is a lookup rather
than anything structural. run_fuse() used to inline the same five-branch
`if ctype == ...` chain twice, once for the initial population and once inside
the basin-hopping search, and the two copies had already drifted apart: the
search copy carried leftover debug prints and, more importantly, forgot to pass
`use_spglib` to CHGNet, so a run with `use_spglib=False` quietly got the
default of True for every structure the search relaxed.

Failure convention, preserved from the original: any exception from a backend
is turned into a rejected structure - the atoms are handed back untouched with
`energy = 1e20` and `converged = False` - rather than aborting the run. A
structure search is expected to propose things that make calculators fall over,
so one bad candidate should cost one candidate, not the whole job. An
unrecognised `ctype` produces the same result, which is also what the original
if-chain did by falling through every branch.
"""
from __future__ import annotations

from fuse202.calculators.run_chgnet import run_chgnet
from fuse202.calculators.run_gulp import run_gulp
from fuse202.calculators.run_multiple_calculators import run_calculators
from fuse202.calculators.run_qe import run_qe
from fuse202.calculators.run_vasp import run_vasp

# Energy assigned to a structure whose calculation failed - large enough that
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

	Returns (atoms, energy, converged). On any failure - including an
	unrecognised ctype - returns the atoms unchanged with FAILED_ENERGY and
	converged=False.
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
		# the search moves on. Exception rather than a bare except, so
		# KeyboardInterrupt and SystemExit still stop the run.
		return atoms, FAILED_ENERGY, False

	# Unrecognised ctype - same outcome as a failure, matching what the original
	# if-chain did when no branch matched.
	return atoms, FAILED_ENERGY, False
