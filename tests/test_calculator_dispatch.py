"""Tests for the single calculator dispatch point.

run_fuse() used to inline the same five-branch `if ctype == ...` chain twice,
once for the initial population and once in the basin-hopping search. The two
copies had drifted: the search copy carried leftover debug prints and did not
forward `use_spglib` to CHGNet. Both now route through
fuse202.calculators.dispatch.run_calculator().
"""
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import fuse202.calculators.dispatch as dispatch
from fuse202.calculators.dispatch import FAILED_ENERGY, SUPPORTED_CALCULATORS, run_calculator

BACKEND_FOR = {
	'gulp': 'run_gulp',
	'vasp': 'run_vasp',
	'qe': 'run_qe',
	'chgnet': 'run_chgnet',
	'mixed': 'run_calculators',
}


@pytest.mark.parametrize("ctype", SUPPORTED_CALCULATORS)
def test_each_ctype_routes_to_its_backend(ctype, monkeypatch):
	calls = []

	def fake(*args, **kwargs):
		calls.append(ctype)
		return "relaxed-atoms", -42.0, True

	monkeypatch.setattr(dispatch, BACKEND_FOR[ctype], fake)
	atoms, energy, converged = run_calculator(ctype, "atoms")

	assert calls == [ctype], f"{ctype} did not reach {BACKEND_FOR[ctype]}"
	assert (atoms, energy, converged) == ("relaxed-atoms", -42.0, True)


@pytest.mark.parametrize("ctype", SUPPORTED_CALCULATORS)
def test_no_other_backend_is_invoked(ctype, monkeypatch):
	"""Routing must be exclusive - picking one backend must not run another."""
	invoked = []
	for name in BACKEND_FOR.values():
		monkeypatch.setattr(
			dispatch, name,
			lambda *a, _n=name, **k: (invoked.append(_n), ("atoms", 0.0, True))[1],
		)
	run_calculator(ctype, "atoms")
	assert invoked == [BACKEND_FOR[ctype]]


@settings(max_examples=25, deadline=None)
@given(
	ctype=st.sampled_from(SUPPORTED_CALCULATORS),
	exc_type=st.sampled_from([ValueError, RuntimeError, OSError, TypeError, KeyError,
	                          IndexError, ZeroDivisionError, AttributeError]),
)
def test_any_backend_failure_becomes_a_rejected_structure(ctype, exc_type):
	"""A structure search is expected to propose things that make calculators
	fall over, so one bad candidate costs one candidate rather than the whole
	run. Whatever the backend raises, the result is the atoms handed back
	untouched with FAILED_ENERGY and converged=False.

	Patched by hand rather than with monkeypatch: that fixture is
	function-scoped, so under @given it would be set up once and leak across
	examples.
	"""
	def raising(*args, **kwargs):
		raise exc_type("backend blew up")

	name = BACKEND_FOR[ctype]
	original = getattr(dispatch, name)
	setattr(dispatch, name, raising)
	try:
		sentinel = object()
		atoms, energy, converged = run_calculator(ctype, sentinel)
	finally:
		setattr(dispatch, name, original)

	assert atoms is sentinel, "the original atoms should be returned untouched"
	assert energy == FAILED_ENERGY
	assert converged is False


@pytest.mark.parametrize("ctype", ["", "GULP", "nonsense", None])
def test_unrecognised_ctype_is_treated_as_a_failure(ctype):
	"""What the original if-chain did by falling through every branch."""
	sentinel = object()
	atoms, energy, converged = run_calculator(ctype, sentinel)
	assert atoms is sentinel
	assert energy == FAILED_ENERGY
	assert converged is False


@pytest.mark.parametrize("interrupt", [KeyboardInterrupt, SystemExit])
def test_interrupts_are_not_swallowed(interrupt, monkeypatch):
	"""The original used a bare `except:`, which caught KeyboardInterrupt and
	SystemExit too - so a long run could not be Ctrl-C'd out of during a
	calculation. Only Exception is caught now."""
	def interrupted(*args, **kwargs):
		raise interrupt

	monkeypatch.setattr(dispatch, "run_gulp", interrupted)
	with pytest.raises(interrupt):
		run_calculator('gulp', "atoms")


def test_use_spglib_reaches_chgnet(monkeypatch):
	"""The search-loop copy of the dispatch chain forgot to forward use_spglib
	to CHGNet, so a run configured with use_spglib=False silently got the
	default of True for every structure the search relaxed. Deduplicating the
	two copies fixed that; this pins it."""
	seen = {}

	def fake_chgnet(atoms, **kwargs):
		seen.update(kwargs)
		return atoms, -1.0, True

	monkeypatch.setattr(dispatch, "run_chgnet", fake_chgnet)
	run_calculator('chgnet', "atoms", use_spglib=False)
	assert seen.get('use_spglib') is False
