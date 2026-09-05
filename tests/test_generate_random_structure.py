"""Tests for the retry loops in structure generation - the ones that used to
make run_fuse() hang.

Manual testing of run_fuse() (with a faked energy calculator, so no external
binary was involved) hung indefinitely at "Generating Initial Population" for
two different realistic-looking compositions (SrTiO3 and Ca3Ti2O7, using
parameters straight out of examples/input_files/fgen_input.py). Two compounding
bugs caused it, both now fixed:

1. get_new_structure() (structure/make_new_structure.py) called
   generate_random_structure() again every time it came back empty, with no cap
   of its own - so whenever a composition could never satisfy
   error_check_structure(), it retried forever.
2. generate_random_structure()'s own budget (max_attempts) only worked when
   candidates were cleanly *rejected*. Two bare `except: continue` blocks
   skipped straight back to the top of the loop *before* `attempts += 1` ran,
   so if either call raised on every try, attempts never advanced and the cap
   never fired.

These tests now assert the fixed behaviour: both loops terminate. They are
written so that a regression re-introduces a *failing test*, not a hanging
suite - each one bounds the work with a call counter and asserts the loop gave
up on its own rather than relying on a timeout.
"""
import pytest

import fuse202.structure.generate_random_structure as grs
import fuse202.structure.make_new_structure as mns


def test_generate_random_structure_gives_up_after_max_attempts(monkeypatch):
	"""Clean-rejection path: every candidate is rejected by error_check_structure,
	so the attempt cap should fire and hand back atoms=None."""
	monkeypatch.setattr(grs, "create_random_string", lambda *a, **k: ("string", "instructions"))
	monkeypatch.setattr(grs, "create_random_instructions", lambda *a, **k: "instructions")
	monkeypatch.setattr(grs, "assemble_structure", lambda *a, **k: ([1, 2, 3], "instructions"))
	monkeypatch.setattr(grs, "error_check_structure", lambda *a, **k: 0)  # always reject

	atoms, string, instructions, accept = grs.generate_random_structure(
		target_atoms=5,
		cubic_solutions={}, tetragonal_solutions={}, hexagonal_solutions={},
		orthorhombic_solutions={}, monoclinic_solutions={},
		atoms_per_fu=5, ideal_density=1.0, density_cutoff=0.4,
		check_bonds=True, btol=0.25, system_type='neutral',
		fu=[38, 22, 8, 8, 8], composition={'Sr': 1, 'Ti': 1, 'O': 3},
		bondtable={}, ap=4.672, check_distances=True, dist_cutoff=1.0,
	)

	assert atoms is None
	assert accept == 0


def test_generate_random_structure_terminates_when_every_attempt_raises(monkeypatch):
	"""Regression test for the swallowed-exception bug: create_random_string()
	raising on every call must still count towards the attempt cap. Before the
	fix `attempts += 1` sat after the `except: continue`, so the counter never
	advanced and this looped forever.

	The call counter both proves the cap fired and bounds the damage if it ever
	regresses - the assertion fails on call 5000 instead of the suite hanging.
	"""
	calls = {"n": 0}

	def always_raises(*args, **kwargs):
		calls["n"] += 1
		if calls["n"] > 5000:
			raise AssertionError(
				"generate_random_structure() did not give up - the attempt cap "
				"is not counting failed attempts again"
			)
		raise ValueError("simulated failure inside create_random_string")

	monkeypatch.setattr(grs, "create_random_string", always_raises)

	atoms, string, instructions, accept = grs.generate_random_structure(
		target_atoms=5,
		cubic_solutions={}, tetragonal_solutions={}, hexagonal_solutions={},
		orthorhombic_solutions={}, monoclinic_solutions={},
		atoms_per_fu=5, ideal_density=1.0, density_cutoff=0.4,
		check_bonds=True, btol=0.25, system_type='neutral',
		fu=[38, 22, 8, 8, 8], composition={'Sr': 1, 'Ti': 1, 'O': 3},
		bondtable={}, ap=4.672, check_distances=True, dist_cutoff=1.0,
	)

	assert atoms is None
	assert calls["n"] == 1000, "should stop at exactly max_attempts (1000) tries"


def test_generate_random_structure_lets_keyboard_interrupt_through(monkeypatch):
	"""The retry blocks used to use a bare `except:`, which swallows
	KeyboardInterrupt and SystemExit too - meaning a stuck run could not even be
	Ctrl-C'd out of. They now catch Exception, so interrupts propagate."""

	def interrupted(*args, **kwargs):
		raise KeyboardInterrupt

	monkeypatch.setattr(grs, "create_random_string", interrupted)

	with pytest.raises(KeyboardInterrupt):
		grs.generate_random_structure(
			target_atoms=5,
			cubic_solutions={}, tetragonal_solutions={}, hexagonal_solutions={},
			orthorhombic_solutions={}, monoclinic_solutions={},
			atoms_per_fu=5, ideal_density=1.0, density_cutoff=0.4,
			check_bonds=True, btol=0.25, system_type='neutral',
			fu=[38, 22, 8, 8, 8], composition={'Sr': 1, 'Ti': 1, 'O': 3},
			bondtable={}, ap=4.672, check_distances=True, dist_cutoff=1.0,
		)


def test_get_new_structure_raises_instead_of_retrying_forever(monkeypatch):
	"""Regression test for the hang itself: when generate_random_structure()
	can never produce a structure, get_new_structure() must give up with an
	actionable error rather than asking again indefinitely."""
	calls = {"n": 0}

	def always_fails(*args, **kwargs):
		calls["n"] += 1
		if calls["n"] > 500:
			raise AssertionError(
				"get_new_structure() did not give up - its retry cap is gone again"
			)
		return None, None, None, 0

	monkeypatch.setattr(mns, "generate_random_structure", always_fails)

	with pytest.raises(RuntimeError, match="Could not generate a valid random structure"):
		mns.get_new_structure(
			composition={'Sr': 1, 'Ti': 1, 'O': 3}, max_atoms=10, imax_atoms=10,
			restart=False, max_ax=10, density_cutoff=0.4, check_bonds=True, btol=0.25,
			check_distances=True, system_type='neutral', dist_cutoff=1.0, vac_ratio=4,
			atoms_per_fu=5, imax_fus=2, max_fus=2, cubic_solutions={}, tetragonal_solutions={},
			hexagonal_solutions={}, orthorhombic_solutions={}, monoclinic_solutions={},
			bondtable={}, ideal_density=1.0, fu=[38, 22, 8, 8, 8], ap=4.672, use_spglib=False,
		)

	assert calls["n"] == 100, "should stop after exactly max_rounds (100) rounds"
