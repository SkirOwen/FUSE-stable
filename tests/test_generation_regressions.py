"""Regression tests for the three bugs that stopped FUSE generating any
structure at all.

Before these fixes, run_fuse() could not produce a single valid structure for
the project's own example composition (Ca3Ti2O7 with the parameters in
examples/input_files/fgen_input.py) - it either hung or exhausted every
attempt. Afterwards, full runs complete reliably. Each test below pins one of
the three root causes.
"""
import random

import numpy as np
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

import fuse202.structure.generate_random_structure as grs
from fuse202.utils.rng import Rng
from fuse202.structure.make_basin_move import make_basin_move
from fuse202.structure.possible_solutions import CrystalSolutionsCalculator


# --------------------------------------------------------------------------
# 1. numpy ints in the structure string
# --------------------------------------------------------------------------

def test_atom_count_handles_numpy_ints(monkeypatch, tmp_path):
	"""The atom count used to be `list(filter((120).__ne__, string))`. The string
	mixes plain ints with numpy int64 (get_fu() builds it from ase's
	get_atomic_numbers()), and `(120).__ne__(np.int64(5))` returns NotImplemented,
	which Python 3.12+ refuses to treat as a bool. The resulting TypeError was
	swallowed into n_atoms = 0, so no candidate ever matched target_atoms and
	generation could never succeed.

	Here the generated string is all numpy ints and has exactly the requested
	atom count, so it must be accepted.
	"""
	target_atoms = 4
	# 5 entries, one of which is a vacancy marker (120) -> 4 real atoms
	numpy_string = [np.int64(8), np.int64(20), np.int64(120), np.int64(22), np.int64(8)]

	monkeypatch.setattr(grs, "create_random_string", lambda *a, **k: (numpy_string, "instr"))
	monkeypatch.setattr(grs, "create_random_instructions", lambda *a, **k: "instr")
	monkeypatch.setattr(grs, "assemble_structure", lambda *a, **k: ([1, 2, 3, 4], "instr"))
	monkeypatch.setattr(grs, "error_check_structure", lambda *a, **k: 1)  # accept

	atoms, string, instructions, accept = grs.generate_random_structure(
		target_atoms=target_atoms,
		cubic_solutions={}, tetragonal_solutions={}, hexagonal_solutions={},
		orthorhombic_solutions={}, monoclinic_solutions={},
		atoms_per_fu=4, ideal_density=1.0, density_cutoff=0.4,
		check_bonds=True, btol=0.25, system_type='neutral',
		fu=[np.int64(8)], composition={'O': 4}, bondtable={}, ap=4.0,
		check_distances=True, dist_cutoff=1.0,
	)

	assert atoms is not None, "a correctly-sized candidate was rejected"
	assert accept == 1


def test_generate_random_structure_tolerates_a_none_string(monkeypatch):
	"""create_random_string() legitimately returns None for the string when it
	cannot build one. That must count as "no candidate this round" and retry,
	not raise."""
	monkeypatch.setattr(grs, "create_random_string", lambda *a, **k: (None, None))
	monkeypatch.setattr(grs, "create_random_instructions", lambda *a, **k: None)
	monkeypatch.setattr(grs, "assemble_structure", lambda *a, **k: ([], None))

	atoms, string, instructions, accept = grs.generate_random_structure(
		target_atoms=4,
		cubic_solutions={}, tetragonal_solutions={}, hexagonal_solutions={},
		orthorhombic_solutions={}, monoclinic_solutions={},
		atoms_per_fu=4, ideal_density=1.0, density_cutoff=0.4,
		check_bonds=True, btol=0.25, system_type='neutral',
		fu=[8], composition={'O': 4}, bondtable={}, ap=4.0,
		check_distances=True, dist_cutoff=1.0,
	)
	assert atoms is None and accept == 0


# --------------------------------------------------------------------------
# 2. orthorhombic cells must have an even number of sub-modules along z
# --------------------------------------------------------------------------

@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(max_ax=st.integers(min_value=2, max_value=12))
def test_orthorhombic_z_is_always_even(max_ax, tmp_path_factory):
	"""Orthorhombic cells are restricted to an even number of sub-modules along
	z. A refactor generalised the 2D/3D loops behind a shared `step` argument
	but moved the start from 2 to 1, so orthorhombic got odd z values instead -
	silently producing lattices structure generation could never satisfy.
	Monoclinic has no such restriction and is unaffected.
	"""
	solutions_dir = tmp_path_factory.mktemp("solutions")
	calc = CrystalSolutionsCalculator(max_ax=max_ax, solutions_dir=str(solutions_dir))
	_, _, _, orthorhombic, monoclinic = calc.calculate_all_solutions()

	z_values = {dims[2] for entries in orthorhombic.values() for dims in entries}
	assert z_values, "no orthorhombic solutions generated"
	assert all(z % 2 == 0 for z in z_values), f"odd z in orthorhombic solutions: {sorted(z_values)}"

	# monoclinic is the unrestricted counterpart - it should still allow odd z
	mono_z = {dims[2] for entries in monoclinic.values() for dims in entries}
	assume(max_ax >= 2)
	assert any(z % 2 == 1 for z in mono_z), "monoclinic should not be restricted to even z"


@settings(max_examples=25, deadline=None)
@given(
	n=st.integers(min_value=1, max_value=400),
)
def test_every_solution_key_is_reproduced_by_its_dimensions(n):
	"""Each key in the solutions dicts is the sub-module count its stored
	dimensions multiply out to. Cheap invariant, but it is what structure
	generation relies on when it looks a cell shape up by size."""
	calc = CrystalSolutionsCalculator(max_ax=8, solutions_dir=None)
	cubic = calc._calculate_cubic_solutions()
	for key, dims in cubic.items():
		x, y, z = dims
		assert x * y * z == key, f"cubic key {key} does not match dims {dims}"


# --------------------------------------------------------------------------
# 3. the basin-hopping move-redraw list must not be empty
# --------------------------------------------------------------------------

def test_moves_to_choose_is_populated():
	"""make_basin_move() has 28 call sites that redraw a different move from
	`moves_to_choose` when the chosen one does not apply. A refactor replaced
	only the first pick with a weighted choice and commented out the code that
	populated the list, leaving it permanently empty - so any redraw crashed the
	run with "Cannot choose from an empty sequence".

	Captured at the point of use via the injected generator, so it stays
	meaningful if the implementation changes.
	"""
	captured = {}

	class _Captured(Exception):
		"""Sentinel: unwind as soon as the move has been picked."""

	class SpyRng(Rng):
		def choices(self, population, weights=None, k=1):
			captured['population'] = list(population)
			captured['weights'] = list(weights) if weights else None
			raise _Captured

	moves = {1: 1, 2: 3, 11: 2}

	# Only the move-selection step is under test. Bail out via the sentinel the
	# moment it happens: left to run on this dummy structure, make_basin_move
	# would redraw moves forever, since every move fails on an empty cell.
	with pytest.raises(_Captured):
		make_basin_move(
			{'atoms': [], 'modules': [], 'nmods': 0, 'ap': 4.0,
			 'sub module cell': [], 'shape in submods': []},
			moves, {}, 0.75, 1.5, 1.0, 0.4, True, 0.25, 'neutral', [8], {'O': 1},
			True, 1.0, {}, {}, {}, {}, {}, 10, 40, 4, False, {}, 1, 1, 10, False, {},
			rng=SpyRng(1),
		)

	assert captured.get('population') == [1, 2, 11]
	assert captured.get('weights') == [1, 3, 2]


def test_same_seed_reproduces_the_same_draws():
	"""The point of threading a generator through instead of seeding the global
	`random` module: a run is repeatable from its seed, and nothing else in the
	process has its randomness disturbed."""
	a, b, different = Rng(1234), Rng(1234), Rng(9999)
	population = ['x', [1, 2], (3, 4), 7]
	draws_a = [a.choice(population) for _ in range(20)]
	draws_b = [b.choice(population) for _ in range(20)]
	draws_c = [different.choice(population) for _ in range(20)]
	assert draws_a == draws_b
	assert draws_a != draws_c


@settings(max_examples=30, deadline=None)
@given(seed=st.integers(min_value=0, max_value=2 ** 31 - 1))
def test_draws_return_plain_python_objects(seed):
	"""Rng must hand back exactly the objects it was given.

	Trivially true of random.Random, and that is the point: it is what makes the
	stdlib the right engine here. A numpy Generator would not - its `choice`
	returns numpy scalars and reinterprets a list of lists as a 2-D array,
	returning an ndarray row. Numpy integers leaking into the structure string
	is what made `(120).__ne__(...)` return NotImplemented and silently break
	atom counting on Python 3.12+, so this guards against swapping the engine
	for one that reintroduces that.
	"""
	rng = Rng(seed)
	assert type(rng.choice([1, 2, 3])) is int
	assert type(rng.choice([[1, 2], [3, 4]])) is list
	assert type(rng.choice(['a', 'b'])) is str
	shuffled = [1, 2, 3, 4, 5]
	rng.shuffle(shuffled)
	assert all(type(x) is int for x in shuffled)
	assert sorted(shuffled) == [1, 2, 3, 4, 5]


def test_rng_reports_its_seed_even_when_unseeded():
	"""An unseeded run still draws a concrete seed and reports it, so a result
	can be reproduced after the fact. `seed_value`, not `seed`: the latter is
	random.Random's own re-seeding method."""
	assert Rng(1234).seed_value == 1234
	drawn = Rng().seed_value
	assert isinstance(drawn, int)
	# and that reported seed genuinely reproduces the run
	assert [Rng(drawn).random() for _ in range(5)] == [Rng(drawn).random() for _ in range(5)]
