"""Tests for the post-calculation structure checks.

run_fuse() applied this identical sequence in two places (initial population
and basin-hopping search); it now lives in structure/validation.py.
"""
import pytest
from ase import Atoms
from hypothesis import given, settings
from hypothesis import strategies as st

from fuse202.structure.validation import (
	FAILED_ENERGY,
	check_relaxed_structure,
	shortest_interatomic_distance,
)


@pytest.fixture
def spread_out_cell():
	"""Two atoms 5 A apart in a 10 A cell - no short contacts, even across the
	periodic boundary."""
	return Atoms(symbols=['Sr', 'Ti'], positions=[(0, 0, 0), (5, 0, 0)],
	             cell=[10, 10, 10], pbc=True)


def test_good_structure_passes_and_energy_becomes_per_atom(spread_out_cell):
	energy, converged = check_relaxed_structure(
		spread_out_cell, -100.0, True, expected_atoms=2, dist_cutoff=1.0, e_prec=1.e-5)
	assert converged is True
	assert energy == -50.0  # -100 eV over 2 atoms


def test_losing_atoms_during_relaxation_is_rejected(spread_out_cell):
	"""A relaxation that changes the atom count produced something that is no
	longer the composition being searched."""
	energy, converged = check_relaxed_structure(
		spread_out_cell, -100.0, True, expected_atoms=3, dist_cutoff=1.0, e_prec=1.e-5)
	assert converged is False
	assert energy == FAILED_ENERGY / 2  # rejection happens before the per-atom divide


def test_atoms_on_top_of_each_other_are_rejected():
	"""An attractive-looking energy for a physically impossible structure is
	exactly what these checks exist to throw away."""
	collided = Atoms(symbols=['Sr', 'Ti'], positions=[(0, 0, 0), (0.2, 0, 0)],
	                 cell=[10, 10, 10], pbc=True)
	energy, converged = check_relaxed_structure(
		collided, -500.0, True, expected_atoms=2, dist_cutoff=1.0, e_prec=1.e-5)
	assert converged is False
	assert energy == FAILED_ENERGY / 2


def test_an_already_failed_calculation_stays_failed(spread_out_cell):
	energy, converged = check_relaxed_structure(
		spread_out_cell, FAILED_ENERGY, False, expected_atoms=2, dist_cutoff=1.0, e_prec=1.e-5)
	assert converged is False


def test_shortest_distance_sees_through_the_periodic_boundary():
	"""Two atoms far apart inside the cell can still be neighbours across it -
	which is why the check builds a 2x2x2 supercell rather than measuring
	within the cell alone."""
	atoms = Atoms(symbols=['Sr', 'Ti'], positions=[(0, 0, 0), (9.5, 0, 0)],
	              cell=[10, 10, 10], pbc=True)
	assert shortest_interatomic_distance(atoms) == pytest.approx(0.5)


@settings(max_examples=30, deadline=None)
@given(
	total_energy=st.floats(min_value=-5000, max_value=-1, allow_nan=False),
	n_atoms=st.integers(min_value=1, max_value=8),
)
def test_accepted_energy_is_always_the_per_atom_average(total_energy, n_atoms):
	"""Whatever the backend reports as a total, downstream code compares
	energies per atom - structures of different sizes are only comparable that
	way."""
	atoms = Atoms(symbols=['Sr'] * n_atoms,
	              positions=[(4.0 * i, 0, 0) for i in range(n_atoms)],
	              cell=[4.0 * n_atoms + 10, 20, 20], pbc=True)
	energy, converged = check_relaxed_structure(
		atoms, total_energy, True, expected_atoms=n_atoms, dist_cutoff=1.0, e_prec=1.e-5)
	assert converged is True
	assert energy == pytest.approx(total_energy / n_atoms, abs=1e-4)
