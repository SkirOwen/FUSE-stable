"""Sanity checks applied to a structure after an energy calculation.

A relaxation can return a physically meaningless structure, with atoms lost,
gained, or driven on top of each other, together with an energy that would
otherwise look attractive to the search. These checks reject such results so
that the search treats them as bad candidates.
"""
from __future__ import annotations

from decimal import Decimal

# Matches the energy a failed calculation is given in calculators.dispatch.
FAILED_ENERGY = 1.e20


def shortest_interatomic_distance(atoms) -> float:
	"""Return the shortest distance between any two atoms.

	A 2x2x2 supercell is used so that contacts through the periodic boundary
	are included. Zero self distances on the diagonal are ignored.

	Parameters
	----------
	atoms : ase.Atoms
		The structure to measure.

	Returns
	-------
	float
		Shortest interatomic distance in Angstroms.
	"""
	repeated = atoms.repeat([2, 2, 2])
	all_distances = repeated.get_all_distances()
	non_zero = [
		distance
		for row in all_distances
		for distance in row
		if distance != 0
	]
	shortest = min(non_zero)
	return shortest


def check_relaxed_structure(
		atoms,
		energy,
		converged,
		*,
		expected_atoms: int,
		dist_cutoff: float,
		e_prec,
) -> tuple[float, bool]:
	"""Validate a relaxed structure and convert its energy to eV per atom.

	A structure that changed atom count during relaxation, or that contains a
	contact at or below `dist_cutoff`, is marked unconverged and given
	FAILED_ENERGY. The energy is divided by the atom count after any rejection,
	so a rejected structure carries FAILED_ENERGY divided by its atom count.

	Parameters
	----------
	atoms : ase.Atoms
		The structure as returned by the calculator.
	energy : float
		Total energy reported by the calculator.
	converged : bool
		Whether the calculator reported convergence.
	expected_atoms : int
		Atom count the structure had before relaxation.
	dist_cutoff : float
		Shortest permitted interatomic contact, in Angstroms.
	e_prec : float
		Precision the energy is rounded to.

	Returns
	-------
	tuple of (float, bool)
		Energy in eV per atom, and the updated convergence flag.

	Notes
	-----
	The energy is converted with `float()` before `Decimal`. Backends do not all
	return a Python float: CHGNet returns numpy.float32, which Decimal rejects
	outright, so omitting the conversion makes every CHGNet run fail on the
	first structure it relaxes.
	"""
	if len(atoms) != expected_atoms:
		converged = False
		energy = FAILED_ENERGY

	if shortest_interatomic_distance(atoms) <= dist_cutoff:
		converged = False
		energy = FAILED_ENERGY

	energy_per_atom = float(energy) / len(atoms)
	rounded = Decimal(energy_per_atom).quantize(Decimal(str(e_prec)))
	energy_per_atom = float(rounded)
	return energy_per_atom, converged
