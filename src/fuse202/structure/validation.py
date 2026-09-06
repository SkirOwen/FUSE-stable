"""Sanity checks applied to a structure after an energy calculation.

A relaxation can hand back something physically meaningless - atoms lost or
gained, or driven on top of each other - and an energy for it that would
otherwise look attractive to the search. These checks reject those results by
marking them unconverged and assigning the same large energy a failed
calculation gets, so the search treats them as bad candidates rather than
discoveries.

run_fuse() applied this identical sequence in two places, once for the initial
population and once in the basin-hopping search.
"""
from __future__ import annotations

from decimal import Decimal

# Matches the energy a failed calculation is given in calculators.dispatch.
FAILED_ENERGY = 1.e20


def shortest_interatomic_distance(atoms) -> float:
	"""Shortest distance between any two atoms, across periodic images.

	Uses a 2x2x2 supercell so contacts through the cell boundary are seen, and
	ignores the zero self-distances on the diagonal.
	"""
	repeated = atoms.repeat([2, 2, 2])
	all_distances = repeated.get_all_distances()
	non_zero = [
		distance
		for row in all_distances
		for distance in row
		if distance != 0
	]
	return min(non_zero)


def check_relaxed_structure(
		atoms,
		energy,
		converged,
		*,
		expected_atoms: int,
		dist_cutoff: float,
		e_prec,
) -> tuple[float, bool]:
	"""Validate a just-relaxed structure and normalise its energy.

	Returns `(energy_per_atom, converged)`. A structure that lost or gained
	atoms during relaxation, or that contains a contact at or below
	`dist_cutoff`, is marked unconverged and given FAILED_ENERGY.

	Note the ordering, which is preserved from the original: the energy is
	divided by the atom count *after* any rejection, so a rejected structure
	carries FAILED_ENERGY/len(atoms) rather than FAILED_ENERGY itself.
	"""
	if len(atoms) != expected_atoms:
		converged = False
		energy = FAILED_ENERGY

	if shortest_interatomic_distance(atoms) <= dist_cutoff:
		converged = False
		energy = FAILED_ENERGY

	energy = energy / len(atoms)
	energy = float(Decimal(energy).quantize(Decimal(str(e_prec))))
	return energy, converged
