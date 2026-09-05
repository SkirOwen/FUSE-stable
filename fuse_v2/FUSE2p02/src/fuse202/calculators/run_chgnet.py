from __future__ import annotations

import spglib

from ase import Atoms
from ase.io import read, write
from chgnet.model import StructOptimizer
from chgnet.model.model import CHGNet
from pymatgen.core import Structure


def check_isolated_atoms(atoms: Atoms, threshold: float = 6.0) -> bool:
	"""
	Check for isolated atoms in the structure.

	Parameters
	----------
	atoms : Atoms
		ASE Atoms object.
	threshold : float
		Minimum allowable distance for non-isolated atoms.

	Returns
	-------
	bool
		True if no isolated atoms are found, False otherwise.
	"""
	for i in range(len(atoms)):
		distances = atoms.get_distances(i, list(range(len(atoms))), mic=True)
		if min(filter(lambda d: d > 0, distances), default=float('inf')) > threshold:
			return False
	return True


def run_chgnet(
		atoms,
		n_opts: int = 2,
		rel=StructOptimizer(),
		relaxer_opts: dict | None = None,
		opt_class: list[str] | None = None,
		mode: str = 'relax',
		opt_device: str = 'cpu',
		use_spglib: bool = True
) -> tuple[Atoms, float, bool]:
	"""
	Run CHGNet relaxation or single-point energy calculation.

	Parameters
	----------
	atoms :
		ASE Atoms object representing the structure.
	n_opts : int
		Number of optimization steps.
	rel : StructOptimizer
		Optimizer instance.
	relaxer_opts : dict | None
		Relaxation options, including fmax, steps, and verbosity.
	opt_class : list[str] | None
		Optimizer class names.
	mode : str
		"relax" for structure relaxation, "single" for single-point calculation.
	opt_device : str
		Device used for optimization (e.g., 'cpu' or 'cuda').
	use_spglib : bool
		Whether to use SPGLIB for standardization.

	Returns
	-------
	tuple[Atoms, float, bool]
		Updated Atoms object, energy, and convergence status.
	"""
	# First do a quick test to make sure that there are no isolated atoms (> 6 angstroms from their nearest neighbour)
	if opt_class is None:
		opt_class = ['FIRE', 'BFGSLineSearch']

	if relaxer_opts is None:
		# TODO: use kwargs here to pass option to relaxer
		relaxer_opts = {
			'fmax': [0.1, 0.05],
			'steps': [250, 750],
			'verbose': [True, True]
		}

	# Check if isolated and if so stops the run
	for i in range(len(atoms)):
		j = list(range(len(atoms)))
		j.remove(i)
		distances = atoms.get_distances(i, j, mic=True)
		if min(distances) > 6.0:
			converged = False
			energy = 1.e20
			print(min(distances))
			return atoms, energy, converged

	chgnet = CHGNet.load()

	write("temp.cif", atoms)
	temp_atoms = Structure.from_file("temp.cif")

	if mode == 'relax':
		try:
			for i in range(n_opts):
				relaxer = StructOptimizer(optimizer_class=opt_class[i], use_device=opt_device)
				prediction = relaxer.relax(
					temp_atoms,
					fmax=relaxer_opts['fmax'][i],
					steps=relaxer_opts['steps'][i],
					verbose=relaxer_opts['verbose'][i],
				)
				temp_atoms = prediction['final_structure']

			result = chgnet.predict_structure(temp_atoms)
			energy = result['e'] * len(temp_atoms)
			forces = result['f']

			temp_atoms.to_file("temp.cif")
			print(abs(forces.max()) <= relaxer_opts['fmax'][-1])

			if abs(forces.max()) <= relaxer_opts['fmax'][-1]:
				converged = True
			else:
				converged = False

			atoms = read("temp.cif")

			if use_spglib:
				try:
					lattice, positions, numbers = spglib.standardize_cell(atoms, symprec=1.e-5)
					temp2 = Atoms(numbers=numbers, pbc=True)
					temp2.cell = lattice
					temp2.set_scaled_positions(positions)
					atoms = temp2.copy()
				# print("I'm using SPGLIB!")

				except:
					# print("I failed at using SPGLIB!!")
					pass

		except:
			converged = False
			energy = 1.e20

	if mode == 'single':
		result = chgnet.predict_structure(temp_atoms)
		energy = result['e'] * len(temp_atoms)
		converged = True

	return atoms, energy, converged

# atoms,energy,converged=run_chgnet(atoms=read("after.cif"))

# print(energy,converged)

class CHGNetRunner:
	pass
