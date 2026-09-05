# fuse imports
import fuse202.bonds.bond_table
import math
import numpy
import pandas
from ase.io import *
from ase.visualize import *
from fuse202.structure.extract_modules import *
from fuse202.structure.generate_random_structure import *
from fuse202.structure.possible_solutions import *
# other imports
from random import choice


def get_new_structure(
		composition='',
		max_atoms='',
		imax_atoms='',
		restart='',
		max_ax='',
		density_cutoff='',
		check_bonds='',
		btol='',
		check_distances='',
		system_type='',
		dist_cutoff='',
		vac_ratio='',
		atoms_per_fu='',
		imax_fus='',
		max_fus='',
		cubic_solutions='',
		tetragonal_solutions='',
		hexagonal_solutions='',
		orthorhombic_solutions='',
		monoclinic_solutions='',
		bondtable='',
		ideal_density='',
		fu='',
		ap='',
		use_spglib='',
):
	gen = False
	# generate_random_structure() gives up after its own internal attempt budget
	# and returns atoms=None. Without a cap here we would just ask it again,
	# forever, whenever a composition/parameter combination can never satisfy
	# error_check_structure() - which is exactly how run_fuse() used to hang
	# during "Generating Initial Population" with no output and no way out.
	max_rounds = 100
	rounds = 0
	while gen == False:
		rounds += 1
		if rounds > max_rounds:
			raise RuntimeError(
				f"Could not generate a valid random structure after {max_rounds} rounds "
				f"of {max_rounds * 1000} total attempts for composition {composition}. "
				"No candidate passed error checking, so the search space is probably "
				"over-constrained: try relaxing density_cutoff, btol or dist_cutoff, "
				"raising max_atoms/imax_atoms, or check that the bond table covers "
				"every element in the composition."
			)
		##############################################################################
		# now call the generate random structure function
		# for each iteration of generating n structures:
		target_fu = choice(list(range(1, imax_fus + 1)))
		target_atoms = target_fu * atoms_per_fu

		atoms, string, instructions, accept = generate_random_structure(
			target_atoms,
			cubic_solutions,
			tetragonal_solutions,
			hexagonal_solutions,
			orthorhombic_solutions,
			monoclinic_solutions,
			atoms_per_fu,
			ideal_density,
			density_cutoff,
			check_bonds,
			btol,
			system_type,
			fu,
			composition,
			bondtable,
			ap,
			check_distances,
			dist_cutoff,
			vac_ratio=vac_ratio,
			max_fus=imax_fus,
			target_number_atoms=target_atoms,
		)

		##############################################################################
		# now call write the atoms object to a temporary file & then call the extract modules function
		if use_spglib:
			import spglib
			from ase import Atoms
			try:
				lattice, positions, numbers = spglib.standardize_cell(atoms, symprec=1.e-5)
				temp2 = Atoms(numbers=numbers, pbc=True)
				temp2.cell = lattice
				temp2.set_scaled_positions(positions)
				temp2.numbers = numbers
				atoms = temp2.copy()

			# print(atoms.cell.cellpar())

			except:
				pass

		if atoms is not None:
			write("temp.cif", atoms)
			input_files = ["temp.cif"]
			structure = extract_module(input_files, bondtable)
			gen = True

	# print(structure['sub module cell'])
	# print(structure['shape in submods'])

	return structure
