import os
import sys
from ase.visualize import *
from fuse202.structure.assemble_structure_generator import *
from fuse202.structure.create_random_instructions import *
from fuse202.structure.create_random_string import *
from fuse202.structure.error_check_structure import *
from fuse202.utils.rng import default_rng


def generate_random_structure(
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
		vac_ratio='',
		max_fus='',
		target_number_atoms='',
		rng=None
):
	if rng is None:
		rng = default_rng()

	accept = 0
	n_atoms = 0
	complete = False
	max_attempts = 1000
	attempts = 0

	while not complete:
		# Counted AND checked here at the top. Both have to be up here: the
		# `except ... continue` paths below jump straight back to the top of the
		# loop, skipping anything at the bottom - which is why a
		# create_random_string()/assemble_structure() call that failed every time
		# used to loop forever, never counting an attempt and never reaching the
		# cap. Give up after max_attempts and report failure the same way a run
		# of clean rejections does, by handing back atoms=None.
		attempts += 1
		if attempts > max_attempts:
			accept = 0
			atoms = None
			string = None
			instructions = None
			break
		# print(attempts)
		# try:
		# print("\nstart")
		# generate initial string component and start of instructions
		# before we try this, print out a list of the variables that we're passing to the random string
		# print("cubic_solutions: ",cubic_solutions)
		# print("tetragonal solitions: ",tetragonal_solutions)
		# print("hexagonal_solutions: ",hexagonal_solutions)
		# print("orthorhombic_solutions: ",orthorhombic_solutions)
		# print("monoclinic_solutions: ",monoclinic_solutions)
		# print("atoms_per_fu: ",atoms_per_fu)
		# print("fu: ",fu)
		# print("vac_ratio: ",vac_ratio)
		# print("max_fus: ", max_fus)
		# print("system_type: ", system_type)
		# print("composition: ",composition)
		# print("ap: ",ap)

		try:
			string, instructions = create_random_string(
				cubic_solutions,
				tetragonal_solutions,
				hexagonal_solutions,
				orthorhombic_solutions,
				monoclinic_solutions,
				atoms_per_fu,
				fu,
				vac_ratio=vac_ratio,
				max_fus=max_fus,
				system_type=system_type,
				composition=composition,
				ap=ap,
				rng=rng,
			)
		except Exception:
			continue
		# check to see if it has the correct number of atoms
		# print("string: \n",string)
		# print("hello2")
		instructions = create_random_instructions(string, ap, cubic_solutions, tetragonal_solutions,
		                                          hexagonal_solutions, orthorhombic_solutions, monoclinic_solutions,
		                                          instructions, rng=rng)
		# print("instructions: \n",instructions)

		# print("atoms: \n",atoms)

		# Count the real atoms in the string; 120 (ord('x')) marks a vacancy.
		# This used to be `list(filter((120).__ne__, string))`, which is unsafe:
		# the string mixes plain ints with numpy int64 (get_fu() builds it from
		# ase's get_atomic_numbers()), and `(120).__ne__(np.int64(...))` returns
		# NotImplemented rather than a bool. Python 3.12+ raises
		# "NotImplemented should not be used in a boolean context" for that,
		# which the old bare `except` swallowed into n_atoms = 0 - so no
		# candidate ever matched target_atoms and generation could never
		# succeed on a modern interpreter. A comparison uses normal reflected
		# dispatch and handles both int types correctly.
		# create_random_string() returns None for the string when it fails to
		# build one; that is a legitimate "no candidate this round" result, so
		# treat it as zero atoms and let the loop retry.
		n_atoms = 0 if string is None else len([atom for atom in string if atom != 120])

		try:
			atoms, instructions = assemble_structure(string, instructions, rng=rng)
		except Exception:
			continue

		if len(atoms) > 0:
			accept = error_check_structure(atoms, ideal_density, density_cutoff, check_bonds, btol, system_type, fu,
			                               composition, bondtable, ap, check_distances, dist_cutoff,
			                               target_number_atoms=target_atoms)
		else:
			accept = 0

		# print(accept)
		# accept=1
		# except:
		#	pass

		if n_atoms == target_atoms:
			if accept == 1:
				complete = True

	# view(atoms)
	# print(n_atoms)
	# print("\n")
	# print("string:\n",string,"\ninstructions: \n",instructions)
	# sys.exit()
	return atoms, string, instructions, accept
