from typing import Literal

import numpy
import sys
from ase import Atoms
from ase.io import write
from ase.visualize import *





# start of information from shannon database,
# format = 1st key: atomic symbol, 2nd key(s): charge states, then list containing min bonds then max bonds


def _test_bonds(
		atoms,
		cations: list,
		anions: list,
		charges,
		ap,
		lib: dict[str, dict[int, list[int, int]]],
		system_type='',
		count_bonds=False
):
	if system_type == 'ionic':
		# print(anions)
		anion_rad = [v[-1] for anion in anions for v in lib[anion].values()]
		# anion_rad = []
		# for anion in anions:
		# 	temp = lib[anion]
		# 	keys = list(temp.keys())
		# 	if len(temp) == 1:
		# 		anion_rad.append(temp[keys[0]][-1])
		# 	if len(temp) > 1:
		# 		for j in range(len(temp)):
		# 			anion_rad.append(temp[keys[j]][-1])
		# count_bonds=True

		anion_rad = max(anion_rad)
		if count_bonds:
			print("########## Starting structure ##########")
			print(f"number of atoms: {len(atoms)}")

		ci = [i for i, atom in enumerate(atoms) if atom.symbol in cations]

		atoms_repeated = atoms.repeat([2, 2, 2])
		ai = [i for i, atom in enumerate(atoms_repeated) if atom.symbol in anions]

		# bond_dist=(data[-1]+anion_rad)*1.25
		bond_dist = ap * 0.785
		if count_bonds:
			print(f"bond_distance {bond_dist}")

		wrong = 0
		for idx in ci:
			cat = atoms_repeated[idx]
			try:
				data = lib[atoms_repeated[idx].symbol][charges[cations.index(atoms_repeated[idx].symbol)]]
				br = list(range(data[0] - 1, data[1] + 2))
			except KeyError:
				min_bonds = []
				max_bonds = []
				radius = []
				temp = []
				for j in range(len(list(lib[atoms_repeated[idx].symbol].keys()))):
					temp = lib[atoms_repeated[idx].symbol][list(lib[atoms_repeated[idx].symbol].keys())[j]]
					min_bonds.append(temp[0])
					max_bonds.append(temp[1])
					radius.append(temp[2])
				data = [min(min_bonds), max(max_bonds), sum(radius) / len(radius)]
				br = list(range(data[0] - 1, data[1] + 2))


			nbonds = 0
			temp1 = list(atoms_repeated.get_distances(idx, ai, mic=True))
			# print(temp1)
			try:
				temp1.remove(0.)
			except:
				pass

			# to restore the bond test used in v. 1.02 change the distance comparison to: if temp1[j] <= ap*0.785
			for j in range(len(temp1)):
				if temp1[j] <= bond_dist:
					nbonds += 1

			try:
				data = lib[atoms_repeated[idx].symbol][charges[cations.index(atoms_repeated[idx].symbol)]]
				br = list(range(data[0] - 1, data[1] + 2))
			except KeyError:
				min_bonds = []
				max_bonds = []
				radius = []
				temp = []
				for j in range(len(list(lib[atoms_repeated[idx].symbol].keys()))):
					temp = lib[atoms_repeated[idx].symbol][list(lib[atoms_repeated[idx].symbol].keys())[j]]
					min_bonds.append(temp[0])
					max_bonds.append(temp[1])
					radius.append(temp[2])
				data = [min(min_bonds), max(max_bonds), sum(radius) / len(radius)]
				br = list(range(data[0] - 1, data[1] + 2))
			# print(data)
			if count_bonds:
				print(cat.symbol, " ", nbonds)
			if nbonds not in br:
				wrong += 1
		# print ("total errors", wrong)
		if count_bonds:
			write("test_structure.cif", atoms_repeated)
		# sys.exit()
		try:
			return float(wrong) / float(len(ci))
		except ZeroDivisionError:
			return 1.

	if system_type == 'neutral':
		ci = []
		ai = []
		# atoms.rattle(0.025)
		for i in range(len(atoms)):
			if atoms[i].symbol in cations:
				ci.append(i)
			if atoms[i].symbol in anions:
				ai.append(i)
		wrong = 0
		for i in range(len(ci)):
			cat = atoms[ci[i]]
			nbonds = 0
			for j in range(len(ai)):
				ani = atoms[ai[j]]
				bond = atoms.get_distance(ci[i], ai[j], mic=True)
				if bond <= ap * 0.785:
					nbonds += 1

			min_bonds = []
			max_bonds = []
			radius = []
			temp = []
			for j in range(len(list(lib[atoms[ci[i]].symbol].keys()))):
				temp = lib[atoms[ci[i]].symbol][list(lib[atoms[ci[i]].symbol].keys())[j]]
				min_bonds.append(temp[0])
				max_bonds.append(temp[1])
				radius.append(temp[2])

			data = [min(min_bonds), max(max_bonds), sum(radius) / len(radius)]
			br = list(range(data[0] - 1, data[1] + 2))
			# print cat.symbol, " ",nbonds
			if nbonds not in br:
				wrong += 1
		# print ("total errors", wrong)
		# write("test_structure.cif",atoms)
		# sys.exit()
		try:
			return float(wrong) / float(len(ci))
		except ZeroDivisionError:
			return 1.

# wrong_fract=test_bonds(atoms=atoms,cations=cations,anions=anions,charges=charges,ap=ap,lib=lib)
# print wrong_fract


def test_bonds(
		atoms: Atoms,
		cations: list[str],
		anions: list[str],
		charges: list[int],
		ap: float,
		lib: dict[str, dict[int, list[int, int]]],
		system_type: Literal["ionic", "neutral"],
		count_bonds: bool = False
) -> float:
	""""""

	bond_dist = ap * 0.785

	if system_type == 'ionic':
		# print(anions)
		anion_rad = [v[-1] for anion in anions for v in lib[anion].values()]

		anion_rad = max(anion_rad)
		if count_bonds:
			print("########## Starting structure ##########")
			print(f"number of atoms: {len(atoms)}")

		ci = [i for i, atom in enumerate(atoms) if atom.symbol in cations]

		atoms_repeated = atoms.repeat([2, 2, 2])
		ai = [i for i, atom in enumerate(atoms_repeated) if atom.symbol in anions]

		# bond_dist=(data[-1]+anion_rad)*1.25
		if count_bonds:
			print(f"bond_distance {bond_dist}")

		wrong = 0
		for idx in ci:
			cat = atoms_repeated[idx]

			distances = atoms_repeated.get_distances(idx, ai, mic=True)
			# remove self-distance (distance 0)
			distances = [d for d in distances if d > 0]

			# to restore the bond test used in v. 1.02 change the distance comparison to: if temp1[j] <= ap*0.785
			nbonds = sum(1 for d in distances if d <= bond_dist)

			symbol = cat.symbol
			try:
				cation_charge = charges[cations.index(symbol)]  # TODO: Index Out Of Range possible here
				data = lib[symbol][cation_charge]
			except IndexError:
				min_bonds = []
				max_bonds = []
				radius = []
				for entry in lib[symbol].values():
					min_bonds.append(entry[0])
					max_bonds.append(entry[1])
					radius.append(entry[2])

				data = [min(min_bonds), max(max_bonds), sum(radius) / len(radius)]

			allowed_range = list(range(data[0] - 1, data[1] + 2))

			if count_bonds:
				print(cat.symbol, " ", nbonds)

			if nbonds not in allowed_range:
				wrong += 1

		if count_bonds:
			write("test_structure.cif", atoms_repeated)
		return float(wrong) / float(len(ci)) if ci else 1.0

	elif system_type == 'neutral':
		ci = [i for i, atom in enumerate(atoms) if atom.symbol in cations]
		ai = [i for i, atom in enumerate(atoms) if atom.symbol in anions]

		wrong = 0
		for idx in ci:
			cat = atoms[idx]
			symbol = cat.symbol

			distances = atoms.get_distances(idx, ai, mic=True)
			# remove self-distance (distance 0)
			distances = [d for d in distances if d > 0]

			# to restore the bond test used in v. 1.02 change the distance comparison to: if temp1[j] <= ap*0.785
			nbonds = sum(1 for d in distances if d <= bond_dist)

			min_bonds = []
			max_bonds = []
			radius = []
			for entry in lib[symbol].values():
				min_bonds.append(entry[0])
				max_bonds.append(entry[1])
				radius.append(entry[2])

			data = [min(min_bonds), max(max_bonds), sum(radius) / len(radius)]
			allowed_range = list(range(data[0] - 1, data[1] + 2))

			if nbonds not in allowed_range:
				wrong += 1

		# write("test_structure.cif",atoms)
		return float(wrong) / float(len(ci)) if ci else 1.0

# wrong_fract=test_bonds(atoms=atoms,cations=cations,anions=anions,charges=charges,ap=ap,lib=lib)
# print wrong_fract


def main():
	# # atoms = read("0.cif")
	# cations = ['Y', 'Ti', 'Ba']
	# anions = ['O']
	# charges = [3, 4, 2]
	# ap = 4.95
	# lib = {
	# 	'Y': {3: [6, 9]},
	# 	'Ti': {4: [4, 8]},
	# 	'Ba': {2: [6, 12]}
	# }
	# Define a test library dictionary.
	# For cations, the format is: charge: [min_bonds, max_bonds, average_radius]
	lib = {
		'Y': {3: [6, 9, 1.0]},
		'O': {-2: [1, 2, 1.4]},
	}

	# Create a test structure for an ionic system.
	# Construct a simple structure with one Y and six O atoms arranged octahedrally around Y.
	positions = [
		(2.5, 2.5, 2.5),  # Y at the centre
		(3.5, 2.5, 2.5),  # O on +x
		(1.5, 2.5, 2.5),  # O on -x
		(2.5, 3.5, 2.5),  # O on +y
		(2.5, 1.5, 2.5),  # O on -y
		(2.5, 2.5, 3.5),  # O on +z
		(2.5, 2.5, 1.5),  # O on -z
	]
	symbols = ['Y'] + ['O'] * 6
	cell = [5, 5, 5]
	atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)

	cations = ['Y']
	anions = ['O']
	charges = [3]  # Charge for Y
	ap = 2.5  # Parameter to compute bond distance threshold

	print("Testing ionic system:")
	error_fraction_ionic = _test_bonds(atoms, cations, anions, charges, ap, lib, system_type='ionic', count_bonds=True)
	print("Ionic test error fraction:", error_fraction_ionic)

	print("\nTesting neutral system:")
	error_fraction_neutral = _test_bonds(atoms, cations, anions, charges, ap, lib, system_type='neutral',
	                                    count_bonds=True)
	print("Neutral test error fraction:", error_fraction_neutral)

	print("======")

	print("Testing ionic system:")
	error_fraction_ionic = test_bonds(atoms, cations, anions, charges, ap, lib, system_type='ionic', count_bonds=True)
	print("Ionic test error fraction:", error_fraction_ionic)

	print("\nTesting neutral system:")
	error_fraction_neutral = test_bonds(atoms, cations, anions, charges, ap, lib, system_type='neutral',
	                                    count_bonds=True)
	print("Neutral test error fraction:", error_fraction_neutral)



if __name__ == "__main__":
	main()
