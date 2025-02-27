import random

from ase import *
from torch.optim.radam import radam

import fuse202.possible_solutions
# from fuse_102.example_input_only.gulp.input_gulp_example import imax_atoms


# Function to generate a random string of atomic numbers which can then be used to
# create a structure within fuse

def create_random_string(
		cubic_solutions,
		tetragonal_solutions,
		hexagonal_solutions,
		orthorhombic_solutions,
		monoclinic_solutions,
		atoms_per_fu,
		fu,
		vac_ratio='',
		max_fus='',
		system_type='',
		composition='',
		ap=''
):
	### now need to rework this so that it is sensibly assembling structures in the
	# modular motifs from the original FUSE, the completely random strings are not
	# working very well!

	## try choosing the lattice type first, means we need to partially choose the
	# instructions before the string ##
	# 0 = cubic
	# 1 = tetragonal
	# 2 = hexagonal
	# 3 = orthorhombic
	# 4 = monoclinic
	# 5 = triclinic
	## 
	#############################################################################
	### first need to work out the maximum number of sub-modules ###
	# parameter set that worked!?
	# 1
	# atoms_per_fu:  11
	# fu:  [39, 39, 22, 22, 8, 8, 8, 8, 8, 8, 8]
	# vac_ratio:  4
	# max_fus:  1
	# system_type:  neutral
	# composition:  {'Y': 2, 'Ti': 2, 'O': 7}
	# ap:  4.276364
	# hello2
	random.seed(25)

	max_atoms = len(max_fus * fu)
	max_vac = vac_ratio * max_atoms
	max_sub = int((max_atoms + max_vac) / 4)

	### now need to choos the number of sub-modules, based on possible cell sizes
	accept = 0
	while accept == 0:
		target = 0
		latt_attemp = 0
		lattice = random.choice([0, 1, 2, 3, 4, 5])
		if latt_attemp >= 5000:
			lattice = random.choice([0, 1, 2, 3, 4, 5])
			# lattice=2
			latt_attemp = 0
		latt_attemp += 1

		# while target/len(fu) < 1:

		instructions = [ap]
		instructions.append(lattice)
		if lattice == 0:
			solutions = cubic_solutions
		if lattice == 1:
			solutions = tetragonal_solutions
		if lattice == 2:
			solutions = hexagonal_solutions
		if lattice == 3:
			solutions = orthorhombic_solutions
		if lattice == 4:
			solutions = monoclinic_solutions
		if lattice == 5:
			solutions = monoclinic_solutions

		nsubs_1 = solutions.keys()
		nsubs_1 = list(nsubs_1)
		nsubs_2 = []
		for i in range(len(nsubs_1)):
			if nsubs_1[i] <= max_sub:
				nsubs_2.append(nsubs_1[i])

		nsub = random.choice(nsubs_2)
		target = nsub * 4

		temp_max = int(target / len(fu))
		### set flag to see if we are completed #####################################
		#############################################################################
		## now we have to change this, such that it always returns a string containing the correct number of sub-modules/positions
		### start by generating initial string which we will use to create atoms object
		# choose number of formula units to use
		temp = list(range(1, temp_max + 1))
		if len(temp) > 0:
			fus = random.choice(temp)
		else:
			break
		# append n fus to an empty string
		string = []

		for i in range(fus):
			for j in range(len(fu)):
				string.append(fu[j])
		##########################################################################

		### Decorate with a random number of 120s (use 120 for blank space), maintinaing a multiple of 4 ##
		max_vac = vac_ratio * atoms_per_fu * fus
		max_vac -= string.count(120)
		diff = target - len(string)
		for i in range(diff):
			string.append(120)
		# print("\n")
		# print(string)
		# print("count 8s: ",string.count(8))

		# system_type='neutral'
		### now error check the strings currently only required for ionic structures 
		if system_type == 'neutral':
			random.shuffle(string)
			if len(string) == target:
				accept += 1

		if system_type == 'ionic':  # ideally rather than just randomising it,
			# current way of doing it to ensure that there's a cation at the origin of the sub-module
			symbols = list(composition.keys())
			cats = []
			for i in range(len(symbols)):
				temp = Atoms(symbols[i])
				if composition[symbols[i]][1] > 0:
					cats.append(list(temp.get_atomic_numbers())[0])

			number_of_modules = int(len(string) / 4)
			As = []  # split out the cations
			Bs = []  # everything else
			for i in range(len(string)):
				if string[i] in cats:
					As.append(string[i])
				else:
					Bs.append(string[i])
			if len(string) % 4 != 0:
				for i in range(4 - (len(string) % 4)):
					Bs.append(120)

			# print("number of modules: ",number_of_modules)
			if len(As) < number_of_modules:  # if we have fewer cations than submodules, need to move some vacancies from the anion string to the cation string
				difference = int(number_of_modules - len(As))
				avalible_vacancies = Bs.count(120)
				# print("differnce in As: ",difference)
				# print("vacancies in Bs: ",avalible_vacancies)
				if avalible_vacancies >= difference:
					for z in range(difference):
						As.append(120)
						Bs.remove(120)
			# print("As :",As)
			# print("Bs :",Bs)
			#
			if len(As) > number_of_modules:  # if we have more cations than submodules, need to move some cations accross to the anion string
				difference = int(len(As) - number_of_modules)
				cats_to_move = []
				for z in range(difference):
					chosen = random.choice(As)
					cats_to_move.append(chosen)
					As.remove(chosen)

			string = []
			random.shuffle(As)
			random.shuffle(Bs)
			for i in range(len(As)):
				string.append(As[-1])
				del (As[-1])
				for j in range(3):
					string.append(Bs[-1:][0])
					del (Bs[-1:])

			# shuffle(string)

			# print(string)
			# print("count 8s: ",string.count(8))
			# print("\n")

			# sys.exit()
			if len(string) == target:
				accept += 1

	if accept == 0:
		string = None

	return string, instructions
#


class RandomStructureGenerator:
	"""Generates a random structure based on lattice types and modular motifs."""

	def __init__(
		self,
		cubic_solutions: dict,
		tetragonal_solutions: dict,
		hexagonal_solutions: dict,
		orthorhombic_solutions: dict,
		monoclinic_solutions: dict,
		atoms_per_fu: int,
		fu: list[int],
		vac_ratio: int = 4,
		max_fus: int = 1,
		system_type: str = "neutral",
		composition: dict | None = None,
		ap: float | None = None
	):
		self.solutions_dict = {
			0: cubic_solutions,
			1: tetragonal_solutions,
			2: hexagonal_solutions,
			3: orthorhombic_solutions,
			4: monoclinic_solutions,
			5: monoclinic_solutions,  # Triclinic uses the same as monoclinic
		}
		self.atoms_per_fu = atoms_per_fu
		self.fu = fu
		self.vac_ratio = vac_ratio
		self.max_fus = max_fus
		self.system_type = system_type
		self.composition = composition or {}
		self.ap = ap

		self.max_atoms = len(self.max_fus * self.fu)
		self.max_vac = self.vac_ratio * self.max_atoms
		self.max_sub = (self.max_atoms + self.max_vac) // 4

		random.seed(25)

	def _select_random_lattice(self) -> int:
		"""Randomly selects a lattice type."""
		return random.choice(list(self.solutions_dict.keys()))

	def _select_number_of_submodules(self, lattice: int, max_sub: int) -> int | None:
		"""
		Selects a number of submodules based on possible cell sizes.
		If no possible cell sizes, returns None
		"""
		solutions = self.solutions_dict[lattice]

		valid_nsubs = [n for n in solutions.keys() if n <= max_sub]
		return random.choice(valid_nsubs) if valid_nsubs else None

	def _generate_initial_string(self, target: int) -> list[int]:
		"""Generates an initial atomic configuration string."""
		temp_max = target // len(self.fu)
		fus = random.choice(list(range(1, temp_max + 1))) if temp_max > 0 else None
		if fus is None:
			return []

		# Construct atomic sequence
		string = [atom for _ in range(fus) for atom in self.fu]

		# Fill vacancies (120 represents a blank space)
		vacancies_needed = target - len(string)
		string.extend([120] * vacancies_needed)

		return string

	def _error_check_ionic(self, string: list[int], target: int) -> list[int] | None:
		"""Performs additional validation for ionic structures."""
		symbols = list(self.composition.keys())
		cations = [Atoms(symbol).get_atomic_numbers()[0] for symbol in symbols if self.composition[symbol][1] > 0]

		number_of_modules = len(string) // 4
		As, Bs = [], []

		for element in string:
			if element in cations:
				As.append(element)
			else:
				Bs.append(element)

		# Ensure enough cations are available for the number of modules
		if len(As) < number_of_modules:
			missing_cations = number_of_modules - len(As)
			available_vacancies = Bs.count(120)
			if available_vacancies >= missing_cations:
				As.extend([120] * missing_cations)
				Bs = [b for b in Bs if b != 120]

		# Ensure we don’t have excess cations
		if len(As) > number_of_modules:
			extra_cations = len(As) - number_of_modules
			for _ in range(extra_cations):
				As.remove(random.choice(As))

		# Rebuild the string with submodule placement
		new_string = []
		random.shuffle(As)
		random.shuffle(Bs)

		for _ in range(len(As)):
			new_string.append(As.pop())
			new_string.extend(Bs[-3:])
			del Bs[-3:]

		return new_string if len(new_string) == target else None

	def create_random_string(self) -> tuple[list[int] | None, list]:
		"""Main function to generate a random structure string and associated instructions."""

		accept = False
		instructions = [self.ap]

		while not accept:
			lattice = self._select_random_lattice()
			instructions.append(lattice)

			nsub = self._select_number_of_submodules(lattice, self.max_sub)
			if not nsub:
				return None, instructions  # No valid solution found

			target = nsub * 4
			string = self._generate_initial_string(target)
			if not string:
				return None, instructions

			# Perform error checks based on system type
			if self.system_type == "neutral":
				random.shuffle(string)
				accept = len(string) == target

			elif self.system_type == "ionic":
				string = self._error_check_ionic(string, target)
				accept = bool(string)

		return string, instructions


import statistics
from decimal import Decimal
import fuse202.possible_solutions
from fuse202.bond_table import BOND_DATA


def main():

	max_ax = 40
	solutions = fuse202.possible_solutions.possible_solutions(max_ax)
	cubic_solutions, tetragonal_solutions, hexagonal_solutions, orthorhombic_solutions, monoclinic_solutions = solutions
	print(cubic_solutions)


	imax_atoms = 50
	max_atoms = 50

	composition = {'Ca': 3, 'Ti': 2, 'O': 7}
	atoms_per_fu = sum(composition.values())
	imax_fus: int = imax_atoms // atoms_per_fu
	max_fus: int = max_atoms // atoms_per_fu

	target_fu = random.choice(list(range(1, imax_fus + 1)))

	fu = []
	for element, count in composition.items():
		atomic_number = Atoms(element).get_atomic_numbers()[0]
		fu.extend([atomic_number] * count)

	print(fu)

	symbol_fu = [Atoms(numbers=[v]).get_chemical_symbols()[0] for v in fu]


	vac_ratio = 4
	system_type = "neutral"
	ap_scale = ""
	bondtable = BOND_DATA
	fixed_ap = ""
	ap = cal_ap(ap_scale, bondtable, fixed_ap, symbol_fu, system_type)

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
		ap=ap
	)
	print(string)
	print(instructions)

	generator = RandomStructureGenerator(
		cubic_solutions, tetragonal_solutions, hexagonal_solutions,
		orthorhombic_solutions, monoclinic_solutions,
		atoms_per_fu, fu, vac_ratio=vac_ratio, max_fus=max_fus, system_type='neutral',
		composition=composition, ap=ap
	)

	string_, instructions_ = generator.create_random_string()
	print("Generated Structure:", string_)
	print("Instructions:", instructions_)

	print(string_ == string)
	print(instructions_ == instructions)


def cal_ap(bondtable: dict, symbol_fu: list, system_type: str, ap_scale, fixed_ap) -> float:

	if system_type != "neutral":
		return 0.0

	ap = statistics.fmean(
		statistics.fmean(
			bound_value[-1] for bound_value in bondtable[symbol].values()
		) for symbol in symbol_fu
	)
	ap = 4 * ap
	if ap_scale:
		ap = ap * ap_scale
	if fixed_ap:
		ap = fixed_ap

	return float(Decimal(ap).quantize(Decimal('1e-4')))


if __name__ == "__main__":
	main()
