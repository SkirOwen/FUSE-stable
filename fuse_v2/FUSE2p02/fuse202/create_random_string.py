import sys
from ase import *
from random import choice
from random import shuffle


# Function to generate a random string of atomic numbers which can then be used to
# create a structure within fuse

def create_random_string(cubic_solutions, tetragonal_solutions, hexagonal_solutions, orthorhombic_solutions,
                         monoclinic_solutions, atoms_per_fu, fu, vac_ratio='', max_fus='', system_type='',
                         composition='', ap=''):
	### now need to rework this so that it is sensibly assembling structures in the
	# modular motifs from the original FUSE, the completly random strings are not
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

	max_atoms = len(max_fus * fu)
	max_vac = vac_ratio * max_atoms
	max_sub = int((max_atoms + max_vac) / 4)
	### now need to choos the number of sub-modules, based on possible cell sizes
	accept = 0
	while accept == 0:
		target = 0
		latt_attemp = 0
		lattice = choice([0, 1, 2, 3, 4, 5])
		if latt_attemp >= 5000:
			lattice = choice([0, 1, 2, 3, 4, 5])
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

		nsub = choice(nsubs_2)
		target = nsub * 4

		temp_max = int(target / len(fu))
		### set flag to see if we are completed #####################################
		#############################################################################
		## now we have to change this, such that it always returns a string containing the correct number of sub-modules/positions
		### start by generating initial string which we will use to create atoms object
		# choose number of formula units to use
		temp = list(range(1, temp_max + 1))
		if len(temp) > 0:
			fus = choice(temp)
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
			shuffle(string)
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
					chosen = choice(As)
					cats_to_move.append(chosen)
					As.remove(chosen)

			string = []
			shuffle(As)
			shuffle(Bs)
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
	# print (string)
	# print(len(string))
	# print(lattice)
	return string, instructions
#############################################################################
