import glob
import os
import pickle
import platform
from collections import defaultdict

from fuse202.config import SOLUTIONS_DIR


def return_factors(n: int) -> set:
	"""Return all factors of a given number n."""
	return {i for i in range(1, int(n ** 0.5) + 1) if n % i == 0 for i in (i, n // i)}


def cube_function(x):
	"""Calculate the number of sub-modules for cubic lattices."""
	return 2 * (x ** 3)


def tetragonal_function(x, z):
	"""Calculate the number of sub-modules for tetragonal lattices."""
	return (x ** 2) * z


def orthorhombic_function(x, w, z):
	"""Calculate the number of sub-modules for orthorhombic/monoclinic/triclinic lattices."""
	return x * w * z

#
# def possible_solutions(max_ax, restart: bool):
# 	if not restart:
# 		if platform.system == 'Windows':
# 			temp = glob.glob("*.p")
# 			for i in range(len(temp)):
# 				os.remove(temp[i])
# 		if platform.system in ['Linux', "darwin"] :
# 			os.system("rm *.p")
#
# 	def load_or_generate(filename: str, generator_func):
# 		"""Helper function to load data from cache or generate and save it."""
# 		try:
# 			return pickle.load(open(filename, 'rb'))
# 		except (FileNotFoundError, EOFError):
# 			data = generator_func()
# 			pickle.dump(data, open(filename, 'wb'))
# 			return data
#
# 	# is cubic possible? #########################################################
# 	# cubic is defined as integer solutions to (n**2 * 2*n)-y = 0
# 	# below, is the total number of sub-modules required to make cubic lattices, for n = 1 - 50
# 	cubic_solutions = load_or_generate(
# 		"cubes.p",
# 		lambda: {cube_function(i): [i, i, 2 * i] for i in range(1, max_ax + 1)}
# 	)
#
#
# 	# possible tetragonal solutions? ############################################
# 	# below is a dictionary of all possible x and z pairs for tetragonal cells for y sub-modules, for x & z upto 50.
# 	# some numbers of sub modules have more than one solution, which is why each list is stored as a pair
# 	try:
# 		tetragonal_solutions = pickle.load(open("tetragonal.p", 'rb'))
# 	except:
# 		tetragonal_solutions = {}
# 		for i in range(1, max_ax + 1):
# 			for j in range(2, max_ax + 1, 2):
# 				y = tetragonal_function(i, j)
# 				try:
# 					tetragonal_solutions[y].append([i, j])
# 				except:
# 					tetragonal_solutions[y] = [[i, j]]
#
# 		pickle.dump(tetragonal_solutions, open("tetragonal.p", 'wb'))
#
# 	#############################################################################
#
# 	# possible hexagonal solutions ##############################################
# 	# produced from the same as above, but without the z = even restriction
# 	# some numbers of sub modules have more than one solution, which is why each list is stored as a pair
# 	try:
# 		hexagonal_solutions = pickle.load(open("hexagonal.p", 'rb'))
# 	except:
# 		hexagonal_solutions = {}
# 		for i in range(1, max_ax + 1):
# 			for j in range(2, max_ax + 1):
# 				y = tetragonal_function(i, j)
# 				try:
# 					hexagonal_solutions[y].append([i, j])
# 				except:
# 					hexagonal_solutions[y] = [[i, j]]
#
# 		pickle.dump(hexagonal_solutions, open("hexagonal.p", 'wb'))
#
# 	#############################################################################
#
# 	# possible orthorhombic solutions ###########################################
# 	# below are the possible orthorhombic solutions, stored as above, with the restriction
# 	# that z = even upto 50 sub-mods in every direction have to generate them on the
# 	# fly as having the full list in the file kills jedit!
# 	try:
# 		orthorhombic_solutions = pickle.load(open("orthorhombic.p", 'rb'))
# 	except:
# 		orthorhombic_solutions = {}
# 		for i in range(1, max_ax + 1):
# 			for j in range(1, max_ax + 1):
# 				for k in range(2, max_ax + 1, 2):
# 					y = orthorhombic_function(i, j, k)
# 					try:
# 						orthorhombic_solutions[y].append([i, j, k])
# 					except:
# 						orthorhombic_solutions[y] = [[i, j, k]]
# 		pickle.dump(orthorhombic_solutions, open("orthorhombic.p", 'wb'))
#
# 	# print(orthorhombic_solutions)
#
# 	#############################################################################
#
# 	# possible monoclinic & triclinic dimensions, as orthorhombic but without the
# 	# z = evan restriction
# 	try:
# 		monoclinic_solutions = pickle.load("monoclinic.p", 'rb')
# 	except:
# 		monoclinic_solutions = {}
# 		for i in range(1, max_ax + 1):
# 			for j in range(1, max_ax + 1):
# 				for k in range(1, max_ax + 1):
# 					y = orthorhombic_function(i, j, k)
# 					try:
# 						monoclinic_solutions[y].append([i, j, k])
# 					except:
# 						monoclinic_solutions[y] = [[i, j, k]]
# 		pickle.dump(monoclinic_solutions, open("monoclinic.p", 'wb'))
#
# 	# print(orthorhombic_solutions)
#
# 	#############################################################################
#
# 	return cubic_solutions, tetragonal_solutions, hexagonal_solutions, orthorhombic_solutions, monoclinic_solutions


class CrystalSolutionsCalculator:
	def __init__(self, max_ax: int, solutions_dir=SOLUTIONS_DIR):
		self.max_ax = max_ax
		self.solutions_dir = solutions_dir

	def _clear_pickle_files(self) -> None:
		"""Remove all pickle files in the solutions directory."""
		for file in glob.glob(os.path.join(self.solutions_dir, "*.p")):
			os.remove(file)

	def _load_or_create_solutions(self, filename: str, calculator_func) -> dict:
		"""Load existing solutions from pickle file or create new ones."""
		filepath = os.path.join(self.solutions_dir, filename)
		try:
			with open(filepath, 'rb') as f:
				return pickle.load(f)
		except (FileNotFoundError, EOFError):
			solutions = calculator_func()
			with open(filepath, 'wb') as f:
				pickle.dump(solutions, f)
			return solutions

	def _calculate_cubic_solutions(self) -> dict[int, list[int]]:
		"""Calculate all possible cubic solutions."""
		solutions = {}
		for i in range(1, self.max_ax + 1):
			solutions[cube_function(i)] = [i, i, 2 * i]
		return solutions

	def _calculate_2d_solutions(self, step: int = 1) -> dict[int, list[list[int]]]:
		"""Calculate solutions for tetragonal/hexagonal systems."""
		solutions = {}
		for i in range(1, self.max_ax + 1):
			for j in range(2, self.max_ax + 1, step):
				y = tetragonal_function(i, j)
				solutions.setdefault(y, []).append([i, j])
		return solutions

	def _calculate_3d_solutions(self, step: int = 1) -> dict[int, list[list[int]]]:
		"""Calculate solutions for orthorhombic/monoclinic systems."""
		solutions = {}
		for i in range(1, self.max_ax + 1):
			for j in range(1, self.max_ax + 1):
				for k in range(1, self.max_ax + 1, step):
					y = orthorhombic_function(i, j, k)
					solutions.setdefault(y, []).append([i, j, k])
		return solutions

	def calculate_all_solutions(self, restart: bool = False) -> tuple[dict, ...]:
		"""
		Calculate or load solutions for all crystal systems.

		Args:
			restart (bool): If True, delete existing solution files and recalculate

		Returns:
			Tuple of dictionaries containing solutions for each crystal system
		"""
		if restart:
			self._clear_pickle_files()

		cubic = self._load_or_create_solutions("cubes.p", self._calculate_cubic_solutions)
		tetragonal = self._load_or_create_solutions("tetragonal.p", lambda: self._calculate_2d_solutions(step=2))
		hexagonal = self._load_or_create_solutions("hexagonal.p", lambda: self._calculate_2d_solutions(step=1))
		orthorhombic = self._load_or_create_solutions("orthorhombic.p", lambda: self._calculate_3d_solutions(step=2))
		monoclinic = self._load_or_create_solutions("monoclinic.p", lambda: self._calculate_3d_solutions(step=1))

		return cubic, tetragonal, hexagonal, orthorhombic, monoclinic


def possible_solutions(max_ax: int, restart: bool = False) -> tuple[dict, ...]:
	"""
	Calculate possible solutions for various crystal systems.

	Args:
		max_ax (int): Maximum axis length to consider
		restart (bool): If True, delete existing solution files and recalculate

	Returns:
		Tuple containing dictionaries of solutions for cubic, tetragonal,
		hexagonal, orthorhombic, and monoclinic crystal systems
	"""
	calculator = CrystalSolutionsCalculator(max_ax)
	return calculator.calculate_all_solutions(restart)
