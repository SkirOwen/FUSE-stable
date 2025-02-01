from __future__ import annotations

import itertools as it
import os
import sys

from typing import Iterable

from fuse202.config import SPP_PATH


def assemble_spp(elements: Iterable, output='lib.lib', spp_path=SPP_PATH):
	# os.environ['SPP_PATH']=spp_path
	# print(os.environ['SPP_PATH'])

	# First enumerate the pairs of potentials:
	pairs = list(it.product(elements, repeat=2))

	with open(output, 'w') as o:

		for elem1, elem2 in pairs:
			try:
				f = open(
					os.path.join(spp_path, f"{elem1.upper()}-{elem2.upper()}", f"{elem1.upper()}-{elem2.upper()}.POT"),
					'r').readlines()
			except:
				f = open(
					os.path.join(spp_path, f"{elem2.upper()}-{elem1.upper()}", f"{elem2.upper()}-{elem1.upper()}.POT"),
					'r').readlines()
			o.write("\n")

			for j in f:
				o.write(j)


def get_spp_reader(elements: Iterable, spp_path: str = SPP_PATH):
	"""
	Creates a dispatcher function that reads the correct SPP (Semiempirical Pseudopotential) file on demand.

	:param elements: Iterable of element symbols.
	:param spp_path: Path where SPP files are stored.
	:return: A function that takes (elem1, elem2) and returns the file contents or None if not found.
	"""
	pairs = {tuple(sorted((e1, e2))): None for e1, e2 in it.product(elements, repeat=2)}

	def fetch_potential(elem1: str, elem2: str) -> str | None:
		"""Fetch the potential file content for the given element pair."""
		key = tuple(sorted((elem1, elem2)))
		if key not in pairs:
			print(f"Warning: Invalid element pair ({elem1}, {elem2})")
			return None

		pair_dir = os.path.join(spp_path, f"{key[0].upper()}-{key[1].upper()}")
		pot_file = os.path.join(pair_dir, f"{key[0].upper()}-{key[1].upper()}.POT")

		if not os.path.isfile(pot_file):
			print(f"Warning: Potential file not found for pair ({elem1}, {elem2})")
			return None

		with open(pot_file, 'r') as f:
			return f.read()

	return fetch_potential


# assemble_spp(elements)

def main():
	elements = ['Sr', 'Ti', 'Y', 'O']
	fetch_potential = get_spp_reader(elements, spp_path='/path/to/spp')

	content = fetch_potential('H', 'O')
	if content:
		print(content)  # Prints the potential file's content for H-O


if __name__ == "__main__":
	main()
