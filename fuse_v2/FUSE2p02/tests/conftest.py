import pytest
from ase import Atoms


@pytest.fixture
def octahedral_yo6() -> Atoms:
	"""A single YO6 octahedron: Y at the centre, one O at +-x/+-y/+-z, each 1.0 A away.

	Lifted from fuse202/test_bonds.py's own hand-written main() smoke test, which
	built this exact structure to sanity-check bond counting by eye before there
	was a real test suite to do it.
	"""
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
	return Atoms(symbols=symbols, positions=positions, cell=[5, 5, 5], pbc=True)
