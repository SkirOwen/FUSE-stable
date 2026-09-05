from fuse202 import bond_table
from fuse202.bond_table import BOND_DATA


def test_bond_data_has_expected_entry_shape():
	# format: {symbol: {oxidation_state: [min_bonds, max_bonds, atomic_radius]}}
	assert BOND_DATA['O'][0] == [2, 8, 0.74]


def test_bond_data_every_entry_has_sane_min_max_radius():
	for symbol, charge_states in BOND_DATA.items():
		for charge, (min_bonds, max_bonds, radius) in charge_states.items():
			assert min_bonds <= max_bonds, f"{symbol}[{charge}] has min > max"
			assert radius > 0, f"{symbol}[{charge}] has non-positive radius"


def test_bond_data_only_ever_has_a_single_oxidation_state():
	"""KNOWN LIMITATION (not a crash, just documenting current data): the module
	docstring describes BOND_DATA as keyed by oxidation state ("2nd key(s):
	charge states"), implying multiple entries per element, but every element
	in the table today has exactly one entry under key 0. Anything that relies
	on BOND_DATA to distinguish real oxidation states (e.g. ionic bond
	checking) is working with single-state data in practice - see
	test_bond_checking.py's ionic KeyError test for the concrete failure mode
	this causes."""
	assert all(list(charge_states.keys()) == [0] for charge_states in BOND_DATA.values())


def test_initialize_and_load_bond_table_roundtrip(tmp_path, monkeypatch):
	npz_path = tmp_path / "bondtable.npz"
	monkeypatch.setattr(bond_table, "BOND_TABLE_FILE", str(npz_path))

	assert not npz_path.exists()
	bond_table.initialize_bond_table()
	assert npz_path.exists()

	loaded = bond_table.load_bond_table()
	assert loaded == BOND_DATA


def test_load_bond_table_creates_file_if_missing(tmp_path, monkeypatch):
	npz_path = tmp_path / "bondtable.npz"
	monkeypatch.setattr(bond_table, "BOND_TABLE_FILE", str(npz_path))

	assert not npz_path.exists()
	loaded = bond_table.load_bond_table()
	assert npz_path.exists()
	assert loaded == BOND_DATA
