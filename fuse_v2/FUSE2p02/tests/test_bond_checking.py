import pytest

from fuse202.bond_table import BOND_DATA
# Aliased on import: fuse202.test_bonds.test_bonds is production bond-checking
# code, not a pytest test, but pytest's default collection still picks up any
# `test_*`-named callable sitting in a test module's namespace - importing it
# under its real name here previously made pytest try to collect and call it
# directly, and fail on missing fixtures. See test_bond_table.py's
# "only ever a single oxidation state" test and this file's KeyError test for
# what the function itself is actually being characterized on.
from fuse202.test_bonds import test_bonds as check_bonds


def test_neutral_octahedral_structure_has_no_bonding_errors(octahedral_yo6):
	# Y-O distance is 1.0 A; ap=2.5 -> bond cutoff 1.9625 A, so all 6 O count as
	# bonded. BOND_DATA['Y'][0] = [4, 12, ...], and 6 is within that range.
	error_fraction = check_bonds(
		octahedral_yo6,
		cations=['Y'],
		anions=['O'],
		charges=[0],
		ap=2.5,
		lib=BOND_DATA,
		system_type='neutral',
	)
	assert error_fraction == 0.0


def test_neutral_flags_undercoordinated_cation(octahedral_yo6):
	# Shrink the bond cutoff so only the closest O counts as bonded - 1 bond is
	# outside BOND_DATA['Y'][0]'s [4, 12] range, so this should be flagged wrong.
	error_fraction = check_bonds(
		octahedral_yo6,
		cations=['Y'],
		anions=['O'],
		charges=[0],
		ap=0.5,  # bond cutoff 0.39 A - below every Y-O distance (1.0 A)
		lib=BOND_DATA,
		system_type='neutral',
	)
	assert error_fraction == 1.0


def test_ionic_lookup_raises_keyerror_for_a_charge_not_in_the_bond_table(octahedral_yo6):
	"""KNOWN BUG (not yet fixed): test_bonds()'s ionic branch does
	`charges[cations.index(symbol)]` then `lib[symbol][cation_charge]`, and only
	catches IndexError around that pair. BOND_DATA has exactly one oxidation
	state (key 0) per element (see test_bond_table.py), so any ionic charge
	other than 0 misses the dict lookup and raises an uncaught KeyError instead
	of falling back gracefully.

	This file's own dead `_test_bonds()` (superseded by `test_bonds()`, but left
	in place - see error_check_structure.py, which imports test_bonds, not
	_test_bonds) caught KeyError here and fell back to averaging the bond range
	across whatever charge states *were* present. That fallback did not carry
	over to test_bonds(). Pinning this as a known bug rather than fixing it now
	- fixing it is a Phase 4 decision once run_fuse()'s ionic path (system_type
	is basically only ever 'neutral' in practice today) is actually exercised
	and it's clear what the right fallback behavior should be.
	"""
	with pytest.raises(KeyError):
		check_bonds(
			octahedral_yo6,
			cations=['Y'],
			anions=['O'],
			charges=[3],  # BOND_DATA['Y'] only has key 0, not 3
			ap=2.5,
			lib=BOND_DATA,
			system_type='ionic',
		)


def test_neutral_ignores_charges_argument(octahedral_yo6):
	"""Documents that `charges` is accepted but never read in the 'neutral'
	branch - only the 'ionic' branch uses it. A nonsense charges list changes
	nothing for a neutral-system call."""
	result_with_real_charge = check_bonds(
		octahedral_yo6, cations=['Y'], anions=['O'], charges=[0], ap=2.5,
		lib=BOND_DATA, system_type='neutral',
	)
	result_with_bogus_charge = check_bonds(
		octahedral_yo6, cations=['Y'], anions=['O'], charges=[999], ap=2.5,
		lib=BOND_DATA, system_type='neutral',
	)
	assert result_with_real_charge == result_with_bogus_charge
