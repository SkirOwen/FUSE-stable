"""Tests for the small, largely pure helper functions defined at the bottom of
run_fuse.py (get_fu, cal_ap, cal_ideal_density). Importing fuse202.run_fuse
pulls in its full dependency chain (torch, chgnet, ase, ...) including two
eager CHGNet model loads at import time via mutable default arguments
(`rel=StructOptimizer()` on both run_fuse() and run_chgnet()) - that's a real
cost of testing anything in this module today, and a concrete argument for
Phase 4's run_fuse() breakup to stop doing that.
"""
from fuse202.bonds.bond_table import BOND_DATA
from fuse202.run_fuse import cal_ap, cal_ideal_density, get_fu


def test_get_fu_expands_composition_to_one_atomic_number_per_atom():
	# SrTiO3 -> Sr (Z=38), Ti (Z=22), O x3 (Z=8)
	fu = get_fu({'Sr': 1, 'Ti': 1, 'O': 3})
	assert sorted(fu) == sorted([38, 22, 8, 8, 8])


def test_get_fu_respects_composition_counts():
	fu = get_fu({'Na': 2, 'Cl': 2})
	assert sorted(fu) == [11, 11, 17, 17]


def test_cal_ap_neutral_system_is_positive():
	ap = cal_ap(BOND_DATA, ['Y', 'O'], system_type='neutral', ap_scale='', fixed_ap='')
	assert ap > 0


def test_cal_ap_non_neutral_system_type_short_circuits_to_zero():
	# cal_ap only computes a real value for system_type == 'neutral' - every
	# other value (including 'ionic') returns 0.0 unconditionally today.
	ap = cal_ap(BOND_DATA, ['Y', 'O'], system_type='ionic', ap_scale='', fixed_ap='')
	assert ap == 0.0


def test_cal_ap_fixed_ap_overrides_computed_value():
	ap = cal_ap(BOND_DATA, ['Y', 'O'], system_type='neutral', ap_scale='', fixed_ap=9.999)
	assert ap == 9.999


def test_cal_ap_scale_multiplies_computed_value():
	base_ap = cal_ap(BOND_DATA, ['Y', 'O'], system_type='neutral', ap_scale='', fixed_ap='')
	scaled_ap = cal_ap(BOND_DATA, ['Y', 'O'], system_type='neutral', ap_scale=2.0, fixed_ap='')
	assert scaled_ap == round(base_ap * 2.0, 4)


def test_cal_ideal_density_positive_for_known_composition():
	fu = get_fu({'Sr': 1, 'Ti': 1, 'O': 3})
	density = cal_ideal_density(BOND_DATA, fu)
	assert density > 0
