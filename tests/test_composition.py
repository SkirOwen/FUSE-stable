"""Tests for the Composition domain type."""
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fuse202.domain.composition import Composition

ELEMENTS = st.sampled_from(["H", "Li", "O", "Na", "Al", "Si", "Ca", "Ti", "Fe", "Sr", "Y", "Ba"])


def test_counts_and_symbols():
	composition = Composition({"Ca": 3, "Ti": 2, "O": 7})
	assert composition.symbols == ["Ca", "Ti", "O"]
	assert composition.atoms_per_formula_unit == 12


def test_atomic_numbers_expands_one_entry_per_atom():
	composition = Composition({"Sr": 1, "Ti": 1, "O": 3})
	assert composition.atomic_numbers == [38, 22, 8, 8, 8]


def test_atomic_numbers_are_plain_python_ints():
	"""ase returns numpy integers, which propagate into the structure
	representation and broke atom counting on Python 3.12 and later."""
	numbers = Composition({"Ca": 3, "Ti": 2, "O": 7}).atomic_numbers
	assert all(type(number) is int for number in numbers)


def test_formula_round_trip():
	assert str(Composition.from_formula("Ca3Ti2O7")) == "Ca3Ti2O7"
	assert Composition.from_formula("Ca3Ti2O7").counts == {"Ca": 3, "Ti": 2, "O": 7}


def test_formula_treats_a_bare_symbol_as_one_atom():
	assert Composition.from_formula("SrTiO3").counts == {"Sr": 1, "Ti": 1, "O": 3}


@pytest.mark.parametrize("bad", ["", "Xx2", "3Ca", "Ca3!"])
def test_unparseable_formulas_are_rejected(bad):
	with pytest.raises(ValueError):
		Composition.from_formula(bad)


@pytest.mark.parametrize("bad", [{}, {"Xx": 1}, {"Ca": 0}, {"Ca": -2}, {"Ca": 1.5}])
def test_invalid_compositions_are_rejected(bad):
	with pytest.raises(ValueError):
		Composition(bad)


def test_normalised_reduces_to_the_smallest_whole_ratio():
	assert Composition({"Ca": 6, "Ti": 4, "O": 14}).normalised().counts == {"Ca": 3, "Ti": 2, "O": 7}


def test_normalised_leaves_an_already_reduced_composition_alone():
	composition = Composition({"Ca": 3, "Ti": 2, "O": 7})
	assert composition.normalised().counts == composition.counts


def test_formula_units_in_counts_whole_units():
	composition = Composition({"Sr": 1, "Ti": 1, "O": 3})
	two_units = ["Sr", "Ti", "O", "O", "O"] * 2
	assert composition.formula_units_in(two_units) == 2.0


def test_formula_units_in_rejects_wrong_proportions():
	composition = Composition({"Sr": 1, "Ti": 1, "O": 3})
	assert composition.formula_units_in(["Sr", "Ti", "O"]) is None
	assert composition.formula_units_in(["Sr", "Sr", "Ti", "O", "O", "O"]) is None


def test_atom_counts_reports_zero_for_absent_elements():
	composition = Composition({"Sr": 1, "O": 3})
	assert composition.atom_counts_in(["O", "O"]) == {"Sr": 0, "O": 2}


@settings(max_examples=40, deadline=None)
@given(
	elements=st.lists(ELEMENTS, min_size=1, max_size=4, unique=True),
	counts=st.lists(st.integers(min_value=1, max_value=9), min_size=1, max_size=4),
)
def test_atomic_numbers_length_always_matches_the_atom_count(elements, counts):
	composition = Composition(dict(zip(elements, counts)))
	assert len(composition.atomic_numbers) == composition.atoms_per_formula_unit


@settings(max_examples=40, deadline=None)
@given(
	elements=st.lists(ELEMENTS, min_size=1, max_size=4, unique=True),
	counts=st.lists(st.integers(min_value=1, max_value=9), min_size=1, max_size=4),
	multiplier=st.integers(min_value=1, max_value=6),
)
def test_scaling_a_composition_does_not_change_its_normalised_form(elements, counts, multiplier):
	"""Ca3Ti2O7 and Ca6Ti4O14 describe the same material."""
	base = dict(zip(elements, counts))
	scaled = {symbol: count * multiplier for symbol, count in base.items()}
	assert Composition(base).normalised().counts == Composition(scaled).normalised().counts


@settings(max_examples=40, deadline=None)
@given(
	elements=st.lists(ELEMENTS, min_size=1, max_size=4, unique=True),
	counts=st.lists(st.integers(min_value=1, max_value=9), min_size=1, max_size=4),
)
def test_formula_string_round_trips_through_the_parser(elements, counts):
	composition = Composition(dict(zip(elements, counts)))
	assert Composition.from_formula(str(composition)).counts == composition.counts
