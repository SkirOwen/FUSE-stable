"""Tests for the Structure domain type, including the dict compatibility bridge."""
import pytest
from ase import Atoms

from fuse202.domain.structure import Structure


@pytest.fixture
def srtio3():
	return Atoms(
		symbols=["Sr", "Ti", "O", "O", "O"],
		scaled_positions=[(0, 0, 0), (.5, .5, .5), (.5, .5, 0), (.5, 0, .5), (0, .5, .5)],
		cell=[3.905, 3.905, 3.905],
		pbc=True,
	)


@pytest.fixture
def structure_dict(srtio3):
	"""The dict form the search passes around, as extract_module returns it."""
	return {
		"modules": [srtio3.copy()],
		"sub module cell": [3.905, 3.905, 1.95, 90, 90, 90],
		"shape in submods": [1, 1, 2, 90, 90, 90],
		"nmods": 2,
		"ap": 3.905,
		"atoms": srtio3,
	}


def test_from_dict_then_to_dict_preserves_the_original(structure_dict):
	restored = Structure.from_dict(structure_dict).to_dict()
	for key, value in structure_dict.items():
		assert restored[key] is value or restored[key] == value


def test_attribute_and_key_access_agree(structure_dict):
	structure = Structure.from_dict(structure_dict)
	assert structure.sub_module_cell == structure["sub module cell"]
	assert structure.shape_in_submods == structure["shape in submods"]
	assert structure.optimised == structure["optimised?"]


def test_assignment_through_either_route_is_visible_from_both(structure_dict):
	structure = Structure.from_dict(structure_dict)
	structure["energy"] = -7.5
	assert structure.energy == -7.5
	structure.optimised = True
	assert structure["optimised?"] is True


def test_the_keys_the_search_writes_after_relaxing_all_work(structure_dict):
	"""run_fuse sets exactly these after a calculator returns."""
	structure = Structure.from_dict(structure_dict)
	structure["optimised?"] = True
	structure["energy"] = -6.9
	structure["converged"] = True
	structure["time"] = 1.4
	assert (structure.optimised, structure.energy, structure.converged, structure.time) == (
		True, -6.9, True, 1.4)


def test_an_unknown_key_raises_KeyError_like_a_dict(structure_dict):
	structure = Structure.from_dict(structure_dict)
	with pytest.raises(KeyError):
		structure["not a real key"]
	with pytest.raises(KeyError):
		structure["not a real key"] = 1


def test_from_dict_ignores_keys_it_does_not_know(structure_dict):
	"""Pre built pool entries carry extra keys; those must not break conversion."""
	structure_dict["used?"] = False
	assert Structure.from_dict(structure_dict).ap == 3.905


def test_unset_optional_fields_stay_out_of_the_dict_form(structure_dict):
	"""A structure that has not been relaxed never had a "time" key, and gaining
	one would change what is written to the restart files."""
	assert "time" not in Structure.from_dict(structure_dict).to_dict()
	assert "source" not in Structure.from_dict(structure_dict).to_dict()


def test_copy_is_independent(structure_dict):
	original = Structure.from_dict(structure_dict)
	duplicate = original.copy()
	duplicate.atoms.translate([1.0, 0, 0])
	duplicate.energy = -1.0
	duplicate.sub_module_cell[0] = 99

	assert original.energy != duplicate.energy
	assert original.sub_module_cell[0] == 3.905
	assert original.atoms.positions[0][0] != duplicate.atoms.positions[0][0]


def test_defaults_match_what_the_search_expects_of_a_fresh_structure():
	"""make_basin_move sets these three before returning a trial structure."""
	fresh = Structure()
	assert fresh.optimised is False
	assert fresh.energy == 0.0
	assert fresh.converged is False


def test_survives_a_pickle_round_trip(structure_dict):
	"""Structures are pickled into the restart files, so a resumed run depends
	on this."""
	import pickle

	original = Structure.from_dict(structure_dict)
	original.energy = -7.1
	restored = pickle.loads(pickle.dumps(original))

	assert restored.energy == -7.1
	assert restored["sub module cell"] == original.sub_module_cell
	assert len(restored.atoms) == len(original.atoms)


def test_code_written_for_dicts_works_unchanged_on_a_Structure(structure_dict):
	"""The whole point of the key mapping. A restart file written by v2.3 holds
	plain dicts and one written now holds Structures, so the search has to
	handle both. This is the shape of every existing call site."""

	def as_the_search_does(record):
		record["optimised?"] = True
		record["energy"] = -6.5
		return len(record["atoms"]), record["ap"], record["optimised?"]

	from_dict_form = as_the_search_does(dict(structure_dict))
	from_structure = as_the_search_does(Structure.from_dict(structure_dict))
	assert from_dict_form == from_structure
