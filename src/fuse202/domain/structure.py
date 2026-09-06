"""A candidate crystal structure as it moves through the search.

FUSE builds structures from sub modules tiled through a unit cell, relaxes them
with an energy calculator, and keeps or discards them. This type carries that
candidate and everything recorded about it along the way.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any

# Maps the string keys the search has always used to attribute names. The keys
# are part of the on disk restart format and of every existing input script's
# expectations, so they are kept exactly as they are.
_KEY_TO_FIELD = {
	"atoms": "atoms",
	"modules": "modules",
	"sub module cell": "sub_module_cell",
	"shape in submods": "shape_in_submods",
	"nmods": "nmods",
	"ap": "ap",
	"optimised?": "optimised",
	"energy": "energy",
	"converged": "converged",
	"time": "time",
	"source": "source",
}
_FIELD_TO_KEY = {value: key for key, value in _KEY_TO_FIELD.items()}


@dataclass
class Structure:
	"""One candidate structure and the state the search keeps about it.

	Parameters
	----------
	atoms : ase.Atoms
		The structure itself.
	modules : list of ase.Atoms
		The sub modules it was assembled from.
	sub_module_cell : list
		Sub module cell as [a, b, c, alpha, beta, gamma].
	shape_in_submods : list
		Unit cell shape counted in sub modules, as
		[x, y, z, alpha, beta, gamma].
	nmods : int
		Number of sub modules in the cell.
	ap : float
		Sub module lattice parameter in Angstroms.
	optimised : bool
		Whether a calculator has relaxed this structure yet.
	energy : float
		Energy in eV per atom once relaxed, otherwise zero.
	converged : bool
		Whether the calculator reported convergence.
	time : float or None
		Seconds the relaxation took, when it has run.
	source : str or None
		Where the structure came from, for pre built structures read from disk.

	Notes
	-----
	Instances also support `structure["sub module cell"]` and item assignment,
	using the same string keys the search has always used. The search currently
	passes plain dicts with those keys through several thousand lines of code
	and writes them into restart files, so the mapping is kept while call sites
	move over to attribute access one at a time. `to_dict` and `from_dict`
	convert between the two forms.

	The class is mutable because the search assigns to fields such as `energy`
	and `optimised` in place, in around thirty places.
	"""

	atoms: Any = None
	modules: list = field(default_factory=list)
	sub_module_cell: list = field(default_factory=list)
	shape_in_submods: list = field(default_factory=list)
	nmods: int = 0
	ap: float = 0.0
	optimised: bool = False
	energy: float = 0.0
	converged: bool = False
	time: float | None = None
	source: str | None = None

	@classmethod
	def from_dict(cls, data: dict) -> "Structure":
		"""Build a Structure from the dict form the search uses.

		Parameters
		----------
		data : dict
			Keys as used throughout the search, such as "sub module cell".
			Unknown keys are ignored; absent keys take their default.

		Returns
		-------
		Structure
			The same structure in typed form.
		"""
		attributes = {}
		for key, value in data.items():
			field_name = _KEY_TO_FIELD.get(key)
			if field_name is not None:
				attributes[field_name] = value
		return cls(**attributes)

	def to_dict(self) -> dict:
		"""Return the dict form the search uses.

		Returns
		-------
		dict
			Only the keys that are set, so a structure that has not been
			relaxed does not gain a "time" entry it never had.
		"""
		data = {}
		for structure_field in fields(self):
			value = getattr(self, structure_field.name)
			if structure_field.name in ("time", "source") and value is None:
				continue
			data[_FIELD_TO_KEY[structure_field.name]] = value
		return data

	def copy(self) -> "Structure":
		"""Return a copy with its own atoms and module list.

		Returns
		-------
		Structure
			A new Structure. `atoms` and `modules` are copied so that moving one
			copy's atoms does not move the other's; the remaining fields are
			plain values.
		"""
		duplicate = Structure(**{f.name: getattr(self, f.name) for f in fields(self)})
		if self.atoms is not None:
			duplicate.atoms = self.atoms.copy()
		duplicate.modules = [module.copy() for module in self.modules]
		duplicate.sub_module_cell = list(self.sub_module_cell)
		duplicate.shape_in_submods = list(self.shape_in_submods)
		return duplicate

	def __getitem__(self, key: str):
		try:
			return getattr(self, _KEY_TO_FIELD[key])
		except KeyError:
			raise KeyError(key) from None

	def __setitem__(self, key: str, value) -> None:
		try:
			field_name = _KEY_TO_FIELD[key]
		except KeyError:
			raise KeyError(key) from None
		setattr(self, field_name, value)

	def __contains__(self, key: str) -> bool:
		return key in _KEY_TO_FIELD

	def keys(self):
		"""Return the string keys of the dict form."""
		return self.to_dict().keys()
