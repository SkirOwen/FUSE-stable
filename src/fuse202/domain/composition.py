"""The empirical formula a search is run for.

A composition is the set of elements and their counts in one formula unit, for
example Ca3Ti2O7. It is the input a user supplies and the constraint every
generated structure is checked against.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from math import gcd

from ase.data import atomic_numbers, chemical_symbols

_FORMULA_TERM = re.compile(r"([A-Z][a-z]?)(\d*)")


@dataclass(frozen=True)
class Composition:
	"""Elements and their counts in one formula unit.

	Parameters
	----------
	counts : dict of {str: int}
		Element symbol to the number of atoms of it in one formula unit.

	Raises
	------
	ValueError
		If the composition is empty, names something that is not a chemical
		element, or gives a count that is not a positive integer.

	Notes
	-----
	Counts are stored as plain Python ints and `atomic_numbers` returns plain
	ints, because ase returns numpy integers from `get_atomic_numbers()`. Those
	propagate into the structure representation, where a numpy int reaching
	`(120).__ne__(...)` returns NotImplemented rather than a bool, which
	Python 3.12 and later refuse to evaluate. Converting here means the rest of
	the code never sees a numpy scalar.
	"""

	counts: dict[str, int]

	def __post_init__(self):
		if not self.counts:
			raise ValueError("composition is empty")
		normalised = {}
		for symbol, count in self.counts.items():
			if symbol not in atomic_numbers:
				raise ValueError(f"{symbol!r} is not a chemical element symbol")
			as_int = int(count)
			if as_int != count or as_int < 1:
				raise ValueError(
					f"count for {symbol} must be a positive whole number, got {count!r}"
				)
			normalised[symbol] = as_int
		object.__setattr__(self, "counts", normalised)

	@classmethod
	def from_formula(cls, formula: str) -> "Composition":
		"""Build a composition from a formula string such as "Ca3Ti2O7".

		Parameters
		----------
		formula : str
			Element symbols each optionally followed by a count. A symbol with
			no count is taken as one atom.

		Returns
		-------
		Composition
			The parsed composition.

		Raises
		------
		ValueError
			If the string contains anything that is not a recognised element
			symbol with an optional count.
		"""
		counts: dict[str, int] = {}
		position = 0
		for match in _FORMULA_TERM.finditer(formula):
			if match.start() != position:
				unparsed = formula[position:match.start()]
				raise ValueError(f"could not parse {unparsed!r} in formula {formula!r}")
			symbol, digits = match.groups()
			counts[symbol] = counts.get(symbol, 0) + (int(digits) if digits else 1)
			position = match.end()
		if position != len(formula):
			raise ValueError(f"could not parse {formula[position:]!r} in formula {formula!r}")
		return cls(counts)

	@property
	def symbols(self) -> list[str]:
		"""Element symbols, in the order they were given."""
		return list(self.counts)

	@property
	def atoms_per_formula_unit(self) -> int:
		"""Total number of atoms in one formula unit."""
		return sum(self.counts.values())

	@property
	def atomic_numbers(self) -> list[int]:
		"""Atomic number of every atom in one formula unit.

		Returns
		-------
		list of int
			One entry per atom, so Ca3Ti2O7 gives twelve entries. Plain Python
			ints, never numpy integers.
		"""
		numbers = []
		for symbol, count in self.counts.items():
			numbers.extend([int(atomic_numbers[symbol])] * count)
		return numbers

	def normalised(self) -> "Composition":
		"""Return the composition reduced to its smallest whole number ratio.

		Returns
		-------
		Composition
			Ca6Ti4O14 becomes Ca3Ti2O7. A composition already in its smallest
			ratio is returned unchanged.
		"""
		divisor = 0
		for count in self.counts.values():
			divisor = gcd(divisor, count)
		reduced = {symbol: count // divisor for symbol, count in self.counts.items()}
		return Composition(reduced)

	def atom_counts_in(self, symbols) -> dict[str, int]:
		"""Count how many atoms of each of this composition's elements appear.

		Parameters
		----------
		symbols : iterable of str
			Chemical symbols, such as `atoms.get_chemical_symbols()`.

		Returns
		-------
		dict of {str: int}
			One entry per element of this composition, including zeros.
		"""
		present = list(symbols)
		return {symbol: present.count(symbol) for symbol in self.counts}

	def formula_units_in(self, symbols) -> float | None:
		"""How many whole formula units the given atoms make up.

		Parameters
		----------
		symbols : iterable of str
			Chemical symbols, such as `atoms.get_chemical_symbols()`.

		Returns
		-------
		float or None
			The number of formula units, or None when the atoms are not a whole
			number of them in the right proportions.
		"""
		counted = self.atom_counts_in(symbols)
		ratios = set()
		for symbol, per_unit in self.counts.items():
			units = counted[symbol] / per_unit
			if units != int(units):
				return None
			ratios.add(int(units))
		if len(ratios) != 1:
			return None
		return float(ratios.pop())

	def __str__(self) -> str:
		parts = []
		for symbol, count in self.counts.items():
			parts.append(symbol if count == 1 else f"{symbol}{count}")
		return "".join(parts)
