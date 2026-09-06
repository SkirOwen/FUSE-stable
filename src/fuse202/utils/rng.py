"""Seedable random number source for FUSE.

FUSE is a stochastic search: the structures it proposes and the basin-hopping
moves it makes are all random draws. Two things follow from that, and this
module exists to provide both.

**Reproducibility.** A run should be repeatable from a single seed recorded in
its output, so a result can be re-derived or a failure re-triggered. That seed
belongs at the entry point (`run_fuse(seed=...)`), threaded down, not baked
into a leaf function.

**No global state.** The generator is an object passed explicitly rather than
the process-wide `random` module. A hardcoded `random.seed(25)` once sat inside
create_random_string(), left over from debugging: because the stdlib seed is
global it pinned the RNG for the *entire program* on every call, so structure
generation returned the same structure every time and basin-hopping move
selection was silently frozen too. Passing a generator makes that class of bug
structurally impossible.

Underneath is numpy's Generator (PCG64), which has better statistical
properties and a more stable cross-platform stream than the stdlib Mersenne
Twister. The methods below deliberately present the stdlib `random` API and
return plain Python objects, because numpy's own equivalents do not:
`Generator.choice` returns numpy scalars, and given a list of lists it treats
the input as a 2-D array and returns a row as an ndarray. Numpy integers
leaking into the structure string is precisely what made `(120).__ne__(...)`
return NotImplemented and silently break structure generation on Python 3.12+,
so these wrappers hand back exactly the objects they were given.
"""
from __future__ import annotations

import numpy as np


class Rng:
	"""Random draws for FUSE, backed by numpy's PCG64 generator."""

	def __init__(self, seed: int | None = None):
		# Draw an explicit seed when none is given, so even an unseeded run can
		# report the seed it used and be reproduced afterwards.
		if seed is None:
			seed = int(np.random.SeedSequence().entropy % (2 ** 32))
		self._seed = int(seed)
		self._generator = np.random.default_rng(self._seed)

	@property
	def seed(self) -> int:
		"""The seed this generator was created with; log it to repeat a run."""
		return self._seed

	def choice(self, seq):
		"""One element of `seq`, returned as-is (never a numpy scalar or view)."""
		if len(seq) == 0:
			raise IndexError("Cannot choose from an empty sequence")
		return seq[int(self._generator.integers(0, len(seq)))]

	def choices(self, population, weights=None, k: int = 1) -> list:
		"""`k` elements drawn with replacement, optionally weighted."""
		if len(population) == 0:
			raise IndexError("Cannot choose from an empty sequence")
		if weights is None:
			indices = self._generator.integers(0, len(population), size=k)
		else:
			total = float(sum(weights))
			if total <= 0:
				raise ValueError("weights must sum to a positive value")
			probabilities = [float(w) / total for w in weights]
			indices = self._generator.choice(len(population), size=k, p=probabilities)
		return [population[int(i)] for i in indices]

	def shuffle(self, seq: list) -> None:
		"""Shuffle a Python list in place, leaving its elements untouched.

		Fisher-Yates rather than Generator.shuffle, which would coerce the list
		to an array and hand back numpy-typed elements.
		"""
		for i in range(len(seq) - 1, 0, -1):
			j = int(self._generator.integers(0, i + 1))
			seq[i], seq[j] = seq[j], seq[i]

	def random(self) -> float:
		"""A float in [0.0, 1.0)."""
		return float(self._generator.random())

	def uniform(self, a: float, b: float) -> float:
		"""A float in [a, b]."""
		return float(self._generator.uniform(a, b))

	def __repr__(self) -> str:
		return f"Rng(seed={self._seed})"


# Fallback for callers that have not been given a generator. Contained to this
# module rather than the stdlib's process-wide state, and any code path that
# cares about reproducibility should be passed an explicit Rng instead.
_default_rng = Rng()


def default_rng() -> Rng:
	"""The shared fallback generator, used when no Rng is passed explicitly."""
	return _default_rng
