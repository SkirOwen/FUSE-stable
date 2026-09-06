"""Seedable random number source for FUSE.

FUSE is a stochastic search, so a run should be repeatable from a single seed
and no part of the code should be able to disturb randomness elsewhere in the
process. Generators are therefore passed explicitly rather than taken from the
process wide `random` module.
"""
from __future__ import annotations

import random


class Rng(random.Random):
	"""A seeded random.Random that remembers the seed it was created with.

	Parameters
	----------
	seed : int or None
		Seed for the generator. When None, a seed is drawn from system entropy
		so that an unseeded run can still report a seed and be repeated later.

	Attributes
	----------
	seed_value : int
		The seed in use. Named `seed_value` because `seed` is random.Random's
		own re seeding method.
	"""

	def __init__(self, seed: int | None = None):
		if seed is None:
			seed = random.SystemRandom().randrange(2 ** 32)
		super().__init__(seed)
		self.seed_value = int(seed)

	def __repr__(self) -> str:
		return f"Rng(seed={self.seed_value})"


# Contained to this module rather than the standard library's process wide
# state. Code that needs reproducibility should be passed an explicit Rng.
_default_rng = Rng()


def default_rng() -> Rng:
	"""Return the shared generator used when no Rng is passed explicitly.

	Returns
	-------
	Rng
		The module level fallback generator.
	"""
	return _default_rng
