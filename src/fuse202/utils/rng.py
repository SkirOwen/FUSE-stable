"""Seedable random number source for FUSE.

FUSE is a stochastic search: the structures it proposes and the basin-hopping
moves it makes are all random draws. Two things follow from that.

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

This is an instance of the standard library's own `random.Random`, which
already provides everything needed - local state, seeding, and `choice`,
`shuffle`, `choices`, `random`. FUSE only ever draws scalars from Python lists
(105 `choice` calls, 13 `shuffle`, a handful of the rest); there is no array
sampling or distribution work, so numpy's Generator would add nothing here
except type friction - its `choice` returns numpy scalars and turns a list of
lists into an ndarray row, and numpy integers leaking into the structure string
is exactly what silently broke atom counting on Python 3.12+.
"""
from __future__ import annotations

import random


class Rng(random.Random):
	"""A seeded `random.Random` that remembers its seed so a run can report it."""

	def __init__(self, seed: int | None = None):
		# Draw an explicit seed when none is given, so even an unseeded run can
		# report the seed it used and be reproduced afterwards.
		if seed is None:
			seed = random.SystemRandom().randrange(2 ** 32)
		super().__init__(seed)
		# Not `self.seed`: that name is random.Random's own re-seeding method.
		self.seed_value = int(seed)

	def __repr__(self) -> str:
		return f"Rng(seed={self.seed_value})"


# Fallback for callers that have not been given a generator. Contained to this
# module rather than the stdlib's process-wide state, and any code path that
# cares about reproducibility should be passed an explicit Rng instead.
_default_rng = Rng()


def default_rng() -> Rng:
	"""The shared fallback generator, used when no Rng is passed explicitly."""
	return _default_rng
