"""Characterizes generate_random_structure()'s retry loop, and the *unbounded*
retry loop one layer up in its caller, get_new_structure() - see
fuse202/generate_random_structure.py and fuse202/make_new_structure.py.

Manual testing of run_fuse() (with a faked energy calculator, so no external
binary was involved) hung indefinitely at "Generating Initial Population" for
two different realistic-looking compositions (SrTiO3 and Ca3Ti2O7, using
parameters straight out of fuse_v2/example_input_files/fgen_input.py). Reading
the source shows why, and it's two compounding bugs:

1. get_new_structure()'s `while gen == False:` loop (make_new_structure.py:43)
   calls generate_random_structure() again every time it comes back empty,
   with no cap of its own - so whenever a composition can never satisfy
   error_check_structure() within generate_random_structure()'s own attempt
   budget, get_new_structure() retries that budget forever.
2. generate_random_structure()'s own budget (max_attempts=1000) only works
   when candidates are cleanly *rejected*. Two bare `except: continue` blocks
   (around create_random_string() and assemble_structure()) skip straight
   back to the top of the loop *before* `attempts += 1` runs - so if either
   call raises on every try, attempts never advances and the cap never fires.

test_generate_random_structure_gives_up_after_max_attempts below executes bug
2's *good* path (clean rejection) and confirms it behaves correctly - that
part isn't broken. The other two tests can't safely be executed: the bare
`except:` in the real code catches literally everything, including an
exception raised from inside a mock specifically to escape the loop, and
multiprocessing termination was ruled out because this process has already
initialized CHGNet on the MPS backend at import time (see
test_run_fuse_helpers.py's module docstring) and forking after a GPU context
is initialized is its own hazard. So bugs 1 and 2's actually-broken paths are
pinned as source-structure assertions instead - honest about the fact that
they characterize the code as written, not an executed reproduction, and
worth re-deriving properly once Phase 4 breaks up this call chain.
"""
import inspect
import re

import fuse202.generate_random_structure as grs
import fuse202.make_new_structure as mns


def test_generate_random_structure_gives_up_after_max_attempts(monkeypatch):
	"""generate_random_structure()'s own `attempts == max_attempts` cap (1000,
	hardcoded) does work when candidates are cleanly rejected (no exception):
	it returns atoms=None rather than looping forever."""
	monkeypatch.setattr(grs, "create_random_string", lambda *a, **k: ("string", "instructions"))
	monkeypatch.setattr(grs, "create_random_instructions", lambda *a, **k: "instructions")
	monkeypatch.setattr(grs, "assemble_structure", lambda *a, **k: ([1, 2, 3], "instructions"))
	monkeypatch.setattr(grs, "error_check_structure", lambda *a, **k: 0)  # always reject

	atoms, string, instructions, accept = grs.generate_random_structure(
		target_atoms=5,
		cubic_solutions={}, tetragonal_solutions={}, hexagonal_solutions={},
		orthorhombic_solutions={}, monoclinic_solutions={},
		atoms_per_fu=5, ideal_density=1.0, density_cutoff=0.4,
		check_bonds=True, btol=0.25, system_type='neutral',
		fu=[38, 22, 8, 8, 8], composition={'Sr': 1, 'Ti': 1, 'O': 3},
		bondtable={}, ap=4.672, check_distances=True, dist_cutoff=1.0,
	)

	assert atoms is None
	assert accept == 0


def test_KNOWN_BUG_swallowed_exception_bypasses_the_attempt_cap():
	"""KNOWN BUG (not yet fixed, not safely executable - see module docstring):
	`attempts += 1` (generate_random_structure.py:114) sits after both
	try/except-continue blocks that wrap create_random_string() and
	assemble_structure(). If either raises on every call, `continue` returns
	to the top of the loop without ever reaching the increment, so
	max_attempts is never reached. Pinned here as a source-order assertion:
	both `except:` blocks must appear before the increment in the function
	body. If this test starts failing because the increment moved earlier (or
	the except clauses became narrower/no longer swallow everything before
	incrementing), that's the bug being fixed - update or delete this test to
	match.
	"""
	source = inspect.getsource(grs.generate_random_structure)
	attempts_increment_pos = source.index("attempts += 1")
	except_positions = [m.start() for m in re.finditer(r"\n\s*except\s*:\s*\n\s*continue", source)]

	assert len(except_positions) == 2, (
		"expected exactly 2 bare `except: continue` blocks in "
		"generate_random_structure() - count changed, re-verify this test's premise"
	)
	assert all(pos < attempts_increment_pos for pos in except_positions), (
		"a bare except-continue now runs AFTER the attempts increment - if so, "
		"the swallowed-exception bug described in this test's docstring is fixed"
	)


def test_KNOWN_BUG_get_new_structure_outer_loop_has_no_attempt_cap():
	"""KNOWN BUG (not yet fixed, not safely executable - see module docstring):
	get_new_structure()'s `while gen == False:` loop (make_new_structure.py:43)
	has no max-attempts guard at all, unlike generate_random_structure()'s own
	`attempts == max_attempts` check. This is the loop that actually hung
	during manual run_fuse() testing. Pinned as a source assertion: no
	attempt-counting construct exists between the `while gen == False:` line
	and the function's end. If this test starts failing because someone added
	a cap, that's the bug being fixed - update or delete this test to match.
	"""
	source = inspect.getsource(mns.get_new_structure)
	loop_start = source.index("while gen == False:")
	loop_body = source[loop_start:]

	assert "attempts" not in loop_body, (
		"get_new_structure() now tracks attempts in its retry loop - if so, "
		"the unbounded-retry bug described in this test's docstring is fixed"
	)
