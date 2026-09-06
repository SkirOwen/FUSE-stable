# FUSE v3, design criteria and direction

Status: draft for discussion. Nothing here is built yet.

## Why v3

FUSE v2 works, but only after nine bugs found in one review session, seven of
which were silent. It could not generate a single crystal structure for its own
documented example, and had not been able to for some time. That is not bad luck.
The bugs cluster tightly around a handful of structural properties, and the same
properties will keep producing them.

Measured on the current package (10,232 lines):

| | count | what it cost us |
|---|---|---|
| bare `except:` | 50 | hid the headline bug for months; made loops un-interruptible |
| `from x import *` | 73 | untraceable name provenance; caused a real circular import |
| `os.chdir` calls | 39 | working directory as hidden global state |
| `print()` calls | 273 | no way to control, capture, or silence output |
| `run_fuse()` parameters | 70 | no caller can hold the whole surface in their head |
| `run_fuse()` lines | 1,798 | untestable as a unit |

And the central data object is a bare `dict` keyed by strings like
`'sub module cell'`, `'shape in submods'` and `'optimised?'`. A missing key is a
runtime `KeyError` at hour six of a cluster job.

**v3 is an engineering release.** The science is not being changed. The goal is
that the same algorithms become code where this class of bug is hard to write.

## Decisions taken

- **Evolve in place**, not rewrite. The science works; it is the scaffolding that
  fails. Every step stays runnable and checkable against the test suite.
- **Same science.** Same algorithms, same search behaviour. Bit-identical output
  is already impossible (the RNG changed in v2.2), but component-level equivalence
  is checkable and will be checked, see *Proving we did not break it*.
- **New API, legacy input still works.** v3 gets the public API it deserves.
  Existing v2 input scripts keep running, via a compatibility layer that maps them
  onto the new API. What happens behind that entry point is free to change.
- **Group now, public later.** Design so a public release needs no rework; do not
  pay the full cost of one yet.

## Design criteria

Each of these exists because something specific went wrong. They are acceptance
criteria, not aspirations.

1. **Failures are loud.** No bare `except:`. Every caught exception names its type
   and says why swallowing it is correct. *Origin: `filter((120).__ne__, ...)`
   raised on every attempt for months, silently caught, turning "cannot generate
   structures" into "runs forever with no output".*

2. **The domain has types.** Structures, compositions and search state are typed
   objects, not dicts keyed by strings with spaces in them. A wrong field name is
   a static error, not an hour-six crash. *Origin: `KeyError: 'optimised?'`.*

3. **No global mutable state.** Randomness, configuration and paths are passed
   explicitly. *Origin: a debug `random.seed(25)` in a leaf function froze the
   entire program's randomness, including move selection.*

4. **Boundaries normalise types.** Values entering the domain from ase, numpy or a
   calculator are converted to Python types at the edge. *Origin: `numpy.int64`
   broke atom counting; `numpy.float32` broke every CHGNet run. Same bug twice.*

5. **Loops terminate provably.** Every retry has a bound and a clear error when it
   is exhausted. *Origin: two unbounded loops; one had no cap at all, the other's
   cap could never fire.*

6. **One implementation of anything.** Where two copies exist today, they have
   already drifted, every single pair we found had. *Origin: the calculator
   dispatch, the restart backup, and the post-calculation checks were each
   duplicated, and each pair differed.*

7. **I/O at the edges.** The science operates on values and returns values.
   Reading and writing files, and choosing where they go, happens in one place and
   never by changing the process working directory. *Origin: 39 `os.chdir` calls.*

8. **Output is structured.** Progress, results and diagnostics go through logging
   and typed result objects, not `print()`. A caller can silence it, capture it, or
   drive a UI from it. *Origin: 273 `print()` calls; a cluster job's only record is
   scraped stdout.*

9. **Reproducible from a seed.** A recorded seed reproduces an entire run,
   including the search. *Origin: currently only the initial population reproduces.*

10. **A newcomer can run it in five minutes.** `uv sync`, one command to check the
    install, one example that works with nothing else installed. Already true as of
    v2.2 and must stay true.

## Paradigm

**Functional core, imperative shell, with protocols at the plug points.**
Not a deep class hierarchy.

- **Functional core.** The science is transformations over values: given a
  structure and parameters, return a new structure. These are pure, easy to test
  in isolation, and cannot corrupt shared state. Most of FUSE's algorithms are
  already shaped this way; they are just tangled up with I/O.

- **Imperative shell.** One orchestration layer owns the loop, the filesystem, the
  restart files and the logging. It is allowed to be stateful because it is the
  only thing that is, and it is thin enough to read.

- **Dataclasses for data and configuration.** `Structure`, `Composition`,
  `SearchSettings`, `CalculatorSettings`. The 70 parameters become a handful of
  grouped, validated objects, which also gives the compatibility layer something
  concrete to map old keyword arguments onto.

- **Protocols, not ABCs, at the plug points.** Three natural seams:
  - `Calculator`, the `(atoms, energy, converged)` contract five backends already
    share. This already works and is already the test suite's swap point.
  - `Move`, the 14 basin-hopping moves currently live in a 2,000-line if/elif
    chain in one function. Each is independently describable and testable.
  - `StructureGenerator`, random generation today, the ML model as an alternative.

  Structural typing keeps the fakes trivial and avoids forcing every implementation
  to inherit from us.

Why not heavy OOP: the failures were not caused by missing polymorphism. Deep
hierarchies would add indirection without addressing a single item on the criteria
list. Why not pure functional: restart files, long-running cluster jobs and
external processes are irreducibly stateful, pretending otherwise costs more than
it saves.

## Code conventions

Binding for v3 and for any new code added to v2.

1. **Tabs for indentation.** Actual tab characters, not four spaces. The existing
   codebase uses tabs and mixing the two is worse than either.
2. **Numpy style docstrings.** Parameters, Returns, and where useful Raises.
3. **Docstrings describe the function, not its history.** What it does, what it
   takes, what it returns. Why a line exists, what bug it guards against, and how
   the code came to look this way belong in ordinary comments next to the code or
   in the commit message. A reader calling the function does not need its
   biography.
4. **No dashes in docstrings.** No em dashes, en dashes, or hyphens used as
   punctuation. Use commas, semicolons, parentheses, or a new sentence. Hyphens
   inside real compound words are fine.
5. **No operations inside return statements.** Compute into a named variable and
   return that. `return results.sorted().where(condition)` hides both the work
   and a useful name; a debugger cannot stop on it and a traceback cannot point
   inside it. A single plain call, `return run_gulp(...)`, is fine.
6. **Follow the Zen of Python.** In this codebase three lines of it carry most of
   the weight. *Explicit is better than implicit*: pass the generator, pass the
   paths, name the exception. *Errors should never pass silently*: 50 bare
   `except:` clauses are how a total failure looked like a slow run. *Flat is
   better than nested*: the deepest bugs were found at six levels of indentation
   inside a 2,000 line function.

## Bugs to fix, and how

Carried from the v2 review. Each is currently documented and, where it was unsafe
to change behaviour, deliberately left alone.

| Bug | Fix |
|---|---|
| `dE_curr` quantised to itself, a no-op, so the Metropolis comparison is never rounded, defeating `e_prec` exactly where it matters | Quantise to `e_prec` like the adjacent `dE_glob` line. **Changes search behaviour**, needs a before/after comparison run, not a silent fix |
| `e_prec` only works for powers of ten (`Decimal.quantize` uses its argument's *exponent*, so `e_prec=0.5` rounds to 1 d.p., `0.25` does nothing) | `round(x / e_prec) * e_prec`, matches the parameter's name, tolerates numpy types |
| numpy scalars reach code assuming Python numbers (caused two separate silent failures) | Normalise at the boundary: cast in `get_fu()` and wherever a calculator's energy enters. Criterion 4 |
| Seed does not reproduce the basin-hopping search | Find the remaining unseeded draw. The fallback generator was ruled out; the source was not identified |
| 50 bare `except:` | Narrow each to the exception actually expected. Several are load-bearing, the calculator dispatch deliberately rejects a structure rather than aborting a run, so each needs reading, not a blanket edit |
| `test_bonds.py` is production code named like a test file | Rename to `bond_checks.py`. Pytest already has to be told not to collect it |
| `_test_bonds()` dead alongside `test_bonds()`; empty `class Fuse`; `assemble_structure_2.py` vs `assemble_structure_generator.py` | Delete the dead, resolve the duplicate by checking which is reachable |
| Structure-generation success depends on rejection sampling that can take thousands of attempts | Not a bug, but worth measuring: with the fixes in place, how often does generation succeed per attempt, and can the constraints be applied earlier? |

## Telling users when results differ

v3 fixes bugs that change search behaviour, so a rerun of an old calculation can
legitimately produce a different answer. Silently returning different numbers for
the same input is the worst possible outcome for a scientific code, so v3 says so.

Three mechanisms, in increasing intrusiveness:

1. **A divergence log**, `docs/v3-divergences.md`. One entry per behavioural
   change: what v2 did, what v3 does, which bug caused it, and whether published
   results could be affected. Written as each fix lands, not reconstructed later.
2. **A line in every run's output header**, next to the seed already recorded
   there, naming the version and pointing at that log. Someone comparing two
   `output.txt` files then sees immediately that the versions differ and where to
   read why.
3. **An explicit warning when a v2 input script is run through the compatibility
   layer**, since that is exactly the person most likely to be repeating old work:

   ```
   Running a v2 input script under FUSE v3.
   v3 fixes 4 bugs that change search behaviour; results will not match v2
   run for run. See docs/v3-divergences.md. Pass legacy_quiet=True to silence.
   ```

The warning is worth the noise. The alternative is a researcher quietly getting
different structures and not knowing whether the science or the software changed.

## Proving we did not break it

"Same science" is only meaningful if it is checkable. Three levels, in order of
strength:

1. **Component-level differential testing.** For each function extracted or
   rewritten, run the v2.2 implementation and the v3 one on identical inputs and
   compare outputs exactly. `release/v2.2` is a clean snapshot, so both versions
   are available; the mechanics were checked rather than assumed, and there is a
   trap worth knowing about:

   - **Leaf modules** (no intra-package imports) can be loaded directly from a
     v2.2 worktree with `importlib.util.spec_from_file_location` and compared
     in-process. Verified working: v2.2 and current `possible_solutions` produce
     identical output this way.
   - **Modules that import other `fuse202` modules cannot.** They load without
     error, but their `from fuse202.x import y` statements resolve against the
     *currently installed* package, so you silently compare v2.2's function
     running on v3's dependencies. It looks like it works. This is the same shape
     of mistake that produced a wrong conclusion during the v2 review, where a
     comparison against upstream used hand-substituted inputs and appeared to
     exonerate code that was in fact fine for different reasons.
   - **The reliable recipe** is to shadow the package: run the comparison in a
     subprocess whose `sys.path` starts with the v2.2 worktree's `src/`, using the
     project venv's interpreter so dependencies are still available. Verified:
     both `fuse202` and its submodules then resolve to the worktree.

   Gate every extraction on this.

2. **Invariants.** Properties that must hold regardless of implementation:
   composition is preserved through every move, energies are per-atom, error
   fractions lie in [0, 1], accepted structures pass the bond and distance checks.
   Hypothesis is already in the suite for this.

3. **Statistical comparison of whole searches.** Run N searches under each version
   and compare the distribution of energies found. Weak, slow, and the last line of
   defence, but it is the only end-to-end check available once the RNG stream
   differs.

A caveat worth stating plainly: v2.2's own behaviour is only "correct" in the sense
that it is what the science was validated against. Where v3 differs *because a v2
bug was fixed*, that is a deliberate divergence and must be recorded as one.

## Classes: where they earn their place

Not a matter of taste here; it follows from what each piece of the code is.

**Use a class when there is real state with invariants to protect, or a set of
interchangeable implementations.** In FUSE that is three places:

- **`Move`.** Fourteen basin hopping moves currently live in a 2,000 line
  if/elif chain. Each has the same shape: take a structure, return a modified
  structure or report that it does not apply. They are the textbook case for
  interchangeable implementations, and splitting them is the single biggest
  readability win available.
- **`Calculator`.** Five backends, one contract, already the point where the
  test suite substitutes fakes. It works today; formalise it.
- **`StructureGenerator`.** Random generation now, the ML model as an
  alternative, chosen at runtime.

For all three, use `typing.Protocol` rather than an abstract base class. Nothing
needs to inherit from FUSE to satisfy them, fakes stay trivial, and the existing
plain functions already conform without modification.

**Use a dataclass when you are carrying data around.** `Structure`,
`Composition`, `SearchSettings`, `CalculatorSettings`. This is the largest single
improvement available: the dict keyed by `'sub module cell'` and `'optimised?'`
becomes a typed object where a wrong field is caught before the job starts.

**Use a plain function for everything else**, which is most of the science. A
bond check is a function of a structure and a table. A density check is a
function of a structure. Wrapping those in classes adds a `self` that holds
nothing and hides the fact that they are pure.

**Do not use a class as a namespace.** The empty `class Fuse: pass` sitting in
`run_fuse.py` is the warning. A module is already a namespace.

So: dataclasses for data, protocols for the three plug points, functions for the
science, one stateful orchestrator for the run. Roughly two or three real classes
of behaviour, not a hierarchy.

## If this were built from scratch

Recorded as a target to steer toward, not a plan to execute. We decided to evolve
in place, and the value of that decision is that the code keeps working. But it
helps to know where the evolution is heading.

The shape follows from one observation: **FUSE is a pipeline over a Structure,
wrapped in a search loop, talking to slow and unreliable external programs.**

```
fuse/
  domain/          what a structure IS. No I/O, no randomness, no calculators.
    structure.py     Structure, Module, SubModule
    composition.py   Composition, formula parsing, atom counts
    cell.py          lattice shapes and the sub-module tilings
  science/         transformations over the domain. Pure functions.
    generate.py      random structure generation
    assemble.py      building an Atoms object from modules
    validate.py      bonds, density, contacts
    moves/           one module per basin hopping move
  calculators/     the outside world: GULP, VASP, QE, CHGNet
    protocol.py      the (atoms, energy, converged) contract
    gulp.py, vasp.py, qe.py, chgnet.py
  search/          the algorithms
    basin_hopping.py
    state.py         SearchState: population, energies, temperature, counters
  runtime/         the imperative shell
    session.py       owns paths, restart files, and the run loop
    restart.py       save and resume, one place
    reporting.py     logging and the output file
  api.py           the public entry point
  compat.py        v2's run_fuse() keyword arguments mapped onto api.py
  cli.py
```

Four properties matter more than the exact folder names:

1. **`domain/` and `science/` import nothing from `calculators/` or `runtime/`.**
   The science becomes testable without a filesystem or a calculator, which is
   most of what made v2 hard to test at all.
2. **The filesystem lives in one module.** Paths are values passed in, never the
   process working directory. This removes 39 `os.chdir` calls and the whole
   class of bug where code works only when run from the right place.
3. **Search state is one object, not forty locals.** That is what makes restart
   honest: saving a run means serialising one thing, and the reason restart is
   fragile today is that the state is scattered across a 1,798 line function.
4. **One move per file.** Fourteen files of roughly a hundred lines each,
   individually testable, instead of one file of two thousand.

The honest caveat: this is a target, and targets drawn before the work are
usually wrong somewhere. The Move protocol in particular will only prove itself
once three or four real moves have been extracted and the awkward ones,
especially those that consult the pre built structure pool, have had their say.

## Sequencing

Roughly, and each step verifiable on its own:

1. Types first, `Structure` and `Composition`, with the boundary normalisation
   that kills the numpy-scalar class of bug. Everything else gets easier once the
   central object has a shape.
2. Configuration objects, and the compatibility layer mapping v2's 70 keyword
   arguments onto them. Existing input scripts must keep passing their tests.
3. Extract the `Move` protocol and pull the 14 moves out of the 2,000-line chain,
   one at a time, differential-tested against v2.2.
4. The imperative shell: one place that owns paths, restart files and logging;
   remove `os.chdir` and replace `print()` with logging.
5. `run_fuse()`'s remaining decomposition, which by then is mostly done.
6. The behaviour-changing bug fixes, each with a before/after comparison.
7. Public API polish and docs.

## Non-goals

- Changing the search algorithm, the move set, or any scientific parameter default.
- Performance work. Nothing here is known to be too slow; correctness first.
- Modernising gnboss (pinned to TensorFlow 2.3 / Python 3.6-3.8, a separate project).
- Breaking existing input scripts.
