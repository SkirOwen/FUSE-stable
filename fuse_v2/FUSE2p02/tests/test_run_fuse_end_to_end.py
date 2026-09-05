"""End-to-end characterization of run_fuse() itself, using FakeCalculator so no
real GULP/VASP/QE/CHGNet install is needed.

run_fuse() is a single ~1750-line function (fuse202/run_fuse.py) that changes
the working directory and writes a large number of pickle/cif/output files as
it runs, and calls sys.exit() unconditionally at the end of a completed run
(see run_fuse.py around line 1844) - so every test here runs inside an
isolated tmp_path (via monkeypatch.chdir) and expects SystemExit rather than
a normal return. That sys.exit() also means run_fuse() can only sensibly be
called once per process - it isn't designed to be invoked repeatedly or
composed with other code today. This test file exists to pin *that* observed
behavior, not to declare it correct - Phase 4 (breaking up run_fuse) is where
that gets revisited.

get_new_structure() (the random initial-population generator) is ALSO faked
here, not just the energy calculator - see test_generate_random_structure.py
for why: it currently has a real, documented, unbounded-retry bug that made
run_fuse() hang indefinitely for realistic compositions in manual testing
before this fake was added. Faking it here means this test exercises
run_fuse()'s own orchestration (calculator dispatch, basin-hopping control
flow, output writing) without also depending on that separately-tracked bug
being fixed first.
"""
import pytest
from ase import Atoms

from fuse202.run_fuse import run_fuse

from fakes import FakeCalculator


def _canned_srtio3_structure():
	"""A real 5-atom cubic perovskite SrTiO3 cell (a=3.905 A), standing in for
	whatever get_new_structure() would otherwise have produced. Shaped to
	match the dict extract_module() returns (see fuse202/structure/extract_modules.py)."""
	atoms = Atoms(
		symbols=['Sr', 'Ti', 'O', 'O', 'O'],
		scaled_positions=[
			(0.0, 0.0, 0.0),
			(0.5, 0.5, 0.5),
			(0.5, 0.5, 0.0),
			(0.5, 0.0, 0.5),
			(0.0, 0.5, 0.5),
		],
		cell=[3.905, 3.905, 3.905],
		pbc=True,
	)
	return {
		'modules': [atoms.copy()],
		'sub module cell': [3.905, 3.905, 3.905, 90, 90, 90],
		'shape in submods': [1, 1, 1, 90, 90, 90],
		'nmods': 1,
		'ap': 3.905,
		'atoms': atoms,
	}


def test_run_fuse_completes_one_gulp_iteration_with_fake_calculator(tmp_path, monkeypatch):
	monkeypatch.chdir(tmp_path)
	fake_gulp = FakeCalculator(fixed_energy=-10.0, converged=True)
	monkeypatch.setattr("fuse202.run_fuse.run_gulp", fake_gulp.gulp)
	monkeypatch.setattr(
		"fuse202.run_fuse.get_new_structure",
		lambda **kwargs: _canned_srtio3_structure(),
	)

	with pytest.raises(SystemExit):
		run_fuse(
			composition={'Sr': 1, 'Ti': 1, 'O': 3},
			max_atoms=5,
			imax_atoms=5,
			initial_gen=1,
			iterations=1,
			ctype='gulp',
			kwds='opti conp',
			gulp_opts=['\ndump temp.res\n'],
			lib='dummy.lib',
			output_graph_at_end=False,
		)

	assert fake_gulp.call_log, "run_fuse() never invoked the (faked) GULP calculator"
