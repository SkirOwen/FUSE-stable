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
	monkeypatch.setattr("fuse202.calculators.dispatch.run_gulp", fake_gulp.gulp)
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


def _canned_trial_structure():
	"""What make_basin_move() hands back: the same structure shape as
	get_new_structure(), plus the three keys make_basin_move sets itself before
	returning (see structure/make_basin_move.py, near the end). Omitting those
	is not a harmless simplification - the search loop reads ['optimised?']
	directly and dies with KeyError without it."""
	structure = _canned_srtio3_structure()
	structure['optimised?'] = False
	structure['energy'] = 0.0
	structure['converged'] = False
	return structure


def test_run_fuse_drives_the_basin_hopping_search_loop(tmp_path, monkeypatch):
	"""Covers run_fuse()'s search region - trial generation, the *second*
	calculator dispatch, and the accept/reject step. The original end-to-end
	test above never reaches any of this: with iterations=1 it completes the
	initial population and exits, calling make_basin_move zero times.

	Energies cycle downhill-then-uphill so the Metropolis step has to actually
	decide rather than trivially accepting everything.
	"""
	monkeypatch.chdir(tmp_path)
	fake_gulp = FakeCalculator(per_atom_energies=(-10.0, -10.5, -9.5), converged=True)
	monkeypatch.setattr("fuse202.calculators.dispatch.run_gulp", fake_gulp.gulp)
	monkeypatch.setattr(
		"fuse202.run_fuse.get_new_structure",
		lambda **kwargs: _canned_srtio3_structure(),
	)

	basin_moves = []

	def fake_make_basin_move(*args, **kwargs):
		basin_moves.append(1)
		return _canned_trial_structure(), 1, {}

	monkeypatch.setattr("fuse202.run_fuse.make_basin_move", fake_make_basin_move)

	with pytest.raises(SystemExit):
		run_fuse(
			composition={'Sr': 1, 'Ti': 1, 'O': 3},
			max_atoms=5,
			imax_atoms=5,
			initial_gen=1,
			iterations=3,
			search=1,          # basin hopping
			search_gen_bh=1,
			rmax=1000,
			ctype='gulp',
			kwds='opti conp',
			gulp_opts=['\ndump temp.res\n'],
			lib='dummy.lib',
			output_graph_at_end=False,
		)

	assert basin_moves, "search loop never ran - make_basin_move was not called"
	# one call for the initial population, plus one per trial structure
	assert len(fake_gulp.call_log) == 1 + len(basin_moves)
	assert all(backend == 'gulp' for backend in fake_gulp.call_log)

	# every relaxed structure is written out for the record
	written = sorted(p.name for p in (tmp_path / "structures").glob("*.cif"))
	assert written, "no per-structure cif files were written"

	# the running best is tracked on disk, and moves are logged by name
	assert (tmp_path / "current_structure.cif").exists()
	output = (tmp_path / "output.txt").read_text()
	assert "Swap two atoms" in output, "the move description was not logged"
