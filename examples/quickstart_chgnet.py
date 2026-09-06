"""FUSE quickstart: a complete crystal structure prediction run in ~seconds.

This is the example to run first, to confirm your installation works.

Unlike the paper-reproduction inputs in examples/input_files/, this one needs
**no external chemistry code**. Those all end in GULP or VASP, which you must
install and configure separately; this uses the CHGNet machine-learned
potential, which is an ordinary Python package installed by:

	uv sync --extra ml

Run it from a directory you don't mind filling with output files:

	mkdir -p /tmp/fuse-quickstart && cd /tmp/fuse-quickstart
	uv run --project /path/to/FUSE-stable python /path/to/examples/quickstart_chgnet.py

It searches for structures of SrTiO3 (strontium titanate) and should finish in
well under a minute. The settings below are deliberately tiny - enough to
exercise every stage of the code, far too small for a real prediction. See
examples/input_files/ for realistic setups.

What you should see, in order:
  - "Generating Initial Population", then each structure relaxed in turn
  - "Generation completed" with a lowest energy around -7 eV/atom
  - the basin-hopping search making moves, each named
	(e.g. "move : 13 : Random new structure upto current size")
  - structures/I-*.cif   the initial population
	structures/S-*.cif   structures the search proposed
	current_structure.cif, global_minimum.cif, output.txt

If it finishes and writes those files, your installation is working. Both
halves of FUSE - random structure generation and the basin-hopping search -
have then run end to end.
"""
from fuse202.run_fuse import run_fuse

run_fuse(
	# ---------------------------------------------------------------- system
	composition={'Sr': 1, 'Ti': 1, 'O': 3},  # the empirical formula to search
	max_atoms=10,      # cap on atoms per structure (tiny: real runs use 50+)
	imax_atoms=10,     # cap for the initial population specifically

	# ------------------------------------------------------------- run size
	initial_gen=2,     # structures in the starting population (real runs: 20-50)
	iterations=5,      # structures to relax in this run (real runs: hundreds)
	rmax=1000,         # stop after this many structures with no new best

	# ------------------------------------------------------- reproducibility
	# Fixing the seed makes the initial population identical between runs, which
	# is what makes this usable as a smoke test. Remove it for a real search.
	# Note: the seed currently pins the initial population but not the
	# basin-hopping search that follows.
	seed=1,

	# ----------------------------------------------------------- calculator
	# CHGNet alone - no GULP, no VASP, nothing to install separately.
	ctype='chgnet',
	n_opts=1,                                    # one relaxation stage
	opt_class=['FIRE'],                          # the ASE optimiser to use
	relaxer_opts={
		'fmax': [0.2],      # force convergence, eV/A - loose, for speed
		'steps': [20],      # optimiser steps - very few, for speed
		'verbose': [False],
	},
	opt_device='cpu',   # 'cuda' if you have an NVIDIA GPU
	mode='relax',       # 'single' for single-point energies instead

	output_graph_at_end=False,  # skip the matplotlib summary plot
)

print("\nQuickstart finished. If you see structures/*.cif and output.txt here, "
	"your FUSE installation is working.")
