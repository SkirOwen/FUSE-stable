# FUSE v2

This is the current stable version of FUSE v.2, and is provided "as is" under the GNU public licence:


Version 2.02 and later are detailed in the following paper:

Faraday Discussions, 2024, DOI: 10.1039/D4FD00094C

Available as an open access article here:

https://pubs.rsc.org/en/content/articlelanding/2024/fd/d4fd00094c#!divAbstract

We have tested FUSE v.2 primarily on Linux desktops and clusters, and it also runs on macOS.

## Quick start

FUSE uses [uv](https://docs.astral.sh/uv/) to manage its environment. From the
repository root:

```
uv sync --extra ml
```

Then check the installation actually works by running the quickstart, which needs
no external chemistry code and finishes in a few seconds:

```
mkdir -p /tmp/fuse-quickstart && cd /tmp/fuse-quickstart
uv run --project /path/to/FUSE-stable python /path/to/FUSE-stable/examples/quickstart_chgnet.py
```

It runs a complete (very small) structure search for SrTiO3 using the CHGNet
machine-learned potential, and writes .cif files plus an `output.txt`. If those
appear, your installation is working.

To see what is installed and which calculators are available without running
anything:

```
uv run fuse202 check
```

It reports the core dependencies, whether the ML potential is installed, whether
GULP or VASP are configured, and which `ctype` values you can therefore use. It
exits non-zero only if something FUSE cannot work without is missing - not
having GULP or VASP is a normal state, not an error.

## Installation

```
uv sync              # core: everything needed for GULP, VASP and Quantum Espresso
uv sync --extra ml   # also installs the CHGNet ML potential (torch, pymatgen)
```

The core dependencies are the atomic simulation environment (ase 3.23 or later),
numpy, pandas, spglib, func-timeout and matplotlib. The `ml` extra adds chgnet,
torch and pymatgen, and is only needed for the CHGNet calculator or the ML
structure generator.

FUSE also uses external chemistry codes to perform the energy calculations, which
you need access to for real structure prediction runs. Currently GULP and VASP are
supported, with Quantum Espresso implemented but not extensively tested. In
principle any calculator supported by ase should be usable
(https://wiki.fysik.dtu.dk/ase/). If there is an unsupported calculator you would
like to use, please get in touch.

The exception is CHGNet, the ML potential published by B. Deng et al.
(https://www.nature.com/articles/s42256-023-00716-3), which is a Python package
rather than an external program - `uv sync --extra ml` is all it needs. That is
what makes the quickstart above runnable with nothing else installed.

## Repository layout

```
src/fuse202/         the package
  calculators/         GULP, VASP, Quantum Espresso, CHGNet, and the dispatch between them
  structure/           structure generation, assembly, basin-hopping moves, validation
  bonds/               bond table and bond checking
  utils/               paths, file operations, random number generation
tests/               pytest suite (runs without any external chemistry code)
examples/
  quickstart_chgnet.py   the smoke test described above
  input_files/           the six paper-reproduction setups
  probe_structures/      reference results for those experiments
resources/           GULP potential libraries, the SPP library, the gnboss archive
old_versions/        FUSE v1.02, v1.04 and v1.06, kept for reference
```

Each released version also has a `release/vX.Y` branch and a GitHub release whose
source archive contains only that version.

## Running a calculation

A FUSE run is an ordinary Python script that calls `run_fuse()` with your settings.
Copy one of the files in `examples/input_files/` (they are all configured for
Ca3Ti2O7), edit the parameters, and run it from the directory where you want the
output:

```
uv run --project /path/to/FUSE-stable python my_input.py
```

See `examples/input_files/README.txt` for what each of the six example setups
corresponds to in the paper.

## Reproducibility

`run_fuse(seed=...)` fixes the random seed. Leave it unset for a different search
each time; set it to an integer to repeat one. Either way the seed used is printed
and written to `output.txt`, so a run can be reproduced after the fact:

```
Random seed: 2468177647   (pass seed=2468177647 to run_fuse() to repeat this run exactly)
```

Note the current limitation: the seed reproduces the initial population, but not
the basin-hopping search that follows it.

## Running the tests

From the repository root:

```
uv run pytest tests/
```

The suite runs without GULP, VASP, Quantum Espresso or CHGNet installed - the
calculators are faked at a single seam (`fuse202.calculators.dispatch`). See
`tests/fakes.py` if you want to point a test at a real calculator instead.

## Installing the gnboss ML model

The ML structure generater is replacated from the model presented in this paper by G. Cheng et al.:

https://www.nature.com/articles/s41467-022-29241-4

Note: gnboss **must** be installed into its own separate environment - it cannot be added to FUSE's own dependencies or installed alongside them. It pins `tensorflow==2.3.2`, which only ships wheels for CPython 3.6-3.8, whereas FUSE v2 requires Python 3.11 or newer; its `numpy`, `pandas` and `pymatgen` pins are likewise years behind what FUSE resolves. This is why FUSE calls gnboss as a subprocess through the `gn_boss_command` parameter (pointing at that separate environment's python) rather than importing it. Its pinned dependency list is readable at `resources/gnboss-requirements.txt` without needing to unpack the archive.

to install, first unpack the "resources/gnboss_template.zip" archive and go into the first directory.

create a fresh anaconda environment:
```
conda create --name gnboss python=3.6.13
```
activate the environment:
```
conda activate gnboss
```
then install the required packages:
```
python -m pip install -r requirements.txt
```
then deactivate the model by typing:
```
conda deactivate.
```
## Statistical proxy potentials (SPPs)

To use the statistical proxy potential library provided, unzip the "resources/SPP_light.zip" and take note of path to the unpacked directory.

The statistical proxy potentials are implemented for use with GULP, and are published in this paper by D. Antypov et al.:

https://chemrxiv.org/engage/chemrxiv/article-details/65292c4a8bab5d20554598dd

## Environment variables

To use all of the features of FUSE v.2 you will need to set some of the following.
Only three are read by FUSE itself - `ASE_GULP_COMMAND` (via ase), `GNBOSS_TEMP`
and `SPP_PATH`, the latter two by the example input scripts. The rest are read by
ase or used by you when constructing paths and commands.

for GULP:
"ASE_GULP_COMMAND" : the path where your gulp executable can be found, followed by "gulp < gulp.gin > gulp.got"
"GULP_LIB": This can be set to just empty quotation marks.

for VASP:
"VASP_PP_PATH": The path to the directory containing vasp pseudo-potentials
"VASP_SCRIPT": The path to the python script for running vasp. For the example input files provieded here, this just needs to be "vasp.py"

for using the gnboss ML model (it runs in its own conda environment - see above):
"CONDA": the path to your conda executable, normally in the "bin" directory of your anaconda install. Used when building the `gn_boss_command` you pass to run_fuse().
"GNBOSS": the path to the conda environment you created for gnboss, typically under the "envs" folder of your anaconda install. Also used when building `gn_boss_command`.
"GNBOSS_TEMP": the path to the unpacked gnboss template. Read directly by the mlgen example input scripts.

for using statistical proxy potentials (this also requires GULP):
"SPP_PATH": The path to the statistical proxy potential directory unpacked above, so that FUSE can asseble any required potential sets.
