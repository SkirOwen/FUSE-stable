"""Test doubles for FUSE's external calculators (GULP/VASP/QE/CHGNet).

None of these are installed on a typical dev machine or in CI, so run_fuse()
and friends need a way to run without them. run_gulp/run_vasp/run_qe/
run_chgnet/run_calculators all share the same (atoms, energy, converged)
return contract - FakeCalculator matches that shape so it can be swapped in
via monkeypatch wherever fuse202 imports those functions by name (they're
brought in with `from fuse202.calculators.run_gulp import *` etc., so patch the name on
the *importing* module, e.g. `monkeypatch.setattr(fuse202.run_fuse, "run_gulp",
fake.gulp)`, not on fuse202.calculators.run_gulp itself).

To exercise a real calculator when one is actually available on a given
machine, point the same monkeypatch target at the real function instead of
a FakeCalculator method - the call signature is unchanged either way, so no
other test code needs to change.
"""
import itertools
from dataclasses import dataclass, field


@dataclass
class FakeCalculator:
	"""Deterministic stand-in for FUSE's calculator backends.

	By default every call returns the same fixed_energy/converged pair
	regardless of the structure passed in. call_log records which backend was
	invoked, so tests can assert on dispatch (e.g. "gulp was invoked exactly
	twice") without caring about the fake energy value itself.

	`per_atom_energies` drives the basin-hopping search loop, which can only
	make an accept/reject decision if successive structures differ in energy.
	Pass a sequence and each call takes the next value, cycling when exhausted.

	Energy units: the real backends return a *total* energy, which run_fuse()
	then divides by len(atoms). Both fields here are specified **per atom** and
	multiplied by the atom count on the way out, so a test asking for -10.5
	gets exactly -10.5 eV/atom recorded downstream rather than something that
	depends on how many atoms the structure happens to have.
	"""
	fixed_energy: float = -10.0
	converged: bool = True
	per_atom_energies: tuple = ()
	call_log: list = field(default_factory=list)
	_cycle: object = field(default=None, init=False, repr=False)

	def _next_per_atom_energy(self):
		if not self.per_atom_energies:
			return self.fixed_energy
		if self._cycle is None:
			self._cycle = itertools.cycle(self.per_atom_energies)
		return next(self._cycle)

	def _respond(self, atoms, label):
		self.call_log.append(label)
		energy = self._next_per_atom_energy()
		try:
			energy = energy * len(atoms)
		except TypeError:  # atoms is a placeholder (default ''), not a real Atoms
			pass
		return atoms, energy, self.converged

	def gulp(self, atoms='', shel=None, kwds='', opts='', lib='', produce_steps='',
	         gulp_command='gulp < gulp.gin > gulp.got', gulp_timeout=''):
		return self._respond(atoms, 'gulp')

	def vasp(self, atoms='', vasp_opts='', kcut='', produce_steps='', dist_cutoff=1.0, ftol=0.5):
		return self._respond(atoms, 'vasp')

	def qe(self, atoms=None, qe_opts=None, kcut=None, produce_steps=''):
		return self._respond(atoms, 'qe')

	def chgnet(self, atoms, n_opts=2, rel=None, relaxer_opts=None, opt_class=None,
	           mode='relax', opt_device='cpu', use_spglib=True):
		return self._respond(atoms, 'chgnet')
