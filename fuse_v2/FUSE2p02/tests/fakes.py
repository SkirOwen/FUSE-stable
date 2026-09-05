"""Test doubles for FUSE's external calculators (GULP/VASP/QE/CHGNet).

None of these are installed on a typical dev machine or in CI, so run_fuse()
and friends need a way to run without them. run_gulp/run_vasp/run_qe/
run_chgnet/run_calculators all share the same (atoms, energy, converged)
return contract - FakeCalculator matches that shape so it can be swapped in
via monkeypatch wherever fuse202 imports those functions by name (they're
brought in with `from fuse202.run_gulp import *` etc., so patch the name on
the *importing* module, e.g. `monkeypatch.setattr(fuse202.run_fuse, "run_gulp",
fake.gulp)`, not on fuse202.run_gulp itself).

To exercise a real calculator when one is actually available on a given
machine, point the same monkeypatch target at the real function instead of
a FakeCalculator method - the call signature is unchanged either way, so no
other test code needs to change.
"""
from dataclasses import dataclass, field


@dataclass
class FakeCalculator:
	"""Deterministic stand-in for FUSE's calculator backends.

	Every call returns the same fixed_energy/converged pair regardless of the
	structure passed in, which is enough to drive run_fuse()'s basin-hopping
	loop without it being able to distinguish trial structures by energy.
	call_log records what was called and with what, so tests can assert on
	dispatch (e.g. "gulp was invoked exactly twice") without caring about the
	fake energy value itself.
	"""
	fixed_energy: float = -10.0
	converged: bool = True
	call_log: list = field(default_factory=list)

	def _respond(self, atoms, label):
		self.call_log.append(label)
		return atoms, self.fixed_energy, self.converged

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
