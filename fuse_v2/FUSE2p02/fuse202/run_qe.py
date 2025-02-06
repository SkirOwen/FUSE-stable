import glob
import shutil
import math
import numpy as np
import os
import re
##########################################################################################
import shlex
import sys
from ase import *
from ase.calculators.espresso import Espresso
from ase.io import read, write

from fuse202.constant import Ry_TO_eV


def cellpar(atoms: Atoms) -> tuple:
	cell = atoms.cell
	a = np.linalg.norm(cell[0])
	b = np.linalg.norm(cell[1])
	c = np.linalg.norm(cell[2])

	alpha = np.arccos(np.dot(cell[1], cell[2]) / (b * c)) * 180. / np.pi
	beta = np.arccos(np.dot(cell[0], cell[2]) / (a * c)) * 180. / np.pi
	gamma = np.arccos(np.dot(cell[0], cell[1]) / (a * b)) * 180. / np.pi

	return a, b, c, alpha, beta, gamma


def check_convergence(output_file: str = "") -> bool:

	# Assuming file exists
	with open(output_file, "r") as f:
		lines = f.readlines()

	converged = False

	if not any("JOB DONE." in line for line in lines):
		converged = False

	if any('nstep' in line for line in lines):
		if 'A final scf calculation at the relaxed structure.' in lines:
			converged = True
	elif not any('nstep' in line for line in lines):
		if any('convergence has been achieved' in line for line in lines):
			converged = True

	return converged


def remove_temp_files(pattern: str = "pwscf*") -> None:
	""""""
	for filename in glob.glob(pattern):
		if os.path.isdir(filename):
			shutil.rmtree(filename)
		else:
			os.remove(filename)


def run_qe(atoms: Atoms, qe_opts: dict, kcut, produce_steps=''):
	new_atoms = atoms.copy()
	for i in range(len(list(qe_opts.keys()))):
		if len(kcut) >= 2:
			temp_kcut = kcut[i]
		else:
			temp_kcut = kcut
		cell = new_atoms.get_cell_lengths_and_angles()
		kp = []
		kp = (
			int(math.ceil(temp_kcut / cell[0])),
			int(math.ceil(temp_kcut / cell[1])),
			int(math.ceil(temp_kcut / cell[2]))
		)
		calc = qe_opts[str(list(qe_opts.keys())[i])]
		calc.set(kpts=kp)
		calc = qe_opts[str(list(qe_opts.keys())[i])]
		new_atoms.set_calculator(calc)
		try:
			energy = new_atoms.get_potential_energy()
		except:
			converged = False
			energy = 1.e20
		if produce_steps == True:
			label = str('atoms' + str(i + 1) + '.cif')
			write(label, new_atoms)

		# remove any output / tempory files before starting"
		remove_temp_files(pattern="pwscf*")

	converged = check_convergence("espresso.pwo")
	#    print("converged:" +str(converged)+str("\n\n"))

	temp_atoms = new_atoms.repeat([2, 2, 2])
	temp1 = temp_atoms.get_all_distances()
	temp2 = []
	for i in range(len(temp1)):
		for j in range(len(temp1[i])):
			if temp1[i][j] != 0:
				temp2.append(temp1[i][j])
	if min(temp2) <= 1.2:
		converged = False

	# convert from Ry to eV
	energy = energy / Ry_TO_eV

	return new_atoms, energy, converged

# qe_opts={
# '1':Espresso(pseudopotentials={'Sr':'Sr.pbe-spn-kjpaw_psl.1.0.0.UPF',
#                               'Nb':'Nb.pbe-spn-kjpaw_psl.1.0.0.UPF',
#                               'O':'O.pbe-n-kjpaw_psl.1.0.0.UPF'},
#             input_data={'control':{'calculation':'vc-relax',
#                                    'pseudo_dir':'./pp/',
#                                    'outdir':'./',
#                                    'etot_conv_thr':0.0001,
#                                    'forc_conv_thr':0.001,
#                                    'disk_io':'low',
#                                    'nstep':50},
#                         'system':{'ecutwfc':29.39945,
#                                   'occupations':'smearing',
#                                   'smearing':'gaussian',
#                                   'degauss':0.0073498},
#                         'electrons':{'electron_maxstep':50,
#                                      'conv_thr':0.000001}},
#             koffset=(0,0,0))
# }
#
# atoms,energy,converged=run_qe(atoms=read('SrNbO3_converged.cif'),qe_opts=qe_opts,kcut=[20])
#
# print(str('energy: ')+str(energy))
# print(str('converged: ')+str(converged))
