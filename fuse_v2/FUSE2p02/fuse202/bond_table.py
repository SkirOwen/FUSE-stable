# bonding environments for neutral systems in FUSE, atomic radii used from VESTA version 3.4.8
# format: 1st keys: element symbold, 2nd keys: oxidation state; values: 0 = min bond, 1 = max bond, 2= atomic radius
# all CN numbers start from 1, upto either max in shannon table or  at least 12 if metallic radii
import numpy as np

BOND_TABLE_FILE = "bondtable.npz"

BOND_DATA = {
	'H': {0: [1, 2, 0.46]},
	'Li': {0: [4, 12, 1.57]},
	'Be': {0: [3, 12, 1.12]},
	'B': {0: [3, 6, 0.81]},
	'C': {0: [3, 6, 0.77]},
	'N': {0: [3, 6, 0.74]},
	'O': {0: [2, 8, 0.74]},
	'F': {0: [2, 6, 0.72]},
	'Na': {0: [4, 12, 1.91]},
	'Mg': {0: [4, 12, 1.60]},
	'Al': {0: [4, 12, 1.43]},
	'Si': {0: [4, 6, 1.18]},
	'P': {0: [4, 6, 1.10]},
	'S': {0: [4, 6, 1.04]},
	'Cl': {0: [3, 6, 0.99]},
	'K': {0: [4, 12, 2.35]},
	'Ca': {0: [5, 12, 1.97]},  # database starts at 6, but 5 co-ord known in (Y,Ca,Sr)3Ga2O6
	'Sc': {0: [6, 12, 1.64]},
	'Ti': {0: [4, 12, 1.47]},
	'V': {0: [4, 12, 1.35]},
	'Cr': {0: [4, 12, 1.29]},
	'Mn': {0: [4, 12, 1.37]},
	'Fe': {0: [4, 12, 1.26]},
	'Co': {0: [4, 12, 1.25]},
	'Ni': {0: [4, 12, 1.25]},
	'Cu': {0: [2, 12, 1.28]},
	'Zn': {0: [4, 12, 1.37]},
	'Ga': {0: [4, 12, 1.53]},
	'Ge': {0: [4, 6, 1.22]},
	'As': {0: [4, 6, 1.21]},
	'Se': {0: [4, 6, 1.04]},
	'Br': {0: [3, 6, 1.14], },
	'Rb': {0: [6, 14, 2.50]},
	'Sr': {0: [6, 12, 2.15]},
	'Y': {0: [6, 12, 1.82]},
	'Zr': {0: [4, 12, 1.60]},
	'Nb': {0: [6, 12, 1.47]},
	'Mo': {0: [4, 12, 1.40]},
	'Tc': {0: [4, 12, 1.35]},
	'Ru': {0: [4, 12, 1.34]},
	'Rh': {0: [6, 12, 1.34]},
	'Pd': {0: [2, 12, 1.37]},
	'Ag': {0: [2, 12, 1.44]},
	'Cd': {0: [4, 12, 1.52]},
	'In': {0: [4, 12, 1.67]},
	'Sn': {0: [4, 12, 1.58]},
	'Sb': {0: [4, 12, 1.41]},
	'Te': {0: [3, 12, 1.37]},
	'I': {0: [3, 6, 1.33]},
	'Xe': {0: [4, 6, 2.18]},
	'Cs': {0: [6, 12, 2.72]},
	'Ba': {0: [6, 12, 2.24]},
	'La': {0: [6, 12, 1.88]},
	'Ce': {0: [6, 12, 1.82]},
	'Pr': {0: [6, 12, 1.82]},
	'Nd': {0: [6, 12, 1.82]},
	'Pm': {0: [6, 12, 1.81]},
	'Sm': {0: [6, 12, 1.81]},
	'Eu': {0: [6, 12, 2.06]},
	'Gd': {0: [6, 12, 1.79]},
	'Tb': {0: [6, 12, 1.77]},
	'Dy': {0: [6, 12, 1.77]},
	'Ho': {0: [6, 12, 1.76]},
	'Er': {0: [6, 12, 1.75]},
	'Tm': {0: [6, 12, 1.00]},
	'Yb': {0: [6, 12, 1.94]},
	'Lu': {0: [6, 12, 1.72]},
	'Hf': {0: [4, 12, 1.59]},
	'Ta': {0: [6, 12, 1.47]},
	'W': {0: [4, 12, 1.41]},
	'Re': {0: [4, 12, 1.37]},
	'Os': {0: [4, 12, 1.35]},
	'Ir': {0: [6, 12, 1.36]},
	'Pt': {0: [4, 12, 1.39]},
	'Au': {0: [4, 12, 1.44]},
	'Hg': {0: [2, 12, 1.55]},
	'Tl': {0: [6, 12, 1.71]},
	'Pb': {0: [4, 12, 1.75]},
	'Bi': {0: [5, 12, 1.82]},
	'Po': {0: [6, 12, 1.77]},
	'At': {0: [6, 6, 0.62]},
	'Fr': {0: [6, 6, 1.00]},
	'Ra': {0: [8, 12, 2.35]},
	'Ac': {0: [6, 12, 2.03]},
	'Th': {0: [6, 12, 1.80]},
	'Pa': {0: [6, 12, 1.63]},
	'U': {0: [2, 12, 1.56]},
	'Np': {0: [6, 12, 1.56]},
	'Pu': {0: [6, 12, 1.64]},
	'Am': {0: [6, 12, 1.73]},
	'Cm': {0: [6, 12, 1.00]},
	'Bk': {0: [6, 12, 1.00]},
	'Cf': {0: [6, 12, 1.00]},
	'No': {0: [6, 12, 1.00]}
}

import os
def initialize_bond_table(force_rebuild: bool = False) -> None:
	"""Creates the bond table file if it does not exist or if force_rebuild is True."""
	if not os.path.exists(BOND_TABLE_FILE) or force_rebuild:
		print(f"Creating {BOND_TABLE_FILE}...")
		np.savez(BOND_TABLE_FILE, bond_table=BOND_DATA)
	else:
		print(f"{BOND_TABLE_FILE} already exists. Use force_rebuild=True to recreate it.")


def load_bond_table() -> dict:
	"""Loads the bond table from the .npz file, initializing it if necessary."""
	if not os.path.exists(BOND_TABLE_FILE):
		print(f"{BOND_TABLE_FILE} not found. Creating it now...")
		initialize_bond_table()

	data = np.load(BOND_TABLE_FILE, allow_pickle=True)
	return data["bond_table"].item()


# Lazy loading: Only load bond table when explicitly called
if __name__ == "__main__":
	initialize_bond_table()  # Runs only if executed directly
