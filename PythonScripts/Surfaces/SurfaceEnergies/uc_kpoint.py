import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.build import molecule
from ase.calculators.vasp import Vasp
mpl.style.use("seaborn-poster")

infile = sys.argv[1]
uc = read(infile)

KPTS = [2, 3, 4, 5, 6, 8, 10, 12]
TE = []

for k in KPTS:
	calc = Vasp(
		directory='uc-kpts-{0}'.format(k),
		xc='PBE',
		kpts=[k, k, k],
		encut=520,
		atoms=uc)

	TE.append(uc.get_potential_energy())
	if None in TE:
		calc.abort()

# consider the change in energy from lowest energy state
TE = np.array(TE)
TE -= TE.min()

print(KPTS)
print(TE)

plt.plot(KPTS, TE)
plt.xlabel('No. k-points (MP)')
plt.ylabel('Total Energy (eV)')
plt.savefig('uc-kpt-convergence.png')
