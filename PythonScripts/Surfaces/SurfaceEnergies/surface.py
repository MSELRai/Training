import sys
from ase import Atom, Atoms
from ase.io import read, write
from ase.build import sort, surface, add_adsorbate
from ase.constraints import FixAtoms

# Unit Cells obtained from materials project
# Alpha Mo2C (mp-1552) most predominat surface is (121)
# Beta Mo2C (mp-1221498) most predominant surface is (1 0 -1 1) in hexagonal coordinates, which translates to (101) in xyz

infile_uc = sys.argv[1]
vacuum = 7.5               #vacuum added to both sides of cell
miller_idx = (1,2,1)
nrepeat_x = 4
nrepeat_y = 4
nlayers = 4

unit_cell = read(infile_uc)
catalyst = surface(unit_cell, miller_idx, nlayers, vacuum).repeat((nrepeat_x, nrepeat_y, 1))

calc = Vasp(xc='pbe',
        command="vasp_gam_vtst",
	directory='optimize',
	encut=520,
        prec='Normal',
	ivdw=11,
        ismear=0,
        sigma=0.05,
        ibrion=2,
	isif=2,
        nsw = 1000,
        ediff=1.0E-5,
        ediffg=-0.03,
        algo='V',
        lreal='Auto',
        ncore=48,
        setups='recommended',
        restart=False,
        kpts=(1, 1, 1))

catalyst.calc = calc

fix_z = 0.5 * catalyst.get_cell()[2][2]

c = FixAtoms(indices=[atom.index for atom in atoms if atom.position[2] <= fix_z])
catalyst.set_constraint(c)

catalyst = sort(catalyst)

en = catalyst.get_potential_energy()
print('Potential energy: {:.3f} eV'.format(en))

