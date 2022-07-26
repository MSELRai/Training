from ase.io import read, write
from ase.build import molecule
from ase.calculators.vasp import Vasp

infile = 'POSCAR-optimized-MFIReO3.vasp'

# Need to add dispersion interactions xc = optb88-vdw or ivdw = 11
# Compare DFT+D3 and optB88-vdW

atoms = read(infile)

calc = Vasp(xc='pbe',
        command="ibrun vasp_gam",
        encut=520,
        prec='Normal',
        ivdw=11,
        ismear=0,
        sigma=0.05,
        ibrion=2,
        isif=3,
        nsw = 1000,
        ediff=1.0E-5,
        ediffg=-0.03,
        algo='V',
        lreal='Auto',
        ncore=48,
        setups='recommended',
        restart=False,
        gamma=True,
        kpts=(1, 1, 1))

atoms.calc = calc
en = atoms.get_potential_energy()
print('Potential energy: {:.3f} eV'.format(en))
