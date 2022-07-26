import sys
import numpy as np
from ase import Atom, Atoms
from ase.io import read, write
from ase.build import sort, surface, add_adsorbate
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp

infile_uc = sys.argv[1]
vacuum = 7.5
miller_idx = [(1,0,0), (1,1,1), (1,2,1), (1,0,1), (0,1,1)]
nrepeat_x = 1
nrepeat_y = 1
nlayers = 4
kpoints = 4

# Create Bulk Material
unit_cell = read(infile_uc)

#Optimize bulk structure
calc_bulk_volume = Vasp(xc='pbe', ivdw=11,
    command="ibrun vasp_std_vtst",
	directory='bulk',
	encut=520,
    prec='Normal',
    ismear=2,
    sigma=0.2,
	ispin=2,
    ibrion=2,
	isif=3,
    nsw = 1000,
    ediff=1.0E-5,
    ediffg=-0.03,
    algo='V',
	lorbit=11,
    setups='recommended',
    restart=False,
    kpts=(kpoints, kpoints, kpoints))

unit_cell.calc = calc_bulk_volume
energy_bulk_volume = unit_cell.get_potential_energy()


calc_bulk = Vasp(xc='pbe', ivdw=11,
    command="ibrun vasp_std_vtst",
	directory='bulk',
    isif=2,
	nsw=1000,
	restart=True)

unit_cell.calc = calc_bulk
energy_bulk = unit_cell.get_potential_energy()


print('Bulk Lattice Vector a: ', unit_cell.get_cell()[0])
print('Bulk Lattice Vector b: ', unit_cell.get_cell()[1])
print('Bulk Lattice Vector c: ', unit_cell.get_cell()[2])
print('Bulk Energy: {:.3f} eV'.format(energy_bulk))
print('')

SURFACE = []
ENERGY = []
MIDX = []

for idx in miller_idx:

	#Create Surface from optimized UC
	#Surface has two surface terminations, need to check both!
	for t in ("T1", "T2"):
		system = "slab_"+t+"/"+str(idx)

		slab = surface(read("bulk/CONTCAR"), idx, nlayers, vacuum).repeat((nrepeat_x, nrepeat_y, 1))
		fix_z =  0.5 * slab.get_cell()[2][2]

		c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] <= fix_z])
		slab.set_constraint(c)
		slab.set_pbc(True)
		slab = sort(slab)
		#Calculate Slab Energies
		calc_slab_scf = Vasp(xc='pbe', ivdw=11,
        command="ibrun vasp_std_vtst",
		directory=system,
		encut=520,
        prec='Normal',
        ismear=2,
        sigma=0.2,
		ispin=2,
        ibrion=-1,
		isif=2,
        nsw = 0,
        ediff=1.0E-5,
        ediffg=-0.03,
		algo='V',
		setups='recommended',
		lorbit=11,
		restart=False,
		kpts=(kpoints, kpoints, 1))


		slab.calc = calc_slab_scf
		energy_slab_unrelaxed = slab.get_potential_energy()


		calc_slab_go = Vasp(xc='pbe', ivdw=11,
		command="ibrun vasp_std_vtst",
			directory=system,
		nsw = 1000,
			ibrion = 2,
			restart=True)

        slab.calc = calc_slab_go
		energy_slab_relaxed = slab.get_potential_energy()


		# Do some Math
		n_atoms_slab = len(slab)
		n_atoms_bulk = len(unit_cell)

		relaxation_energy = energy_slab_relaxed - energy_slab_unrelaxed

		a = slab.get_cell()[0]
		b = slab.get_cell()[1]
		a_cross_b = np.cross(a,b)
		cross_sectional_area = np.linalg.norm(a_cross_b)

		surface_energy_unrelaxed = 0.5 * (energy_slab_unrelaxed - n_atoms_slab * energy_bulk / n_atoms_bulk)
		surface_energy = (surface_energy_unrelaxed + relaxation_energy)/ cross_sectional_area


		#Bookkeeping
		SURFACE.append(system)
		ENERGY.append(surface_energy)
		MIDX.append(idx)

		if None in ENERGY:
			calc_slab_go.abort()


print("Surface | Miller Index |Surface Energy (eV / Ang**2)")
print("-------------------------------------")
for i in range(0, len(ENERGY)):
	print(SURFACE[i]," | ", MIDX[i]," | ", ENERGY[i])

