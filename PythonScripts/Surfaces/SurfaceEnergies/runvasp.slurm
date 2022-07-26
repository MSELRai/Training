#!/bin/bash
#SBATCH --job-name=@@@        
#SBATCH --output=@@@.o%j      
#SBATCH --error=@@@.e%j       
#SBATCH --partition=skx-normal    
#SBATCH --nodes=1                 
#SBATCH --ntasks=48               
#SBATCH --time=01:30:00
#SBATCH --mail-user=wnw36@msstate.edu
#SBATCH --mail-type=ALL
#SBATCH --account=TG-CHE140141    

module swap intel/17.0.4
module load vasp/5.4.4.pl2


script="go.py"
system="@@@"


./$script $system
#ibrun vasp_gam > my_vasp.out
