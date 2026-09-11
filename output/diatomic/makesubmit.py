import os
import mets

p = os.getcwd()
mols = mets.mols()

for mol in mols:
    wdir = f"{p}/{mol}"
    if not os.path.exists(f"{wdir}"):
        os.system(f"mkdir -p {wdir}")
    os.chdir(f"{wdir}")
    
    f = open(f"{mol}.in", "w")
    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write(f"#SBATCH --job-name={mol}\n")
    f.write("#SBATCH -e out_%j\n")
    f.write("#SBATCH -o out_%j\n")
    f.write("#SBATCH --mem-per-cpu=500\n")
    f.write("#SBATCH -t 1-00:00:00\n")
    f.write("#SBATCH -n 8\n")
    f.write("#SBATCH -p normal\n")
    f.write("\n")
    f.write("python3 ../makeinp.py\n")
    f.close()
    
    print(f"submited {mol}")
    os.system(f"sbatch {mol}.in")
