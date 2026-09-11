import os
import mets

p = os.getcwd()
mols = mets.mols()

for mol in mols:
    wdir = f"{p}/{mol}"
    if not os.path.exists(f"{wdir}"):
        os.system(f"mkdir -p {wdir}")
    os.chdir(f"{wdir}")
    
    os.system("python3 ../makeinp.py\n")
