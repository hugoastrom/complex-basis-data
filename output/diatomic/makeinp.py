import os
import numpy as np
import mets

def makesubmit(mols, rvals):
    p = os.getcwd()

    mol = p.split("/")[-1]
    for r in rvals:
        wdir = f"{p}/{r}"
        if not os.path.exists(f"{wdir}"):
            os.system(f"mkdir -p {wdir}")
        os.chdir(wdir)

        xyzfile = open(f"{mol}.xyz", "w")
        xyzfile.write("2\n")
        xyzfile.write("\n")
        xyzfile.write(f"{mols[mol]['atoms'][0]} 0.0 0.0 0.0\n")
        xyzfile.write(f"{mols[mol]['atoms'][1]} 0.0 0.0 {r}\n")
        xyzfile.close()

        for M in mols[mol]["M"]:

            runfile = open(f"{M}.run", "w")
            runfile.write(f"System {mol}.xyz\n")
            runfile.write("Basis aug-cc-pVTZ\n")
            runfile.write("Method HF\n")
            runfile.write("LinearSymmetry true\n")
            runfile.write(f"Multiplicity {M}\n")
            runfile.write("Guess core\n")
            runfile.write("OptLM false\n")
            runfile.write("ComplexBasis true\n")
            runfile.close()

            conv = False
            try:
                with open(f"{M}.stdout") as f:
                    for line in f:
                        if "Converged to" in line:
                            conv = True
            except:
                conv = False
            if not conv:
                os.system(f"erkale_complex_orbs_omp {M}.run &> {M}.stdout")

def main():
    mols = mets.mols()
    rvals = mets.rvals()
    p = os.getcwd()
    
    makesubmit(mols, rvals)

if __name__=="__main__":
    main()
