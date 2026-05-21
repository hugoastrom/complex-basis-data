#!/bin/python

import os

systems = {"B": {5: {"aug-cc-pVTZ": ["0.40", "0.38", "0.36", "0.34", "0.32", "0.30", "0.28", "0.26", "0.24", "0.22", "0.20", "0.18", "0.16", "0.14", "0.12", "0.10"]}}, "P": {8: {"AHGBSP3-9": ["0.10", "0.08", "0.06", "0.04", "0.02", "0.00"]}}, "Cl": {8: {"AHGBSP3-9": ["0.10", "0.08", "0.06", "0.04", "0.02", "0.00"], "aug-cc-pVTZ": ["0.12", "0.10"]}}, "Ar": {8: {"AHGBSP3-9": ["0.08", "0.06", "0.04", "0.02", "0.00"], "aug-cc-pVTZ": ["0.12", "0.10"]}}}

p = os.getcwd()
atom = p.split("/")[-1]
if atom in systems:
    for state in systems[atom]:
        for basis in systems[atom][state]:
            for ib,b in enumerate(systems[atom][state][basis]):
                wdir = f"{p}/{b}"
                os.chdir(wdir)
                occfile = f"../../{atom}_{state}.occs"
                na = 0
                nb = 0
                with open(occfile) as f:
                    for line in f:
                        na += int(line.split()[0])
                        nb += int(line.split()[1])
                M = na - nb  + 1

                if ib == 0:
                    chk = f"SaveChk {atom}.chk\n"
                    e_core = 0.0
                    with open(f"{basis}_{state}.stdout") as f:
                        for line in f:
                            if "Converged" in line:
                                e_core = float(line.split()[-1].replace("!",""))

                    e_sap = 0.0
                    with open(f"sap_{basis}_{state}.stdout") as f:
                        for line in f:
                            if "Converged" in line:
                                e_sap = float(line.split()[-1].replace("!",""))

                    guess = "sapfit" if e_sap < e_core else "core"

                else:
                    guess = "sapfit"
                    bprev = systems[atom][state][basis][ib - 1]
                    chk = f"LoadChk ../{bprev}/{atom}.chk\nSaveChk {atom}.chk\n"
                
                f = open(f"erkale_guess_{basis}_{state}.run", "w")
                f.write("System ../atoms.xyz\n")
                f.write(f"Basis {basis}\n")
                f.write("Decontract *\n")
                f.write("Method HF\n")
                f.write("LinearOccupations -1\n")
                f.write(f"LinearOccupationFile {occfile}\n")
                f.write("LinearSymmetry true\n")
                f.write(f"LinearB {b}\n")
                f.write(f"Multiplicity {M}\n")
                f.write(f"Guess {guess}\n")
                f.write("OptLM false\n")
                f.write("Complexbas true\n")
                f.write(chk)
                f.close()

                os.system(f"erkale_complex_orbs_omp erkale_guess_{basis}_{state}.run  &> erkale_guess_{basis}_{state}.stdout")
                os.system(f"grep 'Converged' erkale_guess_{basis}_{state}.stdout")
