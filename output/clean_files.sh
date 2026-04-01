\rm *.occs
for atom in Al Ar B Be C Cl F H He Li Mg N Na Ne O S Si P; do
    \rm ${atom}/*/*.{run,dat,chk}
done
\rm */atoms.xyz */out* */submit*
