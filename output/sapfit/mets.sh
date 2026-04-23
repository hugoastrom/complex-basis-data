#!/bin/bash

# Magnetic field values
Bvals_gto=( $(echo ""|awk '{for(i=0;i<=60;i+=2) {printf(" %.2f",i/100)}}') )

# Gaussian basis sets
#bases=({,aug-}cc-pV{D,T,Q,5}Z {,A}HGBSP{1,2,3}-{5,7,9})
bases=(aug-cc-pVTZ AHGBSP3-9)
# 6-311++G(3df,3pd) def2-TZVP)
