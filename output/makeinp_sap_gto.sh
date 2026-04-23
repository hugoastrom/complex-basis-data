#!/bin/bash

. ../mets.sh

atom=$(basename $PWD)
if [[ -f ../${atom}_0.occs ]]; then
    if [[ ! -f atoms.xyz ]]; then
        cat > atoms.xyz <<EOF
1
$atom
$atom 0 0 0
EOF
    fi
fi

p=${PWD}
for((iB=0;iB<${#Bvals_gto[@]};iB++)); do
    wdir=${p}/${Bvals_gto[iB]}
    if [[ ! -d $wdir ]]; then
	mkdir -p $wdir
    fi
    cd $wdir
    for((iocc=0;;iocc++)); do
	   occfile="../../${atom}_${iocc}.occs"
           if [[ ! -f ${occfile} ]]; then
               break
           fi

           # Find out multiplicity
           M=$(cat ${occfile} |awk 'BEGIN {na=0; nb=0} {na+=$1; nb+=$2} END {print na-nb+1}')
           mmax=$(cat ${occfile} |awk 'BEGIN {mmax=0} {if($3>mmax) {mmax=$3}} END {print mmax}')

           # GTO calculations
           for basis in ${bases[@]}; do
               cat > sap_${basis}_${iocc}.run <<EOF
System ../atoms.xyz
Basis ${basis}
Decontract *
Method HF
LinearOccupations -1
LinearOccupationFile ${occfile}
LinearSymmetry true
LinearB ${Bvals_gto[iB]}
Multiplicity $M
Guess sapfit
OptLM false
Complexbas true
MaxInitIter 100
EOF
	       conv=$(grep "Converged to" sap_${basis}_${iocc}.stdout)
               if [ -z "$conv" ]; then
                   echo ${atom} ${basis} ${iocc}
                   erkale_complex_orbs_omp sap_${basis}_${iocc}.run &> sap_${basis}_${iocc}.stdout
               fi
           done
       done
       done
