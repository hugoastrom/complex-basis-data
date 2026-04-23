#!/bin/bash

#export PATH=$PATH:/home/work/erkale/openmp/src:/home/work/erkale/openmp/src/basistool:/home/work/HelFEM/objdir/src/:/home/work/erkale/serial/src:/home/work/erkale/serial/src/basistool:/home/work/HelFEM/objdir/src/
#export ERKALE_LIBRARY=/home/work/erkale/basis/

. ../back_mets.sh

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
    let ichk=iB-1
    wdir=${p}/${Bvals_gto[iB]}
    if [[ ! -d $wdir ]]; then
	mkdir -p $wdir
    fi
    cd $wdir
    for((iocc=0;;iocc++)); do
    #for ((iocc=8;;iocc++)); do
	   occfile="../../${atom}_${iocc}.occs"
           if [[ ! -f ${occfile} ]]; then
               #echo ../${atom}_${iocc}.occs does not exist
               break
           fi
           Eref=0.0

           # Find out multiplicity
           M=$(cat ${occfile} |awk 'BEGIN {na=0; nb=0} {na+=$1; nb+=$2} END {print na-nb+1}')
           mmax=$(cat ${occfile} |awk 'BEGIN {mmax=0} {if($3>mmax) {mmax=$3}} END {print mmax}')

           # GTO calculations
           for basis in ${bases[@]}; do
	       if [ $iB -eq 0 ]; then
		   LoadChk=""
	       else
		   LoadChk="LoadChk ../${ichk}_${iocc}_${basis}.chk"
	       fi
	       cat > ${basis}_${iocc}.run <<EOF
System ../atoms.xyz
Basis $basis
Decontract *
Method HF
LinearOccupations -1
LinearOccupationFile ${occfile}
LinearSymmetry true
LinearB ${Bvals_gto[iB]}
Multiplicity $M
Guess sapfit
OptLM false
Logfile ${basis}_${iocc}.log
SaveChk ../${iB}_${iocc}_${basis}.chk
$LoadChk
EOF
               if [[ ! -f ${basis}_${iocc}.log ]]; then
                   echo ${atom} ${basis}_${iocc}
                   erkale_omp ${basis}_${iocc}.run &> ${basis}_${iocc}.stdout
               fi
           done
	   if [ -n $Loadchk ]; then
	       \rm -v ../${ichk}_${iocc}_*.chk
	   fi
       done
       done
