#!/bin/bash

. ../mets.sh

atom=$(basename $PWD)

p=${PWD}
for((iB=0;iB<${#Bvals_fem[@]};iB++)); do
    wdir=${p}/${Bvals_fem[iB]}
    if [[ ! -d $wdir ]]; then
	mkdir -p $wdir
    fi
    cd $wdir
    for((iocc=0;;iocc++)); do
	   occfile="../../${atom}_${iocc}.occs"
           if [[ ! -f ${occfile} ]]; then
               #echo ${occfile} does not exist
               break
           fi
           Eref=0.0

           # Find out multiplicity
           M=$(cat ${occfile} |awk 'BEGIN {na=0; nb=0} {na+=$1; nb+=$2} END {print na-nb+1}')
           mmax=$(cat ${occfile} |awk 'BEGIN {mmax=0} {if($3>mmax) {mmax=$3}; if(-$3>mmax) {mmax=-$3}} END {print mmax}')
           # Reference calculations
           for((lfem=${mmax};;lfem+=2)); do
                  if [[ ! -f fem_${lfem}_${iocc}.log ]]; then
                      echo ${atom} fem occ=${iocc} lfem=${lfem} B=${Bvals_fem[iB]}
                      \cp -p ${occfile} occs.dat

		      # Load wave function from previous calculation
		      let lprevious=lfem-2
		      loadchk="fem_${lprevious}_${iocc}.chk"
		      if [[ -f $loadchk ]]; then
	  		#	  load="--load=$loadchk"
			  load=""
		      else
			  load=""
		      fi
		      echo $iB
		      if(($iB == 0)); then
			  rmax=100
			  nelem=7
		      else
			  rmax=40
			  nelem=5
		      fi
                      atomic --Z=${atom} --method=HF --lmax=${lfem} --mmax=${mmax} --Rmax=${rmax} --nelem=${nelem} --M=$M --Bz=${Bvals_fem[iB]} --readocc=-1 --iguess=0 --save=fem_${lfem}_${iocc}.chk ${load} &> fem_${lfem}_${iocc}.log
                  fi
                  if(( lfem > mmax )); then
                      let llast=lfem-2
                      Ecur=$(grep "Total                 energy:" fem_${lfem}_${iocc}.log |tail -n1 |awk '{print $NF}')
                      Elast=$(grep "Total                 energy:" fem_${llast}_${iocc}.log |tail -n1 |awk '{print $NF}')
                      dE=$(echo $Elast $Ecur | awk '{printf("%e\n",$1-$2)}')
                      convd=$(echo $dE $femthr | awk '{print $1<$2}')
                      echo "Ecur=$Ecur Elast=$Elast dE=$dE femthr=$femthr convd=$convd"
                      if(( convd )); then
                          Eref=${Ecur}
                          break
                      fi
                  fi
              done
	       # Remove checkpoints
	       \rm -v fem_*_${iocc}.chk
              done
              done
	       
