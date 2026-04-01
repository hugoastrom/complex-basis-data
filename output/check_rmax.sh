#!/bin/bash

. ../mets.sh

atom=$(basename $PWD)

p=${PWD}
wdir=${p}/0.00
if [[ ! -d $wdir ]]; then
    mkdir -p $wdir
fi
cd $wdir
for((iocc=0;;iocc++)); do
       occfile="../../${atom}_${iocc}.occs"
       if [[ ! -f ${occfile} ]]; then
           break
       fi
       Eref=0.0

       # Find out multiplicity
       M=$(cat ${occfile} |awk 'BEGIN {na=0; nb=0} {na+=$1; nb+=$2} END {print na-nb+1}')
       mmax=$(cat ${occfile} |awk 'BEGIN {mmax=0} {if($3>mmax) {mmax=$3}; if(-$3>mmax) {mmax=-$3}} END {print mmax}')
       # Reference calculations
       for((lfem=${mmax};;lfem+=2)); do
	      let lnext=lfem+2
              if [[ ! -f fem_${lnext}_${iocc}.log ]]; then
                  echo ${atom} fem occ=${iocc} lfem=${lfem} B=${Bvals_fem[iB]}
                  \cp -p ${occfile} occs.dat
                  atomic --Z=${atom} --method=HF --lmax=${lfem} --mmax=${mmax} --Rmax=40 --nelem=7 --M=$M --readocc=-1 --iguess=0 &> fem_${lfem}_${iocc}_small-rmax.log
		  break
              fi
	  done
           Ecur=$(grep "Total                 energy:" fem_${lfem}_${iocc}_small-rmax.log |tail -n1 |awk '{print $NF}')
           Elast=$(grep "Total                 energy:" fem_${lfem}_${iocc}.log |tail -n1 |awk '{print $NF}')
           dE=$(echo $Elast $Ecur | awk '{printf("%e\n",$1-$2)}')
           convd=$(echo $dE $femthr | awk '{print $1<$2}')
           echo "With new Rmax, Ecur=$Ecur and Elast=$Elast dE=$dE femthr=$femthr convd=$convd"
           if(( convd )); then
               Eref=${Ecur}
               \rm -v fem_${lfem}_${iocc}_small-rmax.log
	   else
	       echo "Computing with new Rmax"
	       olddir=${wdir}/old
	       if [[ ! -d $olddir ]]; then
		   mkdir -p $olddir
	       fi
	       \mv -v fem_*_${iocc}.log $olddir
	       for((lfem=${mmax};;lfem+=2)); do
		      if [[ ! -f fem_${lfem}_${iocc}.log ]]; then
			  echo ${atom} fem occ=${iocc} lfem=${lfem} B=${Bvals_fem[iB]}
			  atomic --Z=${atom} --method=HF --lmax=${lfem} --mmax=${mmax} --Rmax=40 --nelem=7 --M=$M --readocc=-1 --iguess=0 &> fem_${lfem}_${iocc}.log
		      fi
		      if(( lfem > mmax )); then
			 let lprev=lfem-2
			 Ecur=$(grep "Total                 energy:" fem_${lfem}_${iocc}.log |tail -n1 |awk '{print $NF}')
			 Elast=$(grep "Total                 energy:" fem_${lprev}_${iocc}.log |tail -n1 |awk '{print $NF}')
			 dE=$(echo $Elast $Ecur | awk '{printf("%e\n",$1-$2)}')
			 convd=$(echo $dE $femthr | awk '{print $1<$2}')
			 echo "Ecur=$Ecur Elast=$Elast dE=$dE femthr=$femthr convd=$convd"
			 if(( convd )); then
			     Eref=${Ecur}
			     break
			 fi
		      fi
		      # Remove checkpoints
		      \rm -v fem_*_${iocc}.chk
		      \rm -v helfem.chk
		  done
		  fi
		  done

