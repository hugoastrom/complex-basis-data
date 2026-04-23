#!/bin/bash


p=$PWD
for at in B Be F He Li O; do
    if [[ ! -d $p/$at ]]; then
	mkdir -p $p/$at
    fi
    cd $p/$at
    #for calc in fem gto; do
    for calc in gto; do
	cat > submit_${calc}.in <<EOF
#!/bin/bash

#SBATCH --job-name=${at}-mag-fields-${calc}
#SBATCH -e out_%j
#SBATCH -o out_%j
#SBATCH --mem-per-cpu=500
#SBATCH -t 14-00:00:00
#SBATCH -n 8
#SBATCH -p normal

../makeinp_sap_${calc}.sh
#../check_rmax.sh
#../makeinp_real_gto.sh

EOF
	if [[ ! -f .submit_${calc} ]]; then
	    echo submit ${at} ${calc}
	    sbatch submit_${calc}.in | tee .submit_${calc}
	fi
    done
done
