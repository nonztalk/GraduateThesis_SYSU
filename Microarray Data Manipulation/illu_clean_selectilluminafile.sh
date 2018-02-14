#!/usr/bin/bash

input=$1

GSE=(`cat ${input} | awk -F, '{print $1}'`)
GSE=(`echo ${GSE[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '`)

cd ./Data

for gse in ${GSE[@]}; do
	echo "Processing ${gse} ... "
	filename=`ls | grep ^${gse}`
	if [[ -z ${filename} ]]; then
		echo "${gse} missing ... "
		echo ${gse} >> missfile_illu
		continue
	fi
	cp ${filename} ../Data_illu
done





