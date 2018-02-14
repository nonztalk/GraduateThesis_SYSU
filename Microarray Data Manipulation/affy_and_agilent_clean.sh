#!/bin/bash
# Suitable for Affymetrix and Agilent

inputFile=$1 # under ./Data

GSE=(`cat ${inputFile} | awk -F, '{ print $1; }'`)
GSE=(`echo ${GSE[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '`)

for gse in ${GSE[@]}; do

	echo "Processing ${gse}..."
	filename=`ls | grep ^${gse}_RAW`

	if [[ -z ${filename} ]]; then
		echo ${gse} >> missfile
		echo "${gse} missing ... "
		continue
	fi

	mkdir ./processed # creat ./Data/processed
	mkdir ./keep # creat ./Data/keep
	tar -xf ${filename} -C ./processed

    if [[ $? != 0 ]]; then
    	echo ${filename} >> errfile
    	echo "${filename} can not open ... "
    	rm -rf ./processed
    	rm -rf ./keep
    	continue
    fi

	GSM=(`cat ${inputFile} | grep ^${gse} | awk -F, '{ print $2; }'`)
	cd ./processed # current ./Data/processed

	for gsm in ${GSM[@]}; do
		cp -r ./`ls | grep ^${gsm}` ../keep
	done

	cd ../keep # current ./Data/keep
	tar -cf ../../Data_affy/${filename} *.gz 
	cd ../ # current ./Data
	rm -rf ./processed
	rm -rf ./keep 

done






