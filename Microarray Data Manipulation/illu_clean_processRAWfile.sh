#!/usr/bin/bash

# Current direction Data_illu
for gse in `ls | grep _RAW.tar`; do
	mkdir ./process
	echo "Processing ${gse} ... "
	tar -xf ${gse} ./process
	if [[ $? != 0 ]]; then
		echo "${gse} error"
		echo ${gse} >> errfile_illu
		rm -rf ./process
		continue
	fi
	cd ./process
	rm -f `ls | grep -v *.txt.gz`
	tar -cf ../../Data_illu/${gse} *
	cd ../
	rm -rf ./process
done
