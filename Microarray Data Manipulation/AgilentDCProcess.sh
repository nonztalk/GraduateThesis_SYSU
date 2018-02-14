#!/bin/bash

for tar in *.tar; do
	echo "processing ${tar} ... "
	GSE=${tar%%_*}
	GPL=`cat /data2/Dengkaiwen/Data_agil/gpl_agil_re.csv | grep ${GSE} | awk -F, '{print $1;}' | uniq`
	File=${GPL}_${GSE}
	tar -xf ${tar} 
	ls | grep ".txt" | grep "GSM" | cut -d "." -f 1 | cut -d "_" -f 1 | paste <(ls | grep ".txt" | grep "GSM") - > fileInfo.txt
	Rscript runExpAgilentDC.R ${File} > ${File}.log 2>&1
	rm -rf GSM*
	rm fileInfo.txt
done