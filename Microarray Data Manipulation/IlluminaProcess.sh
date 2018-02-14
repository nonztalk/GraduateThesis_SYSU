#!/bin/bash

for file in *.txt.gz; do
	GSE=${file%%_*}
	GPL=`cat /data2/Dengkaiwen/Data_illu/gpl_illu_re.csv | grep ${GSE} | awk -F, '{print $1;}' | uniq`
	Rscript runExpIllumina.R ${file} ${GPL} > ${GPL}_${GSE}.log 2>&1
done
