#!/usr/bin/bash

for tar in *.tar; do
	GSE=${tar%%_*}
	GPL=`cat ../gpl_affy_re.csv | grep ${GSE} | awk -F, '{ print $1; }' | uniq`
	File=${GPL}_${GSE}
	tar -xvf ${tar}
	rm *.CHP.gz *.chp.gz
	ls *.gz | cut -d "." -f 1 | cut -d "_" -f 1 | paste <(ls *.gz) -> fileInfo.txt
	Rscript runExpAffyBasic.R ${File} > ${File}.log 2>&1
	rm -f *.gz
done