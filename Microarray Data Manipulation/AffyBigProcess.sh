#!/bin/bash

for tar in *.tar; do
	echo "processing ${tar}"
	GSE=${tar%%_*}
	GPL=`cat /data2/Dengkaiwen/Data_affy/AffyBig/gpl_affy_re.csv | grep ${GSE} | awk -F, '{print $1;}' | uniq`
	HuName=`cat /data2/Dengkaiwen/Data_affy/AffyBig/GPLHuGene | grep ${GPL} | awk -F, '{print $2;}'`
	gplName=`ls /data2/ArrayReAnnotation/Library | grep ${HuName}`
	File=${GPL}_${GSE}
	LibraryPath=/data2/ArrayReAnnotation/Library/${gplName}
	tar -xf ${tar}
	rm *.CHP.gz *.chp.gz
	gunzip *.gz
	ls *.CEL | cut -d "." -f 1 | cut -d "_" -f 1 | paste <(ls *.CEL) - > fileInfo.txt
	apt-probeset-summarize -a rma-bg,quant-norm,pm-only,plier -c ${LibraryPath}/${gplName}.clf -p ${LibraryPath}/${gplName}.pgf --qc-probesets ${LibraryPath}/${gplName}.qcc -o ${File}_genelevel *.CEL
	grep -v "^#" *_genelevel/rma-bg.quant-norm.pm-only.plier.summary.txt > ${File}_norm.txt
	Rscript runExpAffyBig.R ${File} ${HuName} > ${File}.log 2>&1
	rm -rf *_genelevel
	rm -f GSM*
	rm -f fileInfo.txt
done