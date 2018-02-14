#!/bin/bash

for File in *.txt.gz; do
	gunzip ${File}
	Target=${File%%_*}
	Rscript runExpAffymiss.R ${Target} > ${Target}.log 2>&1 
done
