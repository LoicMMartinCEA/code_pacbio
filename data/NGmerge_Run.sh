#!/bin/bash

inFile_1=$1
inFile_2=$2 
nbCore=$3

logFile=${inFile_1:0:-11}"_ngmerge_details.log"
outFile=${inFile_1:0:-11}"_NGmerge_assembled.fastq"

#compress output (z) or not (y)
compress=z

echo "debut"

time NGmerge -n ${nbCore} -v -${compress} -1 ${inFile_1} -2 ${inFile_2} -o ${outFile} > ${logFile}

echo "fin"
