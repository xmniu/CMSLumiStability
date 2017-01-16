#!/bin/bash
inputfile=$1
inputjson=$2
outputdir=$3
min=$4
max=$5

for ((jsonnum=${min}; jsonnum<${max};jsonnum++)); do

root -l -b << EOF
gSystem->Load("selectRunRun_C.so")
selectRunRun("${inputfile}","${outputdir}/FlatNtuple/test_${jsonnum}.root",${jsonnum})
.q
EOF

echo "${outputdir}/FlatNtuple/test_${jsonnum}.root" > ${outputdir}/List/in_${jsonnum}.txt
python makeJson.py --mask ${inputjson} --out ${outputdir}/JSON/lumi_${jsonnum}.txt --tree Events --list ${outputdir}/List/in_${jsonnum}.txt --lumi-branch lumiSec
pileupCalc.py -i  ${outputdir}/JSON/lumi_${jsonnum}.txt --inputLumiJSON pileup_latest.txt  --calcMode true --minBiasXsec 69200 --maxPileupBin 60 --numPileupBins 60 ${outputdir}/PileupHistogram/pileup_${jsonnum}.root

root -l -b << EOF
gSystem->Load("make_f_rw_80X_C.so")
make_f_rw_80X("${outputdir}", $jsonnum)

gSystem->Load("check_f_rw_80X_C.so")
check_f_rw_80X("${outputdir}", $jsonnum)
.q
EOF

done
