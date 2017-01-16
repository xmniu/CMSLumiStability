#!/bin/bash

outputdir=$1
min=$2
max=$3

#set up brilcalc environment
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH

for ((jsonnum=${min}; jsonnum<${max};jsonnum++)); do
#brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json --hltpath HLT_Ele27_WPTight_Gsf_v* -i JSON/lumi_${jsonnum}.txt -u /pb | tail -n 4 | head -n 1 | awk '{print $12}'
#brilcalc lumi --normtag pccLUM15001 --hltpath HLT_Ele27_WPTight_Gsf_v* -i JSON/lumi_${jsonnum}.txt -u /pb | tail -n 4 | head -n 1 | awk '{print $12}'
#brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json --hltpath HLT_IsoMu24_v* -i JSON/lumi_${jsonnum}.txt -u /pb | tail -n 4 | head -n 1 | awk '{printf $12}'
brilcalc lumi --normtag pccLUM15001 -i ${outputdir}/JSON/lumi_${jsonnum}.txt -u /pb --hltpath HLT_IsoMu24_v* | tail -n 4 | head -n 1 | awk '{printf $12}' 
#brilcalc lumi -i JSON/lumi_${jsonnum}.txt -u /pb --hltpath HLT_IsoMu24_v* | tail -n 4 | head -n 1 | awk '{printf $12}'
#brilcalc lumi --hltpath HLT_IsoMu24_v* -i JSONEle/lumi_${jsonnum}.txt -u /pb | tail -n 4 | head -n 1 | awk '{print $12}'
#brilcalc lumi --hltpath HLT_Ele27_WPTight_Gsf_v* -i JSONEle/lumi_${jsonnum}.txt -u /pb | tail -n 4 | head -n 1 | awk '{print $12}'
printf ", "
done
