#!/bin/bash
inputjson=$1

#set up brilcalc environment
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH

#prepare bcidLumiDCS.csv for next step 
brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i ${inputjson} --byls -o bcidLumiDCS.csv -u /pb
#brilcalc lumi --normtag pccLUM15001 -i lumi.txt --byls -o bcidLumiDCS.csv -u /pb

#before proceeding next step, check if you need to update total luminosity and the granularity you want to have for the stability study 
#if yes, do next line to get total luminosity
#brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json -i Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt -u /pb

#write the output luminosity boundary arrays into a text file to be used in selection of Zmm events 
python prodlumibound.py > lumibound.txt
