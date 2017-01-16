#!/bin/bash

lepton=$1
ibound=$2
datafilelabel=$3
#recotype=$4

workdir="/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"
subdir="Lumi"

#datafilelabel="2016H"
datatritype="tag_noIso_probe_noIso"

CMSSW_BASE="/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/"
TOP="$PWD"

cd $CMSSW_BASE
eval `scramv1 runtime -sh`
cd $TOP

#TOP="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Lumi"

root -l -b << EOF
gSystem->Load("${workdir}/Acceptance/CEffUser1D_cc.so")
gSystem->Load("${workdir}/Acceptance/CEffUser2D_cc.so")
gSystem->Load("${workdir}/Acceptance/computeAccSelZmmBinned_per25invpb_C.so")
computeAccSelZmmBinned_per25invpb("${workdir}/Acceptance/zmm_80X_eos.conf", "${workdir}/Tools", "${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/", "${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/", "${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/Z/", 1, 4, 0, ${ibound})
.q
EOF
