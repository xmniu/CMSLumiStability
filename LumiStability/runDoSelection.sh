#!/bin/bash
datafilelabel=$1

workdir="/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"
subdir="Lumi"

datatritype="tag_noIso_probe_noIso"

#TOP="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Lumi"

CMSSW_BASE="/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/"
TOP="$PWD"

cd $CMSSW_BASE
eval `scramv1 runtime -sh`
cd $TOP

root -l -b << EOF
gSystem->Load("${workdir}/Selection/selectZmm_batch_C.so")
selectZmm("${workdir}/Selection/Conf/zmm_${datafilelabel}_eos.conf", "${TOP}/${datafilelabel}", 0, 0, 5, 0)
.q
EOF

cp ${TOP}/${datafilelabel}/ntuples/data_select.root ${filedir}/${subdir}/2016/Mu/ntuples/data_select_${datafilelabel}_${datatritype}.root
