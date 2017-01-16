#!/bin/bash

lepton=$1
ibound=$2
datafilelabel=$3

workdir="/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV"
filedir="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency"
subdir="Lumi"

mcfiledir="80XMC"
mcfilelabel="80XwoIso"
mctritype="tag_IsoMu22_probe_IsoMu22"

#datafilelabel="SyncATLAS"
#datafilelabel="Fill5199ReReco"
#datatritype="tag_IsoMu24_probe_IsoMu24"
datatritype="tag_noIso_probe_noIso"

#TOP="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Lumi"

CMSSW_BASE="/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/"
TOP="$PWD"

cd $CMSSW_BASE
eval `scramv1 runtime -sh`
cd $TOP

root -l -b << EOF
gSystem->Load("${workdir}/Efficiency/selectProbes${lepton}Eff_Data_25invpb_C.so")
selectProbes${lepton}Eff_Data_25invpb("${filedir}/${subdir}/2016/${lepton}/ntuples/data_select_${datafilelabel}_${datatritype}.root","$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}HLTEff/Data", 0, 0, 0, ${ibound})
selectProbes${lepton}Eff_Data_25invpb("${filedir}/${subdir}/2016/${lepton}/ntuples/data_select_${datafilelabel}_${datatritype}.root","$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}StaEff/Data", 4, 0, 0, ${ibound})
selectProbes${lepton}Eff_Data_25invpb("${filedir}/${subdir}/2016/${lepton}/ntuples/data_select_${datafilelabel}_${datatritype}.root","$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}SITEff/Data", 8, 0, 0, ${ibound})

gSystem->Load("${workdir}/Efficiency/selectProbes${lepton}Eff_MC_25invpb_C.so")
selectProbes${lepton}Eff_MC_25invpb("${filedir}/${subdir}/${mcfiledir}/${lepton}/ntuples/zmm_select_${mcfilelabel}_norw_${mctritype}.root", "${workdir}/Tools", "$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}HLTEff/MC", 0, 1, 1, ${ibound})
selectProbes${lepton}Eff_MC_25invpb("${filedir}/${subdir}/${mcfiledir}/${lepton}/ntuples/zmm_select_${mcfilelabel}_norw_${mctritype}.root", "${workdir}/Tools", "$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}StaEff/MC", 4, 1, 1, ${ibound})
selectProbes${lepton}Eff_MC_25invpb("${filedir}/${subdir}/${mcfiledir}/${lepton}/ntuples/zmm_select_${mcfilelabel}_norw_${mctritype}.root", "${workdir}/Tools", "$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}SITEff/MC", 8, 1, 1, ${ibound})

gSystem->Load("${workdir}/Efficiency/TagAndProbe/plotEff_C.so")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}BEpteta.bins",0,0,0,0,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}HLTEff/MC/probes.root","${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}HLTEff/2CT","png",1,0,0,"${lepton}","HLT",0.7,1.02,10)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}BEpteta.bins",0,0,0,0,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}SITEff/MC/probes.root","${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}SITEff/2CT","png",1,0,0,"${lepton}","SIT",0.7,1.02,10)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}BEpteta.bins",0,0,0,0,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}StaEff/MC/probes.root","${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}StaEff/2CT","png",1,0,0,"${lepton}","Sta",0.7,1.02,10)
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}BEpteta.bins",0,0,0,0,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}HLTEff/Data/probes.root","${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}HLTEff/2MG","png",1,0,0,"${lepton}","HLT",0.7,1.02,10,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}HLTEff/MC/probes.root")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}BEpteta.bins",2,6,2,6,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}StaEff/Data/probes.root","${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}StaEff/2MG","png",1,0,0,"${lepton}","Sta",0.7,1.02,10,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}StaEff/MC/probes.root")
plotEff("${workdir}/Efficiency/TagAndProbe/${lepton}BEpteta.bins",2,1,2,1,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}SITEff/Data/probes.root","${filedir}/${subdir}/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}SITEff/2MG","png",1,0,0,"${lepton}","SIT",0.7,1.02,10,"$TOP/${lepton}_${datafilelabel}_${datatritype}_DataSet_${ibound}/${lepton}SITEff/MC/probes.root")
.q
EOF
