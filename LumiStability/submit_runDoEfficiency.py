#!/usr/bin/env python
import os
import subprocess
import commands

for ibound in range(0,15):
#for ibound in [6]:
#    for iIDvar in range(0,11):
#    for iIDvar in [11]:
	cmd = "/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV/LumiStability/runDoEfficiency.sh Mu "+str(ibound)+" Fill5401"
	print cmd
	bsubs_cmd = "bsub -q 1nh -R 'pool > 4000' -C 0 -o" + \
		"/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV/LumiStability/output "+cmd
	status,output=commands.getstatusoutput(bsubs_cmd)
	print output
