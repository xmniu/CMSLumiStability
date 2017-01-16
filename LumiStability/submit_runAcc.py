#!/usr/bin/env python
import os
import subprocess
import commands

for ibound in range(0,49):
#for ibound in [62]:
	cmd = "/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV/LumiStability/runAcc_fornoIso.sh Mu "+str(ibound)+" Fill5287 ReReco"
	print cmd
	bsubs_cmd = "bsub -q 8nh -R 'pool > 4000' -C 0 -o" + \
		"/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV/LumiStability/output "+cmd
	status,output=commands.getstatusoutput(bsubs_cmd)
	print output
