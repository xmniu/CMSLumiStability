#!/usr/bin/env python
import os
import subprocess
import commands

#fills = ['Fill5401','Fill5405','Fill5406','Fill5416','Fill5418','Fill5421','Fill5423','Fill5424','Fill5427','Fill5433','Fill5437']
#fills = ['Fill5205','Fill5206','Fill5209','Fill5210','Fill5211','Fill5213','Fill5219','Fill5222','Fill5223','Fill5229','Fill5251','Fill5253','Fill5254','Fill5256','Fill5257','Fill5258']
#fills = ['Fill5261','Fill5264','Fill5265','Fill5266','Fill5267','Fill5270','Fill5274','Fill5275','Fill5276','Fill5277','Fill5279','Fill5282','Fill5287','Fill5288']
#fills = ['Fill5331','Fill5338','Fill5339','Fill5340','Fill5345','Fill5352','Fill5355','Fill5386','Fill5393','Fill5394','Fill5395']
fills = ['Fill5340','Fill5393','Fill5394']
for fill in fills:
#for ibound in range(0,8):
#for ibound in [0]:
#    for iIDvar in range(0,11):
#    for iIDvar in [11]:
	cmd = "/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV/LumiStability/runDoSelection.sh "+fill
	print cmd
	bsubs_cmd = "bsub -q 8nh -R 'pool > 4000' -C 0 -o" + \
		"/afs/cern.ch/user/x/xniu/ZCounting/CMSSW_8_0_14/src/MitEwk13TeV/LumiStability/output "+cmd
	status,output=commands.getstatusoutput(bsubs_cmd)
	print output
