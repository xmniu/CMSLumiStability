#!/bin/env python
csvFile=open("bcidLumiDCS.csv")
lines=csvFile.readlines()

LS_bound=[]
Run_bound=[]
Fill_bound=[]
del_lumi=[]
rec_lumi=[]

lumiBinner=0
lumiDelBlock=0

for line in range(2,len(lines)-5):
	elements=lines[line].strip().split(',')
	if int(line) < 3:
		LS_bound.append(int(elements[1].strip().split(':')[0]))
		Run_bound.append(int(elements[0].strip().split(':')[0]))
		Fill_bound.append(int(elements[0].strip().split(':')[1]))

	lumiBinner=lumiBinner+float(elements[6])
	lumiDelBlock=lumiDelBlock+float(elements[5])
	if lumiBinner > 10.0: # update total lumi number in pb (to replace 2318) and dataset you want to have (to replace 100)
		LS_bound.append(int(elements[1].strip().split(':')[0]))
		Run_bound.append(int(elements[0].strip().split(':')[0]))
		Fill_bound.append(int(elements[0].strip().split(':')[1]))
		del_lumi.append(round(lumiDelBlock,3))
		rec_lumi.append(round(lumiBinner,3))
		lumiBinner=0
		lumiDelBlock=0
print "LS_bound["+str(len(LS_bound))+"]="
print LS_bound
print "Run_bound["+str(len(Run_bound))+"]="
print Run_bound
print "Fill_bound["+str(len(Fill_bound))+"]="
print Fill_bound
print "del_Lumi["+str(len(del_lumi))+"]="
print del_lumi
print "rec_Lumi["+str(len(rec_lumi))+"]="
print rec_lumi
