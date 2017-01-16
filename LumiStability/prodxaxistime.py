#!/bin/env python
import datetime
from datetime import timedelta

csvFile=open("Fill5395/bcidLumiDCS.csv")
lines=csvFile.readlines()

LS_bound=[]
Run_bound=[]
Fill_bound=[]
Time_bound=[]
Time_bin_length=[]
Time_bin_array=[]
del_lumi=[]
rec_lumi=[]

lumiBinner=0
lumiDelBlock=0

timeCounter=0

for line in range(2,len(lines)-5):
	elements=lines[line].strip().split(',')

        runnum = int(elements[0].strip().split(':')[0])
        fillno = int(elements[0].strip().split(':')[1])
        lumsec = int(elements[1].strip().split(':')[0])

        timemm = int(elements[2].strip().split('/')[0])
        timedd = int(elements[2].strip().split('/')[1])
        timeyy = int(elements[2].strip().split('/')[2].split(' ')[0])
        timehr = int(elements[2].strip().split('/')[2].split(' ')[1].split(':')[0])
        timemi = int(elements[2].strip().split('/')[2].split(' ')[1].split(':')[1])
        timese = int(elements[2].strip().split('/')[2].split(' ')[1].split(':')[2])

        nowtime = datetime.datetime.combine(datetime.date(timeyy, timemm, timedd), datetime.time(timehr, timemi, timese))

        lumrec = float(elements[6])
        lumdel = float(elements[5])
        deadtime = float(lumrec/lumdel)

	if int(line) < 3:
		LS_bound.append(lumsec)
		Run_bound.append(runnum)
		Fill_bound.append(fillno)
		Time_bound.append(elements[2])
		pretime = nowtime

	lumiBinner=lumiBinner+lumrec
	lumiDelBlock=lumiDelBlock+lumdel

	if lumiBinner > 10.0: # update total lumi number in pb (to replace 2318) and dataset you want to have (to replace 100)
                LS_bound.append(lumsec)
                Run_bound.append(runnum)
                Fill_bound.append(fillno)
                Time_bound.append(elements[2])
		del_lumi.append(round(lumiDelBlock,3))
		rec_lumi.append(round(lumiBinner,3))
		lumiBinner=0
		lumiDelBlock=0

#		timebin = int((nowtime - pretime).seconds)
                timebin = (nowtime - pretime).seconds / 60.
		timeCounter += timebin
		Time_bin_length.append(round(timebin,3))
		Time_bin_array.append(round(timeCounter,3))
		pretime = nowtime

print "LS_bound["+str(len(LS_bound))+"]="
print LS_bound
print "Run_bound["+str(len(Run_bound))+"]="
print Run_bound
print "Fill_bound["+str(len(Fill_bound))+"]="
print Fill_bound
print "del_Lumi["+str(len(del_lumi))+"]="
#print del_lumi
print ", ".join(format(x, "2.3f") for x in del_lumi)
print "rec_Lumi["+str(len(rec_lumi))+"]="
#print rec_lumi
print ", ".join(format(x, "2.3f") for x in rec_lumi)

print "Time_bound["+str(len(Time_bound))+"]="
print Time_bound
print "Time_bin_length["+str(len(Time_bin_length))+"]="
#print Time_bin_length
print ", ".join(format(x, "2.3f") for x in Time_bin_length)
print "Time_bin_array["+str(len(Time_bin_array))+"]="
#print Time_bin_array
print ", ".join(format(x, "2.3f") for x in Time_bin_array)
