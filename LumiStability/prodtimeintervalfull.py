#!/bin/env python
import datetime
from datetime import timedelta

csvFile=open("Fill5395/bcidLumiDCS.csv")
lines=csvFile.readlines()

LS_bound = [1, 33, 67] 
Run_bound = [282842, 282842, 282842] 
Time_window=[]

ibound=0
timeRec=0
for line in range(2,len(lines)-5):
#for line in range(2,5):
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
#	deadtime = 1.

        if int(line) < 3:
                pretime = nowtime
                prerunnum = runnum
                prelumsec = lumsec

	if runnum == Run_bound[ibound] and \
	   lumsec == LS_bound[ibound]:
		timeRec += int((nowtime - pretime).seconds) * deadtime
		print timeRec
		Time_window.append(round(timeRec,3)) 
		timeRec = 0
		ibound += 1
		if ibound == len(LS_bound):
			print runnum, lumsec
			print "Time_window["+str(len(Time_window))+"]="
			print Time_window
			print ", ".join(format(x, "4.3f") for x in Time_window)


        elif runnum == prerunnum and \
             lumsec == prelumsec+1:
                timeRec += int((nowtime - pretime).seconds) * deadtime
                print int((nowtime - pretime).seconds) * deadtime
                print "Normal"

        elif runnum == prerunnum and \
             lumsec > prelumsec+1:
                timeRec += int(23) * deadtime
                print 23 * deadtime
                print "LumiSec Gap"
                print runnum, lumsec

        elif runnum >  prerunnum:
                timeRec += int(23) * deadtime
                print 23 * deadtime
                print "End of Run"
                print runnum, lumsec

        elif runnum == prerunnum and \
             lumsec > prelumsec+1:
                print "Error type A!"

        else:
                print "Error ..."

	prerunnum = runnum
#        prefillno = fillno
	prelumsec = lumsec

#	pretimemm = timemm
#	pretimedd = timedd
#	pretimeyy = timeyy
#	pretimehr = timehr
#	pretimemi = timemi
#	pretimese = timese

	pretime = nowtime
	print nowtime

print "Time_window["+str(len(Time_window))+"]="
print Time_window 
print '[%.3f]' % ', '.join(Time_window)
