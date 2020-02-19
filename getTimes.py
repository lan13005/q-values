import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess

nentries = np.array([5000,10000,20000,30000, 40000]#,20000,400000,80000,160000,320000,640000])
times = []

for i in range(len(nentries)):
	exchangeVar=["sed","-i","s@^nentries=.*@nentries="+str(nentries[i])+"@g","test.py"]
	print " ".join(exchangeVar)
	subprocess.Popen(exchangeVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
	subprocess.Popen("python run.py",shell=True).wait()
	#subprocess.Popen("mv log_split.txt log_split.txt",shell=True).wait()
	#subprocess.Popen("mv log_bcal.txt log_bcal.txt",shell=True).wait()
	#subprocess.Popen("mv log_fcal.txt log_fcal.txt",shell=True).wait()

	timeString_split = os.popen("tail log_split.txt -n1").read()
	timeString_bcal = os.popen("tail log_bcal.txt -n1").read()
	timeString_fcal = os.popen("tail log_fcal.txt -n1").read()
	time_split = int(float(timeString_split.split(" ")[1].rstrip().lstrip()))
	time_bcal = int(float(timeString_bcal.split(" ")[1].rstrip().lstrip()))
	time_fcal = int(float(timeString_fcal.split(" ")[1].rstrip().lstrip()))
	totalTime = time_split+time_bcal+time_fcal
	times.append(totalTime)

times = np.array(times)
fig = plt.figure(figsize=(12,6))
plt.plot(nentries,times)




print(totalTime)
