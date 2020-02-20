#!/apps/bin/python 

import os
import numpy as np
import subprocess
import matplotlib.pyplot as plt

nParams = np.array([4000,8000,16000,32000,64000,128000])
nameTag = "nEntries"
times = []

for i in range(len(nParams)):
	#os.system("rm -f log.txt")
	exchangeVar=["sed","-i","s@^nentries=.*@nentries="+str(nParams[i])+"@g","run.py"]
	print " ".join(exchangeVar)
	subprocess.Popen(exchangeVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
	subprocess.Popen("python run.py > log.txt",shell=True).wait()

	timeString = os.popen("tail log.txt -n1").read()
	time = int(float(timeString.split(" ")[1].rstrip().lstrip()))
	times.append(time)

times = np.array(times)
fig, ax  = plt.subplots(figsize=(12,6))
ax.plot(nParams,times)
ax.set_xlabel("Number of entries considered")
ax.set_ylabel("Q-Factor execution time")
for i in range(len(nParams)):
	coordText = "("+str(nParams[i])+","+str(times[i])+")"
	ax.annotate(coordText,(nParams[i],times[i]))

plt.savefig("qFactorScaling_"+nameTag+".pdf")

