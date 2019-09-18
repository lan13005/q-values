import subprocess
import multiprocessing

# Should we output lots of text to slow things down?
verbose = "false"

# How many nearest neighbors to consider?
kDim=200

# USING RNG TO SELECT HISTOGRAMS TO WRITE OUT
# need a shift in the seed since srand(0) = srand(1)
seedShift=123125
# How much events should I write out histograms for?
numberEventsToSave=5

# How much events should we show?
override_nentries="true"
nentries=1000

# how much processes to spawn?
nProcess=4

subprocess.call(['rm','histograms/*'])
subprocess.call(['rm','logs/*'])


def startProcess(iProcess):
	subprocess.call(["root","-l","-b","-q","run_parallel.C("+str(iProcess)+","+str(kDim)+","+str(numberEventsToSave)+","+str(nProcess)+","+str(seedShift)+","+str(nentries)+","+override_nentries+","+verbose+")"])

if __name__ == "__main__":
	jobs=[]
	for iProcess in range(nProcess):
		p=multiprocessing.Process(target=startProcess, args=(iProcess,))
		jobs.append(p)
		p.start()
