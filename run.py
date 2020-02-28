#!/usr/bin/python 

import subprocess
import os
import time
from itertools import combinations


# Remaking some of the old directories so we dont mix old results with new results
subprocess.Popen("rm -r logs", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
subprocess.Popen("mkdir logs", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

runOverAll=False
if not runOverAll:
	subprocess.Popen("rm -r histograms",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	subprocess.Popen("rm -r diagnosticPlots",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	subprocess.Popen("mkdir -p histograms/fcal",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	subprocess.Popen("mkdir -p histograms/bcal",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	subprocess.Popen("mkdir -p histograms/split",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	subprocess.Popen("mkdir -p diagnosticPlots/fcal",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	subprocess.Popen("mkdir -p diagnosticPlots/bcal",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	subprocess.Popen("mkdir -p diagnosticPlots/split",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()

start_time = time.time()

kDim=300
numberEventsToSavePerProcess=2
nProcess=2
seedShift=1212
nentries=1000
override_nentries=1
verbose=0
detector="split"

subprocess.Popen("rm postQ_"+detector+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
subprocess.Popen("rm postQValHists_"+detector+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
subprocess.Popen("rm qvalResults_"+detector+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

# so we need to add single quotes which will include the double quotes we need when passing it as an argument to the main program. If we include double quotes here it will actually be included in th parsing of the text in the program
varStringBase='cosTheta_X_cms;cosTheta_eta_gjs;phi_eta_gjs;phi_X_relativeToBeamPol'#;phi_X_cms;cosThetaHighestEphotonIneta_gjs;cosThetaHighestEphotonInpi0_cms;vanHove_omegas'
#varVec=np.array(varStringBase.rstrip().split(";"))
varVec=varStringBase.rstrip().split(";")

def runOverCombo(combo,nentries):
	""" 
	This function takes in an arguement like (3,) to use the 3rd variable as the distance metric. (1,2) would be using the first 2
	"""
	tagVec=["0" for i in range(len(varVec))]
	for ele in combo:
	    tagVec[ele]="1"
	tag="".join(tagVec)
	selectedVar = [varVec[ele] for ele in combo]
	varString=";".join(selectedVar)
	numVar=len(varString.split(";"))
	
	# the distance calculation needs dim to know the dimension of phase space. Before we compile it the script needs to know so we have replace before compilation 
	exchangeVar=["sed","-i","s@const int dim=.*;@const int dim="+str(numVar)+";@g","main.h"]
	compileMain=["g++","-o","main","main.C"]
	rootFlags = subprocess.check_output(["root-config","--cflags","--glibs", "--libs"])
	rootFlags = rootFlags.rstrip().split(" ")
	
	rooFitFlags = ["-lRooStats","-lRooFitCore", "-lRooFit"]
	compileMain.extend(rootFlags)
	compileMain.extend(rooFitFlags)
	print "\nStarting new compilation\n----------------------"
	print " ".join(exchangeVar)
	print " ".join(compileMain)
	
	# need to set the detector to use in main.h and makeDiagnosticPlots 
	replaceNumProcess=["sed","-i",'s@detector=".*";@detector="'+detector+'";@g',"main.h"]
	subprocess.Popen(replaceNumProcess, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
	replaceNumProcess=["sed","-i",'s@detector=".*";@detector="'+detector+'";@g',"makeDiagnosticHists.C"]
	subprocess.Popen(replaceNumProcess, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

	
	# most processes shoudl have a wait but for some it doesnt matter. i.e. we have to wait for exchangeVar to run before compileMain
	subprocess.Popen(exchangeVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
	out, err = subprocess.Popen(compileMain, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate(); print out

	# This will replace the dimNum so we can programaticaly scan through variable sets
	#replaceVar=["sed","-i","s@dim="+str(numVar)+"@dim=dimNum@g","main.h"]
	#print " ".join(replaceVar)
	#subprocess.Popen(replaceVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

	subprocess.Popen("rm diagnostic_logs.txt", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
	
	print './main "$kDim" $varString "$numberEventsToSavePerProcess" "$nProcess" "$seedShift" "$nentries" "$override_nentries" "$verbose" &'
	print 'Number of threads: '+str(nProcess)
	allProcess=[]
	executeMain=["./main",str(kDim),varString,str(numberEventsToSavePerProcess),str(nProcess),str(seedShift),str(nentries),str(override_nentries),str(verbose),"&"]
	print " ".join(executeMain)
	subprocess.Popen(executeMain).wait()
	    
	   
	# ------------------------------------
	# run the makeDiagnosticHists program
	# ------------------------------------
	subprocess.Popen("cat logs/process* > diagnostic_logs.txt",shell=True).wait()
	subprocess.Popen("rm -f qvalResults_"+detector+".root",shell=True).wait()
	subprocess.Popen("hadd qvalResults_"+detector+".root logs/results*",shell=True).wait()
	
	if not override_nentries:
	    nentries=int(subprocess.Popen("grep nentries fitResults/etaFitNoAccSub_"+detector+".txt | cut -d' ' -f2", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].rstrip())
	subprocess.Popen("root -l -b -q makeDiagnosticHists.C",shell=True).wait()

	subprocess.Popen("sendmail lng1492@gmail.com < defaultEmail.txt",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
	# ------------------------------------
	# ------------------------------------


numVar=len(varVec)
# We are going pass as arugment a list of lists known as combo. This combo list contains all the lists of combos with numVar elements from the list varVec. If we use the command comboinations(range(3),2) we would get something like [ [1,2], [2,3], [1,3] ]. We can use these as indicies to index a a string of 0's to fill in whether a variable will be in use. i.e. if [1,3] is chosen then the string would be 101 with the second var turnedo off. This is useful when we are doing a scan of which variables we should use. Bruteforce style. 
runOverCombo((0, 1, 2),nentries)#, 3, 4, 5, 6),nentries)
counter=0
for numVar in range(1,len((varVec))+1):
    combos=combinations(range(len(varVec)),numVar)
    for combo in combos:
        counter+=1
        #print combo
        #if counter%4==0:
        #    continue
        #print combo
        #runOverCombo(combo,nentries)
    

print("--- %s seconds ---" % (time.time() - start_time))







