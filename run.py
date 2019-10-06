import subprocess
import os
import time

start_time = time.time()

kDim=200
numberEventsToSavePerProcess=1
nProcess=30
seedShift=125
nentries=200000
override_nentries=1
verbose=0
# so we need to add single quotes which will include the double quotes we need when passing it as an argument to the main program. If we include double quotes here it will actually be included in th parsing of the text in the program
varString='cosTheta_X_cms;phi_X_cms;cosTheta_eta_gjs;phi_eta_gjs'
numVar=len(varString.split(";"))

exchangeVar=["sed","-i","s@dim=dimNum@dim="+str(numVar)+"@g","main.h"]
compileMain=["g++","-o","main","main.C"]
rootFlags = subprocess.check_output(["root-config","--cflags","--glibs"])
rootFlags = rootFlags.rstrip().split(" ")
compileMain.extend(rootFlags)
replaceVar=["sed","-i","s@dim="+str(numVar)+"@dim=dimNum@g","main.h"]
print "\nStarting new compilation\n----------------------"
print " ".join(exchangeVar)
print " ".join(compileMain)
print " ".join(replaceVar)


# most processes shoudl have a wait but for some it doesnt matter. i.e. we have to wait for exchangeVar to run before compileMain
subprocess.Popen(exchangeVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
out, err = subprocess.Popen(compileMain, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate(); print out
subprocess.Popen(replaceVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
subprocess.Popen("rm histograms/*",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
subprocess.Popen("rm logs/*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
subprocess.Popen("rm diagnostic_logs.txt", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

print './main "$iProcess" "$kDim" "$numberEventsToSavePerProcess" "$nProcess" "$seedShift" "$nentries" "$override_nentries" "$verbose" $varString &'
allProcess=[]
for iProcess in range(int(nProcess)):
    executeMain=["./main",str(iProcess),str(kDim),str(numberEventsToSavePerProcess),str(nProcess),str(seedShift),str(nentries),str(override_nentries),str(verbose),varString,"&"]
    print " ".join(executeMain)
    allProcess.append(subprocess.Popen(executeMain)) # this will just output to stdout

for iProcess in range(int(nProcess)):
    allProcess[iProcess].wait()
    
subprocess.Popen("cat logs/log* > diagnostic_logs.txt",shell=True).wait()
subprocess.Popen("rm qvalResults.root",shell=True).wait()
subprocess.Popen("hadd qvalResults.root logs/results*",shell=True).wait()
subprocess.Popen("root -l -b -q 'makeDiagnosticHists.C("+str(nentries)+")'",shell=True).wait()

print("--- %s seconds ---" % (time.time() - start_time))






