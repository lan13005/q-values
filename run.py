import subprocess
import os
import time
from itertools import combinations
import numpy as np

subprocess.Popen("rm etaPlots/*", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
start_time = time.time()

kDim=200
numberEventsToSavePerProcess=5
nProcess=36
seedShift=1261
nentries=200000
override_nentries=0
verbose=0
# so we need to add single quotes which will include the double quotes we need when passing it as an argument to the main program. If we include double quotes here it will actually be included in th parsing of the text in the program
varStringBase='cosTheta_X_cms;phi_X_cms;cosTheta_eta_gjs;phi_eta_gjs;cosThetaHighestEphotonIneta_gjs;cosThetaHighestEphotonInpi0_cms;vanHove_omegas'
varVec=np.array(varStringBase.rstrip().split(";"))


def runOverCombo(combo,nentries):
    """ 
    This function takes in an arguement like (3,) to use the 3rd variable as the distance metric. (1,2) would be using the first 2
    """
    tagVec=["0","0","0","0","0","0","0"]
    for ele in combo:
        tagVec[ele]="1"
    tag="".join(tagVec)
    combo = [ele for ele in combo]
    varString=";".join(varVec[combo])
    #numVar=len(varString.split(";"))
    
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

    if not override_nentries:
        nentries=int(subprocess.Popen("grep nentries fitResults/etaFitNoAccSub.txt | cut -d' ' -f2", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].rstrip())
    subprocess.Popen("root -l -b -q 'makeDiagnosticHists.C("+str(nentries)+",\""+tag+"\")'",shell=True).wait()



#numVar=1
#runOverCombo((3,))
numVar=7
runOverCombo((0, 1, 2, 3, 4, 5, 6),nentries)
counter=0
for numVar in np.arange(len((varVec)))+1:
    combos=combinations(range(len(varVec)),numVar)
    for combo in combos:
        counter+=1
        #print combo
        #if counter%4==0:
        #    continue
        #print combo
        #runOverCombo(combo,nentries)
    

print("--- %s seconds ---" % (time.time() - start_time))







