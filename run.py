#!/usr/bin/python 

import subprocess
import os
import sys
import time
from itertools import combinations
from checkSettings import check # contain function that compares files for keywords like checking discriminating variable bin size and ranges 

start_time = time.time()

kDim=300 # number of neighbors
numberEventsToSavePerProcess=2 # how many histograms (root files) we want to save.
seedShift=1212 # in case we dont want the same q-value histogram we can choose another random seed
nProcess=16 # how many processes to spawn
nentries=5000 # how many combos we want to run over. This should be much larger than kDim or we might get errors
override_nentries=1 # A direct modification for nentries. If = 0 then nentries will not be used. if = 1 then nentries is the number of combos to run over
verbose=1 # how much information we want to output to the logs
weightingScheme="as*bs" # can be {"","as","as*bs"}. for no weights, accidental sub, both accidental and sideband. Accidental weights are passed in through the root trees, sideband weights calculated within
emailWhenFinished="lng1492@gmail.com" # we can send an email when the code is finished, no email sent if empty string
makeGraphs=False # do we want to run makeDiagnosticHists
runFullFit=False # do we want to run the full fit to extract the initialization parameters?

check() # Outputting some checks to make sure getInitParams, main.h, and makeDiagnosticHists agree

# What file we will analyze and what tree to look for
# Also need a tag to save the data to so that we dont overwrite other runs
rootFileLocs=[
        ("/d/grid15/ln16/pi0eta/q-values/degALL_bcal_treeFlat_DSelector.root", "degALL_bcal_tree_flat", "bcal")
        #,("/d/grid15/ln16/pi0eta/q-values/degALL_fcal_treeFlat_DSelector.root", "degALL_fcal_tree_flat", "fcal")
        #,("/d/grid15/ln16/pi0eta/q-values/degALL_split_treeFlat_DSelector.root", "degALL_split_tree_flat", "split")
        ]


# ---- REMEMBER TO CHANGE THE MAIN FUNCTION BELOW TO INCLUDE ALL THE VARIABLES YOU WOULD LIKE TO USE -----
varStringBase="cosTheta_X_cms;cosTheta_eta_gjs;phi_eta_gjs"#;Mpi0s;phi_X_relativeToBeamPols;phi_eta_gjs;phi_X_cms;cosThetaHighestEphotonIneta_gjs;cosThetaHighestEphotonInpi0_cms;vanHove_omegas'
varVec=varStringBase.rstrip().split(";")
# --------------------------------------------------------------------------------------------------------


def execFullFit(rootFileLoc,rootTreeName,fileTag):
    # First we clean the folders 
    os.system("rm -rf fitResults")
    os.system("mkdir fitResults")

    # Change the settings
    sedArgs=["sed","-i",'s@rootFileLoc=".*";@rootFileLoc="'+rootFileLoc+'";@g',"getInitParams.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@rootTreeName=".*";@rootTreeName="'+rootTreeName+'";@g',"getInitParams.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@fileTag=".*";@fileTag="'+fileTag+'";@g',"getInitParams.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@weightingScheme=".*";@weightingScheme="'+weightingScheme+'";@g',"getInitParams.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

    # get the initialization parameters
    subprocess.Popen("root -l -b -q getInitParams.C",shell=True).wait()
      


def runOverCombo(combo,nentries,rootFileLoc,rootTreeName,fileTag):


    # this combo is used when looping through all combinations for phase space variables to determine which set of variables are the best. 
    tagVec=["0" for i in range(len(varVec))]
    for ele in combo:
        tagVec[ele]="1"
    tag="".join(tagVec)
    selectedVar = [varVec[ele] for ele in combo]
    varString=";".join(selectedVar)
    numVar=len(varString.split(";"))
    
    # the distance calculation needs dim to know the dimension of phase space. Before we compile it the script needs to know so we have replace before compilation 
    print("Deleting main and recompiling it")
    subprocess.Popen("rm main", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    exchangeVar=["sed","-i","s@const int dim=.*;@const int dim="+str(numVar)+";@g","main.h"]
    compileMain=["g++","-o","main","main.C"]
    try:
        rootFlags = subprocess.check_output(["root-config","--cflags","--glibs", "--libs"])
    except:
        raise Exception("ROOT not loaded!")
    rootFlags = rootFlags.rstrip().split(" ")
    
    
    rooFitFlags = ["-lRooStats","-lRooFitCore", "-lRooFit"]
    compileMain.extend(rootFlags)
    compileMain.extend(rooFitFlags)
    print("\nStarting new compilation\n----------------------")
    print(" ".join(exchangeVar))
    print(" ".join(compileMain))
    
    # setting rootFile, tree, fileTag and weightingScheme in main.h
    sedArgs=["sed","-i",'s@rootFileLoc=".*";@rootFileLoc="'+rootFileLoc+'";@g',"main.h"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@rootTreeName=".*";@rootTreeName="'+rootTreeName+'";@g',"main.h"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@fileTag=".*";@fileTag="'+fileTag+'";@g',"main.h"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@weightingScheme=".*";@weightingScheme="'+weightingScheme+'";@g',"main.h"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    
    # setting rootFile, tree, fileTag and weightingScheme in makeDiagnosticHists
    sedArgs=["sed","-i",'s@rootFileLoc=".*";@rootFileLoc="'+rootFileLoc+'";@g',"makeDiagnosticHists.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@rootTreeName=".*";@rootTreeName="'+rootTreeName+'";@g',"makeDiagnosticHists.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@fileTag=".*";@fileTag="'+fileTag+'";@g',"makeDiagnosticHists.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@weightingScheme=".*";@weightingScheme="'+weightingScheme+'";@g',"makeDiagnosticHists.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    
    
    # most processes shoudl have a wait but for some it doesnt matter. i.e. we have to wait for exchangeVar to run before compileMain
    subprocess.Popen(exchangeVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
    out, err = subprocess.Popen(compileMain, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    print(out)
    
    # This will replace the dimNum so we can programaticaly scan through variable sets
    #replaceVar=["sed","-i","s@dim="+str(numVar)+"@dim=dimNum@g","main.h"]
    #print " ".join(replaceVar)
    #subprocess.Popen(replaceVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    
    subprocess.Popen("rm diagnostic_logs.txt", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    
    print('./main "$kDim" $varString "$numberEventsToSavePerProcess" "$nProcess" "$seedShift" "$nentries" "$override_nentries" "$verbose" &')
    print('Number of threads: '+str(nProcess))
    executeMain=["./main",str(kDim),varString,str(numberEventsToSavePerProcess),str(nProcess),str(seedShift),str(nentries),str(override_nentries),str(verbose),"&"]
    print(" ".join(executeMain))
    subprocess.Popen(executeMain).wait()
        
       
    # ------------------------------------
    # run the makeDiagnosticHists program
    # ------------------------------------
    if makeGraphs:
    	subprocess.Popen("cat logs/"+fileTag+"/process* > logs/"+fileTag+"/diagnostic_logs.txt",shell=True).wait()
    	subprocess.Popen("hadd diagnosticPlots/"+fileTag+"/qvalResults_"+fileTag+".root logs/"+fileTag+"/results*",shell=True).wait()
    	
    	if not override_nentries:
    	    nentries=int(subprocess.Popen("grep nentries fitResults/etaFitNoAccSub_"+fileTag+".txt | cut -d' ' -f2", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].rstrip())
    	subprocess.Popen("root -l -b -q makeDiagnosticHists.C",shell=True).wait()

    if emailWhenFinished:
        print("Sending program finished email")
        subprocess.Popen("sendmail "+emailWhenFinished+" < defaultEmail.txt",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
    # ------------------------------------
    # ------------------------------------


# Should probably clean the old directories
print("Cleaning up directories from previous runs")
os.system("rm -f main")
os.system("rm -rf logs")
os.system("rm -rf histograms")
os.system("rm -rf diagnosticPlots")


for rootFileLoc, rootTreeName, fileTag in rootFileLocs:
    print("\n\n-------------------------------------")
    print("Starting running of {0}".format(rootFileLoc))
    print("TreeName: {0}".format(rootFileLoc))
    print("FileTag: {0}".format(fileTag))

    # cleaning up some more files
    os.system("mkdir -p logs/"+fileTag)
    os.system("mkdir -p histograms/"+fileTag)
    os.system("mkdir -p diagnosticPlots/"+fileTag)
    os.system("rm -f log_"+fileTag);

    # open a file to output to
    print("Saving output to "+'log_'+fileTag+'.txt')
    sys.stdout = open('log_'+fileTag+'.txt', 'w')

    # every time we run this program we should probably clean the qvalue results stuff since something obviously went wrong to have to run it again
    subprocess.Popen("rm -f diagnosticPlots/"+fileTag+"/postQ_"+fileTag+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    subprocess.Popen("rm -f diagnosticPlots/"+fileTag+"/postQValHists_"+fileTag+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    subprocess.Popen("rm -f diagnosticPlots/"+fileTag+"/qvalResults_"+fileTag+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    
    numVar=len(varVec)
    # We are going pass as arugment a list of lists known as combo. This combo list contains all the lists of combos with numVar elements from the list varVec. If we use the command comboinations(range(3),2) we would get something like [ [1,2], [2,3], [1,3] ]. We can use these as indicies to index a a string of 0's to fill in whether a variable will be in use. i.e. if [1,3] is chosen then the string would be 101 with the second var turnedo off. This is useful when we are doing a scan of which variables we should use. Bruteforce style. 
    if runFullFit:
        execFullFit(rootFileLoc,rootTreeName,fileTag)
    runOverCombo(range(numVar),nentries,rootFileLoc,rootTreeName,fileTag)
    # Use this code block to run over all possible combinations of variables
    #counter=0
    #for numVar in range(1,len((varVec))+1):
    #    combos=combinations(range(len(varVec)),numVar)
    #    for combo in combos:
    #        counter+=1
    #        print combo
    #        if counter%4==0:
    #            continue
    #        print combo
    #        runOverCombo(combo,nentries)
    

print("--- %s seconds ---" % (time.time() - start_time))







