#!/usr/bin/python 

import subprocess
import os
import sys
import time
from itertools import combinations
from checkSettings import check # contain function that compares files for keywords like checking discriminating variable bin size and ranges 

start_time = time.time()

#############################################################################
###################  DEFINING ENVIRONMENT VARIABLES #########################
#############################################################################

_SET_kDim=200 # number of neighbors
_SET_numberEventsToSavePerProcess=2 # how many histograms (root files) we want to save.
_SET_seedShift=1212 # in case we dont want the same q-value histogram we can choose another random seed
_SET_nProcess=16 # how many processes to spawn
_SET_nentries=30000 # how many combos we want to run over. This should be much larger than kDim or we might get errors
_SET_override_nentries=1 # A direct modification for nentries. If = 0 then nentries will not be used. if = 1 then nentries is the number of combos to run over
_SET_verbose=1 # how much information we want to output to the logs
_SET_weightingScheme="as*bs" # can be {"","as","as*bs"}. for no weights, accidental sub, both accidental and sideband. Accidental weights are passed in through the root trees, sideband weights calculated within
_SET_varStringBase="cosTheta_X_cm;cosTheta_eta_gj;phi_eta_gj" # what is the phase space variables to calculate distance in 
_SET_discrimVar="Meta" # discriminating/reference variable
_SET_sideBandVar="Mpi0" # side band subtract on this variable
_SET_accWeight="AccWeight" # the branch to look at to get the accidental weights
_SET_standardizationType="range" # what type of standardization to apply when normalizing the phase space variables 
_SET_emailWhenFinished="lng1492@gmail.com" # we can send an email when the code is finished, no email sent if empty string
_SET_runMakeHists=True # do we want to run makeDiagnosticHists
_SET_runFullFit=False # should we fit the full distribution of the discriminating variable to extract initialization parameters for q-factors?
_SET_runQFactor=False # should we run the q-factor analysis


#check() # Outputting some checks to make sure getInitParams, main.h, and makeDiagnosticHists agree

# What file we will analyze and what tree to look for
# Also need a tag to save the data to so that we dont overwrite other runs
rootFileLocs=[
        ("/d/grid15/ln16/pi0eta/q-values/degALL_bcal_treeFlat_DSelector.root", "degALL_bcal_tree_flat", "bcal")
        #,("/d/grid15/ln16/pi0eta/q-values/degALL_fcal_treeFlat_DSelector.root", "degALL_fcal_tree_flat", "fcal")
        #,("/d/grid15/ln16/pi0eta/q-values/degALL_split_treeFlat_DSelector.root", "degALL_split_tree_flat", "split")
        ]
varVec=_SET_varStringBase.rstrip().split(";")


#############################################################################
###################  FIRST DEFINE SOME FUNCTIONS   #########################
#############################################################################

def reconfigureSettings(fileName, _SET_rootFileLoc, _SET_rootTreeName, Set_fileTag):
    '''
    Setting we need to set in getInitParams, main.h, and makeDiagnosticHists
    '''
    sedArgs=["sed","-i",'s@rootFileLoc=".*";@rootFileLoc="'+_SET_rootFileLoc+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@rootTreeName=".*";@rootTreeName="'+_SET_rootTreeName+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@fileTag=".*";@fileTag="'+_SET_fileTag+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@weightingScheme=".*";@weightingScheme="'+_SET_weightingScheme+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_discrimVar=".*";@s_discrimVar="'+_SET_discrimVar+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_sideBandVar=".*";@s_sideBandVar="'+_SET_sideBandVar+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_accWeight=".*";@s_accWeight="'+_SET_accWeight+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()


def execFullFit(_SET_rootFileLoc, _SET_rootTreeName, Set_fileTag):
    # First we clean the folders 
    print("Cleaning fitResults folder")
    os.system("rm -rf fitResults")
    os.system("mkdir fitResults")
    reconfigureSettings("getInitParams.C", _SET_rootFileLoc, _SET_rootTreeName, Set_fileTag)
    # get the initialization parameters
    proc=subprocess.Popen("root -l -b -q getInitParams.C",shell=True).wait()
      

def runOverCombo(combo,_SET_nentries,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag):
    # clean up the outputs of the "main" program
    print("Cleaning folders that are used by main")
    subprocess.Popen("rm -f diagnosticPlots/"+_SET_fileTag+"/qvalResults_"+_SET_fileTag+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    os.system("rm -f main")
    os.system("rm -rf logs/"+_SET_fileTag)
    os.system("rm -rf histograms/"+_SET_fileTag)
    os.system("mkdir -p logs/"+_SET_fileTag)
    os.system("mkdir -p histograms/"+_SET_fileTag)

    # this combo is used when looping through all combinations for phase space variables to determine which set of variables are the best. 
    tagVec=["0" for i in range(len(varVec))]
    for ele in combo:
        tagVec[ele]="1"
    tag="".join(tagVec)
    selectedVar = [varVec[ele] for ele in combo]
    _SET_varString=";".join(selectedVar)
    numVarInChosenVarStr=len(_SET_varString.split(";"))
    
    # the distance calculation needs dim to know the dimension of phase space. Before we compile it the script needs to know so we have replace before compilation 
    print("Deleting main and recompiling it")
    subprocess.Popen("rm main", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    exchangeVar=["sed","-i","s@const int dim=.*;@const int dim="+str(numVarInChosenVarStr)+";@g","main.h"]
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
    
    # reconfiguring the appropriate files to the current settings
    reconfigureSettings("main.h",_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    sedArgs=["sed","-i",'s@standardizationType=".*";@standardizationType="'+_SET_standardizationType+'";@g',"main.h"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

    # most processes shoudl have a wait but for some it doesnt matter. i.e. we have to wait for exchangeVar to run before compileMain
    subprocess.Popen(exchangeVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
    out, err = subprocess.Popen(compileMain, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    print(out)
    
    
    subprocess.Popen("rm diagnostic_logs.txt", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    
    print('./main "$kDim" $varString $discrimVar $sideBandVar "$numberEventsToSavePerProcess" "$nProcess" "$seedShift" "$nentries" "$override_nentries" "$verbose" &')
    print('Number of threads: '+str(_SET_nProcess))
    executeMain=["./main",str(_SET_kDim),_SET_varString,str(_SET_numberEventsToSavePerProcess),str(_SET_nProcess),str(_SET_seedShift),str(_SET_nentries),str(_SET_override_nentries),str(_SET_verbose),"&"]
    print(" ".join(executeMain))
    subprocess.Popen(executeMain).wait()
        
       
def runMakeGraphs(_SET_fileTag,_SET_emailWhenFinished):
    # clean up the files that are created by makeDiagnosticHists before we rerun it
    print("Cleaning diagnosticPlots folder")
    os.system("rm -rf diagnosticPlots/"+_SET_fileTag)
    os.system("mkdir -p diagnosticPlots/"+_SET_fileTag)
    #subprocess.Popen("rm -f diagnosticPlots/"+_SET_fileTag+"/postQ_"+_SET_fileTag+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    #subprocess.Popen("rm -f diagnosticPlots/"+_SET_fileTag+"/postQValHists_"+_SET_fileTag+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    # ------------------------------------
    # run the makeDiagnosticHists program
    # ------------------------------------
    reconfigureSettings("makeDiagnosticHists.C",_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    subprocess.Popen("cat logs/"+_SET_fileTag+"/process* > logs/"+_SET_fileTag+"/diagnostic_logs.txt",shell=True).wait()
    subprocess.Popen("hadd diagnosticPlots/"+_SET_fileTag+"/qvalResults_"+_SET_fileTag+".root logs/"+_SET_fileTag+"/results*",shell=True).wait()
    subprocess.Popen("root -l -b -q makeDiagnosticHists.C",shell=True).wait()



#############################################################################
###################    BEGIN RUNNING THE PROGRAM    #########################
#############################################################################



for _SET_rootFileLoc, _SET_rootTreeName, _SET_fileTag in rootFileLocs:
    print("\n\n-------------------------------------")
    print("Starting running of {0}".format(_SET_rootFileLoc))
    print("TreeName: {0}".format(_SET_rootFileLoc))
    print("FileTag: {0}".format(_SET_fileTag))

    # open a file to output to
    #print("Saving output to "+'log_'+_SET_fileTag+'.txt')
    #os.system("rm -f log_"+_SET_fileTag+".txt");
    #sys.stdout = open('log_'+_SET_fileTag+'.txt', 'w')


    numVar=len(varVec)
    # We are going pass as arugment a list of lists known as combo. This combo list contains all the lists of combos with numVar elements from the list varVec. If we use the command comboinations(range(3),2) we would get something like [ [1,2], [2,3], [1,3] ]. We can use these as indicies to index a a string of 0's to fill in whether a variable will be in use. i.e. if [1,3] is chosen then the string would be 101 with the second var turnedo off. This is useful when we are doing a scan of which variables we should use. Bruteforce style. 
    if _SET_runFullFit:
        execFullFit(_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    if _SET_runQFactor:
        runOverCombo(range(numVar),_SET_nentries,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    if _SET_runMakeHists:
        runMakeGraphs(_SET_fileTag,_SET_emailWhenFinished)
    if _SET_emailWhenFinished:
        print("Sending program finished email")
        subprocess.Popen("sendmail "+_SET_emailWhenFinished+" < defaultEmail.txt",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
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
    #        runOverCombo(combo,_SET_nentries,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    

print("--- %s seconds ---" % (time.time() - start_time))







