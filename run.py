import subprocess
import os
import sys
import time
from itertools import combinations

start_time = time.time()


#############################################################################
###################  DEFINING ENVIRONMENT VARIABLES #########################
#############################################################################

_SET_nProcess=36 # how many processes to spawn
_SET_kDim=300 # number of neighbors
_SET_nentries=-1 # how many combos we want to run over. Set to -1 to run over all. This should be much larger than kDim or we might get errors .
_SET_numberEventsToSavePerProcess=5 # how many histograms (root files) we want to save.
_SET_seedShift=134131 # in case we dont want the same q-value histogram we can choose another random seed
_SET_nRndRepSubset=0 # size of the random subset of potential neighbors. If nRndRepSubset>nentries when override_nentries=1, the program will not use a random subset.
_SET_standardizationType="range" # what type of standardization to apply when normalizing the phase space variables 
_SET_redistributeBkgSigFits=0 # should we do the 3 different fits where there is 100% bkg, 50/50, 100% signal initilizations. Otherwise we will use the scaled fit parameters from getInitParams
_SET_doKRandomNeighbors=0 # should we use k random neighbors as a test instead of doing k nearest neighbors?
_SET_nBS=0 # number of times we should bootstrap the set of neighbors to calculate q-factors with. Used to extract an error on the q-factors. Set to 0 if you dont want to do BS
_SET_saveBShistsAlso=1 # should we save every bootstrapped histogram also?
_SET_weightingScheme="as" # can be {"","as","as*bs"}. for no weights, accidental sub, both accidental and sideband. Accidental weights and sideband weights are taken from the input tree
_SET_accWeight="AccWeight" # the branch to look at to get the accidental weights
_SET_sbWeight="weightBS" # the branch to look at to get the sideband weight
_SET_uniquenessTracking="UT_noTrackingWeights" # "default" or give a branch name. Default is the default implementation of the DSelector tracking specific particle combos. 
_SET_varStringBase="cosTheta_X_cm;cosTheta_eta_gj;phi_eta_gj;Mpi0g1;Mpi0g2" # what is the phase space variables to calculate distance in 
_SET_discrimVars="Meta;Mpi0" # discriminating/reference variable
_SET_emailWhenFinished=""#lng1492@gmail.com" # we can send an email when the code is finished, no email sent if empty string
_SET_verbose=1 # how much information we want to output to the logs folder
# What file we will analyze and what tree to look for
# Also need a tag to save the data to so that we dont overwrite other runs

rootFileBase="/d/grid15/ln16/pi0eta/q-values/"
#rootFileBase="/home/lawrence/Desktop/gluex/q-values/"
rootFileLocs=[
        # ROOT FILE LOCATION ------ ROOT TREE NAME ------NAME TAG TO SAVE FILES AND FOLDERRS UNDER
        #(rootFileBase+"degALL_bcal_treeFlat_DSelector_UTweights.root", "degALL_bcal_tree_flat", "bcal")
        #,(rootFileBase+"degALL_fcal_treeFlat_DSelector_UTweights.root", "degALL_fcal_tree_flat", "fcal")
        #,(rootFileBase+"degALL_split_treeFlat_DSelector_UTweights.root", "degALL_split_tree_flat", "split")

        (rootFileBase+"degALL_ALL_a0a2_treeFlat_DSelector_UTweights.root", "degALL_ALL_a0a2_tree_flat", "all")
        #(rootFileBase+"degALL_BCAL_a0a2_treeFlat_DSelector_UTweights.root", "degALL_BCAL_a0a2_tree_flat", "bcal")
        #(rootFileBase+"degALL_FCAL_a0a2_treeFlat_DSelector_UTweights.root", "degALL_FCAL_a0a2_tree_flat", "fcal")
        #,(rootFileBase+"degALL_SPLIT_a0a2_treeFlat_DSelector_UTweights.root", "degALL_SPLIT_a0a2_tree_flat", "split")
        ]
#_SET_fitLocationBase="var1Fit_toMain" # the location of the fit file from getInitParms will be searched for in "fitResults/"+fileTag+"/"+"_SET_fitLocationBase+"_"+fileTag+".txt"
_SET_fitLocationBase="var1Vsvar2Fit_toMain" # the location of the fit file from getInitParms will be searched for in "fitResults/"+fileTag+"/"+"_SET_fitLocationBase+"_"+fileTag+".txt"

#############################################################################
###################  DEALING WITH CMDLINE ARGS   #########################
#############################################################################
args=sys.argv
def showHelp():
    print("\n-help\n")
    print("run.py accepts only one argument. A 3 digit binary code is needed:")
    print("1. 0/1 given to run the full fit to extract initialization parameters")
    print("2. 0/1 given to run the q factor analysis")
    print("3. 0/1 given to make all the diagnostic histograms appliyng the q-factors")
    print("i.e. if we want to run the full fit, and q-factor analysis but do not make the graphs then 110 is the argument")


def parseCmdArgs():
    if(len(args)!=2):
        showHelp()
        exit()
    for arg in args[1:]:
        arg=str(arg)
        if(len(arg)==3):
            print("\n--------------")
            if arg[0]=="1":
                _SET_runFullFit=True
                print("Running Full Fit")
            elif arg[0]=="0":
                _SET_runFullFit=False
                print("Skipping Full Fit")
            if arg[1]=="1":
                _SET_runQFactor=True
                print("Running Q-Factor Analysis")
            elif arg[1]=="0":
                _SET_runQFactor=False
                print("Skipping Q Factor Analysis")
            if arg[2]=="1":
                _SET_runMakeHists=True
                print("Running Make Diagnostic Hists")
            elif arg[2]=="0":
                _SET_runMakeHists=False
                print("Skipping Make Diagnostic Hists")
            print("--------------\n")
            return _SET_runFullFit,_SET_runQFactor,_SET_runMakeHists
        else:
            showHelp()
            exit()

_SET_runFullFit,_SET_runQFactor,_SET_runMakeHists = parseCmdArgs()

#############################################################################
###################  FIRST DEFINE SOME FUNCTIONS   #########################
#############################################################################

varVec=_SET_varStringBase.rstrip().split(";")
def reconfigureSettings(fileName, _SET_rootFileLoc, _SET_rootTreeName, Set_fileTag):
    '''
    Settings we need to set in getInitParams, main.h, and makeDiagnosticHists
    '''
    sedArgs=["sed","-i",'s@rootFileLoc=".*";@rootFileLoc="'+_SET_rootFileLoc+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@rootTreeName=".*";@rootTreeName="'+_SET_rootTreeName+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@fileTag=".*";@fileTag="'+_SET_fileTag+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@weightingScheme=".*";@weightingScheme="'+_SET_weightingScheme+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_discrimVars=".*";@s_discrimVar="'+_SET_discrimVars+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_accWeight=".*";@s_accWeight="'+_SET_accWeight+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_sbWeight".*";@s_sbWeight="'+_SET_sbWeight+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_uniquenessTracking=".*";@s_uniquenessTracking="'+_SET_uniquenessTracking+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()


def execFullFit(_SET_rootFileLoc, _SET_rootTreeName, Set_fileTag):
    '''
    This programs cleans appropriate folders and runs getInitParams to extract the initialization parameters for the Q-factor analysis
    '''
    # First we clean the folders 
    print("Cleaning fitResults folder")
    os.system("rm -rf fitResults/"+_SET_fileTag)
    os.system("rm -rf fitResults/"+_SET_fileTag)
    os.system("mkdir -p fitResults/"+_SET_fileTag)
    # get the initialization parameters
    proc=subprocess.Popen("root -l -b -q getInitParams.C",shell=True).wait()
      

def runOverCombo(combo,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag):
    '''
    This is the main driver section that runs runs the "main" program with appropriate settings. 
    Here we can accept a combo flag that selects the combination of varString variables to use. This
    is useful if you are trying to do a scan over all possible combinations of phase space variables
    of all set sizes. You should probably set a nentries or else you will wait a long time....
    '''
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
    rootFlags = subprocess.check_output(["root-config","--cflags","--glibs", "--libs"])
    rootFlags = rootFlags.decode(encoding="utf-8").rstrip().split(" ") #python3 requires decoding first, which python2 doesnt. But it doesn hurt
    
    
    rooFitFlags = ["-lRooStats","-lRooFitCore", "-lRooFit"]
    compileMain.extend(rootFlags)
    compileMain.extend(rooFitFlags)
    print("\nStarting new compilation\n----------------------")
    print(" ".join(exchangeVar))
    print(" ".join(compileMain))
    
    # reconfiguring the appropriate files to the current settings
    #sedArgs=["sed","-i",'s@standardizationType=".*";@standardizationType="'+_SET_standardizationType+'";@g',"main.h"]
    #subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

    # most processes shoudl have a wait but for some it doesnt matter. i.e. we have to wait for exchangeVar to run before compileMain
    subprocess.Popen(exchangeVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
    out, err = subprocess.Popen(compileMain, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    print(out)
    
    
    subprocess.Popen("rm diagnostic_logs.txt", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    
    print('Number of threads: '+str(_SET_nProcess))
    _SET_fitLocation = "fitResults/"+_SET_fileTag+"/"+_SET_fitLocationBase+"_"+_SET_fileTag+".txt"
    print("Looking for initialzation fit in: "+_SET_fitLocation)


    if _SET_nentries == -1:
        _SET_override_nentries=0
    else:
        _SET_override_nentries=1
    _SET_iProcess=0;
    print("Launching processes one second apart...")
    time.sleep(3)
    outLogs=[]
    openProcesses=[]
    for _SET_iProcess in range(_SET_nProcess):
        print("Launching process "+str(_SET_iProcess))
        outLog = open("logs/"+_SET_fileTag+"/output"+str(_SET_iProcess)+".txt","w")
        outLogs.append(outLog)
        executeMain=["./main",str(_SET_kDim),_SET_varString,_SET_standardizationType,_SET_fitLocation,str(_SET_redistributeBkgSigFits), str(_SET_doKRandomNeighbors), \
                str(_SET_numberEventsToSavePerProcess),str(_SET_iProcess),str(_SET_nProcess),str(_SET_seedShift),str(_SET_nentries),str(_SET_nRndRepSubset),str(_SET_nBS),str(_SET_saveBShistsAlso),str(_SET_override_nentries),str(_SET_verbose),"&"]
        print(" ".join(executeMain))
        openProcess = subprocess.Popen(executeMain,stdout=outLog)
        openProcesses.append(openProcess)
    exit_codes = [proc.wait() for proc in openProcesses]

        
       
def runMakeGraphs(_SET_fileTag,_SET_emailWhenFinished):
    '''
    Runs the program that aggregates all the results, from all the threads, and untangles them. Then various plots can be made.
    This untangling is mainly necessary due to the hadd function not ordering things properly 
    '''
    # clean up the files that are created by makeDiagnosticHists before we rerun it
    print("Cleaning diagnosticPlots folder")
    os.system("rm -rf diagnosticPlots/"+_SET_fileTag)
    #os.system("rm -rf diagnosticPlots/"+_SET_fileTag)
    os.system("mkdir -p diagnosticPlots/"+_SET_fileTag)
    #subprocess.Popen("rm -f diagnosticPlots/"+_SET_fileTag+"/postQ_"+_SET_fileTag+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    #subprocess.Popen("rm -f diagnosticPlots/"+_SET_fileTag+"/postQValHists_"+_SET_fileTag+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    # ------------------------------------
    # run the makeDiagnosticHists program
    # ------------------------------------
    subprocess.Popen("cat logs/"+_SET_fileTag+"/process* > logs/"+_SET_fileTag+"/diagnostic_logs.txt",shell=True).wait()
    subprocess.Popen("hadd diagnosticPlots/"+_SET_fileTag+"/qvalResults_"+_SET_fileTag+".root logs/"+_SET_fileTag+"/results*",shell=True).wait()
    subprocess.Popen("root -l -b -q makeDiagnosticHists.C",shell=True).wait()


def combineAllGraphs():
    '''
    After running over all the different dataset we will just add all the histograms together and make one final plot
    '''
    tags = [rootFileLoc[2] for rootFileLoc in rootFileLocs]
    haddHistCmd="hadd diagnosticPlots/postQVal.root"
    haddTreeCmd="hadd diagnosticPlots/postQVal_flatTree.root"
    for tag in tags:
        haddHistCmd = haddHistCmd+" diagnosticPlots/"+tag+"/postQValHists_"+tag+".root"
        #haddTreeCmd = haddTreeCmd+" diagnosticPlots/"+tag+"/postQValTrees_"+tag+".root"
        haddTreeCmd = haddTreeCmd+" diagnosticPlots/"+tag+"/postQ_"+tag+"_flatTree.root"
    print("\n\n")
    print(haddHistCmd)
    print(haddTreeCmd)
    sedArgs=["sed","-i",'s@haddHistCmd=".*";@haddHistCmd="'+haddHistCmd+'";@g',"makeDiagnosticHists_drawSum.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@haddTreeCmd=".*";@haddTreeCmd="'+haddTreeCmd+'";@g',"makeDiagnosticHists_drawSum.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    os.system("root -l -b -q makeDiagnosticHists_drawSum.C")



#############################################################################
###################    BEGIN RUNNING THE PROGRAM    #########################
#############################################################################
for _SET_rootFileLoc, _SET_rootTreeName, _SET_fileTag in rootFileLocs:
    print("\n\n-------------------------------------")
    print("Starting running of {0}".format(_SET_rootFileLoc))
    print("TreeName: {0}".format(_SET_rootFileLoc))
    print("FileTag: {0}".format(_SET_fileTag))

    numVar=len(varVec)
    # We are going pass as arugment a list of lists known as combo. This combo list contains all the lists of combos with numVar elements from the list varVec. If we use the command comboinations(range(3),2) we would get something like [ [1,2], [2,3], [1,3] ]. We can use these as indicies to index a a string of 0's to fill in whether a variable will be in use. i.e. if [1,3] is chosen then the string would be 101 with the second var turnedo off. This is useful when we are doing a scan of which variables we should use. Bruteforce style. 
    reconfigureSettings("helperFuncs.h",_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    if _SET_runFullFit:
        execFullFit(_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    if _SET_runQFactor:
        runOverCombo(range(numVar),_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
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
    
if _SET_runMakeHists:
    combineAllGraphs()

print("--- %s seconds ---" % (time.time() - start_time))







