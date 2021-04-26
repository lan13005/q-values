import subprocess
import os
import sys
import time
from itertools import combinations
from termcolor import colored
from checkCompletion import checkCompletion

start_time = time.time()


#############################################################################
###################  DEFINING ENVIRONMENT VARIABLES #########################
#############################################################################
_SET_nProcess=48 # how many processes to spawn
_SET_kDim=800 # number of neighbors
_SET_nentries=-1 # how many combos we want to run over. Set to -1 to run over all. This should be much significantly larger than kDim or we might get errors .
_SET_seedShift=1341 # in case we dont want to save the same q-value histogram we can choose another random seed
_SET_nRndRepSubset=0 # size of the random subset of potential neighbors. If 0 or > nentries then we will not consider random subsets
_SET_standardizationType="range" # what type of standardization to apply when normalizing the phase space variables 
_SET_redistributeBkgSigFits=0 # should we do the 3 different fits where there is 100% bkg, 50/50, 100% signal initilizations. 
_SET_doKRandomNeighbors=0 # should we use k random neighbors as a test instead of doing k nearest neighbors?
_SET_nBS=0 # number of times we should bootstrap the set of neighbors to calculate q-factors with. Used to extract an error on the q-factors. Set to 0 if you dont want to do BS
_SET_saveBShistsAlso=0 # should we save every bootstrapped histogram also?
_SET_weightingScheme="as" # can be {"","as"}. for no accidental weights, accidental sub.
_SET_accWeight="AccWeight" # the branch to look at to get the accidental weights
_SET_sbWeight="weightBS" # the branch to look at to get the sideband weight
_SET_uniquenessTracking="" # default is "" which will set all event counting weights to 1. Otherwise we can give it a branch to look at
_SET_varStringBase="cosTheta_eta_gj;phi_eta_gj;cosTheta_X_cm" #;cosTheta_X_cm;phi_eta_gj # what is the phase space variables to calculate distance in 
_SET_discrimVars="Mpi0;Meta" # discriminating/reference variable
_SET_mcprocessBranch="mcprocess" # if you use a sum of simulated events with a branch for the reaction processes we can output some more information
_SET_emailWhenFinished="lng1492@gmail.com" # we can send an email when the code is finished, no email sent if empty string
_SET_verbose=1 # how much information we want to output to the logs folder
_SET_runTag="" # 3 folders are outputs of this set of programs {fitResults/diagnosticPlots/histograms}. We can append a runTag to the names allowing us to run multiple q-factors at the same time
_SET_runBatch=0 # (default=0) 0=run on a single computer, 1=submit to condor for batch processing
_SET_saveMemUsage=0 #save a file for the memory usage per process
_SET_numberEventsToSavePerProcess=10 # how many histograms (root files) we want to save. -1 = Save all histograms
_SET_alwaysSaveTheseEvents="229079" # A histogram of the fit will always be saved for these semicolon separated events
_SET_saveBranchOfNeighbors=1 # Should we save a branch containing all the neighbors per event. Size increase ~ (4Bytes per int)*kDim*nentries which is around 3GB


# What file we will analyze and what tree to look for
# Also need a tag to save the data to so that we dont overwrite other runs
rootFileBase="/d/grid13/ln16/q-values-2/"
rootFileLocs=[
        # ROOT FILE LOCATION ------ ROOT TREE NAME ------NAME TAG TO SAVE FILES AND FOLDERRS UNDER
        #(rootFileBase+"degALL_a0a2a2pi1_8288_chi13_treeFlat_DSelector.root ", "degALL_a0a2a2pi1_8288_chi13_tree_flat", "all")
        #(rootFileBase+"degALL_data_2017_mEllipse_8288_treeFlat_DSelector.root", "degALL_data_2017_mEllipse_8288_tree_flat", "all")
        #(rootFileBase+"degALL_malte_mc_treeFlat_DSelector_subsetThird.root", "degALL_malte_mc_tree_flat", "all")
        #(rootFileBase+"degALL_2017_8288_chi13_treeFlat_DSelector.root", "degALL_2017_8288_chi13_tree_flat", "all")
        #(rootFileBase+"degALL_malte_mc_pol035_treeFlat_DSelector.root", "degALL_malte_mc_pol035_tree_flat", "all")
        #(rootFileBase+"degALL_a0a2a2pi1_moreD2_lessD1_treeFlat_DSelector_subsetThird.root", "degALL_a0a2a2pi1_moreD2_lessD1_tree_flat", "all")
        #(rootFileBase+"degALL_a0a2a2pi1_D2neg_treeFlat_DSelector.root", "degALL_a0a2a2pi1_D2neg_tree_flat", "all")

         (rootFileBase+"degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol000-part1.root","degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat","000")
        ,(rootFileBase+"degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol045-part1.root","degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat","045")
        ,(rootFileBase+"degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol090-part1.root","degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat","090")
        ,(rootFileBase+"degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_pol135-part1.root","degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat","135")
        ,(rootFileBase+"degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat_polAMO-part1.root","degALL_flat_2017_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_flat","AMO")
        ]
        # (rootFileBase+"degALL_flat_8288_chi13_2_tLT05_pol000.root", "degALL_flat_8288_chi13_2_tree_flat", "000")
        #,(rootFileBase+"degALL_flat_8288_chi13_2_tLT05_pol045.root", "degALL_flat_8288_chi13_2_tree_flat", "045")
        #,(rootFileBase+"degALL_flat_8288_chi13_2_tLT05_pol090.root", "degALL_flat_8288_chi13_2_tree_flat", "090")
        #,(rootFileBase+"degALL_flat_8288_chi13_2_tLT05_pol135.root", "degALL_flat_8288_chi13_2_tree_flat", "135")
        #,(rootFileBase+"degALL_flat_8288_chi13_2_tLT05_polAMO.root", "degALL_flat_8288_chi13_2_tree_flat", "AMO")
        #,]
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
    if(_SET_kDim >= _SET_nentries and _SET_nentries>0):
        print("We cannot have kDim >= nentries. Since the combo in question cannot be a neighbor to itself we will never have kDim neighbors!")
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
    Settings we need to set in helpFuncs
    '''
    sedArgs=["sed","-i",'s@rootFileLoc=".*";@rootFileLoc="'+_SET_rootFileLoc+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@rootTreeName=".*";@rootTreeName="'+_SET_rootTreeName+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@fileTag=".*";@fileTag="'+_SET_fileTag+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@runTag=".*";@runTag="'+_SET_runTag+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@weightingScheme=".*";@weightingScheme="'+_SET_weightingScheme+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_discrimVar=".*";@s_discrimVar="'+_SET_discrimVars+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_accWeight=".*";@s_accWeight="'+_SET_accWeight+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_sbWeight=".*";@s_sbWeight="'+_SET_sbWeight+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_utBranch=".*";@s_utBranch="'+_SET_uniquenessTracking+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@s_mcprocessBranch=".*";@s_mcprocessBranch="'+_SET_mcprocessBranch+'";@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@kDim=.*;@kDim='+str(_SET_kDim)+';@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()


def execFullFit(_SET_rootFileLoc, _SET_rootTreeName, Set_fileTag):
    '''
    This programs cleans appropriate folders and runs getInitParams to extract the initialization parameters for the Q-factor analysis
    '''
    # First we clean the folders 
    print("Cleaning fitResults"+_SET_runTag+" folder")
    os.system("rm -rf fitResults"+_SET_runTag+"/"+_SET_fileTag)
    os.system("rm -rf fitResults"+_SET_runTag+"/"+_SET_fileTag)
    os.system("mkdir -p fitResults"+_SET_runTag+"/"+_SET_fileTag)
    # get the initialization parameters
    proc=subprocess.Popen("root -l -b -q getInitParams.C",shell=True).wait()
      
def mergeResults():
    '''
    After running the multi process q-factors there will be a bunch of resultsX.root files
    We will hadd them all together and then merge the final result file with the input tree
    so we can have a tree that has everything in it
    '''
    print("\n\n")
    concatRootCmd="hadd -f logs"+_SET_runTag+"/"+_SET_fileTag+"/resultsMerged_"+_SET_fileTag+".root"
    for iProcess in range(_SET_nProcess):
        concatRootCmd=concatRootCmd+" logs"+_SET_runTag+"/"+_SET_fileTag+"/results"+str(iProcess)+".root"
    print(concatRootCmd)
    subprocess.Popen(concatRootCmd,shell=True).wait()
    # Finally if everything successfully exists we can merge out q-factor results with the original input tree
    print("All processes completed without error codes")
    print("Begin merging qfactor results...")
    os.system("root -l -b -q mergeQresults.C")
            

def runOverCombo(combo,_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag):
    '''
    This is the main driver section that runs runs the "main" program with appropriate settings. 
    Here we can accept a combo flag that selects the combination of varString variables to use. This
    is useful if you are trying to do a scan over all possible combinations of phase space variables
    of all set sizes. You should probably set a nentries or else you will wait a long time....
    '''
    # clean up the outputs of the "main" program
    print("Cleaning folders that are used by main")
    os.system("rm -f main")
    os.system("rm -rf logs"+_SET_runTag+"/"+_SET_fileTag)
    os.system("rm -rf histograms"+_SET_runTag+"/"+_SET_fileTag)
    os.system("mkdir -p logs"+_SET_runTag+"/"+_SET_fileTag)
    os.system("cp main.C logs/main.C")
    os.system("cp run.py logs/run.py")
    os.system("mkdir -p histograms"+_SET_runTag+"/"+_SET_fileTag)
    os.system("rm -rf memUsage")
    if _SET_saveMemUsage:
        print("Will be saving memory usuage information of each process")
        os.system("mkdir -p memUsage")
    else:
        print("Not saving memory usuage information for the process")


    # We use this setup to make it easy to loop through all combinations for phase space variables to determine which set of variables are the best. 
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
    
    # most processes shoudl have a wait but for some it doesnt matter. i.e. we have to wait for exchangeVar to run before compileMain
    subprocess.Popen(exchangeVar, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait() # we have to wait for this command to finish before compiling...
    out, err = subprocess.Popen(compileMain, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
    print(out)
    
    print('Number of threads: '+str(_SET_nProcess))
    _SET_fitLocation = "fitResults"+_SET_runTag+"/"+_SET_fileTag+"/"+_SET_fitLocationBase+"_"+_SET_fileTag+".txt"
    print("Looking for initialzation fit in: "+_SET_fitLocation)


    if _SET_nentries == -1:
        _SET_override_nentries=0
    else:
        _SET_override_nentries=1


    _SET_cwd=os.getcwd()
    if _SET_runBatch==1:
        # ----- BATCH 
        os.system("rm -rf condor"+_SET_runTag)
        for i in range(_SET_nProcess):
            os.system("mkdir -p condor"+_SET_runTag+"/job"+str(i))
        sedArgs=["sed","-i","s/queue.*/queue "+str(_SET_nProcess)+"/g","submit.main"]
        subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
        sedArgs=["sed","-i","s@Arguments      = .*@Arguments      = "+str(_SET_kDim)+" "+_SET_varString+" "+_SET_standardizationType+" "+_SET_fitLocation+" "+str(_SET_redistributeBkgSigFits)+" "+str(_SET_doKRandomNeighbors)+" "+str(_SET_numberEventsToSavePerProcess)+" $(Process) "+str(_SET_nProcess)+" "+str(_SET_seedShift)+" "+str(_SET_nentries)+" "+str(_SET_nRndRepSubset)+" "+str(_SET_nBS)+" "+str(_SET_saveBShistsAlso)+" "+str(_SET_override_nentries)+" "+str(_SET_verbose)+" "+_SET_cwd+"@","submit.main"]
        subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
        sedArgs=["sed","-i","s/counts=$((.*-1))/counts=$(("+str(_SET_nProcess)+"-1))/g","submit.sh"]
        subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
        os.system("./submit.sh") # call the submit program to launch our condor_submit program
    else:
        # ----- Dumb Multiprocessing 
        _SET_iProcess=0;
        print("Launching processes one second apart...")
        time.sleep(3)
        outLogs=[]
        errLogs=[]
        openProcesses=[]
        pids=[]
        for _SET_iProcess in range(_SET_nProcess):
            print("Launching process "+str(_SET_iProcess))
            outLog = open("logs"+_SET_runTag+"/"+_SET_fileTag+"/out"+str(_SET_iProcess)+".txt","w")
            errLog = open("logs"+_SET_runTag+"/"+_SET_fileTag+"/err"+str(_SET_iProcess)+".txt","w")
            outLogs.append(outLog)
            errLogs.append(errLog)
            executeMain=["./main",str(_SET_kDim),_SET_varString,_SET_standardizationType,_SET_fitLocation,str(_SET_redistributeBkgSigFits), str(_SET_doKRandomNeighbors), \
                    str(_SET_numberEventsToSavePerProcess),str(_SET_iProcess),str(_SET_nProcess),str(_SET_seedShift),str(_SET_nentries),str(_SET_nRndRepSubset), \
                    str(_SET_nBS),str(_SET_saveBShistsAlso),str(_SET_override_nentries),str(_SET_verbose),_SET_cwd,_SET_alwaysSaveTheseEvents, \
                    str(_SET_saveBranchOfNeighbors),"&"]
            # print out the commands with the appropriate quotations so it can be directly run if needed
            validCommand=executeMain[:-1]
            for i in range(1,len(validCommand)):
               validCommand[i]='"'+validCommand[i]+'"' 
            validCommand=" ".join(validCommand)
            print(validCommand)
            outLog.write("\n--------------------\nCMD:\n"+validCommand+"\n--------------------")
            openProcess = subprocess.Popen(executeMain,stdout=outLog,stderr=errLog)

            if _SET_saveMemUsage:
                # have to specifcy the shell, the first argument
                subprocess.Popen(["sh","./auxilliary/getMemOfProcess.sh",str(_SET_iProcess),str(openProcess.pid),"5000","1","memUsage","&"]) #process id, pid of process, max iterations, sleep time, output Directory 
            pids.append(openProcess.pid)
            openProcesses.append(openProcess)
        exit_codes = [proc.wait() for proc in openProcesses]
        #if checkCompletion(exit_codes)==1:
        #    exit()


       
def runMakeGraphs(_SET_fileTag,_SET_emailWhenFinished):
    '''
    Runs the program that aggregates all the results, from all the threads, and untangles them. Then various plots can be made.
    This untangling is mainly necessary due to the hadd function not ordering things properly 
    '''
    # clean up the files that are created by makeDiagnosticHists before we rerun it
    print("Cleaning diagnosticPlots"+_SET_runTag+" folder")
    os.system("rm -rf diagnosticPlots"+_SET_runTag+"/"+_SET_fileTag)
    os.system("mkdir -p diagnosticPlots"+_SET_runTag+"/"+_SET_fileTag)
    # ------------------------------------
    # run the makeDiagnosticHists program
    # ------------------------------------
    #subprocess.Popen("cat logs"+_SET_runTag+"/"+_SET_fileTag+"/process* > logs"+_SET_runTag+"/"+_SET_fileTag+"/diagnostic_logs.txt",shell=True).wait()
    cmd="root -l -b -q makeDiagnosticHists.C"
    print("running cmd: "+cmd)
    subprocess.Popen(cmd,shell=True).wait()


def combineAllGraphs():
    '''
    After running over all the different dataset we will just add all the histograms together and make one final plot
    postQVal_hists_TAG.root are outputs of makeDiagnosticHists.C
    postQVal_flatTree_TAG.root are outputs of mergeQresults.C
    '''
    tags = [rootFileLoc[2] for rootFileLoc in rootFileLocs]
    haddHistCmd="hadd -f diagnosticPlots"+_SET_runTag+"/postQVal_hists.root"
    haddTreeCmd="hadd -f logs"+_SET_runTag+"/postQVal_flatTree.root"
    for tag in tags:
        haddHistCmd = haddHistCmd+" diagnosticPlots"+_SET_runTag+"/"+tag+"/postQVal_hists_"+tag+".root"
        haddTreeCmd = haddTreeCmd+" logs"+_SET_runTag+"/"+tag+"/postQVal_flatTree_"+tag+".root"
    print("\n\n")
    print(haddHistCmd)
    print(haddTreeCmd)
    os.system("rm -f diagnosticPlots/postQVal_hists.root")
    os.system(haddHistCmd)
    os.system("rm -f diagnosticPlots/postQVal_flatTree.root")
    os.system(haddTreeCmd)
    os.system("root -l -b -q makeDiagnosticHists_drawSum.C")


#############################################################################
###################    BEGIN RUNNING THE PROGRAM    #########################
#############################################################################
for _SET_rootFileLoc, _SET_rootTreeName, _SET_fileTag in rootFileLocs:
    print("\n\n-------------------------------------")
    print("Starting running of {0}".format(_SET_rootFileLoc))
    print("TreeName: {0}".format(_SET_rootTreeName))
    print("FileTag: {0}".format(_SET_fileTag))

    numVar=len(varVec)
    reconfigureSettings("helperFuncs.h",_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    if _SET_runFullFit:
        execFullFit(_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    if _SET_runQFactor:
        runOverCombo(range(numVar),_SET_rootFileLoc,_SET_rootTreeName,_SET_fileTag)
    if _SET_runMakeHists:
        mergeResults()
        runMakeGraphs(_SET_fileTag,_SET_emailWhenFinished)
    if _SET_emailWhenFinished:
        print("Sending program finished email")
        subprocess.Popen("sendmail "+_SET_emailWhenFinished+" < defaultEmail.txt",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()

# Once all the datasets have been run over we can combine all the results. This does suppose the datasets should be combined...
if _SET_runMakeHists:
    combineAllGraphs()


    # Use this code block to run over all possible combinations of variables
    # We are going pass as arugment a list of lists known as combo. This combo list contains all the lists of combos with numVar elements from the list varVec. If we use the command comboinations(range(3),2) we would get something like [ [1,2], [2,3], [1,3] ]. We can use these as indicies to index a a string of 0's to fill in whether a variable will be in use. i.e. if [1,3] is chosen then the string would be 101 with the second var turnedo off. This is useful when we are doing a scan of which variables we should use. Bruteforce style. 
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







