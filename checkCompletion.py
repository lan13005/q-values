import subprocess
from termcolor import colored
import os

def checkCompletion(exit_codes):
    # If we use inside run.py we will have a list of the all the exit codes for the spawned processes
    # We can cross check the exit codes to see if it exited properly
    # If we use it outside run.py we wont have exit codes so we will do a search for the log files where we will skip the exit code check
    if len(exit_codes)!=0:
        list_codes = exit_codes
    else:
        listDir="/d/grid13/ln16/q-values-2/logs/all"
        dirs = os.listdir(listDir)
        list_codes=[]
        for folder in dirs:
            if folder[:3]=="out":
                list_codes.append(-1)

    successful_exits=True
    skip_checks=True # if any of the codes are -1 then we skip all exit code comparisons
    print(colored("--------------------------------------------------","blue"))
    print(colored("process# | exit_code(0=good, -1=skip) | successful?","blue"))
    print(colored("--------------------------------------------------","blue"))
    for icode,list_code in enumerate(list_codes):
        manualCheckStatus=subprocess.check_output("tail -n 1 /d/grid13/ln16/q-values-2/logs/all/out"+str(icode)+".txt",shell=True)[:8]=='nentries'
        if manualCheckStatus:
            manualCheckStatus="success"
        else:
            manualCheckStatus="failed"
        print("process"+str(icode)+" | "+str(list_code)+" | "+manualCheckStatus)
        if (list_code != 0) and (manualCheckStatus != "success"):
            successful_exits = False
        if (list_code == -1) and (manualCheckStatus != "success"):
            skip_checks=False

    if len(exit_codes)!=0:
        if not successful_exits:
            print(colored("Atleast one proccess terminated without success. Exiting program...","red"))
            return 1
        else:
            print(colored("All proccess terminated successfully","green"))
            return 0
    else:
        if skip_checks:
            print(colored("Skipped exit code comparison but all proccess seem to terminate successfully","green"))
            return 0
        else:
            print(colored("Skipped exit code comparison but atleast one proccess terminated without success","red"))
            return 0
