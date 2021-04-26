# this program continually queries the log files to see the progress of the main program

import subprocess
import time
import argparse

with open("run.py","r") as cfg:
    startRootFileLoc=0
    endRootFileLoc=0
    inContext=False
    for iline, line in enumerate(cfg):
        if line[:13]=="_SET_nProcess":
            nProcess=int(line.split("=")[1].split("#")[0].rstrip())
        #if line[:12]=="rootFileLocs":
        #    startRootFileLoc = iline
        #    inContext=True
        #if line[:9]=="        ]":
        #    endRootFileLoc = iline
        #    inContext=False
        #if inContext:
        #    if line.lstrip()[:1]=="(":
        #        tag=line.split(",")[-1].lstrip().rstrip()[:-1][1:-1]


parser=argparse.ArgumentParser("look in log files in subfolder given by tag to determine percent complete")
parser.add_argument("-t",help="folder tag to look for log files")
args=parser.parse_args()
tag=args.t
print("LOOKING IN SUBDIR: "+tag)

folder="logs/"+tag
nentriesPerProc=0
proc_status=[False]*nProcess
for iproc in range(nProcess):
    if not proc_status[iproc]:
        p=subprocess.Popen(["tail", "-n20",folder+"/"+"out"+str(iproc)+".txt"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output=p.communicate()
        allLines=output[0].split("\n") 
        newestLine=""
        while newestLine=="":
            for line in allLines:
                if line[:8]=="Starting":
                    newestLine=line
                if line[:9]=="nentries:":
                    proc_status[iproc]=True 
        time.sleep(0.2)
        current, end = newestLine.split(" ")[2].rstrip().lstrip().split("/")
        current, end = int(current), int(end)
        if iproc==0:
            nentriesPerProc=end
            perc_complete = int(100.0*current/end)
        else:
            perc_complete = int(100.0*(current-(end-nentriesPerProc))/nentriesPerProc)
        if proc_status[iproc]:
            print("process {:d} is complete".format(iproc))
        else:
            print("process {:d} percent completed: ({:d}/{:d})  ---- {:d}%".format(iproc,current,end,perc_complete))
    else:
        print("process {:d} is complete".format(iproc))

