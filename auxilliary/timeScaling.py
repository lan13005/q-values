import subprocess
import time

all_nentries=[20000,40000,80000,160000,320000]
all_kDim=[200,400,800,1600]
all_nProcess=[8,16,32]

with open("timeScaling.txt","w") as timeScalingLog:
    timeScalingLog.write("nentries kDim nProcesses timeElapsed timeElapsePerEntry\n")
    for one_nentries in all_nentries:
        for one_kDim in all_kDim:
            for one_nProcess in all_nProcess:
                sedArgs=["sed","-i","s@_SET_nProcess=.* #@_SET_nProcess="+str(one_nProcess)+" #@g","run.py"]
                subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
                sedArgs=["sed","-i","s@_SET_kDim=.* #@_SET_kDim="+str(one_kDim)+" #@g","run.py"]
                subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
                sedArgs=["sed","-i","s@_SET_nentries=.* #@_SET_nentries="+str(one_nentries)+" #@g","run.py"]
                subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
                start_time = time.time()
                print("Waiting for next iteration of run.py to finish...")
                sedArgs=["python","run.py","111"]
                subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
                print("Execution Time: "+str(time.time()-start_time))
                avgTotalTime=0
                avgTimePerEntry=0
                for iProcess in range(one_nProcess):
                    proc = subprocess.Popen(['tail', '-n', "2", "logs/all/processLog"+str(iProcess)+".txt"], stdout=subprocess.PIPE)
                    lines = proc.stdout.readlines()
                    totalTime = lines[0].split(" ")[2].rstrip().lstrip()
                    timePerEntry = lines[1].split(" ")[3].rstrip().lstrip()
                    avgTotalTime+=int(totalTime)
                    avgTimePerEntry+=int(timePerEntry)
                avgTotalTime/=one_nProcess
                avgTimePerEntry/=one_nProcess
                timeScalingLog.write(str(one_nentries)+" "+str(one_kDim)+" "+str(one_nProcess)+" "+str(avgTotalTime)+" "+str(avgTimePerEntry)+"\n")


