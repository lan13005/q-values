import re
import numpy as np
import matplotlib.pyplot as plt

startTimes = []
with open("/d/grid15/ln16/pi0eta/q-values/logs/bcal/processLog0.txt","r") as inFile:
    for line in inFile.readlines():
        if re.search(r"{}".format("^Start"), line):
            startTimes.append(int(line.split("ms")[0].split(" ")[-1]))

startTimes = np.array(startTimes)
deltaTimes = startTimes[1:]-startTimes[:-1]
