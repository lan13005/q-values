def loadCfg():
    configDict={}
    with open("run.py","r") as cfgFile:
        for line in cfgFile:
            if line[:5]=="_SET_":
                newline=line.rstrip().lstrip()
                newline=newline.split("#")[0].rstrip()
                varName=newline.split("=")[0].rstrip().lstrip()
                varValu=newline.split("=")[1].rstrip().lstrip()
                try:
                    varValu=int(varValu)
                except:
                    pass
                if len(varName.split(","))==1:
                    configDict[varName] = varValu
    return configDict

print(loadCfg())
