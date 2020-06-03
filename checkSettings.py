import re
from termcolor import colored

def searchFile(inputFiles, keywords, findFirst):
    for inputFile in inputFiles:
        foundArray = [False]*len(keywords)
        print("----------------")
        print("Outputting settings in "+inputFile+"\n----------------")
        with open(inputFile, 'r') as f:
            lines = f.readlines()
            for lineNum, line in enumerate(lines):
                line = line.lstrip()
                for keyNum, keyword in enumerate(keywords):
                    if (foundArray[keyNum]):
                        continue
                    if re.search(r"{}".format(keyword), line):
                        print("({0}){1}".format(lineNum,line)),
                        if findFirst:
                            foundArray[keyNum] = True
                        break



def check():
    print(colored("\n\n(LineNumber)LineContent\n","green"))
    print(colored("---------\nGeneral Settings\n---------","red"))
    keywords = ["kDim", "numberEventsToSavePerProcess", "seedShift", "nProcess", "nentries" \
            ,"override_nentries", "verbose", "makeGraphs","emailWhenFinished","weightingScheme"
            ,"runFullFit", "discrimVar", "sideBandVar", "standardizationType" 
            ]
    keywords = ["^"+keyword for keyword in keywords]
    searchFile(["run.py"],keywords,True)
    
    
    print("\n\n")
    print(colored("---------\nBin ranges should be the same when extracting initialization parameters using \
    getInitParams and when using them to get q-values in main.h. Fit range probably should also? \n---------","red"))
    keywords=["binRange","fitRange"]
    keywords = ["^"+keyword for keyword in keywords]
    files = ["main.h", "getInitParams.C"]
    searchFile(files, keywords,False)
    
    print("\n\n")
    shouldIContinue = raw_input("continue? (y/n): ")
    if shouldIContinue == "y":
        print("Settings are OK'd. Continuing...")
    elif shouldIContinue == "n":
        print("Settings are not OK. Exiting...")
        exit(0)
    else:
        print("Not valid answer. Exiting...")
        exit(0)
