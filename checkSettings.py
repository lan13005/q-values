import re
from termcolor import colored

def searchFile(inputFiles, keywords, findFirst, lstripLine):
    '''
    inputFile: The file to search in
    keywords: a list of key words to look for, outputs all occurences
    findFirst: should we only look for the first occurence and dont print any other?
    lstripLine: should we lstrip the line before searching? There could be tabs or white space that we dont care about at the front of the line
    '''
    for inputFile in inputFiles:
        foundArray = [False]*len(keywords)
        print("----------------")
        print("Outputting settings in "+inputFile+"\n----------------")
        with open(inputFile, 'r') as f:
            lines = f.readlines()
            for lineNum, line in enumerate(lines):
                if lstripLine:
                    line = line.lstrip()
                for keyNum, keyword in enumerate(keywords):
                    if (foundArray[keyNum]):
                        continue
                    if re.search(r"{}".format(keyword), line):
                        print("({0}){1}".format(lineNum,line)),
                        if findFirst:
                            foundArray[keyNum] = True
                        break


def askForInput():
    shouldIContinue = raw_input("\n\ncontinue? (y/n): ")
    if shouldIContinue == "y":
        print("Settings are OK'd. Continuing...")
    elif shouldIContinue == "n":
        print("Settings are not OK. Exiting...")
        exit(0)
    else:
        print("Not valid answer. Exiting...")
        exit(0)


def check():
    print(colored("\n\n(LineNumber)LineContent\n","green"))
    print(colored("---------\nGeneral Settings\n---------","red"))
    keywords = ["_SET_"]
    keywords = ["^"+keyword for keyword in keywords]
    searchFile(["run.py"],keywords,False,False)
    askForInput()
    
    print("\n\n")
    print(colored("---------\nBin ranges should be the same when extracting initialization parameters using \
    getInitParams and when using them to get q-values in main.h. Fit range probably should also? \n---------","red"))
    keywords=["double _SET_Discrim"]
    keywords = [keyword for keyword in keywords]
    files = ["main.h", "getInitParams.C"]
    searchFile(files, keywords,False,True)
    askForInput()
    
    print("\n\n")
    keywords=["int numDOF"]
    keywords = [keyword for keyword in keywords]
    keywords = ["^"+keyword for keyword in keywords]
    files = ["helperFuncs.h"]
    searchFile(files, keywords,False,False)
    print("Is the fit functions, parameters to be scaled, and parameter limits set up in helpFuncs.h?")
    askForInput()
