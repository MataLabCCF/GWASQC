import os

def createFolder(folderPath, log = True):
    if not os.path.exists(folderPath):
        os.mkdir(folderPath)
    if log:
        logFile = open(f"{folderPath}/QC.log", "w")
    return logFile

def execute(command, logFile):
    print(f"{command}\n")
    logFile.write(f"{command}\n")
    os.system(command)
    #input("PRESS ENTER")