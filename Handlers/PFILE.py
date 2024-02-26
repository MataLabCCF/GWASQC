def getIDFromPSAM(inputFile):
    listID = []
    PSAM = open(inputFile)

    header = True
    for line in PSAM:
        if header:
            header = False
        else:
            split = line.strip().split()
            listID.append(split[0])

    return listID

def removeSamples(fileToRemove, inputFile, outputName, outputFolder, step, plink2):
    outName = f"{outputFolder}/{outputName}_{step}"
    command = f"{plink2} --pfile {inputFile} --remove {fileToRemove} --make-pfile --out {outName}"
    return command, outName

def removeVariants(fileToRemove, inputFile, outputName, outputFolder, step, plink2):
    outName = f"{outputFolder}/{outputName}_{step}"
    command = f"{plink2} --pfile {inputFile} --exclude {fileToRemove} --make-pfile --out {outName}"
    return command, outName

def extractVariants(fileToExtract, inputFile, outputName, outputFolder, step, plink2):
    outName = f"{outputFolder}/{outputName}_{step}"
    command = f"{plink2} --pfile {inputFile} --extract {fileToExtract} --make-pfile --out {outName}"
    return command, outName

def countPSAM(psamFile):
    file = open(psamFile)

    count = 0
    header = True
    for line in file:
        if header:
            if "#FID" in line or "IID" in line:
                header = False
        else:
            count = count+1
    return count

def countPVAR(pvarFile):
    file = open(pvarFile)

    count = 0
    dictVar = {}
    header = True
    for line in file:
        if header:
            if "#CHROM" in line:
                header = False
        else:
            split = line.split()
            count = count+1
            if split[0] not in dictVar:
                dictVar[split[0]] = 0
            dictVar[split[0]] = dictVar[split[0]]+1

    return dictVar, count