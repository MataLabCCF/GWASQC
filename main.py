import os
import sys
import shutil
import argparse
import numpy as np

from Handlers import VCF
from Handlers import BFILE
from Handlers import PFILE
from Utils.utils import execute, createFolder
from Handlers import COVAR

def separatePerPop(inputFile, outputName, outputFolder, plink2, popFileName, listToRemove, logFile):
    popFile = open(popFileName)
    header = True

    dictCountry = {}
    fileCountry = {}

    for line in popFile:
        if header:
            header = False
        else:
            ID, pop = line.strip().split()

            if pop not in fileCountry:
                fileCountry[pop] = open(f"{outputFolder}/{outputName}_list{pop}.txt", "w")
                listToRemove.append(f"{outputFolder}/{outputName}_list{pop}.txt")
            fileCountry[pop].write(f"{ID}\n")

    for pop in fileCountry:
        fileCountry[pop].close()
        command, newFile = PFILE.extractSamples(f"{outputFolder}/{outputName}_list{pop}.txt", inputFile,
                                                outputName,outputFolder, f"extract_{pop}", plink2)

        execute(command, logFile)


def relationshipControl(inputFile, outputName, outputFolder, plink2, cutoff, NAToRA, python, popFileName, listToRemove, logFile):
    toNAToRA = []

    #Add POP column to calculate relationship per cohort
    if popFileName:
        popFile = open(popFileName)
        header = True

        dictCountry = {}

        for line in popFile:
            if header:
                header = False
            else:
                ID, country = line.strip().split()
                dictCountry[ID] = country

        print(f"Saving the original file {inputFile}.psam to {inputFile}_OLD.psam")
        if os.path.isfile(f"{inputFile}_OLD.psam"):
            os.replace(f"{inputFile}.psam", f"{inputFile}_OLD.psam")
        else:
            os.rename(f"{inputFile}.psam", f"{inputFile}_OLD.psam")
        listToRemove.append(f"{inputFile}_OLD.psam")

        psamOld = open(f"{inputFile}_OLD.psam")
        psamWithCountry = open(f"{inputFile}.psam", "w")

        popIDs = []

        header = True
        for line in psamOld:
            if header:
                split = line.strip().split()
                for data in split:
                    psamWithCountry.write(f"{data}\t")
                psamWithCountry.write("POP\n")
                header = False
            else:
                split = line.strip().split()
                for data in split:
                    psamWithCountry.write(f"{data}\t")
                if split[0] in dictCountry:
                    popID = dictCountry[split[0]]
                else:
                    noCountID = split[0].split("_")
                    popID = dictCountry[noCountID[0]]

                psamWithCountry.write(f"{popID}\n")

                if popID not in popIDs:
                    popIDs.append(dictCountry[split[0]])
        psamWithCountry.close()

        for pop in popIDs:
            command = (f"{plink2} --pfile {inputFile} --make-king-table --out {outputFolder}/{outputName}_{pop}_kinship "
                       f"--keep-if POP = {pop}")
            execute(command, logFile)
            toNAToRA.append(f"{outputFolder}/{outputName}_{pop}_kinship")

    #Relationship with everyone together
    else:
        command = f"{plink2} --pfile {inputFile} --make-king-table --out {outputFolder}/{outputName}_kinship"
        execute(command, logFile)
        toNAToRA.append(f"{outputFolder}/{outputName}_kinship")

    toExclude = open(f"{outputFolder}/{outputName}_NAToRA_Removal.txt", 'w')
    listToRemove.append(f"{outputFolder}/{outputName}_NAToRA_Removal.txt")
    countTotal = 0

    #Convert to NAToRA input (ID ID Kinship) and run NAToRA
    for kingInputFile in toNAToRA:
        if os.path.exists(f"{kingInputFile}.kin0"):
            kingFile = open(f"{kingInputFile}.kin0")
            inputNAToRA = open(f"{outputFolder}/{outputName}_input.txt", "w")
            listToRemove.append(f"{kingInputFile}.kin0")

            header = True
            for line in kingFile:
                if header:
                    header = False
                else:
                    ID1, ID2, NUMSNP, HETHET, IBS0, KINSHIP = line.strip().split()
                    inputNAToRA.write(f"{ID1}\t{ID2}\t{KINSHIP}\n")
            kingFile.close()
            inputNAToRA.close()

            command = (f"{python} {NAToRA} -i {outputFolder}/{outputName}_input.txt "
                       f"-o {outputFolder}/{outputName}_output -c {cutoff}")

            execute(command, logFile)
            outputNatora = open(f"{outputFolder}/{outputName}_output_toRemove.txt")
            count = 0
            for line in outputNatora:
                count = count+1
                countTotal = countTotal + 1
                toExclude.write(f"{line}")
            outputNatora.close()
            print(f"\tWe removed {count} from file {kingInputFile}.kin0")
        else:
            print(f"\tThe file {kingInputFile}.kin0 does not exists. Maybe it is a population with a single sample")

    print(f"f\tWe should {countTotal} from file {inputFile}")
    #os.replace(f"{inputFile}_OLD.psam", f"{inputFile}.psam")

    listToRemove.append(f"{outputFolder}/{outputName}_input.txt")
    listToRemove.append(f"{outputFolder}/{outputName}_output_familyList.txt")
    listToRemove.append(f"{outputFolder}/{outputName}_output_toRemove.txt")

    toExclude.close()
    if countTotal > 0:
        command, inputFile = PFILE.removeSamples(f"{outputFolder}/{outputName}_NAToRA_Removal.txt",
                                                 inputFile, outputName, outputFolder, "Relationship", plink2)

        listToRemove = addPfiles(f"{inputFile}", listToRemove)
        execute(command, logFile)

    addStatsOnLog(f"{inputFile}", "Relationship control with NAToRA", logFile)

    return inputFile, listToRemove


#=====================================================================================================================
def filterHWE(inputFile, outputName, outputFolder, plink2, listToRemove, pvalue, isControl, logFile):
    command = f"{plink2} --pfile {inputFile} --hwe {pvalue} --make-pfile "
    if isControl:
        type = "HWEControl"
        command = command + f"--out {outputFolder}/{outputName}_{type} --keep-if STATUS = control "
        execute(command, logFile)

        HWECase = open(f"{outputFolder}/{outputName}_{type}.pvar")
        varToKeep = open(f"{outputFolder}/{outputName}_toKeep_HWE.txt", "w")
        listToRemove.append(f"{outputFolder}/{outputName}_toKeep_HWE.txt")

        header = True
        for line in HWECase:
            if header:
                if "#CHROM" in line:
                    header = False
            else:
                ID = line.strip().split()[2]
                varToKeep.write(f"{ID}\n")

        HWECase.close()
        varToKeep.close()

        command, returnFile = PFILE.extractVariants(f"{outputFolder}/{outputName}_toKeep_HWE.txt", inputFile,
                                                   outputName, outputFolder, "HWEControl_AllSamples", plink2)
        listToRemove = addPfiles(returnFile, listToRemove
                                 )

    else:
        type = "HWECase"
        command = command + f"--out {outputFolder}/{outputName}_{type}"
        returnFile = f"{outputFolder}/{outputName}_{type}"

    execute(command, logFile)
    listToRemove= addPfiles(f"{outputFolder}/{outputName}_{type}", listToRemove)
    addStatsOnLog(f"{returnFile}", f"{type}", logFile)

    return f"{returnFile}", listToRemove

#=====================================================================================================================
def removeExcessHeterozygosity(inputFile, outputName, outputFolder, plink2, listToRemove, numStd, logFile):
    command = f"{plink2} --pfile {inputFile} --out {outputFolder}/{outputName} --het"
    execute(command, logFile)

    listToRemove.append(f"{outputFolder}/{outputName}.het")
    fileHet = open(f"{outputFolder}/{outputName}.het")
    header = True

    hetVect = []
    hetDict = {}

    for line in fileHet:
        if header:
            header = False
        else:
            ID, observedHom, expectedHom, observedCount, F = line.strip().split()
            het = (int(observedCount)-int(observedHom))/int(observedCount)
            hetDict[ID] = het
            hetVect.append(het)

    fileHet.close()

    sd = np.std(hetVect)
    mean = np.mean(hetVect)

    cutoffMax = mean + (numStd*sd)
    cutoffMin = mean - (numStd*sd)

    count = 0
    toRemoveHet = open(f"{outputFolder}/{outputName}_heterozygosity.txt", "w")
    listToRemove.append(f"{outputFolder}/{outputName}_heterozygosity.txt")

    for ind in hetDict:
        if hetDict[ind] < cutoffMin:
            toRemoveHet.write(f"{ind}\n")
            count = count + 1
        elif hetDict[ind] > cutoffMax:
            toRemoveHet.write(f"{ind}\n")
            count = count + 1
    toRemoveHet.close()

    print(f"We have to remove {count} samples associated to heterozygosity")
    if count > 0:
        command, inputFile = PFILE.removeSamples(f"{outputFolder}/{outputName}_heterozygosity.txt", inputFile, outputName,
                                                 outputFolder, "Heterozygosity", plink2)
        execute(command, logFile)
        listToRemove = addPfiles(f"{inputFile}", listToRemove)

        addStatsOnLog(f"{inputFile}", "Heterozygosity", logFile)

    return inputFile, listToRemove
#=====================================================================================================================
def removeDuplicates(inputFile, outputName, outputFolder, plink2, listToRemove, logFile):
    command = f"{plink2} --pfile {inputFile} --out {outputFolder}/{outputName} --freq counts"
    execute(command, logFile)

    dictCount = {}
    listToRemove.append(f"{outputFolder}/{outputName}.acount")
    countFile = open(f"{outputFolder}/{outputName}.acount")
    header = True
    for line in countFile:
        if header:
            header = False
        else:
            split = line.strip().split()
            dictCount[split[1]] = int(split[-1])
    countFile.close()

    dictUniqueVar = {}
    dictDuplicatePos = {}
    pvar = open(f"{inputFile}.pvar")
    header = True
    for line in pvar:
        if header:
            if "CHROM" in line:
                header = False
        else:
            chrom, pos, ID, ref, alt = line.strip().split()[0:5]
            uniqueID = f"{chrom}:{pos}:{ref}:{alt}"
            uniqueIDAlt = f"{chrom}:{pos}:{alt}:{ref}"

            if uniqueID not in dictUniqueVar and uniqueIDAlt not in dictUniqueVar:
                dictUniqueVar[uniqueID] = ID
            else:
                if uniqueID in dictUniqueVar:
                    oldID = dictUniqueVar[uniqueID]
                    foundID = uniqueID
                else:
                    oldID = dictUniqueVar[uniqueIDAlt]
                    foundID = uniqueIDAlt

                #Build list of duplicates
                if foundID not in dictDuplicatePos:
                    dictDuplicatePos[foundID] = 1
                dictDuplicatePos[foundID] = dictDuplicatePos[foundID] + 1

                if dictCount[ID] > dictCount[oldID]:
                    dictUniqueVar[foundID] = ID
    pvar.close()

    toExtract = open(f"{outputFolder}/{outputName}_toExtractDuplicate.txt", "w")
    listToRemove.append(f"{outputFolder}/{outputName}_toExtractDuplicate.txt")

    for ID in dictUniqueVar:
        toExtract.write(f"{dictUniqueVar[ID]}\n")
    toExtract.close()

    duplicateFile = open(f"{outputFolder}/{outputName}_duplicatePos.txt", "w")
    duplicateFile.write(f"Chrom\tPos\tRef\tAlt\tCount\n")

    for ID in dictDuplicatePos:
        data = ID.split(":")
        duplicateFile.write(f"{data[0]}\t{data[1]}\t{data[2]}\t{data[3]}\t{dictDuplicatePos[ID]}\n")
    duplicateFile.close()

    listToRemove.append(f"{outputFolder}/{outputName}_duplicatePos.txt")
    command, inputFile = PFILE.extractVariants(f"{outputFolder}/{outputName}_toExtractDuplicate.txt",
                                             inputFile,outputName, outputFolder,"duplicate", plink2)

    execute(command, logFile)
    listToRemove = addPfiles(inputFile, listToRemove)

    addStatsOnLog(f"{inputFile}", "Duplicate", logFile)

    return inputFile, listToRemove
#=======================================================================================================================
def removeMissingData(inputFile, outputName, outputFolder, plink2, listToRemove, cutoff, step, logFile):
    command = f"{plink2} --pfile {inputFile} --make-pfile --out {outputFolder}/{outputName}_{step}"
    if step == "geno":
        command = f"{command} --geno {cutoff}"
    elif step == "mind":
        command = f"{command} --mind {cutoff}"

    execute(command, logFile)
    listToRemove = addPfiles(f"{outputFolder}/{outputName}_{step}", listToRemove)

    addStatsOnLog(f"{outputFolder}/{outputName}_{step}", step, logFile)

    return f'{outputFolder}/{outputName}_{step}', listToRemove
#=======================================================================================================================
def removeATCG(inputFile, outputName, outputFolder, plink2, listToRemove, logFile):
    print(f"\tOpen file {inputFile}")
    pvar = open(f"{inputFile}.pvar")
    fileToRemove = open(f"{outputFolder}/{outputName}_listATCG.txt", "w")

    count = 0
    header = True
    for line in pvar:
        if header:
            if "#CHROM" in line:
                header = False
        else:
            varID, ref, alt = line.strip().split()[2:5]
            if ((ref == "A" and alt == "T") or (ref == "T" and alt == "A") or
                    (ref == "C" and alt == "G") or (ref == "G" and alt == "C")):
                fileToRemove.write(f"{varID}\n")
                count = count + 1

    fileToRemove.close()
    print(f"\tWe have {count} redundant variants")

    listToRemove.append(f"{outputFolder}/{outputName}_listATCG.txt")

    if count > 0:
        command, inputFile = PFILE.removeVariants(f"{outputFolder}/{outputName}_listATCG.txt", inputFile, outputName,
                                                 outputFolder, "ATCG", plink2)
        execute(command, logFile)
        listToRemove = addPfiles(inputFile, listToRemove)

    addStatsOnLog(f"{inputFile}", "Redundant variants", logFile)


    return inputFile, listToRemove
#=======================================================================================================================
def sexCheck(inputFile, covarDict, outputName, outputFolder, plink1, plink2, listToRemove, logFile):
    command = f"{plink2} --pfile {inputFile} --make-bed --out {outputFolder}/{outputName}_BFILE --output-chr 26"
    execute(command, logFile)
    listToRemove = addBfiles(f"{outputFolder}/{outputName}_BFILE", listToRemove)

    command = f"{plink1} --bfile {outputFolder}/{outputName}_BFILE  --make-bed --out {outputFolder}/{outputName}_SplitX"
    execute(command, logFile)
    listToRemove = addBfiles(f"{outputFolder}/{outputName}_SplitX", listToRemove)

    #command = (f"{plink1} --bfile {outputFolder}/{outputName}_toSplit --split-x hg38 --make-bed "
    #           f"--out {outputFolder}/{outputName}_SplitX")
    #execute(command, logFile)
    #listToRemove = addBfiles(f"{outputFolder}/{outputName}_SplitX", listToRemove)

    command = (f"{plink1} --bfile {outputFolder}/{outputName}_SplitX --check-sex "
               f"--out {outputFolder}/{outputName}_checkSex")
    execute(command, logFile)
    listToRemove.append(f"{outputFolder}/{outputName}_checkSex.sexcheck")

    sexCheckFile = open(f"{outputFolder}/{outputName}_checkSex.sexcheck")
    samplesWithProblem = open(f"{outputFolder}/{outputName}_sexProblem.txt", 'w')
    listToRemove.append(f"{outputFolder}/{outputName}_sexProblem.txt")
    count = 0

    header = True
    for line in sexCheckFile:
        if header:
            header = False
        else:
            FID, IID, declaredSex, inferredSex, status, F = line.strip().split()
            if status != "OK":
                if declaredSex == "2" and inferredSex == "0":
                    if float(F) > 0.5:
                        samplesWithProblem.write(f"{IID}\n")
                        count = count + 1
                elif declaredSex == "1" and float(F) < 0.8:
                    samplesWithProblem.write(f"{IID}\n")
                    count = count + 1
                else:
                    samplesWithProblem.write(f"{IID}\n")
                    count = count + 1

    sexCheckFile.close()
    samplesWithProblem.close()
    logFile.write(f"During the check-sex we found {count} inconsistencies\n")
    if count > 0:
        command, inputFile= PFILE.removeSamples(f"{outputFolder}/{outputName}_sexProblem.txt", inputFile,
                                                outputName, outputFolder, "Sex-check", plink2)

        execute(command, logFile)
        listToRemove = addPfiles(f"{inputFile}", listToRemove)

        addStatsOnLog(f"{inputFile}", "Check-sex", logFile)

    return inputFile, listToRemove
#=====================================================================================================================
def readStepsFile(stepFile):

    dictSteps = {}
    if stepFile != "":
        file = open(stepFile)
        for line in file:
            split = line.strip().split("\t")
            if len(split) == 1:
                dictSteps[split[0]] = ""
            elif len(split) == 2:
                dictSteps[split[0]] = split[1]
    else:
        dictSteps["sex-check"] = ""
        dictSteps["ATCG"] = ""
        dictSteps["mind"] = 0.05
        dictSteps["geno"] = 0.05
        dictSteps["duplicate"] = ""
        dictSteps["heterozygosity"] = 3
        dictSteps["HWE control"] = "1e-6"
        dictSteps["HWE case"] = "1e-10"
        dictSteps["relationship"] = 0.0884

    return dictSteps
#================================ BEGIN ================================================
def addInformationOnPSAM(inputFile, dictTable, outputName, outputFolder, plink2, listToRemove, logFile):
    print(f"Building a list to remove samples without sex and status")
    indList = PFILE.getIDFromPSAM(f"{inputFile}.psam")

    indWithoutCovar = []
    for ind in indList:
        if ind in dictTable:
            if dictTable[ind]["SEX"] not in ["1", "2"] or dictTable[ind]["STATUS"] not in ["1", "2"]:
                indWithoutCovar.append(ind)
        else:
            indNoCount = ind.split("_")[0]
            if indNoCount in dictTable:
                if dictTable[indNoCount]["SEX"] not in ["1", "2"] or dictTable[indNoCount]["STATUS"] not in ["1", "2"]:
                    indWithoutCovar.append(ind)
            else:
                indWithoutCovar.append(ind)

    print(f"Building the list in {outputFolder}/{outputName}_missingBasicCovar.txt")
    fileToExcludeMissing = open(f"{outputFolder}/{outputName}_missingBasicCovar.txt", "w")
    for ind in indWithoutCovar:
        fileToExcludeMissing.write(f"{ind}\n")
    fileToExcludeMissing.close()

    print(f"Removing the samples and saving at {outputFolder}/{outputName}_nonMissing")
    execute(f"{plink2} --pfile {inputFile} --remove {outputFolder}/{outputName}_missingBasicCovar.txt "
            f"--make-pfile --out {outputFolder}/{outputName}_nonMissing", logFile)

    listToRemove.append(f"{outputFolder}/{outputName}_missingBasicCovar.txt")
    listToRemove = addPfiles(f"{outputFolder}/{outputName}_nonMissing", listToRemove)

    print(f"Saving the original file {outputFolder}/{outputName}_nonMissing.psam to {outputFolder}/{outputName}_nonMissing_OLD.psam")
    if os.path.isfile(f"{outputFolder}/{outputName}_nonMissing_OLD.psam"):
        os.replace(f"{outputFolder}/{outputName}_nonMissing.psam", f"{outputFolder}/{outputName}_nonMissing_OLD.psam")
    else:
        os.rename(f"{outputFolder}/{outputName}_nonMissing.psam", f"{outputFolder}/{outputName}_nonMissing_OLD.psam")
    listToRemove.append(f"{outputFolder}/{outputName}_nonMissing_OLD.psam")

    filePSAMWithCovar = open(f"{outputFolder}/{outputName}_nonMissing.psam", "w")
    filePSAMWithoutCovar = open(f"{outputFolder}/{outputName}_nonMissing_OLD.psam")

    header = True
    covarIDs = COVAR.getCovarList(dictTable)
    for line in filePSAMWithoutCovar:
        if header:
            header = False
            filePSAMWithCovar.write("#IID")
            for covar in covarIDs:
                filePSAMWithCovar.write(f"\t{covar}")
            filePSAMWithCovar.write(f"\n")
        else:
            ind = line.strip().split()[0]
            filePSAMWithCovar.write(f"{ind}")
            if ind not in dictTable:
                ind = ind.split("_")[0]

            for covar in covarIDs:
                if covar not in dictTable[ind]:
                    filePSAMWithCovar.write(f"\tNA")
                elif dictTable[ind][covar] in ["", " "]:
                    filePSAMWithCovar.write(f"\tNA")
                else:
                    filePSAMWithCovar.write(f"\t{dictTable[ind][covar]}")
            filePSAMWithCovar.write(f"\n")

    filePSAMWithCovar.close()
    addStatsOnLog(f"{outputFolder}/{outputName}_nonMissing", "Removing basic covar missing", logFile)

    return f"{outputFolder}/{outputName}_nonMissing", listToRemove
#=====================================================================================================================
def buildSexUpdateFile(indList, covarDict, outputFolder, outputName):
    fileToUpdate = open(f"{outputFolder}/{outputName}_updateSex.txt", "w")

    fileToUpdate.write(f"#IID\tSEX\n")
    for ind in indList:
        if ind in covarDict:
            sex = "0"
            if "SEX" in covarDict[ind]:
                sex = covarDict[ind]["SEX"]
            fileToUpdate.write(f"{ind}\t{sex}\n")
        else:
            indNoCount = ind.split("_")[0]
            if indNoCount in covarDict:
                sex = "0"
                if "SEX" in covarDict[indNoCount]:
                    sex = covarDict[indNoCount]["SEX"]
                fileToUpdate.write(f"{ind}\t{sex}\n")
            else:
                fileToUpdate.write(f"{ind}\t0\n")
    fileToUpdate.close()

    return f"{outputFolder}/{outputName}_updateSex.txt"
#=====================================================================================================================
def convertInputFile(inputFile, outputName, outputFolder, covarDict, plink2, listToRemove, logFile):
    print("\tBuilding the list to be used at \"--update-sex\" (see plink2 documentation)")
    indList = getIndListFromInputFile(inputFile)
    fileToUpdateSex = buildSexUpdateFile(indList, covarDict, outputFolder, outputName)
    listToRemove.append(fileToUpdateSex)

    print("\tConverting to PFILE")
    command, inputFile = getConvertCommand(inputFile, outputFolder, outputName, fileToUpdateSex, plink2)

    if command != "":
        execute(command, logFile)
        listToRemove = addPfiles(f"{inputFile}", listToRemove)

    addStatsOnLog(f"{inputFile}", "Convert to PLINK2 format", logFile)
    return f"{inputFile}", listToRemove
#=====================================================================================================================
def getConvertCommand(inputFile, outputFolder, outputName, fileToUpdateSex, plink2):
    if os.path.exists(f"{inputFile}.vcf"):
        command = f"{plink2} --vcf {inputFile}.vcf --make-pfile --out {outputFolder}/{outputName} --update-sex {fileToUpdateSex} --split-par hg38"
    elif os.path.exists(f"{inputFile}.vcf.gz"):
        command = (f"{plink2} --vcf {inputFile}.vcf.gz --make-pfile --out {outputFolder}/{outputName} "
                   f"--update-sex {fileToUpdateSex} --split-par hg38")
    elif inputFile.endswith(".vcf") or inputFile.endswith(".vcf.gz"):
        command = (f"{plink2} --vcf {inputFile} --make-pfile --out {outputFolder}/{outputName} "
                   f"--update-sex {fileToUpdateSex} --split-par hg38")
    elif os.path.exists(f"{inputFile}.bed"):
        command = (f"{plink2} --bfile {inputFile} --make-pfile --out {outputFolder}/{outputName} "
                   f"--update-sex {fileToUpdateSex} --split-par hg38")
    elif os.path.exists(f"{inputFile}.pgen"):
        return "", ""

    return command, f"{outputFolder}/{outputName}"
#=====================================================================================================================
def getIndListFromInputFile(inputFile):
    if os.path.exists(f"{inputFile}.vcf"):
        return VCF.getIDFromVCF(f"{inputFile}.vcf", False)
    elif os.path.exists(f"{inputFile}.vcf.gz"):
        return VCF.getIDFromVCF(f"{inputFile}.vcf.gz", True)
    elif inputFile.endswith(".vcf"):
        return VCF.getIDFromVCF(f"{inputFile}", False)
    elif inputFile.endswith(".vcf.gz"):
        return VCF.getIDFromVCF(f"{inputFile}", True)
    elif os.path.exists(f"{inputFile}.bed"):
        return BFILE.getIDFromFAM(f"{inputFile}.fam")
    elif os.path.exists(f"{inputFile}.pgen"):
        return PFILE.getIDFromPSAM(f"{inputFile}.psam")
    else:
        logFile.write(f"The input file {inputFile} is not supported (VCF, VCF.GZ, BED/BIM/FAM or PGEN/PVAR/PSAM)\n")
        sys.exit(f"The input file {inputFile} is not supported (VCF, VCF.GZ, BED/BIM/FAM or PGEN/PVAR/PSAM)\n")
#======================================= ADDs =================================================
def addPfiles(filePrefix, listToRemove):
    listToRemove.append(f"{filePrefix}.pgen")
    listToRemove.append(f"{filePrefix}.pvar")
    listToRemove.append(f"{filePrefix}.psam")

    return listToRemove

def addBfiles(filePrefix, listToRemove):
    listToRemove.append(f"{filePrefix}.bed")
    listToRemove.append(f"{filePrefix}.bim")
    listToRemove.append(f"{filePrefix}.fam")

    return listToRemove

def addStatsOnLog(filePrefix, stepName, logFile):
    numInd = PFILE.countPSAM(f"{filePrefix}.psam")
    dictVar, numVar = PFILE.countPVAR(f"{filePrefix}.pvar")

    logFile.write(f"{stepName}:\n")
    logFile.write(f"\t#Samples: {numInd}\n")
    logFile.write(f"\t#Variants: {numVar}\n")

    for chrom in dictVar:
        logFile.write(f"\t\tchr{chrom}: {dictVar[chrom]}\n")

#=====================================================================================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Genotyping QC to GWAS')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-i', '--inputFile', help='Input file without suffix', required=True)
    required.add_argument('-o', '--outputName', help='Output name', required=True)
    required.add_argument('-O', '--outputFolder', help='Output folder', required=True)
    required.add_argument('-I', '--info', help='File with covar information. The pipeline requires at '
                                               'least three columns (ID, Sex, Phenotype)', required=True)

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument('-e', '--erase', help='Set to exclude all temporary files', required=False,
                          default= False, action='store_true')
    optional.add_argument('-s', '--steps', help='File with the order of the steps (one step per line)',
                          required=False, default = "")
    optional.add_argument('-P', '--popFile', required=False,
                          help='File with two columns: Individual ID and Population ID.This approach is important for'
                               'relationship control that should be done per population')
    optional.add_argument('-S', '--savePerPop', required=False,
                          help='Save the QCed data with all pops and per population', action = 'store_true')

    programs = parser.add_argument_group("Programs arguments")
    programs.add_argument('--plink2', required=False, help='Path to Plink 2 (default: plink2)', default="plink2")
    programs.add_argument('--plink1', required=False, help='Path to Plink 1.9 (default: plink)', default="plink")
    programs.add_argument('--NAToRA', required=False, help='Path to NAToRA (default: NAToRA.py)', default="NAToRA.py")
    programs.add_argument('--python', required=False, help='Python > 3 with networkX instaled '
                                                           '(default: python)', default="python")

    args = parser.parse_args()

    logFile = createFolder(args.outputFolder)

    listToRemove = []

    print(f"Reading covar file {args.info}")
    covarDict = COVAR.readCovarFile(args.info)

    print(f"Converting the input file {args.inputFile} to PFILE")
    inputFile, listToRemove = convertInputFile(args.inputFile, args.outputName, args.outputFolder, covarDict, args.plink2,
                                               listToRemove, logFile)

    print(f"Adding covar into the PSAM file ({inputFile})")
    inputFile, listToRemove = addInformationOnPSAM(inputFile, covarDict, args.outputName, args.outputFolder,
                                                              args.plink2, listToRemove, logFile)

    stepsToRun = readStepsFile(args.steps)

    for step in stepsToRun:
        if step == 'sex-check':
            print("Running Sex Check")
            inputFile, listToRemove = sexCheck(inputFile, covarDict, args.outputName, args.outputFolder, args.plink1,
                                               args.plink2, listToRemove, logFile)
        elif step == 'ATCG':
            print("Running A/T C/G")
            inputFile, listToRemove = removeATCG(inputFile, args.outputName, args.outputFolder, args.plink2,
                                                 listToRemove, logFile)
        elif step == "geno" or step == "mind":
            print(f"Running {step}")
            inputFile, listToRemove = removeMissingData(inputFile, args.outputName, args.outputFolder, args.plink2,
                                                 listToRemove, stepsToRun[step], step, logFile)
        elif step == "duplicate":
            print("Running duplicate")
            inputFile, listToRemove = removeDuplicates(inputFile, args.outputName, args.outputFolder, args.plink2,
                                                 listToRemove, logFile)
        elif step == "heterozygosity":
            print("Running excess heterozygosity removal")
            inputFile, listToRemove = removeExcessHeterozygosity(inputFile, args.outputName, args.outputFolder, args.plink2,
                                                 listToRemove, int(stepsToRun[step]), logFile)
        elif step == "HWE control":
            print("Removing variants that failed on HWE exact test in controls")
            inputFile, listToRemove = filterHWE(inputFile, args.outputName, args.outputFolder, args.plink2,listToRemove,
                                                stepsToRun[step], True, logFile)
        elif step == "HWE case":
            print("Removing variants that failed on HWE exact test")
            inputFile, listToRemove = filterHWE(inputFile, args.outputName, args.outputFolder, args.plink2,listToRemove,
                                                stepsToRun[step], False, logFile)
        elif step == "relationship":
            print("Relationship control with NAToRA")
            inputFile, listToRemove = relationshipControl(inputFile, args.outputName, args.outputFolder, args.plink2,
                                                          stepsToRun[step], args.NAToRA, args.python, args.popFile,
                                                          listToRemove, logFile)

    if not os.path.exists(f"{args.outputFolder}/FinalData"):
        os.mkdir(f"{args.outputFolder}/FinalData")
    shutil.copyfile(f"{inputFile}.pgen", f"{args.outputFolder}/FinalData/{args.outputName}_QCed.pgen")
    shutil.copyfile(f"{inputFile}.psam", f"{args.outputFolder}/FinalData/{args.outputName}_QCed.psam")
    shutil.copyfile(f"{inputFile}.pvar", f"{args.outputFolder}/FinalData/{args.outputName}_QCed.pvar")

    if args.savePerPop:
        separatePerPop(inputFile, args.outputName, f"{args.outputFolder}/FinalData/", args.plink2, args.popFile, listToRemove, logFile)


    if args.erase:
        for file in listToRemove:
            print(f"Excluding {file}")
            os.remove(file)

    logFile.close()

