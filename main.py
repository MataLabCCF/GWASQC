import os
import sys
import gzip
import argparse
import Handlers
import numpy as np

from Handlers import VCF
from Handlers import BFILE
from Handlers import PFILE
from Utils.utils import execute, createFolder
from Handlers import COVAR


def removeDuplicates(inputFile, outputName, outputFolder, plink2, listToRemove, logFile):
    command = f"{plink2} --pfile {inputFile} --out {outputFolder}/{outputName} --freq counts"
    execute(command, logFile)

    dictCount = {}

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

                if dictCount[ID] > dictCount[oldID]:
                    dictUniqueVar[foundID] = ID
    pvar.close()

    toExtract = open(f"{outputFolder}/{outputName}_toExtractDuplicate.txt", "w")
    for ID in dictUniqueVar:
        toExtract.write(f"{dictUniqueVar[ID]}\n")
    toExtract.close()

    command, inputFile = PFILE.extractVariants(f"{outputFolder}/{outputName}_toExtractDuplicate.txt",
                                             inputFile,outputName, outputFolder,"duplicate", plink2)

    execute(command, logFile)
    listToRemove = addPfiles(inputFile, listToRemove)

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
        listToRemove = addPfiles(command, listToRemove)

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
        command, inputFile= PFILE.removeSamples(f"{outputFolder}/{outputName}_sexProblem.txt", inputFile, outputName,
                                  outputFolder, "Sex-check", plink2)

        execute(command, logFile)
        listToRemove = addBfiles(f"{inputFile}", listToRemove)

        addStatsOnLog(f"{inputFile}", "Check-sex", logFile)

    return inputFile, listToRemove

def readStepsFile(stepFile):

    dictSteps = {}
    if stepFile != "":
        file = open(stepFile)
        for line in file:
            split = line.strip().split()
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
        dictSteps["HWE control"] = ""
        dictSteps["HWE case"] = ""
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
                else:
                    filePSAMWithCovar.write(f"\t{dictTable[ind][covar]}")
            filePSAMWithCovar.write(f"\n")


    filePSAMWithCovar.close()
    listToRemove = addPfiles(f"{outputFolder}/{outputName}_nonMissing", listToRemove)

    return f"{outputFolder}/{outputName}_nonMissing", listToRemove

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


def convertInputFile(inputFile, outputName, outputFolder, covarDict, plink2, listToRemove, logFile):
    print("\tBuilding the list to be used at \"--update-sex\" (see plink2 documentation)")
    indList = getIndListFromInputFile(inputFile)
    fileToUpdateSex = buildSexUpdateFile(indList, covarDict, outputFolder, outputName)

    print("\tConverting to PFILE")
    command, inputFile = getConvertCommand(inputFile, outputFolder, outputName, fileToUpdateSex, plink2)

    if command != "":
        execute(command, logFile)
        listToRemove = addPfiles(f"{inputFile}", listToRemove)

    addStatsOnLog(f"{inputFile}", "Convert to PLINK2 format", logFile)
    return f"{inputFile}", listToRemove

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


#========================================================================================


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
    optional.add_argument('-p', '--phenoName', help='Name of the phenotype column in the info file '
                                                    '(default: PHENO)', required=False, default="PHENO")
    optional.add_argument('-s', '--steps', help='File with the order of the steps (one step per line)',
                          required=False, default = "")
    optional.add_argument('-P', '--popFile', required=False,
                          help='File with two columns: Individual ID and Population ID.This approach is important for'
                               'relationship control that should be done per population')
    optional.add_argument('-S', '--savePerPop', required=False,
                          help='Save the QCed data with all pops and per population')

    programs = parser.add_argument_group("Programs arguments")
    programs.add_argument('--plink2', required=False, help='Path to Plink 2 (default: plink2)')
    programs.add_argument('--plink1', required=False, help='Path to Plink 1.9 (default: plink)')
    programs.add_argument('--NAToRA', required=False, help='Path to NAToRA (default: NAToRA.py)')
    programs.add_argument('--python', required=False, help='Python > 3 with networkX instaled '
                                                           '(default: python)')

    args = parser.parse_args()

    logFile = createFolder(args.outputFolder, f"{args.outputFolder}/{args.outputName}_QC.log")

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
    stepsToRun = {}
    stepsToRun["geno"] = 0.05
    stepsToRun["duplicate"] = 0.05

    for step in stepsToRun:
        if step == 'sex-check':
            print("Running Sex Check")
            inputFile, listToRemove = sexCheck(inputFile, covarDict, args.outputName, args.outputFolder, args.plink1,
                                               args.plink2, listToRemove, logFile)
        if step == 'ATCG':
            print("Running A/T C/G")
            inputFile, listToRemove = removeATCG(inputFile, args.outputName, args.outputFolder, args.plink2,
                                                 listToRemove, logFile)
        if step == "geno" or step == "mind":
            print(f"Running {step}")
            inputFile, listToRemove = removeMissingData(inputFile, args.outputName, args.outputFolder, args.plink2,
                                                 listToRemove, stepsToRun[step], step, logFile)
        if step == "duplicate":
            print("Running duplicate")
            inputFile, listToRemove = removeDuplicates(inputFile, args.outputName, args.outputFolder, args.plink2,
                                                 listToRemove, logFile)





