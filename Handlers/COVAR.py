import sys


def formatTable(dictTable):
    typeSex = getTypeSex(dictTable)
    typeStatus = getTypeStatus(dictTable)

    for ind in dictTable:
        if typeSex == "01":
            if dictTable[ind]["SEX_Original"] == "0":
                dictTable[ind]["SEX"] = "1"
            elif dictTable[ind]["SEX_Original"] == "1":
                dictTable[ind]["SEX"] = "2"
        elif typeSex == "12":
            dictTable[ind]["SEX"] = dictTable[ind]["SEX_Original"]
        elif typeSex == "MaleFemale":
            if dictTable[ind]["SEX_Original"].lower() == "male":
                dictTable[ind]["SEX"] = "1"
            elif dictTable[ind]["SEX_Original"].lower() == "female":
                dictTable[ind]["SEX"] = "2"
        elif typeSex == "MF":
            if dictTable[ind]["SEX_Original"].lower() == "m":
                dictTable[ind]["SEX"] = "1"
            elif dictTable[ind]["SEX_Original"].lower() == "f":
                dictTable[ind]["SEX"] = "2"

        if typeStatus == "01":
            if dictTable[ind]["STATUS_Original"] == "0":
                dictTable[ind]["STATUS"] = "1"
            elif dictTable[ind]["STATUS_Original"] == "1":
                dictTable[ind]["STATUS"] = "2"
        elif typeStatus == "12":
            dictTable[ind]["STATUS"] = dictTable[ind]["STATUS_Original"]
        elif typeStatus == "Text":
            if isCase(dictTable[ind]["STATUS_Original"]):
                dictTable[ind]["STATUS"] = "2"
            elif isControl(dictTable[ind]["STATUS_Original"]):
                dictTable[ind]["STATUS"] = "1"
            else:
                input(f"Not recognized Status: {dictTable[ind]['STATUS_Original']}")

    return dictTable

def isCase(status):
    status = status.lower()
    if status == "case" or status == "pd case" or status == "affected":
        return True
    return False

def isControl(status):
    status = status.lower()
    if status in ["control", "pd control", "unaffected", "non-affected", "non affected", "not affected"]:
        return True
    return False


def readCovarFile(infoName):
    infoFile = open(infoName)

    dictTable = {}
    dictFields = {}

    lineCount = 0

    header = True
    for line in infoFile:
        if header:
            header = False
            headerLine = line.strip().split("\t")
            for i in range(0, len(headerLine)):
                if headerLine[i].upper() == "ID" or headerLine[i].lower() == "ind" or headerLine[
                    i].lower() == "sample_external_id":
                    dictFields["ID"] = i
                elif headerLine[i].upper() == "SEX" or headerLine[i].upper() == "GENDER":
                    dictFields["SEX"] = i
                elif headerLine[i].upper() == "PHENO" or headerLine[i].upper() == "STATUS" or headerLine[
                    i].upper() == "AFFECTION_STATUS":
                    dictFields["STATUS"] = i
                else:
                    dictFields[headerLine[i].upper()] = i

            exit = False
            if "ID" not in dictFields:
                print("Error: ID not found (expected: ID, ind, sample_external_id). The code ignore caps, so "
                      "\"ind\" and \"Ind\" is the same.")
                exit = True
            if "SEX" not in dictFields:
                print("Error: SEX not found (expected: Sex, gender). The code ignore caps, so \"sex\" and \"sEX\" "
                      "is the same. ")
                exit = True
            if "STATUS" not in dictFields:
                print("Error: STATUS not found (expected: Pheno, status, affection_status). The code ignore caps, "
                      "so \"status\" and \"STATUS\" is the same.")
                exit = True
            if exit:
                sys.exit("")
        else:
            lineCount = lineCount + 1
            data = line.split("\t")
            data[-1] = data[-1].strip()

            if dictFields["ID"] < len(data):
                ID = data[dictFields["ID"]].replace(" ", "_")

                sex = data[dictFields["SEX"]]
                status = data[dictFields["STATUS"]]

            if ID != "":
                if ID in dictTable:
                    sys.exit(f'Error: the ID {ID} is duplicated (lines: {dictTable[ID]["lineCount"]} and {lineCount}\n')

                else:
                    dictTable[ID] = {}
                    dictTable[ID]["SEX_Original"] = sex
                    dictTable[ID]["STATUS_Original"] = status
                    dictTable[ID]["lineCount"] = lineCount

                    for i in range(0, len(data)):
                        if i not in [dictFields["ID"], dictFields["SEX"], dictFields["STATUS"]]:
                            dictTable[ID][headerLine[i].upper()] = data[i]
    infoFile.close()
    dictTable = formatTable(dictTable)
    return dictTable


def getTypeSex(dictTable):
    sexType = ""

    for ind in dictTable:
        if dictTable[ind]["SEX_Original"].isnumeric():
            if dictTable[ind]["SEX_Original"] == "0":
                sexType = "01"
            elif dictTable[ind]["SEX_Original"] == "2":
                sexType = "12"
        else:
            if "male" in dictTable[ind]["SEX_Original"].lower() or "female" in dictTable[ind]["SEX_Original"].lower():
                sexType = "MaleFemale"
            if "m" in dictTable[ind]["SEX_Original"].lower() or "f" in dictTable[ind]["SEX_Original"].lower():
                sexType = "MF"

    if sexType == "":
        sys.exit("The sex was not recognized. We expect 0/1, 1/2, Male/Female, M/F")

    return sexType


def getTypeStatus(dictTable):
    statusType = ""

    for ind in dictTable:
        if dictTable[ind]["STATUS_Original"] != "NA":
            if dictTable[ind]["STATUS_Original"].isnumeric():
                if dictTable[ind]["STATUS_Original"] == "0":
                    statusType = "01"
                elif dictTable[ind]["STATUS_Original"] == "2":
                    statusType = "12"
            else:
                statusType = "Text"

    if statusType == "":
        sys.exit("The status was not recognized. We expect 0/1, 1/2")

    return statusType
def getCovarList(dictCovar):
    returnList = []
    for ind in dictCovar:
        for cov in dictCovar[ind]:
            if cov not in ["STATUS_Original", "SEX_Original", "lineCount"]:
                returnList.append(cov)

        return returnList
