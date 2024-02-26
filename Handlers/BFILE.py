def getIDFromFAM(inputFile):
    listID = []
    FAM = open(inputFile)

    for line in FAM:
        split = line.strip().split()
        listID.append(split[1])

    return listID