import gzip

def getIDFromVCF(inputFile, bgziped):
    listID = []
    if bgziped:
        VCF = gzip.open(inputFile)
    else:
        VCF = open(inputFile)

    for line in VCF:
        if bgziped:
            line = line.decode("utf-8")
        if line.startswith("#CHROM"):
            split = line.strip().split()
            for i in range(9, len(split)):
                listID.append(split[i])
            break
    return listID