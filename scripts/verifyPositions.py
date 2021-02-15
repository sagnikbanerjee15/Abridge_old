
def verifyPositions(positions,chromosome_to_length):
    for eachposition in positions:
        if eachposition.count(":")!=1 or eachposition.count("-")!=1:
            return 0,"Error in format"
        chromosome=eachposition.split(".")[0]
        if chromosome not in chromosome_to_length:
            return 0,"Chromosome not found"
        start,end = eachposition.split(":")[-1].split("-")
        start,end=int(start),int(end)
        if start>chromosome_to_length[chromosome] or end>chromosome_to_length[chromosome] or start>end:
            return 0,"Wrong chromosome positions"
    return 1,""
        