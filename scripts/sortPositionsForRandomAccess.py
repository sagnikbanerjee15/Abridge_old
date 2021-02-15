

def sortPositions(positions):
    sorted_positions = {}
    for eachposition in positions:
        chromosome=eachposition.split(".")[0]
        start,end = eachposition.split(":")[-1].split("-")
        start,end=int(start),int(end)
        sorted_positions[chromosome] = [start,end]
    
    sorted_positions_new={}    
    for chromosome in sorted_positions:
        sorted_positions[chromosome]=sorted(sorted_positions[chromosome],key=lambda x: x[0])
        per_chromosome_positions = []
        i=0
        while i<len(sorted_positions[chromosome]):
            
            i+=1
        
        sorted_positions_new[chromosome] = per_chromosome_positions 