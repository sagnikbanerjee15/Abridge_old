from collections import Counter 
import pprint

def find_ngrams(input_list, n):
    return Counter(input_list[idx : idx + n] for idx in range(len(input_list) - 1)) 
  


input_filename = "/90daydata/maizegdb/sagnik/ABRIDGE/developing_abridge/SRR13711355_0_SE_Aligned.sortedByCoord.out.sam.qual.rle"
fhr=open(input_filename,"r")
for line in fhr:
    pprint.pprint(str(dict(find_ngrams(line.strip(), 3))))
fhr.close()