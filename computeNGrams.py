import pprint

def find_ngrams(input_list, n):
  return list(zip(*[input_list[i:] for i in range(n)]))


input_filename = "/90daydata/maizegdb/sagnik/ABRIDGE/developing_abridge/SRR13711355_0_SE_Aligned.sortedByCoord.out.sam.qual.rle"
fhr=open(input_filename,"r")
for line in fhr:
    pprint.pprint(find_ngrams(list(line.strip()), 3))
fhr.close()