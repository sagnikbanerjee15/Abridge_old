from collections import Counter 
import pprint

def find_ngrams(input_list, n):
    res = Counter(input_list[idx : idx + n] for idx in range(len(input_list) - 1))
    return sorted(res.items(),key=lambda item: (-item[1], item[0])) 
  


input_filename = "/90daydata/maizegdb/sagnik/ABRIDGE/developing_abridge/SRR13711355_0_SE_Aligned.sortedByCoord.out.sam.qual.rle.line3"
fhr=open(input_filename,"r")
for line_number,line in enumerate(fhr):
    print(f"Unigram {line_number+1}")
    print(str(dict(find_ngrams(line.strip(), 1))))
    print("="*100)
    print(f"Bigram {line_number+1}")
    print(str(dict(find_ngrams(line.strip(), 2))))
    print("="*100)
    print(f"Trigram {line_number+1}")
    print(str(dict(find_ngrams(line.strip(), 3))))
    print("="*100)
    print(f"Tetragram {line_number+1}")
    print(str(dict(find_ngrams(line.strip(), 4))))
    print("="*100)
fhr.close()