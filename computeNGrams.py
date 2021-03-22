from collections import Counter 
import pprint

def find_ngrams(input_list, n):
    res = Counter(input_list[idx : idx + n] for idx in range(len(input_list) - 1))
    return sorted(res.items(),key=lambda item: (-item[1], item[0])) 

input_filename = "/90daydata/maizegdb/sagnik/ABRIDGE/developing_abridge/SRR13711355_0_SE_Aligned.sortedByCoord.out.sam.qual.rle.line3"
output_filename = "/90daydata/maizegdb/sagnik/ABRIDGE/developing_abridge/SRR13711355_0_SE_Aligned.sortedByCoord.out.sam.qual.rle.line3.trigramcompression"
fhr=open(input_filename,"r")
for line_number,line in enumerate(fhr):
    bigram = dict(find_ngrams(line.strip(), 2))
    trigram = dict(find_ngrams(line.strip(), 3))
    conversion_table = {}
    starting_ASCII_code_for_single_character_replacement = 74
    for key in trigram:
        #print(key,trigram[key])
        if trigram[key]>1000:
            line = line.replace(key,chr(starting_ASCII_code_for_single_character_replacement))
            conversion_table[chr(starting_ASCII_code_for_single_character_replacement)]=key
            starting_ASCII_code_for_single_character_replacement+=1
    fhw = open(output_filename,"w")
    fhw.write(str(conversion_table)+"\t")
    fhw.write(line)
    fhw.close()
    """
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
    """
fhr.close()