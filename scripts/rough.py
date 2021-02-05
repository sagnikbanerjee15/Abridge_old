
import re
insert_characters = {'A':'!','T':'"','G':"#",'C':"$"}
mismatch_characters = {'A':'(','T':')','G':'*','C':'+'}

def convertRegularCigarToExtendedCigar(regular_cigar,md,seq,actual_extended_cigar):
    
    cigar_as_a_list = ""
    md_as_a_list=""
    extended_cigar = []
    # Remove soft clips
    cigar_values=[int(s) for s in re.findall(r'\d+',regular_cigar)]
    cigar_alphabets=re.findall(r'[A-Z]',regular_cigar)
    if cigar_alphabets[0]=='S':
        seq = seq[cigar_values[0]:]
        cigar_alphabets=cigar_alphabets[1:]
        cigar_values=cigar_values[1:]
    if cigar_alphabets[-1]=='S':
        seq = seq[:-cigar_values[-1]]
        cigar_alphabets=cigar_alphabets[:-1]
        cigar_values=cigar_values[:-1]
        
    for i in range(len(cigar_alphabets)):
        cigar_as_a_list+=cigar_alphabets[i]*cigar_values[i]
    
    md_values = [int(s) for s in re.findall(r'\d+',md)]
    md_alphabets=re.findall(r'[A-Z^]',md)
    
    i=0
    num=0
    while i<len(md):
        if md[i].isdigit():
            num=num*10+int(md[i])
        else:
            if md[i]=="^":
                md_as_a_list += num*"M"
                deleted_nucleotides = ""
                i+=1
                while md[i].isdigit()==False:
                    deleted_nucleotides+=md[i]
                    i+=1
                i-=1
                md_as_a_list += deleted_nucleotides
                #print("deleted_nucleotides: ",deleted_nucleotides)
            else:
                md_as_a_list += num*"M"
                md_as_a_list += md[i]
            #print(md[i],num)
            
            num=0
        i+=1
    md_as_a_list += num*"M"
    print("".join(md_as_a_list))
    cigar_as_a_list= [x for x in cigar_as_a_list]
    md_as_a_list = [x for x in md_as_a_list]

    for i in range(len(cigar_as_a_list)):
        if cigar_as_a_list[i]=='I':
            #insert_Is+="I"
            md_as_a_list.insert(i,'I')
    seq = [x for x in seq]
    
    for i in range(len(md_as_a_list)):
        if md_as_a_list[i] in ['A','T','G','C']:
            if cigar_as_a_list[i]=='D': # a deletion
                seq.insert(i,'-')
            elif cigar_as_a_list[i]=='M':
                seq[i]=mismatch_characters[seq[i]]
        elif md_as_a_list[i] == "I":
            seq[i]=insert_characters[seq[i]]
        
            
    extended_cigar = ""
    counter=0
    i=0
    while i<len(seq):
        if seq[i] in ['A','T','G','C']:
            counter+=1
        elif seq[i]=="-":
            extended_cigar+=str(counter)+"M"
            counter=0
            while seq[i]=='-':
                counter+=1
                i+=1
            extended_cigar+=str(counter)+"D"
            counter=0
            i-=1
        else:
            extended_cigar+=str(counter)+"M"
            extended_cigar+=seq[i]
            counter=0
        i+=1
    extended_cigar+=str(counter)+"M"
    
    # Remove 0M for final cigar
    extended_cigar = extended_cigar.replace("0M","")
    
    print("".join(cigar_as_a_list))
    print("".join(md_as_a_list))
    print("".join(seq))
    print(extended_cigar)
    return extended_cigar



convertRegularCigarToExtendedCigar("10M1609N28M3D29M1I31M","38^ATC0T2T56","TGAAAAAAAAAAAAAAAGTTTGGGTTTTTAAAAGAAGTCTTCTTTTTTTTGGGTTTATGGGCTGAAAAAAATTGTTGTGTTTTTTGTATTTTATTTTTT","")

