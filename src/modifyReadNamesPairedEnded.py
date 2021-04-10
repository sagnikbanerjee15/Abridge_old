#! /usr/bin/env python

import sys

def convertReadNameIndexToReadNameString(read_name_index,read_name_length,alphabets,alphabet_size):
    s=""
    for i in range(read_name_length):
        s+=alphabets[read_name_index[i]]
    return s

def generateNextID(read_name_index,read_name_length,alphabets):
    if read_name_length==0:
        read_name_length+=1
        read_name_index[0]=0
    else:
        # Check if the read is the last element of the read_name_length
        for i in range(read_name_length):
            if(read_name_index!=alphabets[alphabet_size-1]):
                break
        if i==read_name_length: # Last read reached
            read_name_length+=1
            for i in range(read_name_length):
                read_name_index[i]=0
        else:
            i=read_name_length-1
            while i>=0:
                if read_name_index[i] == alphabet_size - 1:
                    read_name_index[i] = 0
                else:
                    read_name_index[i]+=1
                    break
                i-=1
    return convertReadNameIndexToReadNameString(read_name_index,read_name_length,alphabets)
            

def main():
    input_samfilename = sys.argv[0]
    output_samfilename = sys.argv[1]
    alphabets="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,."
    read_name_length=0
    read_name_index = [0 for _ in range(100)]
    alphabet_size=strlen(alphabets)
    
    """fhr=open(input_samfilename,"r")
    for line in fhr:
        if(line[0]=="."):
            continue"""
        
    for i in range(10000000):
        generateNextID(read_name_index, read_name_length, alphabets)
    
    

if __name__ == "__main__":
    main()