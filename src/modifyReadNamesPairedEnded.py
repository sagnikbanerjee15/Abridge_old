#! /usr/bin/env python

import sys

def main():
    input_samfilename = sys.argv[0]
    output_samfilename = sys.argv[1]
    
    fhr=open(input_samfilename,"r")
    for line in fhr:
        if(line[0]=="."):
            continue
    
    

if __name__ == "__main__":
    main()