import logging


def assessMemoryRequirement(fai_index_filename,options):
    chromosome_to_length = {}
    fhr=open(fai_index_filename,"r")
    for line in fhr:
        chromosome,length = line.strip().split()[:2]
        chromosome_to_length[chromosome]=int(length)
    fhr.close()
    
    max_length=max([chromosome_to_length[chromosome] for chromosome in chromosome_to_length])
    if max_length>int(options.max_memory):
        if options.quiet == False:
            logging.info(f"Please restart the program with at least {max_length/(1024*1024*1024)} GB of memory ")
    
    total_length=sum([chromosome_to_length[chromosome] for chromosome in chromosome_to_length])
    length_in_gb = total_length/(1024*1024*1024)
    if length_in_gb > int(options.max_memory):
        if options.quiet == False:
            logging.info("Large genome - each chromosome will be loaded separately for decompression")
    