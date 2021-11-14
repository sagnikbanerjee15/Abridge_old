#! /usr/bin/env python3

########################################################################################################################################################
# The software 'abridge' will compress aligned files to a bare minimum needed for generating assemblies and producing read counts
########################################################################################################################################################

from argparse import RawTextHelpFormatter
import argparse 
import logging
import os
import pprint
import sys
import re
import time
import multiprocessing
import random
import glob
import time
import subprocess
from pprint import pformat


from assessMemoryRequirement import *
from verifyPositions import *
from sortPositionsForRandomAccess import *

def parseCommandLineArguments():
    parser = argparse.ArgumentParser(prog="abridge",description="Compress alignments for storage, decompress from compressed file, view alignments from random locations and generate coverages",formatter_class=RawTextHelpFormatter)
    required_named = parser.add_argument_group('Required arguments')
    optional_named = parser.add_argument_group('Optional arguments')
    
    # Required arguments
    required_named.add_argument("-o","--output_directory",help="Enter the name of the output directory. If nothing is specified then the compressed file will be put in the same location as the input samfile")
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-isam","--inputsamfilenames",help="Enter the name of the alignment file you wish to compress. Alignments in SAM format only is expected. Ensure that the file is sorted by coordinate. Also, files must have the header section with the reference information available. You can compress only one file at a time.")
    input_group.add_argument("-iabr","--inputabrfilenames",help="Enter the name of the compressed alignment files you wish to merge. These files must be compressed using abridge. You can decompress only one file at a time.")
    required_named.add_argument("-g","--genome",help="Enter a single fasta file for the reference",required=True)
    
    compress_decompress_group = parser.add_mutually_exclusive_group(required=True)
    compress_decompress_group.add_argument("-cmp","--compress",help="Set this option if you wish to compress the alignment file",action="store_true")
    compress_decompress_group.add_argument("-dcmp","--decompress",help="Set this option if you wish to decompress the alignment file",action="store_true")
    compress_decompress_group.add_argument("-r","--random",help="Retrieve alignments from random locations",action="store_true")
    compress_decompress_group.add_argument("-H","--header",help="Print only the header of reference sequences during decompression",action="store_true")
    
    # Optional arguments
    optional_named.add_argument("-l","--level",help = "This can accept an integer from the set (1,2,3). If level is set to 1 then abridge will perform the fastest but low compression. abridge will use brotli to compress. Decompression will be fast. Setting level to 2 will prompt abridge to perform the medium level compression using 7z. Compression will take time but decompression will be fast. If level is set to 3 then abridge will perform the best compression using 7paq. Both compression and decompression will take average time to complete", type=int, default = 2)
    optional_named.add_argument("-ss","--ignore_scores",help = "Request abrigde to store the quality scores and the alignment score",action="store_true")
    optional_named.add_argument("-igqual","--ignore_quality_scores",help="Ignore all quality scores",action="store_true")
    optional_named.add_argument("-qual","--quality",help="Enter dummy quality scores while decompressing",default='I')
    optional_named.add_argument("-gsc","--ignore_soft_clippings",help="No soft clippings will be stored. Read will be trimmed down to only the portion which matched to nucleotides in the reference",action="store_true")
    optional_named.add_argument("-gm","--ignore_mismatches",help="All mismatches will be ignored",action="store_true")
    optional_named.add_argument("-gs","--ignore_sequence",help="No nucleotide sequence will be produced during decompression",action="store_true")
    optional_named.add_argument("-gu","--ignore_unmapped_reads",help="Request abridge to discard all reads that are unmapped",action="store_true")
    optional_named.add_argument("-sq","--save_all_quality_scores",help="Request abridge to save all quality scores",action='store_true')
    optional_named.add_argument("-aq","--save_exact_quality_scores",help="Adjust quality scores for matched bases to achieve better encoding. For more details please check ...",action="store_true")
    optional_named.add_argument("-q","--quiet",help="Prevent abridge from printing any log information. By default logging is enables",action = "store_true")
    optional_named.add_argument("-n","--cpu",help="Enter the number of CPU cores to be used. This option will be used during compression or decompression.",default = 1)
    optional_named.add_argument("-run_diagnostics","--run_diagnostics",help="abridge will run diagnostics on the cigar compression and decompression. It will exit on discovery of any discrepancies",action="store_true")
    optional_named.add_argument("-p","--positions",help="Enter the position as chromosome:start-end from which reads will be retrieved")
    optional_named.add_argument("-rp","--read_prefix",help="Enter a read prefix for decompression - valid only for random access")
    optional_named.add_argument("--keep_intermediate_error_files","-kief",help="Set this argument if you wish to preserve the intermediate error files to assess time and memory usage. Default behaviour is to delete those",action="store_true")
    optional_named.add_argument("--error_directory","-edir",help="Enter a directory where all error files will be stored. If nothing is specified then error files will be stored in the output directory")
    
    # Suppressed arguments
    parser.add_argument("--logfilename","-logfilename",help=argparse.SUPPRESS)# Name of the logfile
    parser.add_argument("--files_for_removal","-files_for_removal",help=argparse.SUPPRESS)# Files will be removed later
    parser.add_argument("--softwares","-softwares",help=argparse.SUPPRESS) # Software paths
    parser.add_argument("--single_ended","-single_ended",help=argparse.SUPPRESS)
    parser.add_argument("--reference_to_length","-num_of_reference_sequences",help=argparse.SUPPRESS)
    parser.add_argument("--outputfilename","-outputfilena",help=argparse.SUPPRESS) 
    parser.add_argument("--compile_programs","-compile_programs",action = "store_true", help=argparse.SUPPRESS) # Force abridge to compile the C programs 
    
    # Future enhancements
    compress_decompress_group.add_argument("-ov","--generate_overlapping_coverage",help=argparse.SUPPRESS, action = "store_true") # Future - This option can be used in conjuction with --positions to construct coverage from a specific location # help="Enter the name of the compressed file from which you wish to generate an overlapping coverage of reads ",  
    compress_decompress_group.add_argument("-nov","--generate_non_overlapping_coverage",help=argparse.SUPPRESS, action = "store_true") # help="Enter the name of the compressed file from which you wish to generate a non-overlapping coverage of reads "
    
    # Options for generating coverage
    optional_named.add_argument("-d","--d",help=argparse.SUPPRESS, action = "store_true") # help = "Report the depth at each position in each A feature. Positions reported are one based.  Each position and depth follow the complete A feature.",
    optional_named.add_argument("-bg","--bg",help=argparse.SUPPRESS, action = "store_true")
    optional_named.add_argument("-bga","--bga",help=argparse.SUPPRESS, action = "store_true")
    optional_named.add_argument("-split","--split",help=argparse.SUPPRESS, action = "store_true") # help = "Treat \"split\" BAM or BED12 entries as distinct BED intervals.",
    optional_named.add_argument("-mem","--max_memory",help=argparse.SUPPRESS, default=10) # help="Enter the maximum memory allowed (in GB)"
    optional_named.add_argument("-t","--produce_tags",help=argparse.SUPPRESS,nargs="*") # help="Enter a comma separated list of tags that you want abridge to produce during decompression. By default abridge will generate NH, MD and XS tags."
    return parser.parse_args()

def configureLogger(options):
    if os.path.exists(options.logfilename)==True:
        os.system(f"rm -f {options.logfilename}")
    logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S',level=logging.DEBUG, filename=options.logfilename)
    
def validateCommandLineArguments(options):
    """
    """
    if options.compress == True and options.inputsamfilenames is None and options.inputabrfilenames is not None:
        print("For compression you need to provide a list of space spearated samfiles using -isam ")
        if options.quiet==False:
            logging.info("For compression you need to provide a list of space spearated samfiles using -isam ")
        sys.exit()
    if options.decompress == True and options.inputsamfilenames is not None and options.inputabrfilenames is None:
        print("For decompression you need to provide a list of abridge compressed files using -iabr")
        if options.quiet==False:
            logging.info("For decompression you need to provide a list of abridge compressed files using -iabr")
        sys.exit()
    if options.inputsamfilenames is not None:
        inputfiles = options.inputsamfilenames
    else:
        inputfiles = options.inputabrfilenames
    
    if os.path.exists(inputfiles)==False:
        print(f"The input file {inputfiles} does not exist. Exiting...")
        if options.quiet==False:
            logging.info(f"The input file {inputfiles} does not exist. Exiting...")
        sys.exit()
     
    #Check if the input format is sam
    if options.inputsamfilenames is not None:
        if options.inputsamfilenames[-3:]!="sam" and options.compress==True:
            print(f"The input file {options.inputsamfilenames} needs to be in sam format. Exiting...")
            if options.quiet==False:
                logging.info(f"The input file {options.inputsamfilenames} needs to be in sam format. Exiting...")
            sys.exit()
    
    options.files_for_removal = []
    options.softwares = {}
    abridge_location = subprocess.run(['which', 'abridge'], stdout=subprocess.PIPE).stdout.decode('utf-8').split()[-1]
    if options.header==True:
        options.quiet = True
    
    if options.random==True or options.generate_overlapping_coverage==True or options.generate_non_overlapping_coverage==True:
        options.output_directory = "/".join(options.inputabrfilenames.split("/")[:-1])+"/"+str(int(time.time()))
        os.system(f"mkdir -p {options.output_directory}")
    
    if options.random==True and options.positions==None:
        print("You need to enter at least one location from where reads will be retrieved")
        sys.exit()
        
    if options.random == True and options.inputabrfilenames is None:
        print(f"You need to provide one compressed file for random access")
        sys.exit()
        
    if options.random == True and options.inputsamfilenames is not None:
        print(f"You need to provide one compressed file for random access")
        sys.exit()
        
    if options.error_directory == None:
        options.error_directory = options.output_directory
    else:
        os.system(f"mkdir -p {options.error_directory}")
    
    # save_all_quality_scores has higher precedence over ignore_quality_scores
    if options.save_all_quality_scores == True:
        options.ignore_quality_scores = True
        
    if options.save_all_quality_scores == False and options.save_exact_quality_scores == True:
        options.save_all_quality_scores = True
        options.ignore_quality_scores = True
        
    options.softwares["determineEndedness"] = "/".join(abridge_location.split("/")[:-1])+"/src/determineEndedness"
    options.softwares["compressSamFileSingleEnded"] = "/".join(abridge_location.split("/")[:-1])+"/src/compressSamFileSingleEnded"
    options.softwares["compressSamFileSingleEndedPass2"] = "/".join(abridge_location.split("/")[:-1])+"/src/compressSamFileSingleEndedPass2"
    options.softwares["compressSamFilePairedEnded"] = "/".join(abridge_location.split("/")[:-1])+"/src/compressSamFilePairedEnded"
    options.softwares["splitSamFileIntoEachReferenceSequence"] = "/".join(abridge_location.split("/")[:-1])+"/src/splitSamFileIntoEachReferenceSequence"
    options.softwares["buildABRIDGEIndex"]= "/".join(abridge_location.split("/")[:-1])+"/src/buildAbridgeIndex"
    options.softwares["decompressSamFileSingleEnded"]= "/".join(abridge_location.split("/")[:-1])+"/src/decompressSamFileSingleEnded"
    options.softwares["extractSequencesFromReferences"] = "/".join(abridge_location.split("/")[:-1])+"/scripts/extractSequencesFromReferences.py"
    options.softwares["randomRetrievalSingleEnded"] = "/".join(abridge_location.split("/")[:-1])+"/src/randomRetrievalSingleEnded"
    options.softwares["randomRetrievalPairedEnded"] = "/".join(abridge_location.split("/")[:-1])+"/src/randomRetrievalPairedEnded"
    options.softwares["mergeCompressedFilesSingleEnded"] = "/".join(abridge_location.split("/")[:-1])+"/src/mergeCompressedFilesSingleEnded"
    options.softwares["addTagToSamFile"] = "/".join(abridge_location.split("/")[:-1])+"/src/addTagToSamFile"
    options.softwares["maxReadsMappedToSingleNucleotide"] = "/".join(abridge_location.split("/")[:-1])+"/src/maxReadsMappedToSingleNucleotide"
    options.softwares["compressQualityScoresFile"] = "/".join(abridge_location.split("/")[:-1])+"/src/compressQualityScoresFile"
    options.softwares["deCompressQualityScoresFile"] = "/".join(abridge_location.split("/")[:-1])+"/src/deCompressQualityScoresFile"
    #options.softwares["separateReadsMatePairTogetherAndMatePairAway"] = "/".join(abridge_location.split("/")[:-1])+"/src/separateReadsMatePairTogetherAndMatePairAway"
    #options.softwares["modifyReadNamesPairedEnded"] = "/".join(abridge_location.split("/")[:-1])+"/src/modifyReadNamesPairedEnded"
    options.softwares["compressSamFilePairedEnded"] = "/".join(abridge_location.split("/")[:-1])+"/src/compressSamFilePairedEnded"
    options.softwares["decompressSamFilePairedEnded"] = "/".join(abridge_location.split("/")[:-1])+"/src/decompressSamFilePairedEnded"
    options.softwares["generateCoverage"] = "/".join(abridge_location.split("/")[:-1])+"/src/generateCoverage"
    options.softwares["findMaximumNumberOfReadsInEachLine"] = "/".join(abridge_location.split("/")[:-1])+"/src/findMaximumNumberOfReadsInEachLine"
    options.softwares["maxReadsInEachLine"] = "/".join(abridge_location.split("/")[:-1])+"/src/maxReadsInEachLine"
    
    genome_filename_without_location = options.genome.split("/")[-1]
    cmd = f"perl -pe '/^>/ ? print \"\\n\" : chomp' {options.genome} | tail -n +2 > {options.output_directory}/{genome_filename_without_location}"
    if options.quiet == False:
        logging.info(f"Running command - {cmd}")
    os.system(cmd)
    os.system(f"samtools faidx {options.output_directory}/{genome_filename_without_location}")
    options.genome = f"{options.output_directory}/{genome_filename_without_location}"
    """
    cmd = f"mv {options.genome}.temp {options.genome}"
    os.system(cmd)
    """
    
    if options.decompress==True:
        if options.inputabrfilenames[-8:]!=".abridge":
            print(f"The input file {options.inputabrfilenames} needs to be in abridge format. Exiting...")
            if options.quiet==False:
                logging.info(f"The input file {options.inputabrfilenames} needs to be in abridge format. Exiting...")
            sys.exit()
        
    if options.ignore_soft_clippings == True:
        options.ignore_soft_clippings = 1
    else:
        options.ignore_soft_clippings = 0
        
    if options.ignore_mismatches == True:
        options.ignore_mismatches = 1
    else:
        options.ignore_mismatches = 0
        
    if options.ignore_sequence == True:
        options.ignore_sequence = 1
    else:
        options.ignore_sequence = 0
        
    if options.ignore_unmapped_reads == True:
        options.ignore_unmapped_reads = 1
    else:
        options.ignore_unmapped_reads = 0
        
    if options.ignore_quality_scores == True:
        options.ignore_quality_score_flag = 1
    else:
        options.ignore_quality_score_flag = 0
        
    if options.generate_overlapping_coverage == False and options.generate_non_overlapping_coverage == False:
        if True in [options.d, options.bg, options.bga, options.split]:
            print("Incorrect arguments. Please provide arguments -d -bg -bga or -split only when you wish to generate coverage")
            if options.quiet==False:
                logging.info("Incorrect arguments. Please provide arguments -d -bg -bga or -split only when you wish to generate coverage")
            sys.exit()
            
    if options.generate_overlapping_coverage == True and options.generate_non_overlapping_coverage == True:
        print("You can either generate overlapping or non-overlapping coverage. If you need to generate both please run abridge twice each time with either option")
        if options.quiet == False:
            logging.info("You can either generate overlapping or non-overlapping coverage. If you need to generate both please run abridge twice each time with either option")
        sys.exit()
    
    if options.generate_overlapping_coverage == False or options.generate_non_overlapping_coverage == False:
        if [options.d,options.bg,options.bga].count(True)>1:
            print("You can specify only one among -d, -bg or -bga")
            if options.quiet == False:
                logging.info("You can specify only one among -d, -bg or -bga")
            sys.exit()
    
    if options.inputsamfilenames is not None:
        input_filename_without_location = options.inputsamfilenames.split("/")[-1][:-4]
    else:
        input_filename_without_location = options.inputabrfilenames.split("/")[-1][:-8]
    options.outputfilename = f"{options.output_directory}/{input_filename_without_location}.abridge"

def checkSAMAlignments(options,logging):
    """
    """
    
    if "SO:coordinate" not in open(f"{options.inputsamfilenames}","r").readline():
        print(f"The file {options.inputsamfilenames} is not sorted. Exiting")
        if options.quiet == False:
            logging.info(f"The file {file} is not sorted. Exiting")
        sys.exit()
    
    # Verify that all the headers are exactly same
    headers = []
    
    fhr=open(f"{options.inputsamfilenames}","r")
    header=""
    for line in fhr:
        if line[:3]=='@SQ' or line[:3]=="@HD":
            header+=line
        else:
            break
    headers.append(header)
    if len(set(headers))!=1:
        print("All the headers must be same. Please check your bamfile. Exiting...")
        if options.quiet == False:
            logging.info("All the headers must be same. Please check your bamfile. Exiting...")
            sys.exit()

def cleanUp(options):
    # Remove the genome files
    genome_filename_without_location = options.genome.split("/")[-1]
    options.genome = f"{options.output_directory}/{genome_filename_without_location}"
    os.system(f"rm -rf {options.genome}*")
    
    return
    for file in options.files_for_removal:
        cmd=f"rm -rf {file}"
        os.system(cmd)

def runCommand(eachpinput):
    cmd,dummy = eachpinput
    os.system(cmd)
    
def constructFileNames(input_filename,options):
    name_of_input_file_without_location = input_filename.split("/")[-1] 
    
    name_of_max_input_reads_file = f"{options.output_directory}/{name_of_input_file_without_location}_max_input_reads"
    name_of_file_with_max_commas = f"{options.output_directory}/{name_of_input_file_without_location}_max_commas"
    name_of_file_max_read_length = f"{options.output_directory}/{name_of_input_file_without_location}_max_read_size"
    name_of_total_number_of_alignments_file = f"{options.output_directory}/{name_of_input_file_without_location}_total_number_of_alignments"
    frequency_of_flags_filename = f"{options.output_directory}/{name_of_input_file_without_location}_frequency_of_flags"
    
    _outputfilename = f"{options.output_directory}/{name_of_input_file_without_location}_"
    _index_outputfilename = f"{options.output_directory}/{name_of_input_file_without_location}__index"
    unmapped_outputfilename = f"{options.output_directory}/{name_of_input_file_without_location}_unmapped"
    name_of_file_with_quality_scores = f"{options.output_directory}/{name_of_input_file_without_location}_qual"
    name_of_file_with_quality_scores_rle = f"{options.output_directory}/{name_of_input_file_without_location}_qual_rle"
    name_of_file_dictionary = f"{options.output_directory}/{name_of_input_file_without_location}.dictionary"
    
    error_file_max_input_reads = f"{name_of_max_input_reads_file}.error"
    error_file_compress_ = f"{_outputfilename}_compress.error"
    error_file_prep_abridge_index = f"{_index_outputfilename}.error"
    error_file_qual_rle = f"{name_of_file_with_quality_scores_rle}.error"
    
    error_file_7z_compression = f"{options.error_directory}/{name_of_input_file_without_location}_7z_compression.error"
    error_file_zpaq_compression = f"{options.error_directory}/{name_of_input_file_without_location}_zpaq_compression.error"
    error_file_brotli_compression =  f"{options.error_directory}/{name_of_input_file_without_location}_brotli_compression.error"
    
    # Decompression
    error_file_7z_decompression = f"{options.error_directory}/{name_of_input_file_without_location}_7z_decompression.error"
    error_file_zpaq_decompression = f"{options.error_directory}/{name_of_input_file_without_location}_zpaq_decompression.error"
    error_file_brotli_decompression =  f"{options.error_directory}/{name_of_input_file_without_location}_brotli_decompression.error"
    
    
    delete_these_files = [name_of_max_input_reads_file, name_of_file_with_max_commas, name_of_file_max_read_length, name_of_total_number_of_alignments_file, frequency_of_flags_filename]
    other_files = [_outputfilename, _index_outputfilename, unmapped_outputfilename, name_of_file_with_quality_scores, name_of_file_with_quality_scores_rle, name_of_file_dictionary]
    error_files = [error_file_max_input_reads, error_file_compress_, error_file_qual_rle, error_file_prep_abridge_index, error_file_7z_compression, error_file_zpaq_compression, error_file_brotli_compression, error_file_7z_decompression, error_file_zpaq_decompression, error_file_brotli_decompression]
    #compressed_files = [compressed_abridged_filename_br,compressed_abridged_filename_zpaq]
    return  delete_these_files, other_files, error_files

def compressSamFile(options):
    """
    """
    pool = multiprocessing.Pool(processes=int(options.cpu))
    files_to_be_removed = []
    input_filename = options.inputsamfilenames
    name_of_input_file_without_location = input_filename.split("/")[-1] 
    delete_these_files, other_files, error_files = constructFileNames(input_filename,options) 
    name_of_max_input_reads_file, name_of_file_with_max_commas, name_of_file_max_read_length, name_of_total_number_of_alignments_file, frequency_of_flags_filename = delete_these_files
    _outputfilename, _index_outputfilename, unmapped_outputfilename, name_of_file_with_quality_scores, name_of_file_with_quality_scores_rle, name_of_file_dictionary = other_files
    error_file_max_input_reads, error_file_compress_, error_file_qual_rle, error_file_prep_abridge_index, error_file_7z_compression, error_file_zpaq_compression, error_file_brotli_compression, error_file_7z_decompression, error_file_zpaq_decompression, error_file_brotli_decompression = error_files
    #compressed_abridged_filename_7z, compressed_abridged_filename_br,compressed_abridged_filename_zpaq = compressed_files
    
    if options.single_ended == True:
        # Single ended
        ######################################################################################
        # Compiling programs - will be removed during final version
        ######################################################################################
        compressSamFileSingleEnded = options.softwares["compressSamFileSingleEnded"]
        maxReadsMappedToSingleNucleotide = options.softwares["maxReadsMappedToSingleNucleotide"]
        buildABRIDGEIndex = options.softwares["buildABRIDGEIndex"]
        compressQualityScoresFile = options.softwares["compressQualityScoresFile"]
        deCompressQualityScoresFile = options.softwares["deCompressQualityScoresFile"]
        
        if options.compile_programs == True:
            cmd=f"gcc {compressSamFileSingleEnded}.c -o {compressSamFileSingleEnded} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
            cmd=f"gcc {buildABRIDGEIndex}.c -o {buildABRIDGEIndex} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
            cmd=f"gcc {maxReadsMappedToSingleNucleotide}.c -o {maxReadsMappedToSingleNucleotide} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
            cmd=f"gcc {compressQualityScoresFile}.c -o {compressQualityScoresFile} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
            cmd=f"gcc {deCompressQualityScoresFile}.c -o {deCompressQualityScoresFile} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
        ######################################################################################
        
        ######################################################################################
        # Find maximum number of reads mapped to a single location
        # Find maximum read length
        # Find total number of alignments
        ######################################################################################
        
        all_commands = []
        #for input_filename in options.inputsamfilenames:
        input_filename = options.inputsamfilenames
        if options.quiet == False:
            logging.info(f"Starting compression for {input_filename}")
    
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" {maxReadsMappedToSingleNucleotide} "
        cmd += f" {input_filename} "
        cmd += f" {name_of_max_input_reads_file}"
        cmd += f" {name_of_total_number_of_alignments_file} "
        cmd += f" {name_of_file_max_read_length}"
        cmd += f")"
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_max_input_reads}"
        else:
            cmd += f"2> /dev/null"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        all_commands.append([cmd,"dummy"])
        os.system(cmd)
        #pool.map(runCommand,all_commands)
        
        ######################################################################################
        #  run
        ######################################################################################
        all_commands = []
        #for input_filename in options.inputsamfilenames:
        input_filename = options.inputsamfilenames
        if options.quiet == False:
            logging.info(f"Starting  for {input_filename}")
        max_input_reads_in_a_single_nucl_loc = int(open(name_of_max_input_reads_file,"r").read().strip())
        if options.quiet == False:
            logging.info(f"{name_of_input_file_without_location} max_input_reads_in_a_single_nucl_loc = {max_input_reads_in_a_single_nucl_loc}")
            
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f"{compressSamFileSingleEnded}" # argv[0] - name of the program
        cmd += f" {options.genome}" # argv[1] - genome_filename
        cmd += f" {options.ignore_soft_clippings}" #argv[2]
        cmd += f" {options.ignore_mismatches} " #argv[3]
        cmd += f" {options.ignore_quality_score_flag} " #argv[4]
        cmd += f" {options.ignore_unmapped_reads} " #argv[5]
        cmd += f" {input_filename} " #argv[6] - Input samfile
        cmd += f" {_outputfilename} " #argv[7] - Output filename
        cmd += f" {unmapped_outputfilename} " #argv[8] - Filename for unmapped reads
        if options.run_diagnostics==True:
            cmd += f" 1" # argv[9] - run diagnostics
        else:
            cmd += f" 0" # argv[9] - run diagnostics
        cmd += f" {max_input_reads_in_a_single_nucl_loc} " # argv[10] - max_input_reads_in_a_single_nucl_loc
        cmd += f" {name_of_file_with_max_commas} " # argv[11] - name_of_file_with_max_commas
        if options.save_all_quality_scores == True and options.save_exact_quality_scores == True:
            cmd += f" 1 1 " # argv[12] and argv[13]
        elif options.save_all_quality_scores == True and options.save_exact_quality_scores == False:
            cmd += f" 1 0 " # argv[12] and argv[13]
        elif options.save_all_quality_scores == False and options.save_exact_quality_scores == True: 
            cmd += f" 0 1 " # argv[12] and argv[13]
        elif options.save_all_quality_scores == False and options.save_exact_quality_scores == False:
            cmd += f" 0 0 " # argv[12] and argv[13]
        cmd += f" {name_of_file_with_quality_scores} " # argv[14]
        if options.ignore_scores == True:
            cmd += f" 1 " #argv[15] - save_scores = False
        else:
            cmd += f" 0 " #argv[15] - save_scores = True
        cmd += f") "
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_compress_}"
        else:
            cmd += f"2> /dev/null"
        
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        all_commands.append([cmd,"dummy"])
        os.system(cmd)
        
        ######################################################################################
        # Generate RLE for quality scores
        ######################################################################################
        max_read_length = int(open(name_of_file_max_read_length,"r").read().strip())
        if options.quiet == False:
            logging.info(f"{name_of_input_file_without_location} max_read_length = {max_read_length}")
        if options.save_all_quality_scores == True and options.save_exact_quality_scores==False:
            all_commands = []
            cmd  = f"(/usr/bin/time --verbose "
            cmd += f" {compressQualityScoresFile} "
            cmd += f" {name_of_file_with_quality_scores} "
            cmd += f" {name_of_file_with_quality_scores_rle} "
            cmd += f" {max_read_length} "
            if options.save_exact_quality_scores==True:
                cmd += " 1 "
            else:
                cmd += " 0 "
            cmd +=")"
            cmd += f"1> /dev/null "
            if options.keep_intermediate_error_files == True:
                cmd += f"2> {error_file_qual_rle}"
            else:
                cmd += f"2> /dev/null"
            if options.quiet == False:
                logging.info(f"Running cmd - {cmd}")
            all_commands.append([cmd,"dummy"])
            os.system(cmd)
        
        ######################################################################################
        # Index generation
        ######################################################################################
        max_commas = int(open(name_of_file_with_max_commas,"r").read().strip())
        if options.quiet == False:
            logging.info(f"{name_of_input_file_without_location} max_commas = {max_commas}")
        
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f"{buildABRIDGEIndex} {_outputfilename} {name_of_file_with_quality_scores} {_index_outputfilename} {max_commas} 0 "
        if options.ignore_scores == True:
            cmd += f" 1 " #argv[15] - save_scores = False
        else:
            cmd += f" 0 " #argv[15] - save_scores = True
        cmd += f") "   
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_prep_abridge_index}"
        else:
            cmd += f"2> /dev/null"
        if options.quiet == False:
            logging.info(f"Starting to build index for {name_of_input_file_without_location}")
            logging.info(f"Running cmd - {cmd}")
        all_commands.append([cmd,"dummy"])
        os.system(cmd)
        #pool.map(runCommand,all_commands)
        
    else:
        # Paired ended
        ######################################################################################
        # Compiling programs - will be removed during final version
        ######################################################################################
        maxReadsMappedToSingleNucleotide = options.softwares["maxReadsMappedToSingleNucleotide"]
        compressSamFilePairedEnded = options.softwares["compressSamFilePairedEnded"]
        buildABRIDGEIndex = options.softwares["buildABRIDGEIndex"]
        compressQualityScoresFile = options.softwares["compressQualityScoresFile"]
        deCompressQualityScoresFile = options.softwares["deCompressQualityScoresFile"]
        
        if options.compile_programs == True:
            cmd=f"gcc {buildABRIDGEIndex}.c -o {buildABRIDGEIndex} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
            cmd=f"gcc {maxReadsMappedToSingleNucleotide}.c -o {maxReadsMappedToSingleNucleotide} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
            cmd=f"gcc {compressQualityScoresFile}.c -o {compressQualityScoresFile} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
            cmd=f"gcc {deCompressQualityScoresFile}.c -o {deCompressQualityScoresFile} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
            cmd=f"gcc {compressSamFilePairedEnded}.c -o {compressSamFilePairedEnded} -Ofast"
            if options.quiet == False:
                logging.info(f"Running command - {cmd}")
            os.system(cmd)
        ######################################################################################
        
        
        ######################################################################################
        # Find maximum number of reads mapped to a single location
        # Find maximum read length
        # Find total number of alignments
        ######################################################################################
        all_commands = []
        #for input_filename in options.inputsamfilenames:
        input_filename = options.inputsamfilenames
        if options.quiet == False:
            logging.info(f"Starting  for {input_filename}")
        
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" {maxReadsMappedToSingleNucleotide} "
        cmd += f" {input_filename} "
        cmd += f" {name_of_max_input_reads_file} "
        cmd += f" {name_of_total_number_of_alignments_file} "
        cmd += f" {name_of_file_max_read_length} "
        cmd += f")"
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_max_input_reads}"
        else:
            cmd += f"2> /dev/null"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        all_commands.append([cmd,"dummy"])
        os.system(cmd)
        #pool.map(runCommand,all_commands)
        
        all_commands = []
        #for input_filename in options.inputsamfilenames:
        input_filename = options.inputsamfilenames
        if options.quiet == False:
            logging.info(f"Generating for {input_filename}")
            
        cmd  = f" cat "
        cmd += f" {input_filename} "
        cmd += f" |grep -v ^@ "
        cmd += f" |cut -f2|uniq "
        cmd += f" |sort|uniq > {frequency_of_flags_filename}"
        all_commands.append([cmd,"dummy"])
        os.system(cmd)
        #pool.map(runCommand,all_commands)
        
        all_commands = []
        #for input_filename in options.inputsamfilenames:
        input_filename = options.inputsamfilenames
        if options.quiet == False:
            logging.info(f"Starting  for {input_filename}")
        
        max_input_reads_in_a_single_nucl_loc = int(open(name_of_max_input_reads_file,"r").read().strip())
        if options.quiet == False:
            logging.info(f"{name_of_input_file_without_location} max_input_reads_in_a_single_nucl_loc = {max_input_reads_in_a_single_nucl_loc}")
            
        max_number_of_alignments = int(open(name_of_total_number_of_alignments_file,"r").read().strip())
        if options.quiet == False:
            logging.info(f"{name_of_input_file_without_location} max_number_of_alignments = {max_number_of_alignments}")
        
        max_read_length = int(open(name_of_file_max_read_length,"r").read().strip())
        if options.quiet == False:
            logging.info(f"{name_of_input_file_without_location} max_read_length = {max_read_length}")
        
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" {compressSamFilePairedEnded}" # argv[0] - name of the program
        cmd += f" {options.genome}" # argv[1] - genome_filename
        cmd += f" {options.ignore_soft_clippings}" #argv[2]
        cmd += f" {options.ignore_mismatches} " #argv[3]
        cmd += f" {options.ignore_quality_score_flag} " #argv[4]
        cmd += f" {options.ignore_unmapped_reads} " #argv[5]
        cmd += f" {input_filename} " #argv[6] - Input samfile
        cmd += f" {_outputfilename} " #argv[7] - Output filename
        cmd += f" {unmapped_outputfilename} " #argv[8] - Filename for unmapped reads
        if options.run_diagnostics==True:
            cmd += f" 1" # argv[9] - run diagnostics
        else:
            cmd += f" 0" # argv[9] - run diagnostics
        cmd += f" {max_input_reads_in_a_single_nucl_loc} " # argv[10] - max_input_reads_in_a_single_nucl_loc
        cmd += f" {name_of_file_with_max_commas} " # argv[11] - name_of_file_with_max_commas
        if options.save_all_quality_scores == True and options.save_exact_quality_scores == True:
            cmd += f" 1 1 " # argv[12] and argv[13]
        elif options.save_all_quality_scores == True and options.save_exact_quality_scores == False:
            cmd += f" 1 0 " # argv[12] and argv[13]
        elif options.save_all_quality_scores == False and options.save_exact_quality_scores == True: 
            cmd += f" 0 1 " # argv[12] and argv[13]
        elif options.save_all_quality_scores == False and options.save_exact_quality_scores == False:
            cmd += f" 0 0 " # argv[12] and argv[13]
        cmd += f" {name_of_file_with_quality_scores} " # argv[14]
        cmd += f" {max_number_of_alignments} " #argv[15]
        cmd += f" {max_read_length} " #argv[16]
        cmd += f" {frequency_of_flags_filename} " #argv[17]
        cmd += f" {name_of_file_dictionary} " #argv[18]
        if options.ignore_scores == True:
            cmd += f" 1 " #argv[19] - save_scores = False
        else:
            cmd += f" 0 " #argv[19] - save_scores = True
        cmd += f") "
        #cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_compress_}"
        else:
            cmd += f"2> /dev/null"
        
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        all_commands.append([cmd,"dummy"])
        os.system(cmd)
        #pool.map(runCommand,all_commands)
        
        ######################################################################################
        # Generate RLE for quality scores
        ######################################################################################
        if options.save_all_quality_scores == True and options.save_exact_quality_scores==False:
            cmd  = f"(/usr/bin/time --verbose "
            cmd += f" {compressQualityScoresFile} "
            cmd += f" {name_of_file_with_quality_scores} "
            cmd += f" {name_of_file_with_quality_scores_rle} "
            cmd += f" {max_read_length} "
            if options.save_exact_quality_scores==True:
                cmd += " 1 "
            else:
                cmd += " 0 "
            cmd +=")"
            cmd += f"1> /dev/null "
            if options.keep_intermediate_error_files == True:
                cmd += f"2> {error_file_qual_rle}"
            else:
                cmd += f"2> /dev/null"
            if options.quiet == False:
                logging.info(f"Running cmd - {cmd}")
            all_commands.append([cmd,"dummy"])
            os.system(cmd)

        ######################################################################################
        # Index generation
        ######################################################################################
        max_commas = int(open(name_of_file_with_max_commas,"r").read().strip())
        if options.quiet == False:
            logging.info(f"{name_of_input_file_without_location} max_commas = {max_commas}")
        
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f"{buildABRIDGEIndex} {_outputfilename} {name_of_file_with_quality_scores} {_index_outputfilename} {max_commas} 1 "
        if options.ignore_scores == True:
            cmd += f" 1 " 
        else:
            cmd += f" 0 " 
        cmd += f") "   
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_prep_abridge_index}"
        else:
            cmd += f"2> /dev/null"
        if options.quiet == False:
            logging.info(f"Starting to build index for {name_of_input_file_without_location}")
            logging.info(f"Running cmd - {cmd}")
        all_commands.append([cmd,"dummy"])
        os.system(cmd)
        #pool.map(runCommand,all_commands)
    
    ######################################################################################
    # Generate SHA512 for genome file
    ######################################################################################
    genome_sha512_filename = f"{options.output_directory}/genome_sha512"
    cmd = f"sha512sum {options.genome} > {genome_sha512_filename}.temp"
    os.system(cmd)
    open(genome_sha512_filename,"w").write(open(f"{genome_sha512_filename}.temp","r").read().strip().split()[0])
    os.system(f"rm {genome_sha512_filename}.temp")

    ######################################################################################
    # Compress files
    ######################################################################################
    
    
    ######################################################################################
    # Lowest compression (--1) using brotli
    ######################################################################################
    if options.level == 1:
        if options.quiet == False:
            logging.info("Starting compression with brotli")
        
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" brotli -9 -f -k -n "
        cmd += f" {_outputfilename} {unmapped_outputfilename} {_index_outputfilename} {genome_sha512_filename} "
        if options.single_ended == False:
            cmd += f" {name_of_file_dictionary} "
        if options.save_all_quality_scores == True and options.save_exact_quality_scores == False:
            cmd += f" {name_of_file_with_quality_scores_rle} "
        else:
            cmd += f" {name_of_file_with_quality_scores} " 
        cmd += f") "   
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_brotli_compression}"
        else:
            cmd += f"2> /dev/null"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        
        fhw = open(f"{options.output_directory}/compression_method","w")
        fhw.write("brotli")
        fhw.close()
        
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" 7za a -t7z -m0=lzma2 -mx=1 -mfb=64  -ms=on -mmt={options.cpu}  "
        cmd += f" {options.outputfilename} " # Compressed Output filename
        cmd += f" {_outputfilename}.br {unmapped_outputfilename}.br {_index_outputfilename}.br {genome_sha512_filename}.br {options.output_directory}/compression_method "
        if options.single_ended == False:
            cmd += f" {name_of_file_dictionary}.br "
        if options.save_all_quality_scores == True and options.save_exact_quality_scores == False:
            cmd += f" {name_of_file_with_quality_scores_rle}.br "
        else:
            cmd += f" {name_of_file_with_quality_scores}.br " 
        cmd += f") "   
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_7z_compression}"
        else:
            cmd += f"2> /dev/null"
        os.system(cmd)
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        
    
    ######################################################################################
    # Medium compression (--2) using 7z
    ######################################################################################
    if options.level == 2:
        if options.quiet == False:
            logging.info("Starting compression with 7z")
        fhw = open(f"{options.output_directory}/compression_method","w")
        fhw.write("7z")
        fhw.close()
        
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" 7za a -t7z -m0=lzma2 -mx=9 -mfb=64  -ms=on -mmt={options.cpu}  "
        cmd += f" {options.outputfilename} " # Compressed Output filename
        cmd += f" {_outputfilename} {unmapped_outputfilename} {_index_outputfilename} {genome_sha512_filename} {options.output_directory}/compression_method "
        if options.single_ended == False:
            cmd += f" {name_of_file_dictionary} "
        if options.save_all_quality_scores == True and options.save_exact_quality_scores == False:
            cmd += f" {name_of_file_with_quality_scores_rle} "
        else:
            cmd += f" {name_of_file_with_quality_scores} "
        cmd += f") "   
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_7z_compression}"
        else:
            cmd += f"2> /dev/null"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
    
    ######################################################################################
    # Best compression (--3) using zpaq
    ######################################################################################
    
    # Compress the quality scores with zpaq if save_exact_quality_scores is requested
    if options.level == 3:
        if options.quiet == False:
            logging.info("Starting compression with zpaq")
        fhw = open(f"{options.output_directory}/compression_method","w")
        fhw.write("zpaq")
        fhw.close()
        
        zpaq = f" zpaq "
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" {zpaq} a "
        cmd += f" {options.outputfilename}.zpaq "
        cmd += f" {_outputfilename} {unmapped_outputfilename} {_index_outputfilename} {genome_sha512_filename} "
        if options.single_ended == False:
            cmd += f" {name_of_file_dictionary} "
        if options.save_all_quality_scores == True and options.save_exact_quality_scores == False:
            cmd += f" {name_of_file_with_quality_scores_rle} "
        else:
            cmd += f" {name_of_file_with_quality_scores} "
        cmd += f" -m3 -t{options.cpu} -f "
        cmd += f" -noattributes "
        cmd += f") "
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_zpaq_compression}"
        else:
            cmd += f"2> /dev/null"
        os.system(cmd)
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
            
            
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" 7za a -t7z -m0=lzma2 -mx=1 -mfb=64  -ms=on -mmt=1  "
        cmd += f" {options.outputfilename} " # Compressed Output filename
        cmd += f" {options.outputfilename}.zpaq {options.output_directory}/compression_method "
        cmd += f") "   
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {error_file_7z_compression}"
        else:
            cmd += f"2> /dev/null"
        os.system(cmd)
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
            
    files_to_be_removed.append(f"{options.output_directory}/compression_method")
    files_to_be_removed.append(f"{options.outputfilename}.zpaq")
    
    files_to_be_removed.append(f"{_outputfilename}.br")
    files_to_be_removed.append(f"{unmapped_outputfilename}.br")
    files_to_be_removed.append(f"{_index_outputfilename}.br")
    files_to_be_removed.append(f"{genome_sha512_filename}.br")
    files_to_be_removed.append(f"{name_of_file_dictionary}.br")
    files_to_be_removed.append(f"{name_of_file_with_quality_scores}.br")
    files_to_be_removed.append(f"{name_of_file_with_quality_scores}.rle.br")
    files_to_be_removed.append(f"{options.output_directory}/compression_method")
    files_to_be_removed.append(f"{name_of_file_with_quality_scores_rle}.br")
            
    files_to_be_removed.append(f"{_outputfilename}")
    files_to_be_removed.append(f"{unmapped_outputfilename}")
    files_to_be_removed.append(f"{_index_outputfilename}")
    files_to_be_removed.append(f"{genome_sha512_filename}")
    files_to_be_removed.append(f"{name_of_file_dictionary}")
    files_to_be_removed.append(f"{name_of_file_with_quality_scores}")
    files_to_be_removed.append(f"{name_of_file_with_quality_scores_rle}")    
        
    
    files_to_be_removed.extend(delete_these_files)
    files_to_be_removed.extend(other_files)
    
    for file in files_to_be_removed:
        os.system("rm -f "+file)
    
        
          
def findEndedness(options,logging):
    """
    """
    endedness = []
    fhr=open(f"{options.inputsamfilenames}","r")
    for line in fhr:
        if line[0]!='@':
            samformatflag = int(line.strip().split()[1])
            endedness.append(samformatflag % 2)
            break
    
    if len(set(endedness))!=1:
        print("A mixture of single and paired ended files is not allowed. Exiting...")
        if options.quiet == False:
            logging.info("A mixture of single and paired ended files is not allowed. Exiting...")
            
    if endedness[0] == 0:
        options.single_ended = True
    else:
        options.single_ended = False
    
def collectReferenceSequenceNameAndLength(options,logging):
    """
    """
    reference_to_length = {}
    fhr=open(f"{options.inputfilename[0][:-3]}header","r")
    for line in fhr:
        if line[:3]=="@SQ":
            useless,reference_name,reference_length = line.strip().split("\t")
            reference_name = reference_name.split(":")[-1]
            reference_length = int(reference_length.split(":")[-1])
            reference_to_length[reference_name] = reference_length
    fhr.close()     
    options.reference_to_length = reference_to_length

def findChromosomes(filename):
    chromosomes = []
    fhr=open(filename,"r")
    for line_num,line in enumerate(fhr):
        if line_num==0:continue
        chromosomes.append(line.split()[0])
    fhr.close()
    return list(set(chromosomes))
    
def mergeAbridgeData(options):
    pool = multiprocessing.Pool(processes=int(options.cpu))
    
    cmd=f"samtools faidx {options.genome}"
    if os.path.exists(f"{options.genome}.fai")==False:
        if options.quiet ==False:
            logging.info(f"Running cmd - {cmd}")
        os.system(cmd)
        
    files_to_be_removed = []
    all_commands = []
    for input_filename in options.inputabrfilenames:
        input_filename_without_extension = input_filename.split("/")[-1][:-11]
        _filename = options.output_directory + "/" + input_filename_without_extension + "."
        _index_filename = options.output_directory + "/" + input_filename_without_extension+ "..index"
        unmapped_outputfilename = options.output_directory + "/" + input_filename_without_extension + ".unmapped"
        compressed_abridged_filename = input_filename
        
        cmd=f"7za e {input_filename} -y"
        cmd+=f" -o{options.output_directory}"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        all_commands.append([cmd,"dummy"])
    
        files_to_be_removed.extend([input_filename_without_extension,
                                    _filename,
                                    _index_filename,
                                    unmapped_outputfilename])
    pool.map(runCommand,all_commands)
    
    ######################################################################################
    # Compiling programs - will be removed during final version
    ######################################################################################
    mergeCompressedFilesSingleEnded = options.softwares["mergeCompressedFilesSingleEnded"]
    cmd=f"gcc -o {mergeCompressedFilesSingleEnded} {mergeCompressedFilesSingleEnded}.c -g -Ofast"
    os.system(cmd)
    
    #addTagToSamFile = options.softwares["addTagToSamFile"]
    #cmd=f"gcc -o {addTagToSamFile} {addTagToSamFile}.c -g -Ofast"
    #os.system(cmd)
    ######################################################################################
    
    
    # Compile a list of all chromosomes in all the files to be merged
    all_chromosomes = []
    for input_filename in options.inputabrfilenames:
        input_filename_without_extension = input_filename.split("/")[-1][:-11]
        _index_filename = options.output_directory + "/" + input_filename_without_extension+ "..index"
        all_chromosomes.extend(findChromosomes(_index_filename))
    all_chromosomes=list(set(all_chromosomes))
    all_chromosomes=sorted(all_chromosomes)
    """print(all_chromosomes)
    sys.stdout.flush()"""
    
    all_commands = []
    for chromosome in all_chromosomes:
        merged_output_filename = options.output_directory + "/" + options.output_directory.split("/")[-1] + ".merged." + chromosome 
        cmd=mergeCompressedFilesSingleEnded
        for input_filename in options.inputabrfilenames:
            input_filename_without_extension = input_filename.split("/")[-1][:-11]
            _filename = options.output_directory + "/" + input_filename_without_extension + "."
            _index_filename = options.output_directory + "/" + input_filename_without_extension+ "..index"
            unmapped_outputfilename = options.output_directory + "/" + input_filename_without_extension + ".unmapped"
            compressed_abridged_filename = input_filename
            
            cmd+=" "+_filename
        cmd+=" "+merged_output_filename
        cmd+=" "+chromosome
        all_commands.append([cmd,"dummy"]) 
        if options.quiet ==False:
            logging.info(f"Running cmd - {cmd}")
        ######################################################################################
        # Running in single core mode
        ######################################################################################
        #os.system(cmd)
    pool.map(runCommand,all_commands)
    
    
    for chromosome in all_chromosomes:
        cmd = "cat "
        cmd += options.output_directory + "/" + options.output_directory.split("/")[-1] + ".merged." + chromosome + ">> "
        cmd += " " + options.output_directory + "/" + options.output_directory.split("/")[-1] + ".merged"
        os.system(cmd)
        
        os.system("rm -f "+options.output_directory + "/" + options.output_directory.split("/")[-1] + ".merged." + chromosome) 
    
    for file in files_to_be_removed:
        os.system("rm "+file)
    
def decompressSamFile(options):
    """
    """
    pool = multiprocessing.Pool(processes=int(options.cpu))
    
    input_filename = options.inputabrfilenames[:-8]+".sam"
    name_of_input_file_without_location = input_filename.split("/")[-1] 
    delete_these_files, other_files, error_files = constructFileNames(input_filename,options) 
    name_of_max_input_reads_file, name_of_file_with_max_commas, name_of_file_max_read_length, name_of_total_number_of_alignments_file, frequency_of_flags_filename = delete_these_files
    _outputfilename, _index_outputfilename, unmapped_outputfilename, name_of_file_with_quality_scores, name_of_file_with_quality_scores_rle, name_of_file_dictionary = other_files
    error_file_max_input_reads, error_file_compress_, error_file_qual_rle, error_file_prep_abridge_index, error_file_7z_compression, error_file_zpaq_compression, error_file_brotli_compression, error_file_7z_decompression, error_file_zpaq_decompression, error_file_brotli_decompression = error_files
    
    cmd=f"samtools faidx {options.genome}"
    if os.path.exists(f"{options.genome}.fai")==False:
        if options.quiet ==False:
            logging.info(f"Running cmd - {cmd}")
        os.system(cmd)
    
    decompressSamFileSingleEnded = options.softwares["decompressSamFileSingleEnded"]
    decompressSamFilePairedEnded = options.softwares["decompressSamFilePairedEnded"]
    deCompressQualityScoresFile = options.softwares["deCompressQualityScoresFile"]
    maxReadsInEachLine = options.softwares["maxReadsInEachLine"]
    ######################################################################################
    # Compiling programs - will be removed during final version
    ######################################################################################
    if options.compile_programs == True:
        cmd=f"gcc {decompressSamFileSingleEnded}.c -o {decompressSamFileSingleEnded} -Ofast"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        
        cmd=f"gcc {decompressSamFilePairedEnded}.c -o {decompressSamFilePairedEnded} -Ofast"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        
        cmd=f"gcc {deCompressQualityScoresFile}.c -o {deCompressQualityScoresFile} -Ofast"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        
        cmd=f"gcc {maxReadsInEachLine}.c -o {maxReadsInEachLine} -Ofast"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
    ######################################################################################
    
    # Check the amount memory needed
    if options.ignore_sequence == 0:
        assessMemoryRequirement(f"{options.genome}.fai",options)
    
    files_to_be_removed = []
    
    # Use 7za to unpack the single compressed file
    cmd  = f"7za e {options.inputabrfilenames} -y -mmt={options.cpu} "
    cmd += f" -o{options.output_directory}"
    cmd += f" 1> /dev/null 2>/dev/null"
    if options.quiet == False:
        logging.info(f"Running command - {cmd}")
    os.system(cmd)
    time.sleep(1)
    
    compressor = open(f"{options.output_directory}/compression_method","r").read()
    if compressor == "brotli":
        cmd  = f"brotli -d -k -f "
        cmd += f"{_outputfilename}.br "
        cmd += f"{_index_outputfilename}.br "
        cmd += f"{unmapped_outputfilename}.br "
        if os.path.exists(f"{name_of_file_with_quality_scores_rle}.br") == True:
            cmd += f"{name_of_file_with_quality_scores_rle}.br "
        if os.path.exists(f"{name_of_file_with_quality_scores}.br") == True:
            cmd += f"{name_of_file_with_quality_scores}.br "
        if os.path.exists(f"{name_of_file_dictionary}.br") == True:
            cmd += f"{name_of_file_dictionary}.br "
        files_to_be_removed.extend(f"{_outputfilename}.br {_index_outputfilename}.br {unmapped_outputfilename}.br {name_of_file_with_quality_scores}.rle.br".split())
        os.system(cmd)
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
    elif compressor == "7z":
        # No need to decompress further
        pass
    elif compressor == "zpaq":
        filename = options.inputabrfilenames.split("/")[-1]
        zpaq = " zpaq "
        cmd = f"{zpaq} l {options.output_directory}/{filename}.zpaq > {options.output_directory}/file_list"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        fhr = open(f"{options.output_directory}/file_list","r")
        line_num = 0
        for line in fhr:
            line_num += 1
            if line_num == 4:
                compression_directory = line.strip().split()[-1]
                break
            
        #compression_directory = compression_directory.split("/")[-2]
        fhr.close()
        cmd = f"{zpaq} x {options.output_directory}/{filename}.zpaq {'/'.join(compression_directory.split('/')[:-1])} -to {options.output_directory} -t{options.cpu} "
        os.system(cmd)
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        # Move files
        cmd = f"mv {options.output_directory}/{compression_directory}/* {options.output_directory}"
        os.system(cmd)
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        files_to_be_removed.extend(f"{options.output_directory}/{compression_directory}")
        
    files_to_be_removed.extend([_outputfilename,
                                _index_outputfilename,
                                unmapped_outputfilename,
                                name_of_file_with_quality_scores,
                                f"{name_of_file_with_quality_scores}.rle",
                                f"{options.output_directory}/file_list"])
    
    flag_ignore_mismatches,flag_ignore_soft_clippings, flag_ignore_unmapped_sequences, flag_ignore_quality_score, save_all_quality_scores, save_exact_quality_scores,ignore_scores = list(map(int,open(_outputfilename,"rb").readline().split()))
    if flag_ignore_mismatches == 1:
        options.ignore_mismatches = True
    else:
        options.ignore_mismatches = False
    if flag_ignore_soft_clippings == 1:
        options.ignore_soft_clippings = True
    else:
        options.ignore_soft_clippings = False
    if flag_ignore_unmapped_sequences == 1:
        options.ignore_unmapped_sequences = True
    else:
        options.ignore_unmapped_sequences = False
    if flag_ignore_quality_score == 1:
        options.ignore_quality_score = True
    else:
        options.ignore_quality_score = False
    if save_all_quality_scores == 1:
        options.save_all_quality_scores = True
    else:
        options.save_all_quality_scores = False
    if save_exact_quality_scores == 1:
        options.save_exact_quality_scores = True
    else:
        options.save_exact_quality_scores = False
    if ignore_scores == 1:
        options.ignore_scores = True
    else:
        options.ignore_scores = False
    
    ######################################################################################
    # Convert RLE back to quality scores
    ######################################################################################
    if options.save_all_quality_scores == True and options.save_exact_quality_scores==False:
        cmd  = f"(/usr/bin/time --verbose "
        cmd += f" {deCompressQualityScoresFile} "
        cmd += f" {name_of_file_with_quality_scores_rle} "
        cmd += f" {name_of_file_with_quality_scores} "
        cmd += f") "
        cmd += f"1> /dev/null "
        if options.keep_intermediate_error_files == True:
            cmd += f"2> {options.error_directory}/{name_of_input_file_without_location}_RLE_decompress.error"
        else:
            cmd += f"2> /dev/null"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
    elif options.save_all_quality_scores == False and options.save_exact_quality_scores==False:
        cmd = f"touch {name_of_file_with_quality_scores}"
        os.system(cmd)

    
    max_reads_in_each_line_filename = f"{options.output_directory}/max_reads_in_each_line_filename"
    output_sam_filename = f"{options.output_directory}/{name_of_input_file_without_location[:-4]}.decompressed.sam"
    files_to_be_removed.extend([_index_outputfilename,
                                _outputfilename,
                                unmapped_outputfilename,
                                name_of_file_with_quality_scores,
                                name_of_file_dictionary,
                                max_reads_in_each_line_filename])

    cmd  = f"(/usr/bin/time --verbose "
    if os.path.exists(name_of_file_dictionary)==False:
        cmd += f" {decompressSamFileSingleEnded} " # argv[0] - Name of the program
    else:
        cmd += f" {decompressSamFilePairedEnded} " # argv[0] - Name of the program
    cmd += f" {_index_outputfilename} "# argv[1] - Abridge index, _index_outputfilename
    cmd += f" {options.genome} "# argv[2] - genome file
    cmd += f" {output_sam_filename} "# argv[3] - decompressed SAM filename
    cmd += f" {_outputfilename} "# argv[4] - Name of the  filename
    cmd += f" dummy "# argv[5] - Dummy
    cmd += f" {options.quality} "# argv[6] - Quality of reads
    cmd += f" {options.ignore_sequence}" #argv[7] - Whether or not to produce sequences from genome file
    cmd += f" {unmapped_outputfilename}" #argv[8] - Name of file with unmapped reads
    cmd += f" {name_of_file_with_quality_scores} " #argv[9] - Name of file with quality scores
    if os.path.exists(name_of_file_dictionary)==True:
        cmd += f" {name_of_file_dictionary} " #argv[10]
        
        cmd_1  = f" {maxReadsInEachLine} "
        cmd_1 += f" {_outputfilename} "
        cmd_1 += f" {max_reads_in_each_line_filename} "
        if options.quiet == False:
            logging.info(f"Running command - {cmd_1}")
        os.system(cmd_1)
        
        max_reads_in_each_line = open(max_reads_in_each_line_filename,"r").readline().strip()
        cmd += f" {max_reads_in_each_line} " #argv[11]
    cmd += f") "   
    #cmd += f"1> /dev/null "
    if options.keep_intermediate_error_files == True:
        cmd += f"2> {options.error_directory}/{name_of_input_file_without_location}_decompress.error"
    else:
        cmd += f"2> /dev/null"
    if options.quiet == False:
        logging.info(f"Running command - {cmd}")
    os.system(cmd)
    
    for file in files_to_be_removed:
        if options.quiet == False:
            logging.info(f"Deleting - {file}")
        os.system("rm -f "+file)
    
    
def generateCoverage(_outputfilename,_index_outputfilename,options):
    generateCoverage = options.softwares["generateCoverage"]
    findMaximumNumberOfReadsInEachLine = options.softwares["findMaximumNumberOfReadsInEachLine"]
    
    ######################################################################################
    # Compiling programs - will be removed during final version
    ######################################################################################
    if options.compile_programs == True:
        cmd=f"gcc {generateCoverage}.c -o {generateCoverage} -g -Ofast"
        os.system(cmd)
        
        cmd=f"gcc {findMaximumNumberOfReadsInEachLine}.c -o {findMaximumNumberOfReadsInEachLine} -g -Ofast"
        os.system(cmd)
    ######################################################################################
    if options.d == True:
        options.d = 1
    else:
        options.d = 0
    if options.bg == True:
        options.bg = 1
    else:
        options.bg = 0
    if options.bga == True:
        options.bga = 1
    else:
        options.bga = 0
    if options.split == True:
        options.split = 1
    else:
        options.split = 0
    if options.generate_overlapping_coverage == True:
        options.generate_overlapping_coverage = 1
    else:
        options.generate_overlapping_coverage = 0
    if options.generate_non_overlapping_coverage == True:
        options.generate_non_overlapping_coverage = 1
    else:
        options.generate_non_overlapping_coverage = 0
    
    name_of_input_file_without_location = options.inputabrfilenames.split("/")[-1].split(".abridge")[0]
    dictionary_name = f"{options.output_directory}/{name_of_input_file_without_location}.sam.dictionary"
    if os.path.exists(dictionary_name) == True:
        single = 0
    else:
        single = 1
    
    
    max_reads_in_each_line_filename = f"{options.output_directory}/{name_of_input_file_without_location}.max_read_in_each_line"
    
    cmd  = f"(/usr/bin/time --verbose "
    cmd += f" {generateCoverage} "
    cmd += f" {_outputfilename} " #argv[1]
    cmd += f" {_index_outputfilename} " #argv[2]
    cmd += f" {options.d} " #argv[3]
    cmd += f" {options.bg} " #argv[4]
    cmd += f" {options.bga} " #argv[5]
    cmd += f" {options.split} " #argv[6]
    cmd += f" {options.generate_overlapping_coverage} " #argv[7]
    cmd += f" {options.generate_non_overlapping_coverage} " #argv[8]
    cmd += f" {single} " #argv[9]
    if os.path.exists(dictionary_name)==True:
        cmd += f" {dictionary_name} " #argv[10]
        cmd_1  = f" {findMaximumNumberOfReadsInEachLine} "
        cmd_1 += f" {_outputfilename} "
        cmd_1 += f" 1> {max_reads_in_each_line_filename} "
        cmd_1 += f" 2> /dev/null "
        if options.quiet == False:
            logging.info(f"Running command - {cmd_1}")
        os.system(cmd_1)
        max_reads_in_each_line = open(max_reads_in_each_line_filename,"r").readline().strip()
        cmd += f" {max_reads_in_each_line} " #argv[11]
    cmd += f") "
    if options.keep_intermediate_error_files == True:
        #cmd += f" 1> {options.error_directory}/{name_of_input_file_without_location}_coverage_generation_{options.d}_{options.bg}_{options.bga}_{options.split}_{options.generate_overlapping_coverage}_{options.generate_non_overlapping_coverage}.output "
        cmd += f" 2> {options.error_directory}/{name_of_input_file_without_location}_coverage_generation_{options.d}_{options.bg}_{options.bga}_{options.split}_{options.generate_overlapping_coverage}_{options.generate_non_overlapping_coverage}.error"
    else:
        #cmd += f" 1> /dev/null "
        cmd += f"2> /dev/null"
    if options.quiet == False:
        logging.info(f"Running command - {cmd}")
    print(cmd)
    os.system(cmd)
    
    
    
    
def retrieveAlignmentsRandomly(chromosome,start,end,options):
    
    randomRetrievalSingleEnded = options.softwares["randomRetrievalSingleEnded"]
    randomRetrievalPairedEnded = options.softwares["randomRetrievalPairedEnded"] 
    deCompressQualityScoresFile = options.softwares["deCompressQualityScoresFile"]
    ######################################################################################
    # Compiling programs - will be removed during final version
    ######################################################################################
    if options.compile_programs == True:
        cmd=f"gcc {randomRetrievalSingleEnded}.c -o {randomRetrievalSingleEnded} -g -Ofast"
        os.system(cmd)
        cmd=f"gcc {randomRetrievalPairedEnded}.c -o {randomRetrievalPairedEnded} -g -Ofast"
        os.system(cmd)
        cmd=f"gcc {deCompressQualityScoresFile}.c -o {deCompressQualityScoresFile} -Ofast"
        os.system(cmd)
    ######################################################################################
    
    genome_filename = options.genome
    genome_index_filename = options.genome+".fai"
    
    input_filename = options.inputabrfilenames
    name_of_input_file_without_location = input_filename.split("/")[-1][:-8]
    _outputfilename = options.output_directory + "/" + name_of_input_file_without_location + ".sam_"
    _index_outputfilename = options.output_directory + "/" + name_of_input_file_without_location + ".sam__index"
    unmapped_outputfilename = options.output_directory + "/" + name_of_input_file_without_location + ".sam_unmapped"
    name_of_file_with_quality_scores = f"{options.output_directory}/{name_of_input_file_without_location}.sam_qual"
    compressed_abridged_filename = input_filename
    output_sam_filename = options.output_directory + "/" + name_of_input_file_without_location + ".decompressed.sam"
    outputfilename = options.output_directory + "/" + name_of_input_file_without_location + "_" + options.positions.replace(":","_") + ".sam"
    files_to_be_removed = []
    
    cmd  = f"(/usr/bin/time --verbose "
    cmd += f" {deCompressQualityScoresFile} "
    cmd += f" {name_of_file_with_quality_scores}.rle "
    cmd += f" {name_of_file_with_quality_scores} "
    cmd += f") "
    cmd += f"1> /dev/null "
    if options.keep_intermediate_error_files == True:
        cmd += f"2> {options.error_directory}/{name_of_input_file_without_location}_RLE_decompress.error"
    else:
        cmd += f"2> /dev/null"
    os.system(cmd)
    files_to_be_removed.append(f"{name_of_file_with_quality_scores}")

    dictionary_name = f"{options.output_directory}/{name_of_input_file_without_location}.sam.dictionary"
    #print(dictionary_name)
    if os.path.exists(dictionary_name)==False:
        cmd  = f" {randomRetrievalSingleEnded} " #argv[0]
    else:
        cmd  = f" {randomRetrievalPairedEnded} " #argv[0]
    cmd += f" {_index_outputfilename} " #argv[1]
    cmd += f" {genome_filename} " #argv[2]
    cmd += f" {genome_index_filename} " #argv[3]
    cmd += f" {chromosome} " #argv[4]
    cmd += f" {_outputfilename}  " #argv[5]
    cmd += f" {start} " #argv[6]
    cmd += f" {end} " #argv[7]
    cmd += f" {outputfilename}" #argv[8]
    cmd += f" {options.quality}" #argv[9[
    cmd += f" {options.ignore_sequence}" #argv[10]
    if options.read_prefix is not None:
        cmd += f" {options.read_prefix}" #argv[11]
    else:
        cmd += f" \"\" " #argv[11]
    cmd += f" {name_of_file_with_quality_scores}" #argv[12]
    if os.path.exists(dictionary_name)==True:
        cmd += f" {dictionary_name} " # argv[13]
    #print(cmd)
    #sys.stdout.flush()
    os.system(cmd)
    
    

def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    #print(options)
    #return
    
    if options.inputsamfilenames is not None:
        input_filename_without_location = options.inputsamfilenames.split("/")[-1][:-4]
    else:
        input_filename_without_location = options.inputabrfilenames.split("/")[-1][:-8]
        
    if options.output_directory is None:
        options.output_directory = "/".join(options.inputsamfilenames.split("/")[:-1])
        
    options.logfilename = f"{options.output_directory}/{input_filename_without_location}_progress.log"
    os.system(f"mkdir -p {options.output_directory}")
    
    if options.quiet == False:
        configureLogger(options)
        
    for line in pprint.pformat(options).split('\n'):
        logging.info(line)
        
    validateCommandLineArguments(options)
    
    for line in pprint.pformat(options).split('\n'):
        logging.info(line)

    if options.quiet == False:
        logging.info("Logger has been configured")
    if options.quiet == False:
        logging.info("validateCommandLineArguments() execution is complete")
    
    if options.compress==True:
        checkSAMAlignments(options,logging)
        if options.quiet == False:
            logging.info(f"convertInputToSAM() execution is complete")

    if options.compress==True:
        findEndedness(options,logging)
        if options.quiet == False:
            logging.info(f"findEndedness() execution is complete")
    
    if options.compress==True:
        compressSamFile(options)
        if options.quiet == False:
            logging.info(f"compressSamFile() execution is complete")
    
    if options.decompress==True:
        decompressSamFile(options)
        if options.quiet == False:
            logging.info(f"decompressSamFile() execution is complete")
    
    if options.header==True:
        """
        cmd=f"7za e {options.inputfilename[0]} -y"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        """
        cmd=f"cat {options.inputfilename[0][:-3]}|grep ^@ > {options.inputfilename[0][:-3]}.headers"
        os.system(cmd)
        #print(open(f"{options.inputfilename[0][:-3]}.headers","r").read())
        os.system(f"rm -f {options.inputfilename[0][:-3]}.headers")
        
    if options.random == True:
        input_filename = options.inputabrfilenames
        name_of_input_file_without_location = input_filename.split("/")[-1][:-8]
        _outputfilename = options.output_directory + "/" + name_of_input_file_without_location + ".sam_"
        _index_outputfilename = options.output_directory + "/" + name_of_input_file_without_location + ".sam__index"
        unmapped_outputfilename = options.output_directory + "/" + name_of_input_file_without_location + ".sam_unmapped"
        name_of_file_with_quality_scores = f"{options.output_directory}/{name_of_input_file_without_location}.sam_qual"
        compressed_abridged_filename = input_filename
        output_sam_filename = options.output_directory + "/" + name_of_input_file_without_location + ".decompressed.sam"
    
    
        cmd  = f"7za e {input_filename} -y "
        cmd += f" -o{options.output_directory}"
        cmd += f" 1> /dev/null 2>/dev/null"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        
        if options.inputabrfilenames[-3:]==".br":
            cmd = f"brotli -d -k -f {_outputfilename}.br {_index_outputfilename}.br {unmapped_outputfilename}.br {name_of_file_with_quality_scores}.rle.br"
            os.system(cmd)
        
        cmd=f"cat {_outputfilename}|grep ^@ > {options.output_directory}/{name_of_input_file_without_location}.headers"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        chromosome_to_length = {}
        fhr=open(f"{options.output_directory}/{name_of_input_file_without_location}.headers","r")
        for line in fhr:
            useless,chromosome,length = line.strip().split("\t")
            chromosome=chromosome.split(":")[-1]
            length=length.split(":")[-1]
            #print(useless,chromosome,length)
            chromosome_to_length[chromosome] = int(length)
        fhr.close()
        if options.positions.count(":")!=1 or options.positions.count("-")!=1:
            print("Error in position. Please specify in chromosome:start-end format")
            sys.exit()
        chromosome = options.positions.split(":")[0]
        start,end = options.positions.split(":")[-1].split("-")
        start,end=int(start),int(end)
        if start>chromosome_to_length[chromosome] or end>chromosome_to_length[chromosome] or start>end:
            print("Error in chromosome lengths")
            sys.exit()
        retrieveAlignmentsRandomly(chromosome,start,end,options)
        
        """
        FUTURE VERSION
        correct,message=verifyPositions(options.positions,chromosome_to_length)
        if correct==0:
            print(f"Incorrect positions - {message}")
        else:
            sortPositions(options.positions)
        """
        os.system(f"rm -rf {options.output_directory}")
    
    if options.generate_overlapping_coverage == True or options.generate_non_overlapping_coverage ==True:
        input_filename = options.inputabrfilenames
        name_of_input_file_without_location = input_filename.split("/")[-1][:-8]
        _outputfilename = options.output_directory + "/" + name_of_input_file_without_location + ".sam_"
        _index_outputfilename = options.output_directory + "/" + name_of_input_file_without_location + ".sam__index"
        
        cmd  = f"7za e {input_filename} -y "
        cmd += f" -o{options.output_directory}"
        cmd += f" 1> /dev/null 2>/dev/null"
        if options.quiet == False:
            logging.info(f"Running command - {cmd}")
        os.system(cmd)
        
        if options.inputabrfilenames[-3:]==".br":
            cmd = f"brotli -d -k -f {_outputfilename}.br {_index_outputfilename}.br"
            os.system(cmd)
        
        generateCoverage(_outputfilename,_index_outputfilename,options)
        os.system(f"rm -rf {options.output_directory}")
    
    cleanUp(options)
    if options.quiet == False:
        logging.info(f"cleanUp() execution is complete")
        
    

if __name__ == "__main__":
    main()