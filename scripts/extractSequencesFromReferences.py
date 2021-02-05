
import os
import glob
import time

def readPass2Index(filename):
    whole_index = {}
    fhr=open(filename,"r")
    for line in fhr:
        chromosome,start,end,start_byte,end_byte=line.strip().split("\t")
        if chromosome not in whole_index:
            whole_index[chromosome]=[]
        whole_index[chromosome].append([start,end])
    fhr.close()
    return whole_index


def extractSequenceFromReferenceFile(reference_filename,chromosome_name,whole_index,output_filename,options,logging):
    cmd=f"rm {output_filename}"
    os.system(cmd)
    
    temp_file_num=1
    cmd=f"samtools faidx {reference_filename} "
    for regions in whole_index[chromosome_name]:
        cmd+=f"{chromosome_name}:{regions[0]}-{regions[1]} "
        if len(cmd)>2**14:
            cmd+=f" > {output_filename}.{temp_file_num}"
            """while (os.path.exists(f"{output_filename}.{temp_file_num}")==False):
                time.sleep(5)"""
            if options.quiet ==False:
                logging.info(f"Running cmd - {len(cmd)} {cmd}")
            os.system(cmd)
            
            cmd=f"cat {output_filename}.{temp_file_num} >> {output_filename}"
            os.system(cmd)
            cmd=f"rm {output_filename}.{temp_file_num}"
            os.system(cmd)
            temp_file_num+=1
            
            cmd=f"samtools faidx {reference_filename} "
    
    if cmd!=f"samtools faidx {reference_filename} ":
        cmd+=f" > {output_filename}.{temp_file_num}"
        if options.quiet ==False:
            logging.info(f"Running cmd - {len(cmd)} {cmd}")
        os.system(cmd)
        """while (os.path.exists(f"{output_filename}.{temp_file_num}")==False):
            time.sleep(5)"""
        cmd=f"cat {output_filename}.{temp_file_num} >> {output_filename}"
        os.system(cmd)
        cmd=f"rm {output_filename}.{temp_file_num}"
        os.system(cmd)
    
    cmd=f"perl -pe '/^>/ ? print \"\\n\" : chomp' {output_filename} | tail -n +2 > {output_filename}.temp"
    os.system(cmd)
    cmd=f"mv {output_filename}.temp {output_filename}"
    os.system(cmd)