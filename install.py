#! /usr/bin/env python

import os
import glob

"""
Check the bashrc file to check if finder was previously installed
"""
temp_file=os.path.expanduser('~')+"/"+ '.bashrc.temp'
actual_file=os.path.expanduser('~')+"/"+ '.bashrc'

abridge_installation_locations=[]
fhr=open(os.path.expanduser('~')+"/"+ '.bashrc',"r")
for line in fhr:
    if "export" in line and "Abridge" in line:
        abridge_installation_locations.append(line.strip().split("$PATH:")[-1])
fhr.close()

#print(abridge_installation_locations)
remove_these_indices=[]
for i,each_installation in enumerate(abridge_installation_locations):
    check_this_file="software_identity"
    verify_contents="abridge: compress alignments to a reference for improved storage"
    #print(each_installation+"/"+check_this_file
    if os.path.exists(each_installation+"/"+check_this_file)==True:
        #print(open(each_installation+"/"+check_this_file).read().split("\n")[0])
        #print(verify_contents)
        if verify_contents != open(each_installation+"/"+check_this_file).read().split("\n")[0]:
            remove_these_indices.append(i)
        else:
            print("Match found")
    else:
        print("File not found")
        remove_these_indices.append(i)

#print(remove_these_indices)
for i in remove_these_indices[::-1]:
    abridge_installation_locations.pop(i)
    

#print(abridge_installation_locations)
fhw = open(os.path.expanduser('~')+"/"+ '.bashrc.temp',"w")
fhr=open(os.path.expanduser('~')+"/"+ '.bashrc',"r")
for line in fhr:
    if "export" in line and "Abridge" in line:
        if line.strip().split("$PATH:")[-1] in abridge_installation_locations:
            loc=line.strip().split("$PATH:")[-1] 
            print(f"{loc} will be removed from $PATH")
            continue
        else:
            fhw.write(line)
    else:
        fhw.write(line)

pwd = os.getcwd()
fhw.close()

cmd=f"mv {temp_file} {actual_file}"
os.system(cmd)

os.chdir(pwd+"/src/")
os.system("make")
os.chdir(pwd)








