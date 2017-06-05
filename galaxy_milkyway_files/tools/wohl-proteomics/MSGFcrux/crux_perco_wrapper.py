version=".75A"
date="12/21/2015"


import sys, os, subprocess
import pandas
from os import listdir
from os.path import isfile, join

arguments=sys.argv
print arguments
arguments=arguments[1:]
folder=arguments[0]
arguments=arguments[1:]
exp_group_file=arguments[0]
arguments=arguments[1:]


onlyfiles = [ f for f in listdir(folder) if isfile(join(folder,f)) ]

startingdir=os.getcwd()

#cmd = " ".join(arguments)
#print " ".join(arguments)
processes=[]

os.mkdir("output")
#print onlyfiles,"onlyfiles"
for each in onlyfiles:
    #print each
    cmd=arguments[:]
    cmd.append("--fileroot")
    cmd.append(each.split(".")[0])
    cmd.append(each)
    os.mkdir("output/"+each+"_out")
    os.rename(folder+"/"+each,"output/"+each+"_out/"+each)
    os.chdir("output/"+each+"_out")
    print "---------------------------"
    print cmd,"this will run"
    print "---------------------------"
    processes.append(subprocess.Popen(cmd))
    os.chdir(startingdir)
for each in processes:
    each.wait()



group_information = pandas.read_csv(exp_group_file,sep='\t')

run_dict={} # Key is file_idx, value is file_name.mzML
group_to_run_dict={} # Key is group, value is [1, 2, 3, 4] list of file_idx belonging to runs in the group...
group_to_file_name={} # key is group, value is ["xxxx.mzML", "xxxx.mzML"]
for index,row in group_information.iterrows():
    run_dict[str(row['Crux File Integer'])]=row['Original File Name']+".mzML"
    if row['Fractionation Group ID String'] in group_to_run_dict:
        group_to_run_dict[row['Fractionation Group ID String']].append(str(row['Crux File Integer']))
    else:
        group_to_run_dict[row['Fractionation Group ID String']] = [str(row['Crux File Integer'])]
    if row['Fractionation Group ID String'] in group_to_file_name:
        group_to_file_name[row['Fractionation Group ID String']].append(str(row['Original File Name'])+".mzML")
    else:
        group_to_file_name[row['Fractionation Group ID String']] = [str(row['Original File Name'])+".mzML"]




################ To prepare for percolator file output, we're going to read in the percolator pin files
################ Once they're read in, we'll have to generate a unique ID of file_idx+"_"+scan to merge the pin sequence on
################ We're only going to bring in the sequence from the percolator input file.
################ ===== As of crux version 3.1, we will also begin to fix the charge column...
os.chdir(startingdir)
pin_list=[]
for each in onlyfiles:
    os.chdir("output/"+each+"_out")
    with open("temp.tsv",'wb') as filewriter:
        with open(each,'rb') as filereader:
            linectr=0
            peptide_index=0
            index_ctr=0
            specID_index=0
            scan_index=0
            label_index=0
            charge_indicies={}
            pin_lines=[]
            for eachline in filereader:
                if linectr==0:
                    header=eachline.split("\t")
                    for each_item in header:
                        if each_item == "Peptide":
                            peptide_index=index_ctr
                        elif each_item == "SpecId":
                            specID_index=index_ctr
                        elif each_item == "ScanNr":
                            scan_index=index_ctr
                        elif each_item == "Label":
                            label_index=index_ctr
                        elif "Charge" in each_item:
                            charge_indicies[each_item]=index_ctr
                        index_ctr+=1
                    filewriter.write("SpecId\tlabel\tScanNr\tPeptide\tcharge\n")
                elif linectr==1:
                    pass
                else:
                    thisline=eachline.split("\t")
                    peptide = thisline[peptide_index].split(".",1)[1]
                    peptide = peptide.rsplit(".",1)[0]
                    charge=0
                    for each_charge in charge_indicies:
                        if thisline[charge_indicies[each_charge]]=="1":
                            charge=each_charge[-1]
                    if charge == 0:
                        print "We couldnt find a charge for "+str(thisline[specID_index])+" so its defaulting to zero!"
                    filewriter.write(thisline[specID_index]+"\t"+thisline[label_index]+"\t"+thisline[scan_index]+"\t"+peptide+"\t"+str(charge)+"\n")
                linectr+=1

    new_in=pandas.read_csv("temp.tsv",sep='\t')
    os.remove("temp.tsv")
    #print new_in,"new csv..."
    #new_in['label']=new_in['label'].astype(int)
    new_in=new_in[new_in['label']==1]
    pin_list.append(new_in)
    os.chdir(startingdir)

#target_1_58_58_2
def makeUniqueID(x):
    return x['SpecId'].split("_")[1]+"_"+str(x['ScanNr'])

def makeUniqueIDpsms(x):
    return str(x['file_idx'])+"_"+str(x['scan'])

pin_megaframe=pandas.concat(pin_list)
unique_names=pin_megaframe.apply(makeUniqueID,axis=1)
pin_megaframe['unique_name']=unique_names
pin_megaframe.drop('SpecId',axis=1,inplace=True)
pin_megaframe.drop('ScanNr',axis=1,inplace=True)
pin_megaframe.drop('label',axis=1,inplace=True)


os.chdir(startingdir)
for each in onlyfiles:
    os.chdir("output/"+each+"_out")
    target_psms=[ f for f in listdir("crux-output/") if (isfile(join("crux-output/",f)) and "target.psms.txt" in f)]
    decoy_psms=[ f for f in listdir("crux-output/") if (isfile(join("crux-output/",f)) and "decoy.psms.txt" in f)]
    #target_peptides=[ f for f in listdir("crux-output/") if (isfile(join("crux-output/",f)) and "target.peptides.txt" in f)]
    #decoy_peptides=[ f for f in listdir("crux-output/") if (isfile(join("crux-output/",f)) and "decoy.peptides.txt" in f)]
    #os.chdir("crux-output")
    corrected=[]
    for eachfile in target_psms:
        print "About to read in ",eachfile
        df=pandas.read_csv("crux-output/"+eachfile,sep='\t')
        df['file_idx']=df['file_idx'].astype(int)
        #print df,"this is the dataframe! <--------------------------"
        fileset=set(df['file_idx'])
        for eachone in fileset:
            eachint=int(eachone)
            mask=df[df['file_idx'] == eachone]
            mask['file']=run_dict[str(eachone)]
            #print run_dict[str(eachone)]
            corrected.append(mask)
            #print mask,"this is mask for",str(eachone)
    #print "I was handling ",each
    #print corrected,"the frame itself"
    corrected_df=pandas.concat(corrected)
    corrected_df.drop('charge',axis=1,inplace=True)
    unique_names_psms=corrected_df.apply(makeUniqueIDpsms,axis=1)
    corrected_df['unique_name']=unique_names_psms
    merged_df=corrected_df.merge(pin_megaframe,how='inner',on='unique_name')#
    corrected_df=merged_df
    corrected_df.drop('unique_name',axis=1,inplace=True)
    corrected_df.drop('sequence',axis=1,inplace=True)
    corrected_df.rename(columns={'Peptide':'sequence'},inplace=True)
    #HERE WE CAN MAKE UNIQUE IDS
    #AND THEN DO A LEFT JOIN TO TAKE IN THE OLD (good) SEQUENCES
    #AND THEN GET RID OF THE BAD SEQUENCES AND REPLACE THEM WITH THE GOOD
    #SEQUENCES.
    corrected_df.to_csv("crux-output/"+eachfile,sep='\t',index=False)
    corrected_decoys=[]
    for eachfile in decoy_psms:
        df=pandas.read_csv("crux-output/"+eachfile,sep='\t')
        df['file_idx']=df['file_idx'].astype(int)
        #print df,"this is the dataframe!"
        fileset=set(df['file_idx'])
        for eachone in fileset:
            eachint=int(eachone)
            mask=df[df['file_idx'] == eachone]
            mask['file']=run_dict[str(eachone)]
            corrected_decoys.append(mask)
            #print mask,"this is mask for",str(eachone)
    corrected_df=pandas.concat(corrected_decoys)
    corrected_df.to_csv("crux-output/"+eachfile,sep='\t',index=False)
    #this_run=



    #onlyfiles = [ f for f in listdir(folder) if (isfile(join(folder,f)) and ".txt" in f)]
    os.chdir(startingdir)
