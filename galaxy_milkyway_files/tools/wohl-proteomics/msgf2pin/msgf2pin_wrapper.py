######################################
## MSGF2pin Wapper, Wohlschlegel Lab
## Written by William Barshop
## For use in Galaxy-Proteomics Workflows
## where fractionated data is common.
##
version="0.99H"
##Date : 2016-01-13
##
######################################
from optparse import OptionParser
import subprocess
import copy
import pandas
import sys
import pyteomics
from pyteomics import mzid
import os, gc
import pandas
import math
import numpy as np
#from tqdm import *
#from multiprocessing import Pool, cpu_count
#import joblib
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import cpu_count
#os.system("pip install pyteomics")
#sys.exit(2)


#proton_mass=1.00727646677
#proton_mass=1.007825035
#proton_w_electron_mass=1.007825035+0.00054858026 #This is proton + electron mass.
#neutron_mass=1.0086654
isotope_mass_difference=1.00335483

##Option parsing:
##
##

parser = OptionParser()
parser.add_option("-i","--input_files",dest="inputfiles")
parser.add_option("-d","--decoy_files",dest="decoyfiles")
parser.add_option("-m","--num_matches",dest="nummatches")
parser.add_option("-a","--aa-freq",action="store_true",dest="aafreq")
parser.add_option("-p","--PTM",action="store_true",dest="ptm")
parser.add_option("-e","--enzyme",dest="enzyme")
parser.add_option("-n","--PNGaseF",action="store_true",dest="pngasef")
parser.add_option("-s","--isotope",action="store_true",dest="isotope")
parser.add_option("-b","--databases",dest="databases")
parser.add_option("-c","--cleavages",dest="cleavages")
parser.add_option("-o","--outputfolder",dest="outputfolder")
parser.add_option("-f","--fractions",dest="fractions")
parser.add_option("-q","--pattern",dest="pattern")
parser.add_option("-z","--ppm",action="store_true",dest="ppm")
parser.add_option("-x","--fixppm",dest="fix_ppm",type="float") #For this, you pass in an e-value which will be used to fix the ppm by taking the median ppm value of passing PSMs...
options,args=parser.parse_args()

#if not options.fix_ppm is None:
#    print "IT'S NOT NONE!!!!"


inputfiles=options.inputfiles
decoyfiles=options.decoyfiles
nummatches=options.nummatches
aafreq=options.aafreq
ptm=options.ptm
enzyme=options.enzyme
pngasef=options.pngasef
isotope=options.isotope
databases=options.databases
outputfolder=options.outputfolder
fractions=options.fractions
pattern=options.pattern
cleavages=options.cleavages
fix_ppm=options.fix_ppm

#print options.fix_ppm,"FIX PPM",type(options.fix_ppm)
#pause=raw_input("yo")
##Let's read the fractions file to pull out runs and groups.
##
##

#def applyParallel(dfGrouped, func):
#    with Pool(cpu_count()) as p:
#        ret_list = p.map(func, [group for name, group in dfGrouped])
#    return pandas.concat(ret_list)

def applyParallel(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)


def fileNames(group):
    for index,eachrow in group.iterrows():
        #print index
        #group.loc[index,'file']=eachrow['SpecId'].split("_")[0]
        group.loc[index,'file']=eachrow['SpecId'].rsplit("_",7)[0]
        #eachrow['file']=eachrow['SpecId'].split("_")[0]
    return group
#newfilenames=new_df.apply(lambda x: x['SpecId'].split("_")[0],axis=1)




run_groups={} # a dict where key=group name, value=[list,of,runs,in,group]

bio_cond_groups={}

test_control_groups={}

with open(fractions,'rb') as fraction_file:
    for eachline in fraction_file:
        myline=eachline.strip().split('___')
	if len(myline)>1:
            if myline[1] not in run_groups:
                run_groups[myline[1]]=myline[0]
            else:
                print "YOU HAVE DEFINED TWO GROUPS WITH THE SAME IDENTIFYING STRING. THAT'S BAD. DON'T DO THAT UNLESS YOU REALLY MEAN IT FOR SOME STRANGE REASON."
            if len(myline)>2:
                if myline[1] in bio_cond_groups:
                    bio_cond_groups[myline[1]].append(myline[2])
                else:
                    bio_cond_groups[myline[1]]=[myline[2]]
            if len(myline)>3:
                if myline[1] in test_control_groups:
                    print "SHOULDN'T HAPPEN!!!!!!!!!"
                    print "At this point, you have the same group with two target/control definitions."
                    test_control_groups[myline[1]].append(myline[3])
                else:
                    test_control_groups[myline[1]]=[myline[3]]

#print bio_cond_groups,"These are the biological conditions."

target_run_groups={}
decoy_run_groups={}
for each_key in run_groups:
    target_run_groups[each_key]=[]
    decoy_run_groups[each_key]=[]


#decoy_run_groups=copy.deepcopy(run_groups)

with open(inputfiles,'rb') as inputfile_list:
    print "Reading input file list..."
    for eachline in inputfile_list:
        file_name=eachline.strip()
        #print file_name,"THIS IS FILE NAME"
        for eachgroup in target_run_groups:
            #print eachgroup,file_name,"eachgroup and filename"
            if run_groups[eachgroup] in file_name:
                target_run_groups[eachgroup].append(file_name)
                break

with open(decoyfiles,'rb') as decoyfile_list:
    print "Reading decoy file list..."
    for eachline in decoyfile_list:
        file_name=eachline.strip()
        for eachgroup in decoy_run_groups:
            if run_groups[eachgroup] in file_name:
                
                decoy_run_groups[eachgroup].append(file_name)
                break
    
#print decoy_run_groups, "THESE ARE RUN GROUPS<---------------"

print "building commands..."
processes = []
msgf2pin_cmd="msgf2pin "+"-m "+nummatches+" --enzyme "+enzyme+" --pattern "+pattern

if ptm:
    msgf2pin_cmd+=" --PTM"
if aafreq:
    msgf2pin_cmd+=" --aa-freq"
if pngasef:
    msgf2pin_cmd+=" --PNGaseF"
if isotope:
    msgf2pin_cmd+=" --isotope"
if databases:
    msgf2pin_cmd+=" --databases "+databases

#print "Here are out run groups",run_groups


output_list=[]
print "Taking care of a few file creation issues...."
print "Here's our groups...",run_groups
for each_group in run_groups:
    with open(run_groups[each_group]+".runs",'wb') as writer:
        #print target_run_groups[each_group],"<------------------------------"
        for each_run in target_run_groups[each_group]:
            writer.write(each_run+"\n")
    with open(run_groups[each_group]+".decoys",'wb') as writer:
        for each_run in decoy_run_groups[each_group]:
            writer.write(each_run+"\n")
            #print each_run,"<-----------------------"
    

    #DEVNULL = open(os.devnull, 'wb')
    

    

    tmp_msgf2pin_cmd=msgf2pin_cmd
    tmp_msgf2pin_cmd+=" --outputTab "+outputfolder+run_groups[each_group]+".tsv"
    output_list.append(outputfolder+run_groups[each_group]+".tsv")
    tmp_msgf2pin_cmd+=" "+run_groups[each_group]+".runs"
    tmp_msgf2pin_cmd+=" "+run_groups[each_group]+".decoys"
    tmp_msgf2pin_cmd+=" 2>&1"
    print "RUNNING COMMAND",tmp_msgf2pin_cmd
    #processes.append(subprocess.Popen(tmp_msgf2pin_cmd,shell=True))
    processes.append(subprocess.Popen(tmp_msgf2pin_cmd,shell=True))

print "Starting MSGF2PIN...."
for p in processes:
    print p.wait(),"finished..."

if options.ppm:
    if options.fix_ppm is not None:
        correction_df_list=[]
    print "We're going to add ppm mass accuracy measurements to each identification."
    mzid_dict={}
    ppm_correction_list=[]

    for eachfile in output_list:
        del mzid_dict
        gc.collect()
        mzid_dict={}

        #file_line_ctr=0
        #with open(eachfile,'rb') as filereader:
        #    for eachline in filereader:
        #        file_line_ctr+=1

        with open(eachfile+"_ppm",'wb') as filewriter:
            with open(eachfile,'rb') as filereader:
                linectr=0
                peptide_index=0
                expmass_index=0
                mass_index=0
                index_ctr=0
                specID_index=0
                scan_index=0
                label_index=0
                isotope_error_index=0
                charge1_index=-1
                charge2_index=-1
                charge3_index=-1
                charge4_index=-1
                charge5_index=-1
                charge6_index=-1
                charge7_index=-1
                charge8_index=-1
                new_reconstructed_lines=[]
                print "Now for each line..."
                for eachline in filereader:
                    if linectr==0:
                        header=eachline.split("\t")
                        for each_item in header:
                            if each_item == "Peptide":
                                peptide_index=index_ctr
                            elif each_item == "ExpMass":
                                expmass_index=index_ctr
                            elif each_item == "CalcMass":
                                mass_index=index_ctr
                            elif each_item == "Charge1":
                                charge1_index=index_ctr
                            elif each_item == "Charge2":
                                charge2_index=index_ctr
                            elif each_item == "Charge3":
                                charge3_index=index_ctr
                            elif each_item == "Charge4":
                                charge4_index=index_ctr
                            elif each_item == "Charge5":
                                charge5_index=index_ctr
                            elif each_item == "Charge6":
                                charge6_index=index_ctr
                            elif each_item == "Charge7":
                                charge7_index=index_ctr
                            elif each_item == "Charge8":
                                charge8_index=index_ctr
                            elif each_item == "SpecId":
                                specID_index=index_ctr
                            elif each_item == "ScanNr":
                                scan_index=index_ctr
                            elif each_item == "Label":
                                label_index=index_ctr
                            elif each_item == "IsotopeError":
                                isotope_error_index=index_ctr

                            index_ctr+=1
                            
                        #peptide_index+=1
                        #print "this should be peptide...",header[peptide_index]
                        
                        peptide_start=header[:peptide_index]
                        peptide_start.append("ppm")
                        peptide_start.append("absppm")
                        for each_item_pep in header[peptide_index:]:
                            peptide_start.append(each_item_pep)
                        #print "THE NEW HEADER\n",peptide_start
                        new_reconstructed_lines.append(peptide_start)

                    elif linectr ==1:
                        directions=eachline.rstrip().split("\t")
                        start_directions=directions[:peptide_index]
                        #final_dir=start_directions[-1:]
                        #start_directions=start_directions[:-1]
                        #start_directions.append(final_dir[0])
                        start_directions.append("0")
                        start_directions.append("-1")
                        for each_item_dir in directions[peptide_index:]:
                            start_directions.append(each_item_dir)
                        start_directions.append("\n")
                        new_reconstructed_lines.append(start_directions)
                    else:
                        #Okay, so this is where the magic is actually going to happen.
                        #1. for each line, we'll pull the ExpMass and TheoMass (ie Mass) and figure out the CHARGE.
                        #2. we'll calculate ExpM/Z and TheoM/Z
                        #3. we'll finally use those to calculate ppm and abs(ppm) and cram them into the new row
                        #4. finally, we'll add that new row to the new_reconstructed lines
                        split_line = eachline.split("\t")
                        specID=split_line[specID_index]
                        label=split_line[label_index]
                        isotope_error=int(split_line[isotope_error_index])
                        if label == "1":
                            file_name_ppm=specID.rsplit("_",7)[0]+"_target.mzid"
                        else:
                            file_name_ppm=specID.rsplit("_",7)[0]+"_decoy.mzid"
                        scan_num=split_line[scan_index]

                        #print "Label was ",str(label),"so we will be opening the file",file_name_ppm,"from ",specID
                        if not file_name_ppm in mzid_dict:
                            reader=pyteomics.mzid.read(file_name_ppm)
                            mzid_dict[file_name_ppm]={}
                            #print "We'll have to iterate over the mass spec data to get high res masses..."
                            print "Reading in file",file_name_ppm,"..."
                            for each_scan in reader:
                               #print each_scan
                               #mzid_dict[file_name_ppm][str(int(each_scan['scan number(s)']))]=each_scan                                      ##### OLD: Stored entire mzid scan information in dictionaries... filled memory...
                               mzid_dict[file_name_ppm][str(int(each_scan['scan number(s)']))]=each_scan['SpectrumIdentificationItem'][0]
                            del reader
                            gc.collect()
                        this_id=mzid_dict[file_name_ppm][scan_num]
                        theo_mz=this_id['calculatedMassToCharge']
                        exp_mz=this_id['experimentalMassToCharge']
                        #theo_mz=this_id['SpectrumIdentificationItem'][0]['calculatedMassToCharge']                                         ##### Instead, now we only store the BEST ID from each scan...
                        #exp_mz=this_id['SpectrumIdentificationItem'][0]['experimentalMassToCharge']                                        ##### If needed, we can change that later to be only the exp mz and theo mz...
                        msgf2pin_exp_mz=split_line[expmass_index]


                        #scan_num=int(specID.split("_")[3])
                        #### DETERMINE CHARGE STATE HERE
                        ##
                        ##
                        if charge1_index != -1:
                            charge1=(int(split_line[charge1_index])==1)
                        else:
                            charge1=False
                        if charge2_index != -1:
                            charge2=(int(split_line[charge2_index])==1)
                        else:
                            charge2=False
                        if charge3_index != -1:
                            charge3=(int(split_line[charge3_index])==1)
                        else:
                            charge3=False
                        if charge4_index != -1:
                            charge4=(int(split_line[charge4_index])==1)
                        else:
                            charge4=False
                        if charge5_index != -1:
                            charge5=(int(split_line[charge5_index])==1)
                        else:
                            charge5=False
                        if charge6_index != -1:
                            charge6=(int(split_line[charge6_index])==1)
                        else:
                            charge6=False
                        if charge7_index != -1:
                            charge7=(int(split_line[charge7_index])==1)
                        else:
                            charge7=False
                        if charge8_index != -1:
                            charge8=(int(split_line[charge8_index])==1)
                        else:
                            charge8=False
                        charge_array=[charge1,charge2,charge3,charge4,charge5,charge6,charge7,charge8]
                        mycharge=1
                        for each_charge in charge_array:
                            if each_charge:
                                break
                            else:
                                mycharge+=1
                        #exp_mz=exp_mass/float(mycharge)
                        #theo_mz=theo_mass/float(mycharge)
                        #print "exp massindex ",str(expmass_index)
                        #print "massindex",str(mass_index)
                        #print "exp mass",str(exp_mass)
                        #print "theo mass",str(theo_mass)
                        #print "MSGF2PIN exp mz",str(msgf2pin_exp_mz)
                        #print "MZID exp mz",str(exp_mz)
                        #print "theo mz",str(theo_mz)
                        #ppm = (((exp_mz-((isotope_error*proton_mass)/mycharge)-theo_mz)/theo_mz)*1000000)
                        ppm = (((exp_mz-((isotope_error*isotope_mass_difference)/mycharge)-theo_mz)/theo_mz)*1000000)
                        #ppm = (((exp_mz-theo_mz)/theo_mz)*1000000)
                        absppm=abs(ppm)
                        #print "ppm and absppm",str(ppm),str(absppm)
                        new_line=split_line[:peptide_index]
                        new_line.append(str(ppm))
                        new_line.append(str(absppm))
                        for each_item_line in split_line[peptide_index:]:
                            new_line.append(each_item_line)
                        new_reconstructed_lines.append(new_line)
                    linectr+=1

                if options.fix_ppm is None:
                    for each_newline in new_reconstructed_lines:
                        filewriter.write("\t".join(each_newline))

                else:#We're doing the PPM correction!
                    for each_newline in new_reconstructed_lines:
                        filewriter.write("\t".join(each_newline[:peptide_index+2])+"\n") #We're leaving out the proteins, etc...
                    #for_pandas=new_reconstructed_lines[:]
                    #directions=for_pandas[1]
                    #del for_pandas[1]
                    #for_pandas_fixed=[f.split("\t") for f in for_pandas]
                    #new_df = pandas.DataFrame(for_pandas)
                    new_df=pandas.read_csv(eachfile+"_ppm",sep="\t",skiprows=[1])
                    #print new_df
                    #print new_df.columns,"cols..."
                    new_df['file']=""
                    new_df=new_df[np.logical_and(new_df['lnEValue'] >= (-1.0*math.log(fix_ppm)),new_df['Label']==1)]
                    #print new_df,"filtered df!"
                    #newfilenames=new_df.apply(lambda x: x['SpecId'].split("_")[0],axis=1)
                    grouped_df=new_df.groupby(np.arange(len(new_df))//cpu_count())
                    #print grouped_df,"grouped!"
                    newfnames=applyParallel(grouped_df,fileNames)
                    #!!!!!!!!!!!!!!!!!!!

                    #new_df['file']=newfilenames
                    grouped_df=newfnames.groupby(newfnames.file)
                    corrected_ppm=grouped_df['ppm'].median()
                    ppm_correction_list.append(corrected_ppm)
                    print corrected_ppm,"aggggg"
                    print corrected_ppm.index,"index!"
                    #for index,each in aggregated:
                    #    print each,"each one..."
                    #print type(corrected_ppm),"this is aggregated type..."
                    #print corrected_ppm.iloc[0],"this is position iloc 0..."
                    #correction_df_list.append(corrected_ppm)#We can concat these to make a correction map....
                    #grouped_df.to_csv("grouped.csv")
                    os.remove(eachfile+"_ppm")
                    #new_df.to_csv("newdf.tsv",sep="\t")

                    line_ctr=0
                    with open(eachfile+"_ppm_fixed",'wb') as ppmfixed_file:
                        for each_newline in new_reconstructed_lines:
                            working_line=each_newline[:]
                            #before_file.write("\t".join(working_line))    
                            if line_ctr>=2:
                                #print "before...",working_line
                                #print peptide_index,"this is peptide index..."
                                #print absppm_index,"this is absppm index..."
                                #print ppm_index,"this is ppm index..."
                                print working_line[ppm_index],"ppm before..."
                                working_line[ppm_index]=str(float(working_line[ppm_index])-float( corrected_ppm.loc[working_line[specID_index].rsplit("_",7)[0]] ))
                                #print working_line[ppm_index],"ppm after..."
                                #print "after...",working_line
                                working_line[absppm_index]=str(abs(float(working_line[ppm_index])))
                            elif line_ctr==0:
                                header=working_line
                                index_ctr=0
                                for each_item in header:
                                    if each_item == "Peptide":
                                        peptide_index=index_ctr
                                    elif each_item == "ExpMass":
                                        expmass_index=index_ctr
                                    elif each_item == "CalcMass":
                                        mass_index=index_ctr
                                    elif each_item == "Charge1":
                                        charge1_index=index_ctr
                                    elif each_item == "Charge2":
                                        charge2_index=index_ctr
                                    elif each_item == "Charge3":
                                        charge3_index=index_ctr
                                    elif each_item == "Charge4":
                                        charge4_index=index_ctr
                                    elif each_item == "Charge5":
                                        charge5_index=index_ctr
                                    elif each_item == "Charge6":
                                        charge6_index=index_ctr
                                    elif each_item == "Charge7":
                                        charge7_index=index_ctr
                                    elif each_item == "Charge8":
                                        charge8_index=index_ctr
                                    elif each_item == "SpecId":
                                        specID_index=index_ctr
                                    elif each_item == "ScanNr":
                                        scan_index=index_ctr
                                    elif each_item == "Label":
                                        label_index=index_ctr
                                    elif each_item == "IsotopeError":
                                        isotope_error_index=index_ctr
                                    elif each_item == "ppm":
                                        ppm_index=index_ctr
                                    elif each_item == "absppm":
                                        absppm_index=index_ctr

                                    index_ctr+=1

                            ppmfixed_file.write("\t".join(working_line))
                            line_ctr+=1

                    #print aggregated
                    #print "I\'M HERE!!!"
                    #new_df.to_csv(eachfile+"_ppm")#CHANGE THIS TO WRITE THE FIXED FILE...


    for eachfile in output_list:
        os.remove(eachfile)
        if options.fix_ppm is None:
            print "renaming....",eachfile
            os.rename(eachfile+"_ppm",eachfile)
        else:
            print "renaming....",eachfile
            os.rename(eachfile+"_ppm_fixed",eachfile)
        #with open(eachfile+"_ppm",'wb') as filewriter:
    if options.fix_ppm:
        all_corrected_ppm=pandas.concat(ppm_correction_list)
        print "Writing out the ppm corrections per file for later use..."
        #print all_corrected_ppm
        all_corrected_ppm.index=all_corrected_ppm.index+".mzML"
        all_corrected_ppm.to_csv("ppm_shifts.csv",sep=',')


print "We have finished running MSGF2PIN a bunch of times, the output should be in "+outputfolder
