######################################
## MSGF2pin+SILAC Wapper, Wohlschlegel Lab
## Written by William Barshop
## For use in Galaxy-Proteomics Workflows
## where fractionated data is common.
##
version="1.06"
##Date : 2017-07-22
##
######################################
from optparse import OptionParser
import subprocess
import copy
import pandas
import sys
import pyteomics
from pyteomics import mzid
import os, gc, shutil
import math
import numpy as np
import glob
#from tqdm import *
#from multiprocessing import Pool, cpu_count
#import joblib
import multiprocessing
from joblib import Parallel, delayed
from multiprocessing import cpu_count
from itertools import islice

#os.system("pip install pyteomics")
#sys.exit(2)


proton_mass=1.00727646677
#proton_mass=1.007825035
#proton_w_electron_mass=1.007825035+0.00054858026 #This is proton + electron mass.
#neutron_mass=1.0086654
isotope_mass_difference=1.00335483  #This is for C13!
max_concurrency = 8 # MAXIMUM NUMBER OF PARALLEL MSGF2PIN RUNS.
##Option parsing:
##
##

parser = OptionParser()
parser.add_option("-t","--target_folder",dest="targetfolder")
parser.add_option("-d","--decoy_folder",dest="decoyfolder")
parser.add_option("-m","--num_matches",dest="nummatches")
parser.add_option("-a","--aa-freq",action="store_true",dest="aafreq")
parser.add_option("-p","--PTM",action="store_true",dest="ptm")
parser.add_option("-l","--SILAC",action="store_true",dest="silac")
parser.add_option("-e","--enzyme",dest="enzyme")
parser.add_option("-n","--PNGaseF",action="store_true",dest="pngasef")
parser.add_option("-s","--isotope",action="store_true",dest="isotope")
parser.add_option("-b","--databases",dest="databases")
parser.add_option("-c","--cleavages",dest="cleavages")
parser.add_option("-o","--outputfolder",dest="outputfolder")
parser.add_option("-f","--fractions",dest="fractions")
parser.add_option("-q","--pattern",dest="pattern")
parser.add_option("-z","--ppm",action="store_true",dest="ppm",default=False)
parser.add_option("-u","--diaumpire",action="store_true",dest="diaumpire") #Data comes from DIAUmpire, and we will need to correct the SpecID and Scan # to represent the new mzML indicies.
parser.add_option("-x","--fixppm",dest="fix_ppm",type="float") #For this, you pass in an e-value which will be used to fix the ppm by taking the median ppm value of passing PSMs...
options,args=parser.parse_args()

#if not options.fix_ppm is None:
#    print "IT'S NOT NONE!!!!"

silac=options.silac
target_folder=options.targetfolder
decoyfolder=options.decoyfolder
if not silac:
    nummatches=options.nummatches
else:
    nummatches="15"
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

def applyParallelQuarter(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count()/4)(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)

def writeSILACpin(df,target_location,column_names,header):
    with open(target_location,'wb') as file_writer:
        file_writer.write('\t'.join(column_names))
        file_writer.write('\t'.join(header))
        for index,eachrow in df.iterrows():
            thisline=eachrow['data']
            thisline[0]=eachrow['specid']
            file_writer.write('\t'.join(thisline))
    

def renumber_scans(group):
    i=1
    for index,eachrow in group.iterrows():
        original_file=eachrow['specid'].rsplit("_",6)[0]+".mzid"
        if "_decoy" in original_file:
            original_file="decoys/"+original_file
        else:
            original_file="targets/"+original_file
        group.loc[index,'original_file']=original_file
        group.loc[index,'specid']=eachrow['specid'].rsplit("_",1)[0].split("_",1)[1]+"_"+str(i)
        i+=1
    return group


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
groups_to_mzml={}
mzML_files=glob.glob("*.mzML")
mzid_files=[x.replace(".mzML",".mzid") for x in mzML_files]

with open(fractions,'rb') as fraction_file:
    for eachline in fraction_file:
        myline=eachline.strip().split('___')
        print myline
	if len(myline)>1:
            if myline[0] not in groups_to_mzml:
                print "Looking for",myline[0]
                groups_to_mzml[myline[0]]=glob.glob("*"+myline[0]+"*.mzML")

            if myline[1] not in run_groups:
                run_groups[myline[1]]=[myline[0]]
                #print "joining group ",myline[1]
            else:
                run_groups[myline[1]].append(myline[0])
                #print "YOU HAVE DEFINED TWO GROUPS WITH THE SAME IDENTIFYING STRING. THAT'S BAD. DON'T DO THAT UNLESS YOU REALLY MEAN IT FOR SOME STRANGE REASON."
                print "You should only see this message in the case of fractionation..."
            if len(myline)>2:
                if myline[1] in bio_cond_groups:
                    if myline[2] not in bio_cond_groups[myline[1]]:
                        bio_cond_groups[myline[1]].append(myline[2])
                    else:
                        print "I didn't add",myline[2],"because it was already in the bio cond group for",myline[1]
                else:
                    bio_cond_groups[myline[1]]=[myline[2]]
            if len(myline)>3:
                if myline[1] in test_control_groups:
                    print "SHOULDN'T HAPPEN!!!!!!!!!"
                    print "At this point, you have the same group with two target/control definitions."
                    test_control_groups[myline[1]].append(myline[3])
                else:
                    test_control_groups[myline[1]]=[myline[3]]

print bio_cond_groups,"These are the biological conditions."
print run_groups,"THESE ARE THE RUN GROUPS..."
print groups_to_mzml,"these are the groups<--->mzML"

target_run_groups={}
decoy_run_groups={}
for each_key in run_groups:
    target_run_groups[each_key]=[]
    decoy_run_groups[each_key]=[]


#decoy_run_groups=copy.deepcopy(run_groups)


target_files = glob.glob(os.path.join("targets/",'*.mzid'))
decoy_files = glob.glob(os.path.join("decoys/",'*.mzid'))


#with open(inputfiles,'rb') as inputfile_list:
#    print "Reading input file list..."

output_list=[]
output_files={}

if silac:
    nummatches="15" #This would allow top 3 from 5 labels ... Should work for a while!

msgf2pin_cmd="msgf2pin "+"-m "+nummatches+" --enzyme "+enzyme+" --pattern "+pattern
processes = []
#msgf2pin_cmd="msgf2pin "+"-m "+nummatches+" --enzyme "+enzyme+" --pattern "+pattern
print "building commands..."

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

if not silac:
    for eachfile in target_files:
        #print "we're looking for",eachfile
        for eachgroup in target_run_groups:
            this_group_list=run_groups[eachgroup]
            #print this_group_list,"the group list<-=------------------"
            for each_item in this_group_list:
                #print "examining",each_item
                #print "Is it in",eachfile,"?"
                #if run_groups[eachgroup] in eachfile:
                if each_item in eachfile:
                    target_run_groups[eachgroup].append(eachfile)
                    break
    print "Final Groupings:",target_run_groups

    for eachfile in decoy_files:
        for eachgroup in target_run_groups:
            this_group_list=run_groups[eachgroup]
            for each_item in this_group_list:
                if each_item in eachfile:
                    decoy_run_groups[eachgroup].append(eachfile)
                    break


            #if run_groups[eachgroup] in eachfile:
            #    decoy_run_groups[eachgroup].append(eachfile)
            #    break

    print target_run_groups,"TARGET RUN GROUPS... the ones you need to pay attention to!"

    print decoy_run_groups,"DECOY RUN GROUPS... and this one too!"
    #sys.exit(0)




if silac:

    silac_labels=[]
    file_groups={}
    individual_output_list=[]
    print target_files,"Target files..."
    print decoy_files,"Decoy files..."


    for thisfile in target_files:
        eachfile=thisfile.split("/",1)[1]
        print eachfile,"handling..."
        if eachfile.split("_",1)[0] not in silac_labels:
            silac_labels.append(eachfile.split("_",1)[0])
        if eachfile.split("_",1)[1].replace("_target","") not in file_groups:
            file_groups[eachfile.split("_",1)[1].replace("_target","")]={'targets':[thisfile],'decoys':[]}
        else:
            file_groups[eachfile.split("_",1)[1].replace("_target","")]['targets'].append(thisfile)
    for thisfile in decoy_files:
        eachfile=thisfile.split("/",1)[1]
        print eachfile,"handling..."
        if eachfile.split("_",1)[0] not in silac_labels:
            silac_labels.append(eachfile.split("_",1)[0])
        if eachfile.split("_",1)[1].replace("_decoy","") not in file_groups:
            file_groups[eachfile.split("_",1)[1].replace("_decoy","")]={'decoys':[eachfile],'targets':[]}
        else:
            file_groups[eachfile.split("_",1)[1].replace("_decoy","")]['decoys'].append(thisfile)
    
    for eachgroup in file_groups:
        print "LOOKING IN GROUP",eachgroup
        os.mkdir(eachgroup.replace(".mzid",""))
        targets=file_groups[eachgroup]['targets']
        print targets,"targets..."
        decoys=file_groups[eachgroup]['decoys']
        print decoys,"decoys..."
        with open(eachgroup.replace(".mzid","")+"_targets",'wb') as writer:
            for eachfile in targets:
                writer.write(eachfile+"\n")
        with open(eachgroup.replace(".mzid","")+"_decoys",'wb') as writer:
            for eachfile in decoys:
                writer.write(eachfile+"\n")

        tmp_msgf2pin_cmd=msgf2pin_cmd
        tmp_msgf2pin_cmd+=" --outputTab "+eachgroup.replace(".mzid","")+"/"+eachgroup.replace(".mzid","")+".tsv"
        output_files[eachgroup]=eachgroup.replace(".mzid","")+"/"+eachgroup.replace(".mzid","")+".tsv"
        individual_output_list.append(eachgroup.replace(".mzid","")+"/"+eachgroup.replace(".mzid","")+".tsv")
        tmp_msgf2pin_cmd+=" "+eachgroup.replace(".mzid","")+"_targets"
        tmp_msgf2pin_cmd+=" "+eachgroup.replace(".mzid","")+"_decoys"
        tmp_msgf2pin_cmd+=" 2>&1"
        print "RUNNING COMMAND",tmp_msgf2pin_cmd 
        processes.append(subprocess.Popen(tmp_msgf2pin_cmd,shell=True))
        if len(processes) >= max_concurrency:
            for p in processes:
                print p.wait(),"finished..."
            processes=[]

    #print "MSGF2PIN is running...."
    for p in processes:#Make sure everything is done...
        print p.wait(),"finished..."


    #output_files[eachgroup]=eachgroup.replace(".mzid","")+"/"+eachgroup.replace(".mzid","")+".tsv"
    unfixed_pin_files=[]
    file_groups=[]
    file_group_to_ufpins={}
    file_group_to_dfufpins={}
    file_group_to_dfufstrings={}
    for each_file_input in output_files:
        print "STARTING WORKING ON ",each_file_input.replace(".mzid","")
        file_groups.append(each_file_input.replace(".mzid",""))
        unfixed_pin_files.append(output_files[each_file_input])
        if each_file_input.replace(".mzid","") not in file_group_to_ufpins:
            file_group_to_ufpins[each_file_input.replace(".mzid","")]=output_files[each_file_input]
        else:
            print "This group was already in the dictionary...",each_file_input.replace(".mzid",""),"SHOULD NEVER HAPPEN!<---------------"
            file_group_to_ufpins[each_file_input.replace(".mzid","")]=[file_group_to_ufpins[each_file_input.replace(".mzid","")]]
            file_group_to_ufpins[each_file_input.replace(".mzid","")].append(output_files[each_file_input])

        this_pin_dictlist=[]

        file_group_to_dfufstrings[each_file_input.replace(".mzid","")]=[]
        i=0
        with open(each_file_input.replace(".mzid","")+"/"+each_file_input.replace(".mzid","")+".tsv",'r') as ufpin_reader:
            for eachline in ufpin_reader:
                thisline=eachline.split('\t')
                if i == 0:
                    thisheader=thisline
                elif i == 1:
                    directions=thisline
                else:
                    this_pin_dictlist.append({'specid':thisline[0],'data':thisline,'scan':thisline[0].split("_",1)[1].rsplit("_",1)[0],'label':thisline[0].split("_")[0]})
                #file_group_to_dfufstrings[each_file_input.replace(".mzid","")].append(thisline)
                i+=1
        this_pin_df=pandas.DataFrame(this_pin_dictlist)
        grouped_df=this_pin_df.groupby('scan')
        #print grouped_df,"grouped!"
        newspecids=applyParallelQuarter(grouped_df,renumber_scans)
        file_group_to_dfufpins[each_file_input.replace(".mzid","")]=newspecids

        #print newfnames
        #print thisheader
        #print directions

        #file_group_to_dfufpins[each_file_input.replace(".mzid","")]=pandas.read_csv(each_file_input.replace(".mzid","")+"/"+each_file_input.replace(".mzid","")+".tsv",sep='\t',skiprows=[1],engine='python')
        #with open(each_file_input.replace(".mzid","")+"/"+each_file_input.replace(".mzid","")+".tsv",'r') as header_reader:
        #    for line in islice(header_reader,1,1):
        #        header=line

    #OKAY, NOW WE NEED TO TAKE ALL THE OUTPUTS WE HAVE IN OUR file_group_to_dfufpins DICTIONARY AND COMBINE THEM BY FRACTIONATION GROUP!
    #print output_files
    combined_pin_groups={}
    original_mzid_mapping={}
    for each_group in groups_to_mzml:
        this_grp=[]
        this_mapping={}
        #print groups_to_mzml[each_group]
        this_mzml_set=groups_to_mzml[each_group]
        for each_mzml in this_mzml_set:
            this_grp.append(file_group_to_dfufpins[each_mzml.replace(".mzML","")])
        combined_pin_groups[each_group]=pandas.concat(this_grp)
        this_df=combined_pin_groups[each_group][['specid','original_file']].copy(deep=True)
        this_df=this_df.set_index('specid')
        this_dict=this_df.to_dict()
        #print this_dict,"This is the dict <--------------------"
        for each_mzml in this_mzml_set:
            print "adding dict to...",each_mzml.replace(".mzML","")
            original_mzid_mapping[each_mzml.replace(".mzML","")]=this_dict
        #print combined_pin_groups[each_group]
        #for each_run in target_run_groups~~~~~~~~~~~~~~~~~~~~
    for each_group in combined_pin_groups:
        writeSILACpin(combined_pin_groups[each_group],"output/"+each_group+".tsv",thisheader,directions)
        output_list.append("output/"+each_group+".tsv")
        print "wrote file ","output/"+each_group+".tsv"

    #sys.exit(0)
    
    #print "THIS IS THE FINAL HEADER",header
    #print file_group_to_dfufpins
    #for each_fgroup_pin in 


    #sys.exit(0)
    #Alright, at this point I need to read in all the pin files (can be found in output_files dictionary)
    #Then the first column will be processed.  Spectral rank numbers will be changed so that they are non overlapping
    # .str.split and .str.contains are super fast.... yay pandas...

else: #silac isn't true....
    #target_run_groups={}
    #decoy_run_groups={}
    #for each_key in run_groups:
    #    target_run_groups[each_key]=[]
    #    decoy_run_groups[each_key]=[]


    #decoy_run_groups=copy.deepcopy(run_groups)

    #print decoy_run_groups, "THESE ARE RUN GROUPS<---------------"
    ##print "building commands..."
    ##processes = []
    ##msgf2pin_cmd="msgf2pin "+"-m "+nummatches+" --enzyme "+enzyme+" --pattern "+pattern

    ##if ptm:
    ##    msgf2pin_cmd+=" --PTM"
    ##if aafreq:
    ##    msgf2pin_cmd+=" --aa-freq"
    ##if pngasef:
    ##    msgf2pin_cmd+=" --PNGaseF"
    ##if isotope:
    ##    msgf2pin_cmd+=" --isotope"
    ##if databases:
    ##    msgf2pin_cmd+=" --databases "+databases

    #print "Here are out run groups",run_groups


    output_list=[]
    print "Taking care of a few file creation issues...."
    print "Here's our groups...",run_groups
    for each_group in run_groups:
        with open(each_group+".runs",'wb') as writer:
            print target_run_groups[each_group],"<------------------------------"
            for each_run in target_run_groups[each_group]:
                writer.write(each_run+"\n")
        with open(each_group+".decoys",'wb') as writer:
            for each_run in decoy_run_groups[each_group]:
                writer.write(each_run+"\n")
                #print each_run,"<-----------------------"


    #DEVNULL = open(os.devnull, 'wb')




        tmp_msgf2pin_cmd=msgf2pin_cmd
        tmp_msgf2pin_cmd+=" --outputTab "+outputfolder+each_group+".tsv"
        output_list.append(outputfolder+each_group+".tsv")
        tmp_msgf2pin_cmd+=" "+each_group+".runs"
        tmp_msgf2pin_cmd+=" "+each_group+".decoys"
        tmp_msgf2pin_cmd+=" 2>&1"
        #print "RUNNING COMMAND",tmp_msgf2pin_cmd
        #processes.append(subprocess.Popen(tmp_msgf2pin_cmd,shell=True))
        #processes.append(subprocess.Popen(tmp_msgf2pin_cmd,shell=True))
        print "RUNNING COMMAND",tmp_msgf2pin_cmd
        processes.append(subprocess.Popen(tmp_msgf2pin_cmd,shell=True))
        if len(processes) >= max_concurrency:
            for p in processes:
                print p.wait(),"finished..."
            processes=[]


    print "Starting MSGF2PIN...."
    for p in processes:
        print p.wait(),"finished..."


    pass




#for each_group in run_groups:
#    with open(run_groups[each_group]+".runs",'wb') as writer:
#        #print target_run_groups[each_group],"<------------------------------"
#        for each_run in target_run_groups[each_group]:
#            writer.write(each_run+"\n")
#    with open(run_groups[each_group]+".decoys",'wb') as writer:
#        for each_run in decoy_run_groups[each_group]:
#            writer.write(each_run+"\n")

#    tmp_msgf2pin_cmd=msgf2pin_cmd
#    tmp_msgf2pin_cmd+=" --outputTab "+outputfolder+run_groups[each_group]+".tsv"
#    output_list.append(outputfolder+run_groups[each_group]+".tsv")
#    tmp_msgf2pin_cmd+=" "+run_groups[each_group]+".runs"
#    tmp_msgf2pin_cmd+=" "+run_groups[each_group]+".decoys"
#    tmp_msgf2pin_cmd+=" 2>&1"
#    print "RUNNING COMMAND",tmp_msgf2pin_cmd
#    #processes.append(subprocess.Popen(tmp_msgf2pin_cmd,shell=True))
#    processes.append(subprocess.Popen(tmp_msgf2pin_cmd,shell=True))


if options.diaumpire:
    print "We're going to go ahead and alter the pin files scan numbers and SpecID's to mirror those used in the DIAUmpire converted mzML file!"
    for eachfile in output_list:
        with open(eachfile+"_dia",'wb') as filewriter:
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
                            
                        new_reconstructed_lines.append(header)

                    elif linectr ==1:
                        directions=eachline.split("\t")
                        #directions[len(directions)-1]+="\n"
                        new_reconstructed_lines.append(directions)
                    else:
                        #Okay, so this is where the magic is actually going to happen.
                        #1. for each line, we'll pull the SpecID and the Scan Nr
                        #2. We're going to decrement ScanNr so that it agrees with with DIAUmpire converted mzML files.
                        #3. finally, we'll add that new row to the new_reconstructed lines
                        split_line = eachline.split("\t")
                        specID=split_line[specID_index]
                        ScanNr = split_line[scan_index]

                        specID_corrected = specID.rsplit("_",5)
                        specID_corrected[3] = str( specID_corrected[1] )
                        ScanNr_corrected = str(specID_corrected[3])
                        specID_corrected= "_".join(specID_corrected)
                        split_line[specID_index]=specID_corrected
                        split_line[scan_index]=ScanNr_corrected
                        
                        new_reconstructed_lines.append(split_line)
                    linectr+=1
                print "Writing pin file with DIA fixed scans ",eachfile
                for each_newline in new_reconstructed_lines:
                    filewriter.write("\t".join(each_newline))
        os.rename(eachfile,eachfile+"_orig")
        shutil.copy(eachfile+"_dia",eachfile)


if options.ppm:
    if not options.fix_ppm:#options.fix_ppm is not None:
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
        print "Opening file ",eachfile,"..."
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
                                print label_index,"THIS IS LABEL INDEX"
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
                        #print specID,"specID"
                        #print specID.rsplit("_",7)[0],"key to find..."
                        #print original_mzid_mapping[specID.rsplit("_",7)[0]].keys()
                        if silac:
                            file_name_ppm=original_mzid_mapping[specID.rsplit("_",7)[0]]['original_file'][specID]
                        else:
                            if label == "1":
                                file_name_ppm="targets/"+specID.rsplit("_",7)[0]+"_target.mzid"
                            else:
                                file_name_ppm="decoys/"+specID.rsplit("_",7)[0]+"_decoy.mzid"

                        #print file_name_ppm,"lololol"#[specID]
                        scan_num=split_line[scan_index]

                        #print "Label was ",str(label),"so we will be opening the file",file_name_ppm,"from ",specID
                        if not file_name_ppm in mzid_dict:
                            reader=pyteomics.mzid.read(file_name_ppm)#,retrieve_refs=True)
                            mzid_dict[file_name_ppm]={}
                            #print "We'll have to iterate over the mass spec data to get high res masses..."
                            print "Reading in file",file_name_ppm,"..."
                            for each_scan in reader:
                               #mzid_dict[file_name_ppm][str(int(each_scan['scan number(s)']))]=each_scan                                      ##### OLD: Stored entire mzid scan information in dictionaries... filled memory...
                               if 'scan number(s)' in each_scan:
                                   mzid_dict[file_name_ppm][str(int(each_scan['scan number(s)']))]=each_scan['SpectrumIdentificationItem'][0]
                               elif options.diaumpire:
                                   #print "Didn't find traditional scan number headers.... Assuming DIA-Umpire converted files as inputs, using \"peak list scans\" tag instead."
                                   #mzid_dict[file_name_ppm][str(int(each_scan['spectrumID'].split("=")[1]))]=each_scan['SpectrumIdentificationItem'][0]
                                   #mzid_dict[file_name_ppm][str(int(each_scan['spectrumID'].split("=")[1])-1)]=each_scan['SpectrumIdentificationItem'][0]
                                   mzid_dict[file_name_ppm][str(int(each_scan['spectrumID'].split("=")[1])+1)]=each_scan['SpectrumIdentificationItem'][0]
                               else:
                                   print "WARNING: ASSUMING DIAUMPIRE INPUTS!"
                                   mzid_dict[file_name_ppm][str(int(each_scan['spectrumID'].split("=")[1])+1)]=each_scan['SpectrumIdentificationItem'][0] #We're assuming it's DIAUmpire...
                                   #print "I can't find scan numbers in this... BREAKING!"
                                   #sys.exit("failure")
                            del reader
                            gc.collect()
                        #print "-------------"
                        #print scan_num,type(scan_num)
                        #print file_name_ppm
                        #print str(each_scan)
                        #print each_scan['SpectrumIdentificationItem'][0]
                        #k=0
                        #for each_key in mzid_dict[file_name_ppm].keys():
                        #    if k<20:
                        #        print each_key
                        #        k+=1
                        #    else:
                        #        break
                        #if int(each_scan['spectrumID'].split("=")[1])==15:
                        #    print file_name_ppm,"HAD A 15 SCAN NUMBER...."
                        #    print each_scan
                        #print "----------------------------------"
                        #print split_line
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
                        #print split_line[peptide_index],"this is pep index"
                        #print split_line[peptide_index:],"this is after pep index"
                        #print "I would append",','.join(split_line[peptide_index:])
                        #sys.exit(2)
                        if options.fix_ppm:
                            new_line.append(split_line[peptide_index])
                            new_line.append(','.join(split_line[peptide_index+1:]))
                        else:
                            for each_item_line in split_line[peptide_index:]:
                                new_line.append(each_item_line)

                        #And we'll fix the ExpMass and CalcMass columns from msgf2pin.... These are the singly charged masses.
                        new_line[expmass_index]=str((exp_mz*mycharge)-((mycharge-1)*proton_mass))
                        new_line[mass_index]=str((theo_mz*mycharge)-((mycharge-1)*proton_mass))

                        new_reconstructed_lines.append(new_line)

                    linectr+=1


                print "Writing mass corrections to file and console:"
                for each_newline in new_reconstructed_lines:
                    #print each_newline
                    filewriter.write("\t".join(each_newline))

                if not options.fix_ppm:# is None:
                    pass

                else:#We're doing the PPM correction!

                    for each_newline in new_reconstructed_lines:
                        filewriter.write("\t".join(each_newline[:peptide_index+2])+"\n") #We're leaving out the proteins, etc...
                    #for_pandas=new_reconstructed_lines[:]
                    #directions=for_pandas[1]
                    #del for_pandas[1]
                    #for_pandas_fixed=[f.split("\t") for f in for_pandas]
                    #new_df = pandas.DataFrame(for_pandas)
                    print "Reading in ",eachfile+"_ppm"

                    new_df=pandas.read_csv(eachfile+"_ppm",sep="\t",skiprows=[1],low_memory=False)
                    new_df=new_df[new_df['SpecId'].str.contains("DefaultDirection")==False]
                    new_df=new_df[new_df['Label'].str.contains("Label")==False]
                    new_df['Label']=new_df['Label'].astype(int)
                    new_df['ppm']=new_df['ppm'].astype(float)
                    new_df=new_df[np.logical_and(new_df['lnEValue'] >= (-1.0*math.log(fix_ppm)),new_df['Label']==1)]
                    if len(new_df) == 0:
                        print "There are no confident IDs at the input ppm correction E-Value of "+str(fix_ppm)+" which is equivalent to "+str(-1.0*math.log(fix_ppm))+" as -log10(e)"
                        print "You must relax the threshold to find an E-value to filter on with a reasonable population size..."
                        sys.exit(1)
                    #print new_df['SpecId']
                    new_df['file']=new_df['SpecId'].str.rsplit("_",7).str.get(0)
                    #for index,eachrow in new_df.iterrows():
                    #    #print eachrow['SpecId']
                    #    new_df.loc[index,'file']=eachrow['SpecId'].rsplit("_",7)[0]
                    #print new_df
                    #sys.exit(1)
                    #print new_df.file.unique(),"file... BEFORE"
                    #print new_df
                    #new_df['file']=""
                    #print new_df.file.unique(),"file... AFTER"
                    ##grouped_df=new_df.groupby(np.arange(len(new_df))//cpu_count())
                    ##newfnames=applyParallel(grouped_df,fileNames)
                    #!!!!!!!!!!!!!!!!!!!~~

                    ##new_df['file']=newfilenames
                    #grouped_df=newfnames.groupby(newfnames.file)
                    grouped_df=new_df.groupby('file')
                    #print grouped_df,"this is the grouped...."
                    corrected_ppm=grouped_df['ppm'].median()
                    ppm_correction_list.append(corrected_ppm)
                    #print corrected_ppm,"aggggg"
                    #print corrected_ppm.index,"index!"
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
                                #print working_line[ppm_index],"ppm before..."
                                working_line[ppm_index]=str(float(working_line[ppm_index])-float( corrected_ppm.loc[working_line[specID_index].rsplit("_",7)[0]] ))
                                #print working_line[ppm_index],"ppm after..."
                                #print "after...",working_line
                                working_line[absppm_index]=str(abs(float(working_line[ppm_index])))
                                #print working_line,"before"
                                #sys.exit(1)
                                working_line[:peptide_index].append(working_line[peptide_index+1].replace(",","\t"))
                                #print working_line,"after"
                                #sys.exit(1)




                                #working_line[specID_index]=working_line[specID_index].rsplit("_",1)[0]  #### ADDED IN TO SEE IF PERCOLATOR WONT FUCK THIS UP
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

                            #print working_line
                            if working_line is not None:
                                ppmfixed_file.write("\t".join(working_line))
                            line_ctr+=1

                    #print aggregated
                    #print "I\'M HERE!!!"
                    #new_df.to_csv(eachfile+"_ppm")#CHANGE THIS TO WRITE THE FIXED FILE...


    for eachfile in output_list:
        os.remove(eachfile)
        if not options.fix_ppm: #options.fix_ppm is None:
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
