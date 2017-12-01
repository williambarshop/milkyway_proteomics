import os, sys, re
import optparse
import shutil
import pandas
import numpy
import subprocess
import fnmatch
import csv
import itertools
from joblib import Parallel, delayed
import multiprocessing
from pyteomics import mzid
from concurrent.futures import ThreadPoolExecutor # pip install futures
from subprocess import STDOUT, call


#####################################
#This is the script to add FIDO protein level q-values back into our percolator psms outputs...
#We'll add q-values into the .target.psms files, as well as use pyteomics to read some mzid files...
#
#VERSION 0.8.1
version="0.8.1"
#DATE: 11/21/2017
date="11/21/2017"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the FIDO q-value insertion tool for Galaxy, Wohlschlegel Lab UCLA"
print "Written by William Barshop"
print "Version: ",version
print "Date: ",date



def applyParallel(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)

def applyParallelQuarter(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count()/4)(delayed(func)(group) for name, group in dfGrouped)  # CHANGED TO 1 FOR NOW....
    #retLst = Parallel(n_jobs=1)(delayed(func)(group) for name, group in dfGrouped)  # CHANGED TO 1 FOR NOW....
    return pandas.concat(retLst)

def applyParallelHalf(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count()/2)(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)

def applyParallelOne(dfGrouped, func):
    retLst = Parallel(n_jobs=1)(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)

def makeUnmodSeq(group):
    for index,eachrow in group.iterrows():
        group.loc[index,'unmodified sequence']=re.sub("[\(\[].*?[\)\]]", "", eachrow['sequence'])
    return group

def add_Prot_Info(group):
    #protein_q_dict[protein]=qvalue
    #grouped_df <--- this is the DF grouped by "mzid_df.groupby('file_scan_id')" so we can .get_group(file_scan_id) to give us the rows back!
    for index,eachrow in group.iterrows():
        #if isinstance(eachrow['flanking aa'],int):
        #    print eachrow
        try:
            protein_list = eachrow['protein id'].split(",")
            aa_list = eachrow['flanking aa'].split(",")
        except:
            print eachrow['file'],eachrow['scan']
            print eachrow['flanking aa'],type(eachrow['flanking aa'])
            print eachrow
            print "____________________________________________________"
        #print "protein list is ",protein_list
        q_list=[]
        emp_q_list=[]
        start_stop_pos=[]
        fixed_aa_list=[]
        fill_with=None
        if len(aa_list)==1 and len(protein_list)>1:
            fill_with=aa_list[0]
        for each in itertools.izip_longest(protein_list,aa_list,fillvalue=fill_with):
            if each[0] is None:
                print "The protein list SHOULD NOT be smaller than the flanking amino acid list! BREKAING."
                sys.exit(2)
            prot=each[0].strip()
            aa=each[1]
            if aa is None:
                print "The flanking sequence aa should not be None... BREAKING!"
            #print "working on ",each
            try:
                q_list.append(str(protein_q_dict[prot]))
            except:
                #print "couldnt find ",prot
                q_list.append("1.0")
            try:
                emp_q_list.append(str(empirical_q_dict[prot]))
            except:
                #print "empirical couldnt find ",prot
                emp_q_list.append("1.0")
            #for each in grouped_df.groups:
            #    print each
            #if diaumpire:
            #    corrected_group=eachrow['file_scan_id'].split("_")
            #    corrected_group[1]=str(int(corrected_group[1])) #no correction, or -1? +1?
            #    corrected_group="_".join(corrected_group)
            #    thisgrp=grouped_df.get_group(corrected_group)
            #else:# eachrow['file_scan_id'] in grouped_df.groups:
            #    #for eachone in grouped_df.groups:
            #    #    print eachone
            #    #print "looking for",eachrow['file_scan_id']
            thisgrp=grouped_df.get_group(eachrow['file_scan_id'])
            #sys.exit("FAILURE!")
            #else:
            #    print eachrow['file_scan_id'],"FAILURE...."
            #    print "WARNING: THERE ARE DISCREPANCIES BETWEEN THE SEARCH MZID FILES AND THE PERCOLATOR OUTPUTS"
            #    #sys.exit(1)
            #at this point, we need to find the matching unmodified sequence from thisgrp
            #print "ooking for",eachrow['unmodified sequence']
            #print thisgrp,"my group..."
            #filtergrp=thisgrp[numpy.logical_and(thisgrp['Sequence'].str.contains(eachrow['unmodified sequence']),thisgrp['Sequence'].str.len==len(eachrow['unmodified sequence']))]
            filtergrp=thisgrp[thisgrp['unmodified sequence'] == eachrow['unmodified sequence']]
            #print filtergrp,"that was filtered group..."
            if len(filtergrp) <1:
                pass
                print "CRAP! no length..."
            else:
                for innerindex,innerrow in filtergrp.iterrows():
                    tmp_dict={}
                    #print innerrow,"this is the row"
                    prot_startstop_list=innerrow['proteinacc_start_stop_pre_post_;'].split(";")
                    for eachprot in prot_startstop_list:
                        this_split=eachprot.rsplit("_",4)
                        tmp_dict[this_split[0]]=str(this_split[1])+"_"+str(this_split[2])
                    break
            #print eachrow,"the row!"
            #print tmp_dict,"tmp dict..."
            try:
            #if each in tmp_dict:
                start_stop_pos.append(tmp_dict[prot])
                fixed_aa_list.append(aa)
            #else:
            except:
                if clean:
                    group.loc[index,'clean_up']=0
                    break  # WILL THIS BREAK MESS THINGS UP? <----------------------------------------
                    #print "WE'LL REMOVE",eachrow['file_scan_id'],"for not having",each
                else:
                    pass
                #group.loc[
            ######~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        statistical_q=",".join(q_list)
        empirical_q=",".join(emp_q_list)
        if useEmpirical:
            group.loc[index,'protein q-values']=statistical_q
        else:
            group.loc[index,'protein q-values']=empirical_q
        group.loc[index,'empirical protein q-values']=empirical_q
        group.loc[index,'statistical protein q-values']=statistical_q
        group.loc[index,'peptide prot-indicies']=",".join(start_stop_pos)
        try:
            group.loc[index,'protein inference group']=protein_group_dict[prot] #will be the same for any protein in the group.
            group.loc[index,'protein inference group ID']=protein_to_groupID_dict[prot]
        except:
            group.loc[index,'protein inference group']=eachrow['protein id']
            group.loc[index,'protein inference group ID']=0
        #print aa_list,"b4"
        group.loc[index,'flanking aa']=",".join(fixed_aa_list)
        #print fixed_aa_list,"after",q_list
        #print "----------------------------------------------"
    return group


####################################
#Argument parsing! So much fun!
#We'll use OptParse even though some
#people really rave about argparse...
#
#
# NB: With Optparse, if an option is
#     not specified, it will take a
#     value of None
####################################


parser = optparse.OptionParser()
parser.add_option("--fidoq",action="store",type="string",dest="fidoq")
parser.add_option("--index",action="store_true", dest="index")
parser.add_option("--useEmpirical",action="store_true", dest="useEmpirical")
parser.add_option("--clean",action="store_true", dest="clean") #This is going to be to clean up after crux percolator's problem with --top-match filters!  ONLY ENABLE IF TOP-MATCH IS SUPPOSED TO BE = 1!
parser.add_option("--diaumpire",action="store_true", dest="diaumpire")

(options,args) = parser.parse_args()
global clean
global useEmpirical
useEmpirical=False
clean=False
if options.clean:
    clean=True
if options.useEmpirical:
    useEmpirical=True
if options.index:
    add_index=True
global diaumpire
if options.diaumpire:
    diaumpire=True
else:
    diaumpire=False


#We'll start by building a list of mzid files...
mzid_matches = []
for root, dirnames, filenames in os.walk('targets'):
    for filename in fnmatch.filter(filenames, '*.mzid'):
        mzid_matches.append(os.path.join(root, filename))

basedir=os.getcwd()
#os.chdir('targets')

#processes=[]
#cmd_list=[]
#first=True
global pep_mappings
pep_mappings={}
for each_match in mzid_matches:
    print "reading in file ...",each_match
    reader=mzid.read(each_match,indexed_tags={b'peptideEvidence'},use_index=True,retrieve_refs=True,iterative=True)
    #print "Reading it again...."
    #by_id=mzid.read(each_match,indexed_tags={b'peptideEvidence',b'dBSequence'},use_index=True,retrieve_refs=True) ################## <-------- ARE THESE INDEXED TAGS WORKING?
    print "Beginning to iterate over the PSMs..."
    for each in reader:
        #pepSeq=each['SpectrumIdentificationItem'][0]['PeptideEvidenceRef'][0]['peptideEvidence_ref'].split("_")[2]
        pepSeq=each['SpectrumIdentificationItem'][0]['PeptideSequence']
        if not pepSeq in pep_mappings:
            map_list=[]
            for eachProt in each['SpectrumIdentificationItem'][0]['PeptideEvidenceRef']:
                #eviRef=by_id.get_by_id(each_pep_evi_ref['peptideEvidence_ref'])
                #print eviRef
                start=eachProt['start']
                end=eachProt['end']
                pre=eachProt['pre']
                post=eachProt['post']
                #print eviRef['peptide_ref']
                #Protein Accession...
                accession=eachProt['accession']
                #accession=by_id.get_by_id(by_id.get_by_id(each_pep_evi_ref['peptideEvidence_ref'])['dBSequence_ref'])['accession']
                map_list.append("_".join([accession,str(start),str(end),pre,post]))
            pep_mappings[pepSeq]=";".join(map_list)
    print "Done with file ",each_match
    reader.clear_id_cache()
    #by_id.clear_id_cache()
    reader.clear_tree()
    #by_id.clear_tree()
    del reader
    #del by_id

#print pep_mappings
#sys.exit(0)
#os.chdir(basedir)



#mzid_df = pandas.DataFrame.from_dict(pep_mappings,orient="index")
#mzid_df.reset_index(level=0, inplace=True)
#mzid_df.columns=["Sequence","proteinacc_start_stop_pre_post_;"]




#now we'll get rid of the old stuff...
#del new_df_list

#Now we're going to have to read in all the protein Q-Values into a global dictionary so the multithreading can deal w/ it...
global protein_q_dict #again, for joblib to copy into each instance...
global empirical_q_dict
global protein_group_dict
global protein_groupID_dict

protein_q_dict={}
empirical_q_dict={}
protein_group_dict={}
protein_groupID_dict={}
protein_to_groupID_dict={}

protein_q_df=pandas.read_csv(options.fidoq,sep='\t')
i=1

for index,row in protein_q_df.iterrows():
    protein_group_list=[x.strip() for x in row['protein group'].split("\t")]
    for each_protein in protein_group_list:
        protein_q_dict[each_protein]=row['q-value']
        empirical_q_dict[each_protein]=row['empirical q-value']
        protein_group_dict[each_protein]=";".join(row['protein group'])
        if str(row['protein group']) not in protein_groupID_dict:
            #group_list.append(str(each_group))
            protein_groupID_dict[str(row['protein group'])]=i
            protein_to_groupID_dict[each_protein]=i
            i+=1
        else:
            protein_to_groupID_dict[each_protein]=protein_groupID_dict[str(row['protein group'])]

        protein_group_dict[each_protein]=protein_groupID_dict[str(row['protein group'])]

#Okay, now that we have that, we'll need to read in the percolator psms txt files (we'll do this one at a time to make things easy for us.)

targetpsms_matches = []
for root, dirnames, filenames in os.walk('output'):
    for filename in fnmatch.filter(filenames, '*.target.psms.txt'):
        targetpsms_matches.append(os.path.join(root, filename))


for each_psms_file in targetpsms_matches:
    #okay, we'll first take this file and read it into a pd df...
    print each_psms_file,"loading..."
    this_df = pandas.read_csv(each_psms_file,sep='\t')
    this_df['scan']=this_df['scan'].astype(str)
    this_df['flanking aa'].fillna("NA",inplace=True)
    #this_df_tmp_group=this_df.groupby(numpy.arange(len(this_df))//multiprocessing.cpu_count())

    #re.sub("[\(\[].*?[\)\]]", "", x)

    #if 'unmodified sequence' not in this_df.columns:
    #    remove_unmod=True

    this_df_group=this_df.groupby(numpy.arange(len(this_df))//multiprocessing.cpu_count())
    this_df=applyParallelHalf(this_df_group,makeUnmodSeq)
    this_df['proteinacc_start_stop_pre_post_;']=this_df['unmodified sequence'].map(pep_mappings)

    this_df['file_scan_id']=this_df['file']+"_"+this_df['scan']+"_"+this_df['unmodified sequence']
    this_df['scan']=this_df['scan'].astype(int)
    #this_df['statistical protein q-values']=""
    #this_df['protein inference group']=""
    if options.clean:
        this_df['clean_up']=1
        
    #now we'll iterate over these and grab the groups by groupby
    global grouped_df
    grouped_df=this_df.groupby('file_scan_id')
    eachpsms_grouped_df=this_df.groupby('file_scan_id')
    #this_df.to_csv("TESTING.csv")
    #print this_df,"This is before runnning...."
    #break #~!~~~~~~~~~~~~
    addedProtInfo=applyParallelHalf(eachpsms_grouped_df,add_Prot_Info) # CHANGE BACK TO QUARTER...
    del grouped_df
    print "Done adding protein information...",addedProtInfo
    if options.clean:
        addedProtInfo=addedProtInfo[addedProtInfo['clean_up']==1]
        addedProtInfo.drop('clean_up',axis=1,inplace=True)
        print "After cleaning...",addedProtInfo

    #print addedProtInfo,"thats it for",each_psms_file
    addedProtInfo.to_csv(each_psms_file,sep='\t',index=False)

#global grouped_df #This will copy the DF to all joblib instances running in parallel!
#grouped_df = mzid_df.groupby('file_scan_id')


print "All done!"
