import os, sys, re
import shutil
import optparse
import shutil
import pandas
import numpy
import subprocess
import fnmatch
from joblib import Parallel, delayed
import multiprocessing
from Bio import SeqIO
import glob

#####################################
#This is the script to produce SSL file outputs for skyline spectral library construction
#
#Fraction parsing is taken from from after the final "-" in the file name.  For example, "2015-10-05-wb-HEK293-BioRep1-F1.mzML" 
#would belong to fraction "F1"
#
#VERSION 1.8.0
version="1.8.0"
#DATE: 9/18/2017
date="9/18/2017"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the SSL file converter for Galaxy, Wohlschlegel Lab UCLA"
print "Written by William Barshop"
print "Version: ",version
print "Date: ",date


def applyParallel(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)

def applyParallelQuarter(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count()/4)(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)

def applyParallelHalf(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count()/2)(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)


def fido_filter(unfiltered):#dataframe,float arguments...
    
    #filtered=[]
    #print unfiltered,"This is it GOING IN...."
    #unfiltered.apply(fido_filter_row,axis=1)
    for index,eachrow in unfiltered.iterrows():
        prot_q_vals = eachrow['protein q-values'].split(',')
        protein_ids = eachrow['protein id'].split(',')
        indicies = eachrow['peptide prot-indicies'].split(',')
        new_q_vals=[]
        new_prot_ids=[]
        new_indicies=[]

        for each_tuple in zip(protein_ids,prot_q_vals,indicies):
            if float(each_tuple[1])<=gfido_q_threshold: #KEEP IT
                new_q_vals.append(each_tuple[1])
                new_prot_ids.append(each_tuple[0])
                new_indicies.append(each_tuple[2])
            else:
                pass # THROW IT OUT

        if len(new_prot_ids) >= 1:
            unfiltered.loc[index,'protein q-values']=",".join(new_q_vals)
            unfiltered.loc[index,'protein id']=",".join(new_prot_ids)
            unfiltered.loc[index,'indicies']=",".join(new_indicies)
            unfiltered.loc[index,'fido_filter']=1 #This means we'll keep it when we do our df filter... We'll end up dropping this column later.
        else:
            unfiltered.loc[index,'fido_filter']=0

        #print changed,type(changed),"this is changed..."
    return unfiltered


def fido_filter_row(eachrow):
    #print eachrow,"this is my row!",type(eachrow)
    #print str(gfido_q_threshold),"this is the threshold..."
    prot_q_vals = eachrow['protein q-values'].split(',')
    protein_ids = eachrow['protein id'].split(',')
    indicies = eachrow['peptide prot-indicies'].split(',')
    new_q_vals=[]
    new_prot_ids=[]
    new_indicies=[]

    for each_tuple in zip(protein_ids,prot_q_vals,indicies):
        if float(each_tuple[1])<=gfido_q_threshold: #KEEP IT
            new_q_vals.append(each_tuple[1])
            new_prot_ids.append(each_tuple[0])
            new_indicies.append(each_tuple[2])
        else:
            pass # THROW IT OUT

    if len(new_prot_ids) >= 1:
        eachrow['protein q-values']=",".join(new_q_vals)
        eachrow['protein id']=",".join(new_prot_ids)
        eachrow['indicies']=",".join(new_indicies)
        eachrow['fido_filter']=1 #This means we'll keep it when we do our df filter... We'll end up dropping this column later.
    else:
        eachrow['fido_filter']=0


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

fractions=False

#print sys.argv,"THESE ARE THE ARGS"

parser = optparse.OptionParser()
parser.add_option("--pout",action="store",type="string",dest="operation_folder")
parser.add_option("--qthresh",action="store", type="float", dest="q_threshold")
parser.add_option("--fidoq",action="store", type="float", dest="fido_q_threshold")
parser.add_option("--fido_q_raw",action="store", type="string", dest="fido_q_raw")
parser.add_option("--pRS_prob_thresh",action="store", type="float", dest="prs_prob_threshold")
parser.add_option("--ptmRS_prob_thresh",action="store", type="float", dest="ptmrs_prob_threshold")
parser.add_option("--FLRthresh",action="store", type="float", dest="flr_threshold")
parser.add_option("--LFLR",action="store_true", dest="local_flr")
parser.add_option("--mzml",action="store",type="string",dest="mzml_files")
parser.add_option("--mass_corrected_mzml",action="store",type="string",dest="mc_mzml_files")
parser.add_option("--blib",action="store_true",dest="blib")
parser.add_option("--fasta",action="store",type="string",dest="fasta")
parser.add_option("--ffasta",action="store",type="string",dest="ffasta") # THIS IS THE FILTERED FIDO OUTPUT
parser.add_option("--expgroups",action="store",type="string",dest="exp_group_file")
parser.add_option("--ssl",action="store",type="string",dest="ssl_output_folder")
parser.add_option("--fractions",action="store_true",dest="fractions")
parser.add_option("--OnePPS",action="store_true",dest="one_pps")
parser.add_option("--only_mod",action="store_true",dest="only_mod")
parser.add_option("--highest_intensity",action="store_true",dest="highest_intensity")
parser.add_option("--diaumpire",action="store_true",dest="diaumpire")
parser.add_option("--no_mzml",action="store_true",dest="no_mzml")
#parser.add_option("--saint",action="store",type="string",dest="saint_outputs")


(options,args) = parser.parse_args()

if options.fractions is True:
    fractions=True
else:
    print "We'll give outputs by acquisition (runs), not by experimental set."





#saint=False
#if options.saint_outputs is not None:
#    saint=True
#    saint_interact,saint_bait,saint_prey=options.saint_outputs.split(",")
    


#### Check for FIDO q-filter
fido_q=False
if options.fido_q_threshold is not None:
    fido_q=True
    global gfido_q_threshold #For paralell access via JobLib/Multiprocessing...
    gfido_q_threshold = options.fido_q_threshold
    print "We're going to filter by Fido Q-Value of ",options.fido_q_threshold
        






#### Check for LuciPHOr
luciphor=False
if options.flr_threshold is not None:
    luciphor=True
    print "We will filter by LuciPHOr FLR of ",options.flr_threshold

#### Check for PhosphoRS
phosphoRS=False
if options.prs_prob_threshold is not None:
    phosphoRS=True
    print "We will filter by a PhosphoRS probability of ",options.prs_prob_threshold

#### Check for ptmRS
ptmRS=False
if options.ptmrs_prob_threshold is not None:
    ptmRS=True
    options.ptmrs_prob_threshold=options.ptmrs_prob_threshold*100.0
    print "We will filter by a ptmRS probability of ",options.ptmrs_prob_threshold


psms_files=[]
for root, subFolders, files in os.walk(options.operation_folder):
    for eachfile in files:
        if 'target.psms.txt' in eachfile and 'uncorrected' not in eachfile:
            psms_files.append(str(os.path.join(root,eachfile)))

dataframe_vector=[]
for eachfile in psms_files:
    newdf=pandas.read_csv(eachfile,sep='\t',index_col=False)
    dataframe_vector.append(newdf)
combined_results=pandas.concat(dataframe_vector)
#print combined_results.columns
#print combined_results['file']
#sys.exit(2)
del dataframe_vector

group_information = pandas.read_csv(options.exp_group_file,sep='\t',dtype={'Fractionation Group ID String': object,'Fractionation Group Name':object,'Biological Condition':object})

run_dict={} # Key is file_idx, value is file_name.mzML
rev_run_dict={} #Key is file_nzme.mzML, value is file_idx
group_to_run_dict={} # Key is group, value is [1, 2, 3, 4] list of file_idx belonging to runs in the group...
run_to_group_dict={} # Key is file_idx, value is group...
group_to_file_name={} # key is group, value is ["xxxx.mzML", "xxxx.mzML"]
if fractions:
    fractions_to_run_dict={}

for index,row in group_information.iterrows():
    print row
    run_dict[str(row['Crux File Integer'])]=row['Original File Name']+".mzML"
    rev_run_dict[row['Original File Name']+".mzML"]=str(row['Crux File Integer'])
    if row['Fractionation Group ID String'] in group_to_run_dict:
        group_to_run_dict[row['Fractionation Group ID String']].append(str(row['Crux File Integer']))
    else:
        group_to_run_dict[row['Fractionation Group ID String']] = [str(row['Crux File Integer'])]
    if str(row['Crux File Integer']) in run_to_group_dict:
        #run_to_group_dict[str(row['Crux File Integer'])].append(row['Fractionation Group ID String'])
        print "YOU HAVE MULTIPLE COPIES OF A RUN IN THE EXPERIMENTAL INFORMATION FILE... WARNING!"
    else:
        run_to_group_dict[str(row['Crux File Integer'])]=row['Fractionation Group ID String']

    if row['Fractionation Group ID String'] in group_to_file_name:
        group_to_file_name[row['Fractionation Group ID String']].append(str(row['Original File Name'])+".mzML")
    else:
        group_to_file_name[row['Fractionation Group ID String']] = [str(row['Original File Name'])+".mzML"]
    if fractions:
        fraction_tag=str(row['Original File Name'].rsplit("-",1)[1])
        if fraction_tag in fractions_to_run_dict:
            fractions_to_run_dict[fraction_tag].append(str(row['Crux File Integer']))#str(row['Original File Name'])+".mzML")
        else:
            fractions_to_run_dict[fraction_tag]=[str(row['Crux File Integer'])]


combined_results['file_idx']=combined_results['file_idx'].astype(str)
#print combined_results['file']
combined_results['file']=combined_results['file'].astype(str)
#sys.exit(2)
####################### We'll handle putting in the file names just to be sure this has been handled!
fix_combined=[]
#print run_dict,"This is run dict"
#print rev_run_dict,"this is rev dict"
for each_idx in run_dict:
    mask = combined_results[(combined_results.file.str.contains(run_dict[each_idx]))] # MODIFIED....
    mask['file']=run_dict[each_idx]
    fix_combined.append(mask)
#print fix_combined
combined_results=pandas.concat(fix_combined)
#sys.exit(2)
if options.diaumpire:
    print "DIAUmpire inputs: Decrementing scans to match mzML file index."
    combined_results['scan']=combined_results['scan']-1




if luciphor and 'luci_numPPS' not in combined_results.columns.values.tolist():
    print "LuciPHOr erroneously called for filtering... turning it off!"
    luciphor=False
if phosphoRS and 'pRS_peptideLocalizationProbability' not in combined_results.columns.values.tolist():
    print "phosphoRS erroneously called for filtering... turning it off!"
    phosphoRS=False
if ptmRS and 'ptmRS_peptideLocalizationProbability' not in combined_results.columns.values.tolist():
    print "ptmRS erroneously called for filtering... turning it off!"
    ptmRS=False




if luciphor:
    #print combined_results
    combined_results['luci_numPPS']=combined_results['luci_numPPS'].fillna(0)
    combined_results['luci_numPPS']=combined_results['luci_numPPS'].astype(int).fillna(0)
    combined_results['luci_globalFLR']=combined_results['luci_globalFLR'].astype(float).fillna(0.0)
    combined_results['luci_localFLR']=combined_results['luci_localFLR'].astype(float).fillna(0.0)
if phosphoRS:
    combined_results['pRS_peptideLocalizationProbability']=combined_results['pRS_peptideLocalizationProbability'].astype(float).fillna(1.00)
    combined_results['pRS_numPPS']=combined_results['pRS_numPPS'].astype(int).fillna(0)
if ptmRS:
    combined_results['ptmRS_peptideLocalizationProbability']=combined_results['ptmRS_peptideLocalizationProbability'].astype(float).fillna(100.00)
    combined_results['ptmRS_totalNumPPS']=combined_results['ptmRS_totalNumPPS'].fillna(0).astype(int)
#for i,each in combined_results.iterrows():



new_results={} # KEY IS GROUP NAME
results_per_run={} # KEY IS FILE_IDX
#print "tehse are items in group to run dict",group_to_run_dict






#for each_group in group_to_run_dict:
#    bool_filter=combined_results.copy(deep=True)
#    bool_filter['file_idx']=bool_filter['file_idx'].astype(str)
#    #mask8=combined_results[(bool_filter.file_idx == "8" )]
#    mask = combined_results[(combined_results.file_idx.str.contains('|'.join(group_to_run_dict[each_group])))]         #############Was this inappropriate set to file instead of file_idx? #bad idea to use "in"/contains here...
#    new_results[each_group]=mask # results by group
#    for each_file in set(mask['file_idx']):
#        each_file_mask=mask[(mask.file_idx == each_file)]
#        results_per_run[str(each_file)]=each_file_mask # results by run
#print group_to_file_name,"these are the groups of new interest"
for each_group in group_to_file_name:
    #print "-------------------------------------------------------------"
    #print "working on",each_group
    #bool_filter=combined_results.copy(deep=True)
    #bool_filter['file_idx']=bool_filter['file_idx'].astype(str)
    #mask8=combined_results[(bool_filter.file_idx == "8" )]
    #print "masking for",group_to_file_name[each_group]
    mask = combined_results[(combined_results.file.str.contains('|'.join(group_to_file_name[each_group])))]
    new_results[each_group]=mask # results by group
    #print mask
    print "CONSTRUCTING DICTIONARIES!"
    for each_file in set(mask['file']):
        each_file_mask=mask[(mask.file.str.contains(each_file))]
        results_per_run[str(rev_run_dict[each_file])]=each_file_mask # results by run
        #print each_file,str(rev_run_dict[each_file])
        #print len(each_file_mask)
        #print "---------------------"
#print combined_results
#sys.exit(2)
####################################################



if options.one_pps is True:
    print "Allowing Unambiguous Localization Peptides Through..."

basedir=os.getcwd()
if not fractions: #This is for all times when we have 1-D runs to compare.
    for eachgroup in set(group_information['Fractionation Group ID String']):
        os.chdir(basedir)
        os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
        #os.mkdir(basedir+"/"+options.ssl_output_folder+"/"+eachgroup)

        #if not options.no_mzml:
        for eachfile in group_to_file_name[eachgroup]:
            #shutil.copy(os.path.join(basedir,eachfile),basedir+"/"+options.ssl_output_folder+"/"+eachgroup+"/"+eachfile)
            shutil.copy(os.path.join(basedir,eachfile),basedir+"/"+options.ssl_output_folder+"/"+eachfile)
        if not fractions:
            for eachrun in group_to_run_dict[eachgroup]:
                #print "-----------------------------------------------"
                #print results_per_run.keys()
                this_run_results=results_per_run[eachrun]
                this_run_results['protein q-values']=this_run_results['protein q-values'].astype(str)

                if fido_q:
                    this_run_results['fido_filter']=1
                    this_run_grouped=this_run_results.groupby(numpy.arange(len(this_run_results))//multiprocessing.cpu_count())
                    this_run_results=applyParallelHalf(this_run_grouped,fido_filter)
                    #print type(this_run_results),"and type"
                    this_run_results=this_run_results[this_run_results['fido_filter']==1]
                    this_run_results=this_run_results.drop('fido_filter',axis=1)
                    del this_run_grouped

                if luciphor:
                    if options.local_flr is True:
                        if options.one_pps is True:
                            mask=this_run_results[numpy.logical_or(( this_run_results.luci_localFLR <= options.flr_threshold ) , (this_run_results.luci_numPPS - this_run_results.luci_numRPS == 0))]
                            type2_mask=this_run_results[numpy.logical_and(( this_run_results.luci_localFLR > options.flr_threshold ), (this_run_results.luci_numPPS - this_run_results.luci_numRPS > 0))] #Unambiguous and poorly localized.
                    
                        else:
                            mask=this_run_results[( this_run_results.luci_localFLR <= options.flr_threshold)]
                            type2_mask=this_run_results[( this_run_results.luci_localFLR > options.flr_threshold )] #Unambiguous and poorly localized.

                    else:
                        if options.one_pps is True:
                            mask=this_run_results[numpy.logical_or(( this_run_results.luci_globalFLR <= options.flr_threshold ) , (this_run_results.luci_numPPS - this_run_results.luci_numRPS == 0))]
                            type2_mask=this_run_results[numpy.logical_and(( this_run_results.luci_globalFLR > options.flr_threshold ), (this_run_results.luci_numPPS - this_run_results.luci_numRPS > 0))] #Unambiguous and poorly localized.
                    
                        else:
                            mask=this_run_results[( this_run_results.luci_globalFLR <= options.flr_threshold)]
                            type2_mask=this_run_results[( this_run_results.luci_globalFLR > options.flr_threshold )] #Unambiguous and poorly localized.
                elif phosphoRS:
                    if options.one_pps is True:
                        mask=this_run_results[numpy.logical_or(( this_run_results.pRS_peptideLocalizationProbability >= options.prs_prob_threshold ) , (this_run_results.pRS_numPPS - this_run_results.numModSites == 0))]
                        type2_mask=this_run_results[numpy.logical_and(( this_run_results.pRS_peptideLocalizationProbability >= options.prs_prob_threshold ), (this_run_results.pRS_numPPS - this_run_results.numModSites > 0))] #Unambiguous and poorly localized.
                    else:
                        mask=this_run_results[( this_run_results.pRS_peptideLocalizationProbability >= options.prs_prob_threshold)]
                        type2_mask=this_run_results[( this_run_results.pRS_peptideLocalizationProbability < options.prs_prob_threshold )] #Unambiguous and poorly localized.
                elif ptmRS:
                    #if options.one_pps is True:
                    #    mask=this_run_results[numpy.logical_or(( this_run_results.ptmRS_peptideLocalizationProbability >= options.ptmrs_prob_threshold ) , (this_run_results.ptmRS_numPPS - this_run_results.numModSites == 0))]
                    #    type2_mask=this_run_results[numpy.logical_and(( this_run_results.ptmRS_peptideLocalizationProbability >= options.ptmrs_prob_threshold ), (this_run_results.ptmRS_numPPS - this_run_results.numModSites > 0))] #Unambiguous and poorly localized.
                    #else:
                    #print len(this_run_results),"before filtering..."
                    #print "filtering at ",options.ptmrs_prob_threshold
                    mask=this_run_results[( this_run_results.ptmRS_peptideLocalizationProbability >= options.ptmrs_prob_threshold)]
                    #print mask['ptmRS_peptideLocalizationProbability'].min(),"minimum value..."
                    #for key,eachrow in mask.iterrows():
                    #    if eachrow['ptmRS_peptideLocalizationProbability']<=options.ptmrs_prob_threshold:
                    #        print eachrow
                    #print len(mask),"after filtering..."
                    type2_mask=this_run_results[( this_run_results.ptmRS_peptideLocalizationProbability < options.ptmrs_prob_threshold )] #Unambiguous and poorly localized.
                    #print type2_mask
                else:
                    mask=this_run_results
                if options.only_mod is True:
                    if luciphor:
                        mask=mask[mask['luci_numRPS'] >= 1].copy(deep=True)
                    elif phosphoRS:
                        mask=mask[mask['numModSites'] >= 1].copy(deep=True)
                    elif ptmRS:
                        mask=mask[mask['ptmRS_numMods'] >= 1].copy(deep=True)
                    
                if luciphor:
                    ssl_df=mask[['file','scan','charge','luci_sequence','percolator q-value']]
                    ssl_df.rename(columns={'luci_sequence':'sequence'}, inplace=True)
                elif phosphoRS:
                    ssl_df=mask[['file','scan','charge','sequence','percolator q-value']]
                    #ssl_df.rename(columns={'pRS_sequence':'sequence'}, inplace=True)                    
                elif ptmRS:
                    ssl_df=mask[['file','scan','charge','sequence','percolator q-value']]
                    #ssl_df.rename(columns={'ptmRS_sequence':'sequence'}, inplace=True)                    
                else:
                    ssl_df=mask[['file','scan','charge','sequence','percolator q-value']]
                ssl_df.rename(columns={'percolator q-value':'percolator qvalue'}, inplace=True)
                ssl_df_filtered=ssl_df[(ssl_df['percolator qvalue']<=options.q_threshold)]
                #with open(basedir+"/"+options.ssl_output_folder+"/"+eachgroup+"/"+run_dict[eachrun].replace(".mzML","")+".ssl",'wb') as ssl_writer:
                with open(basedir+"/"+options.ssl_output_folder+"/"+run_dict[eachrun].replace(".mzML","")+".ssl",'wb') as ssl_writer:
                    ssl_df_filtered.rename(columns={'percolator qvalue':'score'}, inplace=True)
                    ssl_df_filtered['score-type']="PERCOLATOR QVALUE"
                    ssl_df_filtered.to_csv(path_or_buf=ssl_writer,sep="\t",index=False,header=True)
                if luciphor:
                    with open(basedir+"/"+options.ssl_output_folder+"/"+run_dict[eachrun].replace(".mzML","")+"_type2.ssl",'wb') as type2_ssl_writer:
                        ssl_df_type2=type2_mask[['file','scan','charge','luci_sequence','percolator q-value']]
                        ssl_df_type2.rename(columns={'percolator q-value':'percolator qvalue'}, inplace=True)
                        ssl_df_type2.rename(columns={'luci_sequence':'sequence'}, inplace=True)
                        ssl_df_type2_filtered=ssl_df_type2[(ssl_df_type2['percolator qvalue']<=options.q_threshold)]
                        ssl_df_type2_filtered.rename(columns={'percolator qvalue':'score'}, inplace=True)
                        ssl_df_type2_filtered['score-type']="PERCOLATOR QVALUE"
                        ssl_df_type2_filtered.to_csv(path_or_buf=type2_ssl_writer,sep="\t",index=False,header=True)
                elif phosphoRS:
                    with open(basedir+"/"+options.ssl_output_folder+"/"+run_dict[eachrun].replace(".mzML","")+"_type2.ssl",'wb') as type2_ssl_writer:
                        #ssl_df_type2=type2_mask[['file','scan','charge','pRS_sequence','percolator q-value']]
                        ssl_df_type2=type2_mask[['file','scan','charge','sequence','percolator q-value']]
                        ssl_df_type2.rename(columns={'percolator q-value':'percolator qvalue'}, inplace=True)
                        #ssl_df_type2.rename(columns={'pRS_sequence':'sequence'}, inplace=True)
                        ssl_df_type2_filtered=ssl_df_type2[(ssl_df_type2['percolator qvalue']<=options.q_threshold)]
                        ssl_df_type2_filtered.rename(columns={'percolator qvalue':'score'}, inplace=True)
                        ssl_df_type2_filtered['score-type']="PERCOLATOR QVALUE"
                        ssl_df_type2_filtered.to_csv(path_or_buf=type2_ssl_writer,sep="\t",index=False,header=True)
                elif ptmRS:
                    with open(basedir+"/"+options.ssl_output_folder+"/"+run_dict[eachrun].replace(".mzML","")+"_type2.ssl",'wb') as type2_ssl_writer:
                        #ssl_df_type2=type2_mask[['file','scan','charge','ptmRS_sequence','percolator q-value']]
                        ssl_df_type2=type2_mask[['file','scan','charge','sequence','percolator q-value']]
                        ssl_df_type2.rename(columns={'percolator q-value':'percolator qvalue'}, inplace=True)
                        #ssl_df_type2.rename(columns={'ptmRS_sequence':'sequence'}, inplace=True)
                        ssl_df_type2_filtered=ssl_df_type2[(ssl_df_type2['percolator qvalue']<=options.q_threshold)]
                        ssl_df_type2_filtered.rename(columns={'percolator qvalue':'score'}, inplace=True)
                        ssl_df_type2_filtered['score-type']="PERCOLATOR QVALUE"
                        ssl_df_type2_filtered.to_csv(path_or_buf=type2_ssl_writer,sep="\t",index=False,header=True)

    if options.highest_intensity:
        #We'll need the pin files for this...
        os.chdir(basedir)
        pin_files=[]
        for root,dirs,files in os.walk(basedir):
            for file in files:
                if file.endswith(".pin"):
                    pin_files.append(os.path.join(root,file))
        pins_df=[]
        for each_pin in pin_files:
            cols = pandas.read_csv(each_pin, nrows=1,sep="\t",index_col=False).columns.tolist()[:-1]
            pins_df.append(pandas.read_csv(each_pin,usecols=cols,index_col=False,sep="\t",skiprows=[1],engine='python'))
        merged_pins=pandas.concat(pins_df)
        merged_pins['file_idx']=merged_pins['SpecId'].str.split("_").str.get(1)
        merged_pins['file']=merged_pins['file_idx'].map(run_dict)
        os.chdir(basedir+"/"+options.ssl_output_folder)
        for file in glob.glob("*.ssl"):
            this_ssl=pandas.read_csv(file,index_col=False,sep='\t')
            this_ssl_with_MS2=pandas.merge(this_ssl,merged_pins[['file','ScanNr','lnMS2IonCurrent']],left_on=['file','scan'],right_on=['file','ScanNr'])
            this_ssl_with_MS2=this_ssl_with_MS2.drop_duplicates()
            this_ssl_with_MS2_filtered=this_ssl_with_MS2.loc[this_ssl_with_MS2.groupby(['file','sequence','charge'])["lnMS2IonCurrent"].idxmax()]
            this_ssl_with_MS2_filtered.drop(labels=['lnMS2IonCurrent','ScanNr'],axis=1,inplace=True)
            os.rename(file,file.rsplit(".",1)[0]+".redundant_library")
            this_ssl_with_MS2_filtered.to_csv(path_or_buf=file,sep="\t",index=False,header=True)
        print "We're stripped down the ssl files to only the highest explained MS2 Ion Current per peptide ID per run."


    if options.blib:
        print "Let's build some blib files..."
        os.chdir(basedir)
        os.chdir(basedir+"/"+options.ssl_output_folder)
        #cmd = '/galaxy-central/tools/wohl-proteomics/ssl_converter/blibbuild/bin/BlibBuild -a milkyway-galaxy *.ssl combined_spectral_lib.blib'
        #cmd = 'ls *.ssl | /galaxy-central/tools/wohl-proteomics/ssl_converter/blibbuild/bin/BlibBuild -a milkyway-galaxy -m 1000M combined_spectral_lib.blib'
        #cmd = 'wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibBuild.exe -a milkyway-galaxy -m 1000M *.ssl redundant_spectral_lib.blib'
        #print "running command ",cmd
        #subprocess.call(cmd,shell=True)
        #filtercmd='wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibFilter.exe redundant_spectral_lib.blib filtered_spectral_lib.blib'
        #subprocess.call(filtercmd,shell=True)
        for file in glob.glob("*type2.ssl"):
            os.rename(file,file.split(".")[0]+".not_ssl")
            print "protected",file,"from inclusion in spectral lib..."
        cmd = 'wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibBuild.exe -a milkyway-galaxy -m 1000M *.ssl combined_spectral_lib.blib'
        print "running command ",cmd
        subprocess.call(cmd,shell=True)
        for file in glob.glob("*.not_ssl"):
            os.rename(file,file.split(".")[0]+".ssl")
        filtercmd='wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibFilter.exe combined_spectral_lib.blib filtered.blib'
        subprocess.call(filtercmd,shell=True)

        if options.mc_mzml_files is not None: #and not options.no_mzml:
            for file in glob.glob("*.mzML"):
                os.remove(file)
                print "removing ",file
            os.chdir(basedir)
            for file in glob.glob("mc_*.mzML"):
                print "replacing ",file
                shutil.copy(file,basedir+"/"+options.ssl_output_folder+"/"+file.split("_",1)[1])

else: #This is for when fractions is true... so we'll organize the output into fractions
      #fractions_to_run_dict[fraction_tag].append(str(row['Original File Name'])+".mzML")
    blib_cmds=[]
    filter_cmds=[]
    
    for eachfraction in fractions_to_run_dict:
        fraction_set=fractions_to_run_dict[eachfraction]
        for each_chrom_fraction in fraction_set:
            if not os.path.isdir(basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/"):
                os.mkdir(basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/")
            #if not options.no_mzml:
            shutil.copy(os.path.join(basedir,run_dict[each_chrom_fraction]),basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/"+run_dict[each_chrom_fraction])

    for eachfraction in fractions_to_run_dict:
        os.chdir(basedir)
        fraction_set=fractions_to_run_dict[eachfraction]
        for eachrun in fraction_set:
            this_run_results=results_per_run[eachrun]
            if fido_q:
                this_run_results['fido_filter']=0
                this_run_grouped=this_run_results.groupby(numpy.arange(len(this_run_results))//multiprocessing.cpu_count())
                this_run_results=applyParallelHalf(this_run_grouped,fido_filter)
                print this_run_results.columns,"this is columns..."
                this_run_results=this_run_results[this_run_results['fido_filter']==1]
                this_run_results=this_run_results.drop('fido_filter',axis=1)
                del this_run_grouped

            #print "------------------------------"
            #print "each run is ",eachrun
            #print set(this_run_results['file_idx']),"idx set..."
            #print set(this_run_results['file']),"file set..."
            #print "------------------------------"

            os.chdir(basedir)
            os.chdir(options.operation_folder+str(run_to_group_dict[eachrun])+".pin_out/crux-output/")

            if luciphor:
                if options.local_flr is True:
                    if options.one_pps is True:
                        mask=this_run_results[numpy.logical_or(( this_run_results.luci_localFLR <= options.flr_threshold ) , (this_run_results.luci_numPPS - this_run_results.luci_numRPS == 0))]
                        type2_mask=this_run_results[numpy.logical_and(( this_run_results.luci_localFLR >= options.flr_threshold ), (this_run_results.luci_numPPS - this_run_results.luci_numRPS > 0))] #Unambiguous and poorly localized.
                    
                    else:
                        mask=this_run_results[( this_run_results.luci_localFLR <= options.flr_threshold)]
                        type2_mask=this_run_results[( this_run_results.luci_localFLR >= options.flr_threshold )] #Unambiguous and poorly localized.
                else:
                    if options.one_pps is True:
                        mask=this_run_results[numpy.logical_or(( this_run_results.luci_globalFLR <= options.flr_threshold ) , (this_run_results.luci_numPPS - this_run_results.luci_numRPS == 0))]
                        type2_mask=this_run_results[numpy.logical_and(( this_run_results.luci_globalFLR >= options.flr_threshold ), (this_run_results.luci_numPPS - this_run_results.luci_numRPS > 0))] #Unambiguous and poorly localized.
                    
                    else:
                        mask=this_run_results[( this_run_results.luci_globalFLR <= options.flr_threshold)]
                        type2_mask=this_run_results[( this_run_results.luci_globalFLR >= options.flr_threshold )] #Unambiguous and poorly localized.
            elif phosphoRS:
                if options.one_pps is True:
                    mask=this_run_results[numpy.logical_or(( this_run_results.pRS_peptideLocalizationProbability >= options.prs_prob_threshold ) , (this_run_results.pRS_numPPS - this_run_results.numModSites == 0))]
                    type2_mask=this_run_results[numpy.logical_and(( this_run_results.pRS_peptideLocalizationProbability >= options.prs_prob_threshold ), (this_run_results.pRS_numPPS - this_run_results.numModSites > 0))] #Unambiguous and poorly localized.
                else:
                    mask=this_run_results[( this_run_results.pRS_peptideLocalizationProbability >= options.prs_prob_threshold)]
                    type2_mask=this_run_results[( this_run_results.pRS_peptideLocalizationProbability < options.prs_prob_threshold )] #Unambiguous and poorly localized.
            elif ptmRS:
                #if options.one_pps is True:
                #    mask=this_run_results[numpy.logical_or(( this_run_results.ptmRS_peptideLocalizationProbability >= options.ptmrs_prob_threshold ) , (this_run_results.ptmRS_numPPS - this_run_results.ptmRS_numMods == 0))]
                #    type2_mask=this_run_results[numpy.logical_and(( this_run_results.ptmRS_peptideLocalizationProbability >= options.ptmrs_prob_threshold ), (this_run_results.ptmRS_numPPS - this_run_results.ptmRS_numMods > 0))] #Unambiguous and poorly localized.
                
                mask=this_run_results[( this_run_results.ptmRS_peptideLocalizationProbability >= options.ptmrs_prob_threshold)]
                type2_mask=this_run_results[( this_run_results.ptmRS_peptideLocalizationProbability < options.ptmrs_prob_threshold )] #Unambiguous and poorly localized.

            else:
                mask=this_run_results
            if options.only_mod is True:
                if luciphor:
                    mask=mask[mask['luci_numRPS'] >= 1].copy(deep=True)
                elif phosphoRS:
                    mask=mask[mask['numModSites'] >= 1].copy(deep=True)
                elif ptmRS:
                    mask=mask[mask['ptmRS_numMods'] >= 1].copy(deep=True)
            if luciphor:
                ssl_df=mask[['file','scan','charge','luci_sequence','percolator q-value']]
                ssl_df.rename(columns={'luci_sequence':'sequence'}, inplace=True)
            elif phosphoRS:
                ssl_df=mask[['file','scan','charge','sequence','percolator q-value']]
                #ssl_df.rename(columns={'pRS_sequence':'sequence'}, inplace=True)
            elif ptmRS:
                ssl_df=mask[['file','scan','charge','sequence','percolator q-value']]
                #ssl_df.rename(columns={'ptmRS_sequence':'sequence'}, inplace=True)
            else:
                ssl_df=mask[['file','scan','charge','sequence','percolator q-value']]
            ssl_df.rename(columns={'percolator q-value':'percolator qvalue'}, inplace=True)
            ssl_df_filtered=ssl_df[(ssl_df['percolator qvalue']<=options.q_threshold)]
            with open(basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/"+run_dict[eachrun].replace(".mzML","")+".ssl",'wb') as ssl_writer:
                ssl_df_filtered.rename(columns={'percolator qvalue':'score'}, inplace=True)
                ssl_df_filtered['score-type']="PERCOLATOR QVALUE"
                ssl_df_filtered.to_csv(path_or_buf=ssl_writer,sep="\t",index=False,header=True)
                if luciphor:
                    with open(basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/"+run_dict[eachrun].replace(".mzML","")+"_type2.ssl",'wb') as type2_ssl_writer:
                        ssl_df_type2=type2_mask[['file','scan','charge','luci_sequence','percolator q-value']]
                        ssl_df_type2.rename(columns={'percolator q-value':'percolator qvalue'}, inplace=True)
                        ssl_df_type2.rename(columns={'luci_sequence':'sequence'}, inplace=True)
                        ssl_df_type2_filtered=ssl_df_type2[(ssl_df_type2['percolator qvalue']<=options.q_threshold)]
                        ssl_df_type2_filtered.rename(columns={'percolator qvalue':'score'}, inplace=True)
                        ssl_df_type2_filtered['score-type']="PERCOLATOR QVALUE"
                        ssl_df_type2_filtered.to_csv(path_or_buf=type2_ssl_writer,sep="\t",index=False,header=True)
                elif phosphoRS:
                    with open(basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/"+run_dict[eachrun].replace(".mzML","")+"_type2.ssl",'wb') as type2_ssl_writer:
                        #ssl_df_type2=type2_mask[['file','scan','charge','pRS_sequence','percolator q-value']]
                        ssl_df_type2=type2_mask[['file','scan','charge','sequence','percolator q-value']]
                        ssl_df_type2.rename(columns={'percolator q-value':'percolator qvalue'}, inplace=True)
                        #ssl_df_type2.rename(columns={'pRS_sequence':'sequence'}, inplace=True)
                        ssl_df_type2_filtered=ssl_df_type2[(ssl_df_type2['percolator qvalue']<=options.q_threshold)]
                        ssl_df_type2_filtered.rename(columns={'percolator qvalue':'score'}, inplace=True)
                        ssl_df_type2_filtered['score-type']="PERCOLATOR QVALUE"
                        ssl_df_type2_filtered.to_csv(path_or_buf=type2_ssl_writer,sep="\t",index=False,header=True)
                elif ptmRS:
                    with open(basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/"+run_dict[eachrun].replace(".mzML","")+"_type2.ssl",'wb') as type2_ssl_writer:
                        #ssl_df_type2=type2_mask[['file','scan','charge','ptmRS_sequence','percolator q-value']]
                        ssl_df_type2=type2_mask[['file','scan','charge','sequence','percolator q-value']]
                        ssl_df_type2.rename(columns={'percolator q-value':'percolator qvalue'}, inplace=True)
                        #ssl_df_type2.rename(columns={'ptmRS_sequence':'sequence'}, inplace=True)
                        ssl_df_type2_filtered=ssl_df_type2[(ssl_df_type2['percolator qvalue']<=options.q_threshold)]
                        ssl_df_type2_filtered.rename(columns={'percolator qvalue':'score'}, inplace=True)
                        ssl_df_type2_filtered['score-type']="PERCOLATOR QVALUE"
                        ssl_df_type2_filtered.to_csv(path_or_buf=type2_ssl_writer,sep="\t",index=False,header=True)
        if options.blib:
            print "We're going to build blib files!"
            os.chdir(basedir)
            #os.chdir(basedir+"/"+options.ssl_output_folder)
            os.chdir(basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/")
            #cmd = '/galaxy-central/tools/wohl-proteomics/ssl_converter/blibbuild/bin/BlibBuild -a milkyway-galaxy *.ssl '+eachfraction.replace("-","")+".ssl" #combined_spectral_lib.blib
            for file in glob.glob("*type2.ssl"):
                os.rename(file,file.split(".")[0]+".not_ssl")
                print "protected",file,"from inclusion in spectral lib..."

            command_folder=basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/"
            cmd = 'wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibBuild.exe -a milkyway-galaxy -m 1000M {0}*.ssl {0}combined_spectral_lib.blib'.format(command_folder)
            blib_cmds.append(cmd)
            print "storing command to run later",cmd
            #subprocess.call(cmd,shell=True)
            filtercmd='wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibFilter.exe {0} {1}'.format(command_folder+"combined_spectral_lib.blib",command_folder+"filtered_spectral_lib.blib")
            filter_cmds.append(filtercmd)
            print "storing command for filter later",filtercmd
            #subprocess.call(filtercmd,shell=True)

            #cmd = 'ls *.ssl | /galaxy-central/tools/wohl-proteomics/ssl_converter/blibbuild/bin/BlibBuild -a milkyway-galaxy -m 1000M combined_spectral_lib.blib'
            #print "Running command ...",cmd
            #subprocess.call(cmd,shell=True)

if fractions and options.blib:
    chunk_size=8 # Max concurrency...
    job_list=[blib_cmds[i:i + chunk_size] for i in range(0, len(blib_cmds), chunk_size)]
    for each_jobset in job_list:
        processes=[]
        for each_job in each_jobset:
            print "Running ...",each_job
            processes.append(subprocess.Popen(each_job,shell=True))
        for each_proc in processes:
            each_proc.wait()

    job_list=[filter_cmds[i:i + chunk_size] for i in range(0, len(filter_cmds), chunk_size)]
    for each_jobset in job_list:
        processes=[]
        for each_job in each_jobset:
            print "Running Filter...",each_job
            processes.append(subprocess.Popen(each_job,shell=True))
        for each_proc in processes:
            each_proc.wait()

    for eachfraction in fractions_to_run_dict:
        os.chdir(basedir)
        fraction_set=fractions_to_run_dict[eachfraction]
        for eachrun in fraction_set:
            os.chdir(basedir)
            os.chdir(basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/")

            if options.mc_mzml_files is not None:
                file_list=[file for file in glob.glob("*.mzML")]
                for file in file_list:
                    os.remove(file)
                os.chdir(basedir)
                for file in glob.glob("mc_*.mzML"):
                    if file.split("_",1)[1] in file_list:
                        shutil.copy(file,basedir+"/"+options.ssl_output_folder+"/"+eachfraction.replace("-","")+"/"+file.split("_",1)[1])
            for file in glob.glob("*.not_ssl"):
                os.rename(file,file.split(".")[0]+".ssl")

            os.chdir(basedir)


if options.fido_q_threshold is not None:
    os.chdir(basedir)
    #While we're here, let's go ahead and handle database filtering!
    #fido_q_df = pandas.read_csv(options.fido_q_raw,sep=",")
    fido_q_df = pandas.read_csv(options.fido_q_raw,sep="\t")
    fido_q_df=fido_q_df[fido_q_df['q-value']<=options.fido_q_threshold] #filter down...
    proteins_to_keep=fido_q_df['protein group'].unique().tolist()#this is the list of proteins we want to keep.   

    fasta_ids = list()
    #ensure no duplicates...
    with open("temp_unique.fasta",'a') as out_fasta:
        with open(options.fasta,'rb') as fasta_file:
            for each_seq in SeqIO.parse(fasta_file,"fasta"):
                if each_seq.id not in fasta_ids:
                    fasta_ids.append(each_seq.id)
                    SeqIO.write(each_seq,out_fasta,'fasta')
    #read in unique fasta as dict...
    with open("temp_unique.fasta",'rb') as fasta_file:
        fasta_dict=SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta"))

    new_fasta=[]
    for eachprotein in proteins_to_keep:
        if eachprotein in fasta_dict:
            new_fasta.append(fasta_dict[eachprotein])
    
        #print new_fasta,len(new_fasta),len(proteins_to_keep)
    shutil.copy(options.fasta,basedir+"/"+options.ssl_output_folder+"/"+"UNFILTERED_"+options.fasta)
    with open(basedir+"/"+options.ssl_output_folder+"/"+"FIDO_FILTERED_"+options.fasta,'wb') as fasta_writer:
        SeqIO.write(new_fasta,fasta_writer,"fasta")
    print "FIDO filtered FASTA written!"
    if options.ffasta is not None:
        shutil.copy(basedir+"/"+options.ssl_output_folder+"/"+"FIDO_FILTERED_"+options.fasta,options.ffasta)
        print "Copied the filtered fasta to the output location specified..."
else:
    os.chdir(basedir)
    shutil.copy(options.fasta,basedir+"/"+options.ssl_output_folder+"/"+"UNFILTERED_"+options.fasta)

os.chdir(basedir)
if options.no_mzml:
    mzml_files=[]
    for root, dirs, files in os.walk("."):
        for file_name in fnmatch.filter(files, '*.mzML'):
            mzml_files.append(os.path.join(root,file_name))
            #print "will have to remove ",file_name

    for each_file in mzml_files:
        os.remove(each_file)



print "All done!"
sys.exit(0)
