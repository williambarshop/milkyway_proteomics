import os, sys, re, gc
import optparse
import shutil
import subprocess
import pandas
import h5py
import numpy
import operator
from pyteomics import mzml
#from pyopenms import *
from Bio import SeqIO
from subprocess import Popen
import uniprot as uni
#####################################
#This is a script to combine the outputs of many tools into
#an Rdata image!
#
#VERSION 0.98A
version="0.98A"
#DATE: 4/05/2017
date="4/05/2017"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the Rdata merger for Galaxy, Wohlschlegel Lab UCLA"
print "Written by William Barshop"
print "Version: ",version
print "Date: ",date

basedir=os.getcwd()


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
#parser.add_option("--folder",action="store",type="string",dest="operation_folder",default=".")
parser.add_option("--skyline_chromatograms",action="store",type="string",dest="skyline_chromatograms")
parser.add_option("--skyline_boundaries",action="store",type="string",dest="skyline_boundaries")
parser.add_option("--nsaf_table",action="store",type="string",dest="nsaf_table")
parser.add_option("--spc_table",action="store",type="string",dest="spc_table")
parser.add_option("--peptide_nsaf_table",action="store",type="string",dest="peptide_nsaf_table")
parser.add_option("--peptide_spc_table",action="store",type="string",dest="peptide_spc_table")
parser.add_option("--saint_table",action="store",type="string",dest="saint_table")
parser.add_option("--fido_roc",action="store",type="string",dest="fido_roc")
parser.add_option("--fasta",action="store",type="string",dest="fasta")
parser.add_option("--job_id",action="store",type="string",dest="job_id")
parser.add_option("--mzml",action="store",type="string",dest="mzml")
parser.add_option("--msstats_comparison",action="store",type="string",dest="msstats_comparison")
parser.add_option("--msstats_condition",action="store",type="string",dest="msstats_condition")
parser.add_option("--msstats_quantification",action="store",type="string",dest="msstats_quantification")
parser.add_option("--msstats_skyline_input",action="store",type="string",dest="msstats_skyline_input")
parser.add_option("--RScriptOutput",action="store",type="string",dest="RScriptOutput")
parser.add_option("--RDataOutput",action="store",type="string",dest="RDataOutput")
parser.add_option("--MS2Output",action="store",type="string",dest="MS2Output")
parser.add_option("--MS2Output_id",action="store",type="string",dest="MS2Output_id")
parser.add_option("--experiment_file",action="store",type="string",dest="experiment_file")
parser.add_option("--psm_table",action="store",type="string",dest="psm_table")
parser.add_option("--filter_decoys",action="store_true",dest="filter_decoys")
parser.add_option("--zlib",action="store_true",dest="zlib")
parser.add_option("--renameProteinType",action="store",type="string",default="accToGene",dest="renameProteinType") #Rename proteins from accToGene, accToProtein, geneToProtein (for TGGT/ToxoDB)
parser.add_option("--analysis_type",action="store",type="string",dest="analysis_type")

(options,args) = parser.parse_args()

print "Now we're going to prepare the R script for Rdata image generation"

group_information = pandas.read_csv(options.experiment_file,sep='\t')


arguments=['python','/galaxy-central/tools/wohl-proteomics/wohl_skyline/job_history_runner.py','--job_id',str(options.job_id),"--output_file","job_history.json"]
print "Running job:"
print " ".join(arguments)
job_process=subprocess.Popen(arguments,stdout=subprocess.PIPE, stderr=subprocess.PIPE)

#### We'll read in the Skyline chromatograms data, group on the File column, and remove decoys
#decoys are duplicates of the following columns : fileName ModifiedSequence PrecursorCharge FragmentIon ProductCharge IsotopeLabelType
#let's get to it.
#
# NB NB NB NB: This code assumes that 1. Decoys come after targets. AND 2. Decoys are 1:1 with targets.  There are no targets without decoys, and no target transitions without decoy transitions.
if options.analysis_type=="lfq":
    chromatogram_df=pandas.read_csv(options.skyline_chromatograms,sep="\t",low_memory=False)

    if options.filter_decoys:
        #groups=chromatogram_df.groupby(by='FileName')
        #groups=chromatogram_df.groupby(by=['FileName','PeptideModifiedSequence','PrecursorCharge','FragmentIon','ProductCharge','IsotopeLabelType'])
        groups=chromatogram_df.groupby(by=['FileName','PeptideModifiedSequence','IsotopeLabelType'])
        filtered_df=[]
        first=True
        for groupKey,group in groups:
            filtered_df.append(group.head(n=(len(group)/2)))
        filtered_dataframe=pandas.concat(filtered_df)
    else:
        filtered_dataframe=chromatogram_df
    filtered_dataframe.to_csv("filtered_chromatograms.tsv",sep="\t",index=False)
    print "Finished outputting the filtered chromatograms..."






with open(options.fasta) as fasta_file:
    sequences = []
    identifiers = []
    descriptions = []
    lengths=[]
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        identifiers.append(seq_record.id)
        sequences.append(str(seq_record.seq))
        descriptions.append(seq_record.description.replace("\t",","))
        lengths.append(len(seq_record.seq))
fasta_df=pandas.DataFrame({'id':identifiers,'description':descriptions,'Sequence':sequences})

fasta_df.to_csv("fasta_df.tsv",sep="\t",index=False)

#We're going to read in the MS2 scans via Pyteomics, and then store them in a pandas df.
#Then we'll write them to a TSV file <----- New version implements HDF5 instead to keep it out of memory!
mzml_files=options.mzml.split("___")
mzml_file_scans=[]
ms2fail=False
first=True
h5file=h5py.File(options.MS2Output,'w')

#group_information    
for each_file in mzml_files:
    if first:
        first=False
    else:
        h5file.close()
        gc.collect()
        h5file=h5py.File(options.MS2Output,'a')
    #this_group=h5file.create_group(str(group_information.loc[group_information['Original File Name']==each_file.rsplit(".",1)[0],'Crux File Integer'].iloc[0]))  #This makes groups by crux integer
    this_group=h5file.create_group(each_file)
    try:
        print "Reading in ",each_file
        mzml_reader=mzml.read(each_file,iterative=True)
        for each_scan in mzml_reader:
            if each_scan['ms level']==2:
                this_dataset=this_group.create_dataset(str(each_scan['index']),compression="gzip",compression_opts=9,dtype="float32",data=numpy.column_stack((each_scan['m/z array'],each_scan['intensity array'])).T)
                this_dataset.attrs["scan_index"]=each_scan['index']
        print "done with it........\n"
        del mzml_reader
        del each_scan
        h5file.flush()
    except:
        print "Failed to read in file",each_file
        ms2fail=True
#if not ms2fail:
    #scans_df=pandas.DataFrame.from_dict(mzml_file_scans)
    #scans_df.to_csv("ms2_scans.tsv",sep="\t",index=False)
    #del scans_df

h5file.close()
#job_process.wait()
#out, err = job_process.communicate()
#print out,err
#print "Finished the job processing."


job_process.wait()
out, err = job_process.communicate()
print out,err
print "Finished the job processing."


print "Handling the renaming!" #This is done only for both the PSM_table and the MSstats outputs...

psm_table=pandas.read_csv(options.psm_table,sep='\t',low_memory=False)


### The code below is a mixture of things that we need in LFQ and Qual analysis...
### It is meant to query uniprot, and generate name mappings for the MSstats csv table, and the PSM table.
if options.analysis_type=="lfq":
    combined_results=pandas.read_csv(options.msstats_comparison,sep=',',low_memory=False)
    #combined_results=pandas.read_csv(options.msstats_comparison,sep=',',index_col=False)
    print combined_results
    combined_results['backup']=combined_results['Protein']

protein_series=psm_table['protein id']
uni_mapping_dict={}
uniprot_pattern=r'([A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]|[OPQ][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9])'
protein_extracted=protein_series.str.extractall(uniprot_pattern)[0].unique()#for the PSM table...
protein_full_name_series_list=protein_series.str.split(",",expand=False) #This will be used later to properly adjust positioning of AA's...


if "accToGene" in options.renameProteinType or options.renameProteinType is None: #This means we'll take the acc's and replace with GENE NAMES
    if options.analysis_type=="lfq":
        combined_results['uniprot_acc']=combined_results['Protein'].str.extract(uniprot_pattern,expand=False)
    #uni_mapping_dict_combined=uni.map(list(combined_results['uniprot_acc'].unique()),f='ACC',t='GENENAME')
    uni_mapping_dict_combined=uni.map(list(protein_extracted),f='ACC',t='GENENAME')
    for each_acc in uni_mapping_dict_combined:
        catch=uni_mapping_dict_combined[each_acc]
        if len(catch)>0:
            uni_mapping_dict[str(each_acc)]=str(catch.pop())
    second_query=[]
    for each_acc in protein_extracted:
        if each_acc not in uni_mapping_dict_combined.keys():
            second_query.append(each_acc)
    if len(second_query)>0:
        uni_mapping_dict_two=uni.map(list(second_query),f='ACC',t='ID')
        for each_acc in uni_mapping_dict_two:
            try:
                catch=uni_mapping_dict_combined[each_acc]
            except:
                uni_mapping_dict[str(each_acc)]=each_acc #This is a failed conversion!
                pass

            if len(catch)>0:
                #print catch,type(catch)
                uni_mapping_dict[str(each_acc)]=str(catch.pop())
    if options.analysis_type=="lfq":    
        combined_results['Gene ID']=combined_results['uniprot_acc'].map(uni_mapping_dict)
        combined_results.drop('uniprot_acc',axis=1,inplace=True)
        missing_protein_names=combined_results[combined_results['Protein'].isnull()]
        for each_key,each_row in missing_protein_names.iterrows():
            combined_results.loc[each_key,'Gene ID']=combined_results.loc[each_key,'backup']



elif "accToProtein" in options.renameProteinType:  # THIS MEANS WE'LL TAKE AND REPLACE WITH ProteinS
    combined_results['uniprot_acc']=combined_results['Protein'].str.extract(uniprot_pattern)
    #uni_mapping_dict_combined=uni.map(list(combined_results['uniprot_acc'].unique()),f='ACC',t='ID')
    uni_mapping_dict_combined=uni.map(list(protein_extracted),f='ACC',t='ID')
    for each_acc in uni_mapping_dict_combined:
        catch=uni_mapping_dict_combined[each_acc]
        if len(catch)>0:
            uni_mapping_dict[str(each_acc)]=str(catch.pop())
    second_query=[]
    for each_acc in protein_extracted:
        if each_acc not in uni_mapping_dict_combined.keys():
            second_query.append(each_acc)
    if len(second_query)>0:
        uni_mapping_dict_two=uni.map(list(second_query),f='ACC',t='GENENAME')
        for each_acc in uni_mapping_dict_two:
            try:
                catch=uni_mapping_dict_combined[each_acc]
            except:
                uni_mapping_dict[str(each_acc)]=each_acc #This is a failed conversion!
                pass
            if len(catch)>0:
                uni_mapping_dict[str(each_acc)]=str(catch.pop())

    if options.analysis_type=="lfq":    
        combined_results['Gene ID']=combined_results['uniprot_acc'].map(uni_mapping_dict)
        combined_results.drop('uniprot_acc',axis=1,inplace=True)
        missing_protein_names=combined_results[combined_results['Protein'].isnull()]
        for each_key,each_row in missing_protein_names.iterrows():
            combined_results.loc[each_key,'Gene ID']=combined_results.loc[each_key,'backup']


elif "geneToProtein" in options.renameProteinType:  # THIS MEANS WE'LL TAKE AND REPLACE WITH PROTEIN NAMES
    #protein_full_name_series_list.apply(pd.Series).stack().reset_index(drop=True)
    #protein_full_name_series_collapsed=protein_full_name_series_list.apply(pandas.Series).stack().reset_index(drop=True)
    combined_results['uniprot_acc']=combined_results['Protein'].str.extract(uniprot_pattern)
#    uni_mapping_dict_combined=uni.map(list(combined_results['uniprot_acc'].unique()),f='GENENAME',t='ID')
    uni_mapping_dict_combined=uni.map(list(protein_extracted),f='GENENAME',t='ID')
    for each_acc in uni_mapping_dict_combined:
        catch=uni_mapping_dict_combined[each_acc]
        if len(catch)>0:
            #print catch,type(catch)
            uni_mapping_dict[str(each_acc)]=str(catch.pop())
    second_query=[]
    for each_acc in protein_extracted:
        if each_acc not in uni_mapping_dict_combined.keys():
            second_query.append(each_acc)
    if len(second_query)>0:
        uni_mapping_dict_two=uni.map(list(second_query),f='ACC',t='ID')
        for each_acc in uni_mapping_dict_two:
            try:
                catch=uni_mapping_dict_combined[each_acc]
            except:
                uni_mapping_dict[str(each_acc)]=each_acc #This is a failed conversion!
                pass
            if len(catch)>0:
                uni_mapping_dict[str(each_acc)]=str(catch.pop())

    if options.analysis_type=="lfq":    
        combined_results['Gene ID']=combined_results['uniprot_acc'].map(uni_mapping_dict)
        combined_results.drop('uniprot_acc',axis=1,inplace=True)
        missing_protein_names=combined_results[combined_results['Protein'].isnull()]
        for each_key,each_row in missing_protein_names.iterrows():
            combined_results.loc[each_key,'Gene ID']=combined_results.loc[each_key,'backup']


elif "norename" in options.renameProteinType:
    print "Not renaming proteins... leaving as is!"

if options.analysis_type=="lfq":    
    combined_results.drop('backup',axis=1,inplace=True)
    combined_results.to_csv('msstats_comparison.csv',sep=',',index=False)

#We're going to handle the phospho mods, if they exist, and if they do, make a column for KSEA analysis.

#if 'phospho_numMods' in psm_table.columns:
for each_column_name in [x for x in psm_table.columns if "_numMods" in x and "ptmRS_numMods" not in x]:
    this_mod_name=each_column_name.rsplit("_",1)[0] # This will be the mod name.
    psm_table[this_mod_name+'_numMods']=psm_table[this_mod_name+'_numMods'].astype(float).fillna(0.0)
    #if this_mod_name="
    #psm_table['KSEA sites']=""
    working_set=psm_table[psm_table[this_mod_name+'_numMods']>0.0]
    for index,each_row in working_set.iterrows():
        site_prob=each_row[this_mod_name+'_siteProbabilities']

        this_mod_site_scores_split=each_row[this_mod_name+"_siteProbabilities"].split(';')
        this_mod_site_dict={}
        #this_mod_number_of_sites=len(this_mod_site_scores_split)
        this_mod_number_of_sites=int(each_row[this_mod_name+'_numMods'])
        for eachsite in this_mod_site_scores_split:
            eachsite=eachsite.strip()
            score=float(eachsite.split()[1])
            factor=int(eachsite.split()[0].rsplit('x',1)[1].replace(":","").strip()) #We don't need this for the phospho case, but we'll leave it anyway and instead place the AA in, instead.
            site_index=str(int(eachsite.split()[0].split('(')[1].split(')')[0])-1)+"_"+str(eachsite[0]) #AA stored first here.   This is formatted as "7_Y" which would be the index and the AA modified separated by "_"
            this_mod_site_dict[site_index]=score
        this_mod_sorted_sites = sorted(this_mod_site_dict.items(), key=operator.itemgetter(1), reverse=True)
        this_mod_sorted_sites=this_mod_sorted_sites[:this_mod_number_of_sites]
        #print this_mod_sorted_sites
        #sys.exit(2)
        this_mod_sorted_sites=[(x[0].split("_")[0],x[0].split("_")[1]) for x in this_mod_sorted_sites]
        mod_str=""
        mod_logo=""
        logo_length=6# AA's on each side of the modification site...
        for eachsite in this_mod_sorted_sites:
            AA_position=eachsite[0]
            mod_AA=eachsite[1]

            for each_protein in each_row['proteinacc_start_stop_pre_post_;'].split(";"):
                all_uniprots=set(re.findall(uniprot_pattern,each_protein.rsplit("_",4)[0]))
                all_mappings=[]
                for each_uniprot in all_uniprots:
                    id_filter=fasta_df[fasta_df['id']==each_uniprot]
                    if len(id_filter)==0:
                        #print "WARNING: Couldn't find ",each_uniprot,"in the FASTA db"
                        id_filter=fasta_df[fasta_df['id']==each_protein.rsplit("_",4)[0]]
                        #if len(id_filter)==0:
                        #    print "WARNING: Still can't find",each_protein.rsplit("_",4)[0],"in the FASTA..."
                    for id_index,id_each_row in id_filter.iterrows():
                        prot_index=int(AA_position)+int(each_protein.rsplit("_",4)[1])-1
                        prot_length=len(id_each_row['Sequence'])
                        if prot_index<logo_length-1:
                            n_term_offset="-"*(logo_length-prot_index)
                            mod_logo+=each_uniprot+":"+n_term_offset+id_each_row['Sequence'][prot_index-(logo_length-(logo_length-prot_index)):prot_index+logo_length+1]+";"
                        elif (prot_index+logo_length)>=prot_length:
                            offset_buffer=logo_length-(prot_length-1-prot_index)
                            c_term_offset="-"*offset_buffer
                            mod_logo+=each_uniprot+":"+id_each_row['Sequence'][prot_index-logo_length:prot_index+(logo_length-(logo_length-((prot_length)-prot_index)))]+c_term_offset+";"
                            #print mod_logo
                            #print prot_length,prot_index
                        else:
                            #mod_logo+=each_uniprot+":"+id_each_row['Sequence'][prot_index-logo_length:prot_index+logo_length]+c_term_offset+";"
                            mod_logo+=each_uniprot+":"+id_each_row['Sequence'][prot_index-logo_length:prot_index+logo_length+1]+";"
                        #id_each_row['Sequence']
                    #print mod_logo
                    #these_uniprots=set(re.findall(uniprot_pattern,each_protein.rsplit("_",4)[0]))
                    try:
                        #for each_these_uniprot in these_uniprots:
                        uni_mappings=uni_mapping_dict[each_uniprot]
                    except:
                        #for each_these_uniprot in these_uniprots:
                        #uni_mappings=each_uniprot
                        uni_mappings=each_protein.rsplit("_",4)[0]
                    all_mappings.append(uni_mappings)
                for each_mapping in set(all_mappings): # This will always only have one, unless we change the append above to extend, and make the uni_mapping_dict into a list storage...
                    mod_str+=each_mapping+"_"+str(mod_AA)+str(int(AA_position)+int(each_protein.rsplit("_",4)[1]))+";"
        mod_str=mod_str[:-1]
        mod_str=';'.join(set([x for x in mod_str.split(";")])) #cleanup!
        mod_logo=mod_logo[:-1]
        psm_table.set_value(index,col=this_mod_name+' prot-sites',value=mod_str)
        psm_table.set_value(index,col=this_mod_name+' motifs',value=mod_logo)


psm_table.to_csv('psm_table.tsv',sep='\t',index=False)

with open("Rdata_Script.R",'wb') as script_writer:
    #script_writer.write("comparison_csv<-read.csv(\"{0}\",check.names=FALSE)\n".format(options.msstats_comparison))
    script_writer.write("experiment_design<-read.csv(\"{0}\",sep=\"\\t\",check.names=FALSE)\n".format(options.experiment_file))
    script_writer.write("psm_table<-read.csv(\"{0}\",sep=\"\\t\",check.names=FALSE)\n".format('psm_table.tsv'))
    #script_writer.write("psm_table<-read.csv(\"{0}\",sep=\"\\t\",check.names=FALSE)\n".format(options.psm_table))
    #script_writer.write("chromatograms<-read.csv(\"{0}\",sep=\"\\t\",check.names=FALSE)\n".format(options.skyline_chromatograms))
    script_writer.write("spc_table<-read.csv(\"{0}\",check.names=FALSE,sep=\"\\t\")\n".format(options.spc_table))
    script_writer.write("nsaf_table<-read.csv(\"{0}\",check.names=FALSE,sep=\"\\t\")\n".format(options.nsaf_table))
    script_writer.write("peptide_spc_table<-read.csv(\"{0}\",check.names=FALSE,sep=\"\\t\")\n".format(options.peptide_spc_table))
    script_writer.write("peptide_nsaf_table<-read.csv(\"{0}\",check.names=FALSE,sep=\"\\t\")\n".format(options.peptide_nsaf_table))
    script_writer.write("fido_roc<-read.csv(\"{0}\",check.names=FALSE)\n".format(options.fido_roc))
    script_writer.write("ms2_dataset_id<-\"{0}\"\n".format(options.MS2Output_id))
    #script_writer.write("library(Biostrings)\n")
    script_writer.write("fasta_df <- read.csv(\"{0}\",check.names=FALSE,sep=\"\\t\")\n".format("fasta_df.tsv"))
    script_writer.write("library(rjson)\n")
    script_writer.write("job_history <- fromJSON(file=file(\"job_history.json\"))\n")


    #need to handle using qual for the string on qual compare... the SAINT table is suppressed on the MilkyWay Shiny side, so we treat them as a single type of analysis for the viewer.
    if options.analysis_type=="qualCompare":
        script_writer.write("analysis_type<-\"{0}\"\n".format('qual'))
    else:
        script_writer.write("analysis_type<-\"{0}\"\n".format(options.analysis_type))




    if options.analysis_type=="lfq":
        script_writer.write("chromatograms<-read.csv(\"{0}\",sep=\"\\t\",check.names=FALSE)\n".format("filtered_chromatograms.tsv"))
        script_writer.write("comparison_csv<-read.csv(\"{0}\",check.names=FALSE)\n".format('msstats_comparison.csv'))
        script_writer.write("condition_plot_csv<-read.csv(\"{0}\",check.names=TRUE)\n".format(options.msstats_condition))
        script_writer.write("boundaries_csv<-read.csv(\"{0}\",check.names=FALSE)\n".format(options.skyline_boundaries))
        script_writer.write("quantification_csv<-read.csv(\"{0}\",check.names=FALSE)\n".format(options.msstats_quantification))
        script_writer.write("msstats_skyline_input<-read.csv(\"{0}\",check.names=FALSE)\n".format(options.msstats_skyline_input))
        script_writer.write("saint_table<-read.csv(\"{0}\",check.names=FALSE,sep=\"\\t\")\n".format(options.saint_table))

    elif options.analysis_type=="qualCompare":
        script_writer.write("saint_table<-read.csv(\"{0}\",check.names=FALSE,sep=\"\\t\")\n".format(options.saint_table))

    elif options.analysis_type=="qual":
        pass

    #script_writer.write("job_history <- read.csv(\"job_history.tsv\",sep=\"\\t\",check.names=FALSE)\n")
    #script_writer.write("job_history <- stream_in(file(\"job_history.json\"))\n")
    #script_writer.write("library(jsonlite)\n")
    #script_writer.write("job_history <- stream_in(file(\"job_history.json\"))\n")
    #if not ms2fail:
    #    script_writer.write("ms2_scans <- read.csv(\"ms2_scans.tsv\",sep=\"\\t\",check.names=FALSE)\n")
    #else:
    #    print "CAUTION::: MS2 SCANS ARE NOT BEING INCLUDED DUE TO FAILURES IN READING."
    #script_writer.write("library(Biostrings)\n")


    script_writer.write("save.image(file=\"{0}\")\n".format(options.RDataOutput))
    #script_writer.write("save.image(file=\"{0}\")\n".format('image.RData'))


#OKAY.... The R Script has been written!
#We're going to execute the R script now!
print "Copying RScript back to Galaxy..."
shutil.copy('Rdata_Script.R',options.RScriptOutput)
subprocess.check_call(['Rscript', 'Rdata_Script.R'],shell=False,stderr=sys.stdout.fileno())


#once that's done, we'll load the Rdata and cram in the chromatograms...
#from rpy2.robjects import pandas2ri
#pandas2ri.activate()
#import rpy2.robjects as robjects
#import pandas.rpy.common as com
#from rpy2.robjects.packages import importr

#robjects.r['load']("image.RData")
#filtered_chromatogram_df
#chromatogram_df = com.convert_to_r_dataframe(filtered_chromatogram_df)
#chromatograms = pandas2ri.py2ri(filtered_chromatogram_df)
#base=importr('base')
#base.save_image("final_rdata.RData")


#print "Moving Rdata image to final output location...."
#shutil.copy('image.RData',options.RDataOutput)
#shutil.copy('final_rdata.RData',options.RDataOutput)

print "All done!"
