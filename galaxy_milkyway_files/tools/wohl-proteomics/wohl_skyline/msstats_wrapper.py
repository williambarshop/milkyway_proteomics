import os, sys, re
import optparse
import shutil
import pandas
import numpy
import gc
import subprocess
import uniprot as uni
from natsort import natsorted, ns
import warnings

#From stack overflow, to redirect pandas stderr warnings to stdout for Galaxy's sake....#
def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))

warnings.showwarning = customwarn
#####################################

#####################################
#This is a script to combine the output reports from
#Skyline, in preparation for MSstats!  Let's get started.
#
#VERSION 0.91D
version="0.91D"
#DATE: 10/01/2017
date="10/01/2017"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the MSstats wrapper for Galaxy, Wohlschlegel Lab UCLA"
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
parser.add_option("--folder",action="store",type="string",dest="operation_folder",default=".")
parser.add_option("--galaxy-csv",action="store",type="string",dest="galaxy_csv")
parser.add_option("--experiment_file",action="store",type="string",dest="experiment_file")
parser.add_option("--remove_decoys",action="store_true",dest="remove_decoys")
parser.add_option("--remove_precursors",action="store_true",dest="remove_precursors")
parser.add_option("--rename",action="store_true",dest="rename") #Rename columns.
parser.add_option("--renameProteinType",action="store",type="string",default="accToGene",dest="renameProteinType") #Rename proteins from accToGene, accToProtein, geneToProtein (for TGGT/ToxoDB)
parser.add_option("--remove_empty",action="store_true",dest="remove_empty") 
parser.add_option("--merge_isotopes",action="store_true",dest="merge_isotopes") #Merge isotopes.
parser.add_option("--peptide_level",action="store_true",dest="peptide_level") #no proteins, just peptides.
parser.add_option("--zeros",action="store_true",dest="zeros") # 0 ---> np.nan
parser.add_option("--zero_to_one",action="store_true",dest="zero_to_one") # 0 ---> 1
parser.add_option("--remove_truncated_peaks",action="store_true",dest="remove_truncated_peaks")
parser.add_option("--remove_single_run",action="store_true",dest="remove_single_run") #Makes it so that you have to have measurement from >1 run.
parser.add_option("--mprophet_q",action="store",type="float",dest="mprophet_q") #We'll throw out things above this q-value threshold.
parser.add_option("--fractionated",action="store_true",dest="fractionated")  # Using the final -??? as the fractionation identifier
parser.add_option("--remove_repeated_peptides",action="store_true",dest="remove_repeated_peptides")
parser.add_option("--removeProteinsByText",action="store",type="string",dest="remove_proteins_by_text")

parser.add_option("--featureSubset",action="store",type="string",dest="feature_subset")
parser.add_option("--featureSubsetN",action="store",type="int",dest="feature_subset_N")   # For use with "topN"
parser.add_option("--remove_proteins_with_interference",action="store_true",dest="remove_interfered_proteins")
parser.add_option("--censoredInt",action="store",type="string",dest="censoredInt",default="0")
parser.add_option("--cutoffCensored",action="store",type="string",dest="cutoffCensored",default="minFeature") # Can be "minFeature", "minFeatureNRun", "minRun"
parser.add_option("--fillIncompleteRows",action="store_true",dest="fillIncompleteRows")  # If data is missing, we'll put it back in with "NA" values
parser.add_option("--MBimpute",action="store_true",dest="MBimpute")
parser.add_option("--noimpute",action="store_true",dest="noimpute")
parser.add_option("--remove50missing",action="store_true",dest="remove50missing")
parser.add_option("--normalization",action="store",type="string",dest="normalization")   # Can be "Quantile", "equalizeMedians", "globalStandards"
parser.add_option("--normalization_protein",action="store",type="string",dest="normalize_protein") #Can be a comma separated list, but must match the FASTA header...
parser.add_option("--maxQuantileforCensored",action="store",type="string",dest="maxQuantileforCensored")   # Can be used to emulate 3.4 behavior
parser.add_option("--fillMissingFeatures",action="store_true",dest="fillMissingFeatures")

################# OUTPUTS ################################
parser.add_option("--processedDataOutput",action="store",type="string",dest="processedOutput")
parser.add_option("--comparisonOutput",action="store",type="string",dest="comparisonOutput")
parser.add_option("--RDataOutput",action="store",type="string",dest="RDataOutput")
parser.add_option("--QCplotOutput",action="store",type="string",dest="QCplotOutput")
parser.add_option("--ProfilePlotOutput",action="store",type="string",dest="profilePlotOutput")
parser.add_option("--ConditionPlotOutput",action="store",type="string",dest="conditionPlotOutput")
parser.add_option("--RScriptOutput",action="store",type="string",dest="RScriptOutput")
parser.add_option("--quantificationOutput",action="store",type="string",dest="quantificationOutput")
parser.add_option("--quantificationConditionOutput",action="store",type="string",dest="quantificationConditionOutput")
parser.add_option("--conditionPlotCSVOutput",action="store",type="string",dest="conditionPlotCSVOutput")


################## BELOW THIS ARE PLOTTING OPTIONS ##############################   These are actually all going to be moved into a separate tool
#parser.add_option("--significance",action="store",type="float",dest="significance") # For the volcano plots...
#parser.add_option("--FCthreshold",action="store",type="float",dest="FCthreshold") # FC threshold For the volcano plots...
#parser.add_option("--ylimUp",action="store",type="float",dest="ylimUp") # ylimUp threshold for the plots
#parser.add_option("--ylimDown",action="store",type="float",dest="ylimDown") # ylimDown threshold for plots
#parser.add_option("--xlimUp",action="store",type="float",dest="xlimUp") # xlimUp threshold for Volcano plots
#parser.add_option("--numProtein",action="store",type="int",dest="numProtein",default=180) # Number of proteins per heatmap... Max is 180
#parser.add_option("--clustering",action="store",type="string",dest="clustering",default="protein") # clustering type for heatmap... Can be "protein", "comparison", "both"
#parser.add_option("--proteinName",action="store_true",dest="proteinName") # On volcano plot, draw protein names?

#parser.add_option("--dotSize",action="store",type="int",dest="dotSize",default=3)
#parser.add_option("--textSize",action="store",type="int",dest="textSize",default=4)
#parser.add_option("--legendSize",action="store",type="int",dest="legendSize",default=7)
#parser.add_option("--width",action="store",type="int",dest="width",default=10)
#parser.add_option("--height",action="store",type="int",dest="height",default=10)




(options,args) = parser.parse_args()

def help():
    print "Syntax is as follows:"
    print "python msstats_merger.py --folder=<input folder> --remove_decoys --rename"
    print "--folder is a required input, it is optional to remove decoys from the MSstats input file..."

if options.operation_folder is None:
    print "Please provide a --folder=<input folder> argument!"
    print "It is necessary so that I can run..."
    help()
    sys.exit(1)

def appendFraction(x):
    try:
        #appended_name=str(x['Peptide Modified Sequence']+"_"+str(x['File Name'].split(".")[0].rsplit("-",1)[1]))# MAY HAVE TO CHANGE x['FILE NAME'] TO STR(x['FILE NAME']).... !!!!!!!!!!!! # 3/3/2016 -- BACK TO THIS... See: https://groups.google.com/forum/#!searchin/msstats/multiple$20methods/msstats/ZzP3Q8hGXBY/oTYo60cfovMJ
        appended_name=str(x['Peptide Modified Sequence']+"-Frac"+str(x['File Name'].split(".")[0].rsplit("-",1)[1]))# MAY HAVE TO CHANGE x['FILE NAME'] TO STR(x['FILE NAME']).... !!!!!!!!!!!!   # 3/3/2016 changed to - to make sure we aren't screwing with MSstats
        return appended_name
    except:
        print "FAILED ON CORRECTING FRACTION NAMES ON",x
        sys.exit(0)

def fixFileName(x):
    return str(x['File Name'].split('.')[0].rsplit("-",1)[0])#[:-3])

def peptide_level_fixer(x):
    return x.split("_")[0]
    #return x['Peptide Modified Sequence'].split("_")[0]

os.chdir(options.operation_folder)
if options.galaxy_csv is not None:
    shutil.copy(options.galaxy_csv,"MSstats_input.csv")
if os.path.isfile("MSstats_combined_input.csv"):
    os.remove("MSstats_combined_input.csv")
    print "Deleted MSstats_combined_input.csv to make room for a new one..."
os.chdir(basedir)


    
infiles = []
for file in os.listdir(options.operation_folder):
    if file.endswith(".csv"):
        if "MSstats" in file:
            infiles.append(file)
        else:
            print "ignoring file ",str(file)," because it does not contain \"MSstats\" in the file name"

#for root, subFolders, files in os.walk(options.operation_folder):
#    for eachfile in files:
#        if '.csv' in eachfile[-4:]:#We'll only look for files which end with ".csv", which is the expected file format for the MSstats input files!
#            infiles.append(str(os.path.join(root,eachfile)))
print "about to read csv..."
dataframe_vector=[]
os.chdir(options.operation_folder)
for eachfile in infiles:
    newdf=pandas.read_csv(eachfile,sep=',',index_col=False)#Used to have low_memory=False)
    dataframe_vector.append(newdf)
    del newdf
    gc.collect()
os.chdir(basedir)
print "done reading csvs..."
#Peptide Modified Sequence <----- This will be the column to edit... and it'll be based on "File Name".split(".")[0].rsplit("-",1)[1] which should be like "F1" or "F2", etc...
combined_results=pandas.concat(dataframe_vector)
if options.remove_decoys:
    combined_results.rename(columns={'Protein Name':'ProteinName'},inplace=True)
    mask=combined_results[numpy.invert(combined_results.ProteinName.str.contains("Decoys"))]
    combined_results=mask
    combined_results.rename(columns={'ProteinName':'Protein Name'},inplace=True)
if options.remove_precursors:
    combined_results.rename(columns={'Fragment Ion':'FragmentIon'},inplace=True)
    mask=combined_results[numpy.invert(combined_results.FragmentIon.str.contains("precursor"))]
    combined_results=mask
    combined_results.rename(columns={'FragmentIon':'Fragment Ion'},inplace=True)
if options.fractionated:
    fixed_peptide_names=combined_results.apply(appendFraction,axis=1)    
    fixed_file_names=combined_results.apply(fixFileName,axis=1)   
    combined_results['Peptide Modified Sequence']=fixed_peptide_names
    combined_results['File Name']=fixed_file_names
    combined_results=combined_results.copy(deep=True)

if options.normalization=="globalStandards":
    protein_norm_dict={}
    for each_item in options.normalize_protein.split(","):
        if "::" in each_item:
            peptide_list=[item for item in each_item.split("::")[1:]]
            protein_norm_dict[each_item.split("::")[0]]=peptide_list
        else:
            protein_norm_dict[each_item]=[]
    for each_protein in protein_norm_dict:
        this_pep_list=protein_norm_dict[each_protein]
        if len(this_pep_list)==0:
            #combined_results##############################################
            combined_results.loc[combined_results['Protein Name'].str.contains(each_protein),"Standard Type"]="globalStandard"
        else:
            for each_peptide in this_pep_list:
                combined_results.loc[combined_results['Protein Name'].str.contains(each_protein) & combined_results['Peptide Modified Sequence']==each_peptide,"Standard Type"]="globalStandard"
print "about to remove proteins by text..."
if options.remove_proteins_by_text is not None and options.remove_proteins_by_text is not "":
    proteins_to_remove=options.remove_proteins_by_text.split(",")
    for each_protein in proteins_to_remove:
        each_protein=each_protein.strip()
        print "Removing protein",each_protein,"from the analysis."
        combined_results=combined_results[numpy.invert(combined_results['Protein Name'].str.contains(each_protein))]


#combined_results.sort(columns='Protein Name',inplace=True)
combined_results.sort_values(by='Protein Name',inplace=True)


############## We're going to try to extract uniprot acc's from each protein name!
#######
#combined_results['gene_name']=
#uni.map([x for x in combined_results['uniprot_acc'] if x is not numpy.nan])

#combined_results.loc[combined_results['uniprot_acc'] == numpy.nan]['uniprot_acc']=combined_results['Protein Name']

#combined_results['Protein Name']=combined_results.apply(lambda x: 
combined_results['backup']=combined_results['Protein Name']
uni_mapping_dict={}
uniprot_pattern=r'([A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]|[OPQ][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9])'
if "accToGene" in options.renameProteinType or options.renameProteinType is None: #This means we'll take the acc's and replace with GENE NAMES
    combined_results['uniprot_acc']=combined_results['Protein Name'].str.extract(uniprot_pattern)
    #print combined_results['uniprot_acc']
    uni_mapping_dict_combined=uni.map(list(combined_results['uniprot_acc'].unique()),f='ACC',t='GENENAME')
    for each_acc in uni_mapping_dict_combined:
        catch=uni_mapping_dict_combined[each_acc]
        #catch=uni.map(each_acc,f='ACC',t='GENENAME')
        if len(catch)>0:
            #print catch,type(catch)
            uni_mapping_dict[str(each_acc)]=str(catch.pop())
        #else:
        #    catch=uni.map(each_acc,f='ACC',t='ID')
        #    if len(catch)>0:
        #        uni_mapping_dict[str(each_acc)]=str(catch[each_acc].pop())
        #    else:
        #        uni_mapping_dict[str(each_acc)]=str(each_acc)
    second_query=[]
    for each_acc in combined_results['uniprot_acc'].unique():
        if each_acc not in uni_mapping_dict_combined.keys():
            second_query.append(each_acc)
    if len(second_query)>0:
        uni_mapping_dict_two=uni.map(list(second_query),f='ACC',t='ID')
        for each_acc in uni_mapping_dict_two:
            try:
                catch=uni_mapping_dict_combined[each_acc]
            except:
                pass
                #uni_mapping_dict[str(each_acc)]=str(each_acc)
            #catch=uni.map(each_acc,f='ACC',t='GENENAME')
            if len(catch)>0:
                #print catch,type(catch)
                uni_mapping_dict[str(each_acc)]=str(catch.pop())

    combined_results['Protein Name']=combined_results['uniprot_acc'].map(uni_mapping_dict)
    combined_results.drop('uniprot_acc',axis=1,inplace=True)
    missing_protein_names=combined_results[combined_results['Protein Name'].isnull()]
    for each_key,each_row in missing_protein_names.iterrows():
        combined_results.loc[each_key,'Protein Name']=combined_results.loc[each_key,'backup']



elif "accToProtein" in options.renameProteinType:  # THIS MEANS WE'LL TAKE AND REPLACE WITH PROTEIN NAMES
    combined_results['uniprot_acc']=combined_results['Protein Name'].str.extract(uniprot_pattern)
    #print combined_results['uniprot_acc']
    uni_mapping_dict_combined=uni.map(list(combined_results['uniprot_acc'].unique()),f='ACC',t='ID')
    for each_acc in uni_mapping_dict_combined:
        catch=uni_mapping_dict_combined[each_acc]
        #catch=uni.map(each_acc,f='ACC',t='GENENAME')
        if len(catch)>0:
            #print catch,type(catch)
            uni_mapping_dict[str(each_acc)]=str(catch.pop())
        #else:
        #    catch=uni.map(each_acc,f='ACC',t='ID')
        #    if len(catch)>0:
        #        uni_mapping_dict[str(each_acc)]=str(catch[each_acc].pop())
        #    else:
        #        uni_mapping_dict[str(each_acc)]=str(each_acc)
    second_query=[]
    for each_acc in combined_results['uniprot_acc'].unique():
        if each_acc not in uni_mapping_dict_combined.keys():
            second_query.append(each_acc)
    if len(second_query)>0:
        uni_mapping_dict_two=uni.map(list(second_query),f='ACC',t='GENENAME')
        for each_acc in uni_mapping_dict_two:
            try:
                catch=uni_mapping_dict_combined[each_acc]
            except:
                pass
                #uni_mapping_dict[str(each_acc)]=str(each_acc)
            #catch=uni.map(each_acc,f='ACC',t='GENENAME')
            if len(catch)>0:
                #print catch,type(catch)
                uni_mapping_dict[str(each_acc)]=str(catch.pop())

    combined_results['Protein Name']=combined_results['uniprot_acc'].map(uni_mapping_dict)
    combined_results.drop('uniprot_acc',axis=1,inplace=True)
    missing_protein_names=combined_results[combined_results['Protein Name'].isnull()]
    for each_key,each_row in missing_protein_names.iterrows():
        combined_results.loc[each_key,'Protein Name']=combined_results.loc[each_key,'backup']


elif "geneToProtein" in options.renameProteinType:  # THIS MEANS WE'LL TAKE AND REPLACE WITH PROTEIN NAMES
    combined_results['uniprot_acc']=combined_results['Protein Name'].str.extract(uniprot_pattern)
    #print combined_results['uniprot_acc']
    uni_mapping_dict_combined=uni.map(list(combined_results['uniprot_acc'].unique()),f='GENENAME',t='ID')
    for each_acc in uni_mapping_dict_combined:
        catch=uni_mapping_dict_combined[each_acc]
        #catch=uni.map(each_acc,f='ACC',t='GENENAME')
        if len(catch)>0:
            #print catch,type(catch)
            uni_mapping_dict[str(each_acc)]=str(catch.pop())
        #else:
        #    catch=uni.map(each_acc,f='ACC',t='ID')
        #    if len(catch)>0:
        #        uni_mapping_dict[str(each_acc)]=str(catch[each_acc].pop())
        #    else:
        #        uni_mapping_dict[str(each_acc)]=str(each_acc)
    second_query=[]
    for each_acc in combined_results['uniprot_acc'].unique():
        if each_acc not in uni_mapping_dict_combined.keys():
            second_query.append(each_acc)
    if len(second_query)>0:
        uni_mapping_dict_two=uni.map(list(second_query),f='ACC',t='ID')
        for each_acc in uni_mapping_dict_two:
            try:
                catch=uni_mapping_dict_combined[each_acc]
            except:
                pass
                #uni_mapping_dict[str(each_acc)]=str(each_acc)
            #catch=uni.map(each_acc,f='ACC',t='GENENAME')
            if len(catch)>0:
                #print catch,type(catch)
                uni_mapping_dict[str(each_acc)]=str(catch.pop())

    combined_results['Protein Name']=combined_results['uniprot_acc'].map(uni_mapping_dict)
    combined_results.drop('uniprot_acc',axis=1,inplace=True)
    missing_protein_names=combined_results[combined_results['Protein Name'].isnull()]
    for each_key,each_row in missing_protein_names.iterrows():
        combined_results.loc[each_key,'Protein Name']=combined_results.loc[each_key,'backup']


elif "norename" in options.renameProteinType:
    print "Not renaming proteins... leaving as is!"

combined_results.drop('backup',axis=1,inplace=True)

#combined_results.to_csv("testing_output_BEFORE.csv")

#print uni_mapping_dict
#try:
#    pass
#except:
#    pass

#for index,row in combined_results.iterrows():
#    if combined_results.loc[index,'uniprot_acc'] is numpy.nan:
#        pass
#    else:
#       combined_results.loc[index,'Protein Name']=uni.map(row['uniprot_acc'],f='ACC',t='GENENAME')


#print uni_mapping_dict
#print combined_results[combined_results['Peptide Modified Sequence']=="INTQWLLTSGTTEANAWK"]

#combined_results.to_csv("testing_output.csv")

#sys.exit(0)






if options.mprophet_q is not None:
    #We'll take the results below the Q-value cutoff and set them to zero
    #combined_results=combined_results[combined_results['annotation_QValue']<=options.mprophet_q]
    #combined_results=combined_results.copy(deep=True)  # BAD
    #combined_results[combined_results['annotation_QValue']>options.mprophet_q]['Area']=numpy.nan #ADVISED TO CHANGE TO ZERO VALUES IN SKYLINE BY M. CHOI
    #combined_results[combined_results['annotation_QValue']>options.mprophet_q]['Area']=0
    combined_results.loc[combined_results['annotation_QValue']>options.mprophet_q,'Area']=0
    #combined_results=combined_results.copy(deep=True)  # BAD PRACTICE...
    #combined_results.loc[combined_results['annotation_QValue']>options.mprophet_q,'Area']=numpy.nan


if options.remove_truncated_peaks:
    combined_results['Area']



if options.merge_isotopes:
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    combined_results['unique_name']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']+"_"+combined_results['File Name']+"_"+combined_results['Protein Name']
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(int)
    groups=combined_results.groupby(by=['unique_name'],as_index=False).agg({'Area':numpy.sum})
    #ADD MASK HERE TO SUM PRECURSOR AREAS IN A FUTURE UPDATE TO ALLOW FOR DIA DATA
    combined_results.drop('Area',1,inplace=True)
    merged_results=pandas.merge(combined_results,groups,on=['unique_name'])
    merged_results.drop_duplicates(subset='unique_name',inplace=True)
    merged_results['Fragment Ion']="sum"
    combined_results=merged_results
    column_list=combined_results.columns.tolist()
    column_list.append(column_list.pop(column_list.index('Standard Type')))
    column_list.append(column_list.pop(column_list.index('Truncated')))
    combined_results.reindex(columns=column_list)
    combined_results.drop('unique_name',1,inplace=True)
    #combined_results.to_csv("TEMP_TEST.csv",sep=",",index=False)


if options.peptide_level:
    fixed_names=combined_results['Peptide Modified Sequence'].apply(peptide_level_fixer)
    #fixed_names=[row[0] for row in fixed_names_list]
    combined_results['Protein Name']=fixed_names
    combined_results=combined_results.copy(deep=True)
    combined_results.sort_values(by='Protein Name',inplace=True)
    

if options.remove_repeated_peptides:
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    combined_results['unique_name']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']+"_"+combined_results['Fragment Ion']+"_"+combined_results['File Name']
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(int)
    print "Removing duplicate peptides...",len(combined_results)
    combined_results.drop_duplicates(subset='unique_name',keep=False,inplace=True)
    print "Done!",len(combined_results)
    combined_results.drop('unique_name',1,inplace=True)

    
if options.remove_empty:
    if options.merge_isotopes:
        #mask = combined_results['Area'].apply(lambda x: numpy.isnan(x))
        #combined_results=combined_results[mask]
        combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
        combined_results['unique_name']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']+"_"+combined_results['Protein Name']
        combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(int)
        groups=combined_results.groupby(by=['unique_name'],as_index=False).agg({'Area':numpy.nansum}) # CHANGED TO NANSUM
        groups.rename(columns={'Area':'TransitionAreaSum'},inplace=True)
        merged_results=pandas.merge(combined_results,groups,on=['unique_name'])
        #crap=merged_results[numpy.invert(merged_results['TransitionAreaSum']>0)]
        #print list(merged_results.columns.values)
        #crap=combined_results[numpy.invert(combined_results['Area_y']>0)]
        #print "FILTERING CRAP:",crap
        #pause=raw_input("paused....")
        combined_results=merged_results[merged_results['TransitionAreaSum']>0]
        combined_results.drop('unique_name',1,inplace=True)
        combined_results.drop('TransitionAreaSum',1,inplace=True)   # DROP TRANSITION AREA SUM, TOO.... BUT LEFT FOR NOW FOR TESTING
        ##transition_groups=combined_results.groupby(by='unique_name',as_index=True)['Area'].sum()
        ##transition_groups=transition_groups.apply(lambda x: numpy.isnan(x) or x<=0.0)
        ##transition_groups=transition_groups[transition_groups==True]
        ##index_list=list(transition_groups.index.values)
        ##combined_results=combined_results[numpy.invert(combined_results.unique_name.str.contains("|".join(index_list)))]
        
        
        #print transition_groups,type(transition_groups)
        #print combined_results,"BEFORE"
        #pause=raw_input("PAUSE...")
        
        #print combined_results,"AFTER"
        #sys.exit(0)
        #empty_transitions=[]
        #empty_transitions=transition_groups.
        #for index,each in transition_groups.iteritems():
        #    if numpy.isnan(each) or each<=0.0:
        #        empty_transitions.append(str(index))

        #for eachtransition in empty_transitions:
        #    combined_results=combined_results[numpy.invert(combined_results.unique_name.str.contains(eachtransition))]
        


    protein_groups=combined_results.groupby(by='Protein Name',as_index=False).agg({'Area':numpy.nansum})#['Area'].sum()
    protein_groups.rename(columns={'Area':'ProteinAreaSum'},inplace=True)
    merged_results=pandas.merge(combined_results,protein_groups,on=['Protein Name'])
    combined_results=merged_results[merged_results['ProteinAreaSum']>0]
    combined_results.drop('ProteinAreaSum',1,inplace=True)   # DROP TRANSITION AREA SUM, TOO.... BUT LEFT FOR NOW FOR TESTING
    
    combined_results['Area'].replace(0,numpy.nan)
    del protein_groups
    
    
    #empty_proteins=[]
    #for index,each in protein_groups.iteritems():
    #    #print index
    #    #if "AIPT[+80]VNHSGTFSPQAPVPTTVPVVDVR" in str(index):
    ##        #print each,"THIS IS IT FOR AIPT[+80]VNHSGTFSPQAPVPTTVPVVDVR"
    # #       #sys.exit(0)
    #    if numpy.isnan(each) or each<=0.0:
    #        if "AIPT[+80]VNHSGTFSPQAPVPTTVPVVDVR" in str(index):
    #            print each,"THIS IS IT FOR AIPT[+80]VNHSGTFSPQAPVPTTVPVVDVR"
    #        empty_proteins.append(str(index))
    #combined_results.rename(columns={'Protein Name':'ProteinName'},inplace=True)#
    #
    #for eachprotein in empty_proteins:
    #    combined_results=combined_results[numpy.invert(combined_results.ProteinName.str.contains(eachprotein))]
    #   combined_results.rename(columns={'ProteinName':'Protein Name'},inplace=True)

if options.zeros:
    print "We're going to replace all zero values with np.nan!"
    combined_results['Area'].replace(0,numpy.nan)

elif options.zero_to_one:
    print "We're going to replace all zero values with 1!\n This is usually a BAD idea... Like, almost ALWAYS."
    combined_results['Area'].replace(0,1)
    if options.zeros:
        print "WARNING: CANNOT REPLACE ZEROS WITH ONE-- AND NP NAN.  REPLACING WITH NAN TAKES PRECEDENT."

#else:
#    print "CAUTION: THERE MAY BE ZEROS IN DATASET WHICH MAY INTERFERE WITH STATISTICAL ANALYSIS."

if options.remove_single_run:
    feature_groups=combined_results.groupby(by=['Protein Name','Peptide Modified Sequence'],as_index=False).count()
    #if options.merge_isotopes:
    #    feature_groups=feature_groups[feature_groups['Area']>1] #>1 means at least 2 runs have data (summed isotopes)
    #else:
    #    feature_groups=feature_groups[feature_groups['Area']>3] #>3 means at least 2 runs have data with 3 isotopes
    feature_groups['unique_name']=feature_groups['Protein Name']+feature_groups['Peptide Modified Sequence']
    feature_groups=feature_groups[['Area','unique_name']]
    feature_groups.rename(columns={'Area':'AreaCount'},inplace=True)
    #feature_groups['keep']=1
    
    combined_results['unique_name']=combined_results['Protein Name']+combined_results['Peptide Modified Sequence']
    
    #merged_results=pandas.merge(combined_results,groups,on=['unique_name'])
    combined_results=pandas.merge(combined_results,feature_groups,on=['unique_name'])#,indicator=True)#,how="right")
    if options.merge_isotopes:
        combined_results=combined_results[combined_results['AreaCount']>1]
    else:
        combined_results=combined_results[combined_results['AreaCount']>3]
    
    
    combined_results.drop('unique_name',1,inplace=True)
    combined_results.drop('AreaCount',1,inplace=True)
    #print combined_results


    
#Okay, now we'll go ahead and print out the final MSstats input file!
if options.rename:
    combined_results.rename(columns={'Protein Name':'ProteinName','Peptide Modified Sequence':'PeptideSequence','Precursor Charge':'PrecursorCharge','Fragment Ion':'FragmentIon','Product Charge':'ProductCharge','Isotope Label Type':'IsotopeLabelType','Standard Type':'StandardType','File Name':'Run','Area':'Intensity'},inplace=True)


if options.mprophet_q is not None:
    combined_results = combined_results.drop('annotation_QValue', 1)
    combined_results = combined_results.drop('annotation_Score',1)
	
#if options.remove_decoys:
#    if not options.rename:
#        combined_results.rename(columns={'Protein Name':'ProteinName'},inplace=True)
#    mask=combined_results[numpy.invert(combined_results.ProteinName.str.contains("Decoys"))]
#    if not options.rename:
#        mask.rename(columns={'ProteinName':'Protein Name'},inplace=True)
#    mask.to_csv(str(os.path.join(basedir,options.operation_folder))+"/MSstats_combined_input.csv",sep=",",index=False)
#else:

if options.fillMissingFeatures:
    bioreplicate_dict={}
    condition_dict={}

    for each_run in combined_results['Run'].unique():
        temp=combined_results[combined_results['Run']==each_run]
        bioreplicate_dict[each_run]=temp['BioReplicate'].unique()[0]
        condition_dict[each_run]=temp['Condition'].unique()[0]

    grouped_df=combined_results.groupby(["PeptideSequence","PrecursorCharge","FragmentIon","ProductCharge"])
    concat_list=[]
    correct_length=len(bioreplicate_dict.keys())
    for name,eachgroup in grouped_df:
        if len(eachgroup)!=correct_length:
            for each_name in bioreplicate_dict.keys():#name_list:
                if each_name not in eachgroup['Run'].unique():
                    new_row=eachgroup.head(n=1).copy(deep=True)
                    new_row['Run']=each_name
                    new_row['Intensity']=numpy.nan
                    new_row['Condition']=condition_dict[each_name]
                    new_row['BioReplicate']=bioreplicate_dict[each_name]
                    concat_list.append(new_row)
    concat_list.append(combined_results)
    combined_results=pandas.concat(concat_list)

    #combined_results=pandas.concat([combined_results,new_row])


combined_results=combined_results.drop('Mass Error PPM', 1)
combined_results=combined_results.drop('Average Mass Error PPM', 1)

combined_results.to_csv(str(os.path.join(basedir,options.operation_folder))+"/MSstats_combined_input.csv",sep=",",index=False)

print "We have now written the combined MSstats input file to "+str(os.path.join(basedir,options.operation_folder))+"\\MSstats_combined_input.csv"
#print "Now finished.  Quitting!"

print "Now we're going to prepare the R script for MSstats"

#Let's start by reading in the experiment structure.
group_information = pandas.read_csv(options.experiment_file,sep='\t')


with open("MSstats_Script.R",'wb') as script_writer:
    script_writer.write("library(MSstats)\n")
    script_writer.write("library(preprocessCore)\n")
    if os.name=="nt":
        script_writer.write(str("setwd(\""+str(os.path.join(basedir,options.operation_folder))+"\")\n").replace("\\","\\\\"))   #We're going to set the current directory...
        script_writer.write(str("raw<-read.csv(\""+str(os.path.join(basedir,options.operation_folder))+"/MSstats_combined_input.csv"+"\")\n").replace("\\","\\\\"))   #We will load in the input CSV file! (In this case by absolute path, though that's not necessary...)
    else:
        script_writer.write("setwd(\""+str(os.path.join(basedir,options.operation_folder))+"\")\n")   #We're going to set the current directory...
        script_writer.write("raw<-read.csv(\""+str(os.path.join(basedir,options.operation_folder))+"/MSstats_combined_input.csv"+"\")\n")   #We will load in the input CSV file! (In this case by absolute path, though that's not necessary...)

    ##### MSstats dataProcess #####
    script_writer.write("TMP_result<-dataProcess(raw=raw,")
    if options.fillIncompleteRows:
        script_writer.write("fillIncompleteRows=TRUE,")
    else:
        script_writer.write("fillIncompleteRows=FALSE,")

    if "NULL" in options.maxQuantileforCensored:
        script_writer.write("maxQuantileforCensored=NULL,")
        print "MAX QUANTILE FOR CENSORED HAS BEEN TURNED OFF."
    else:
        #print "MAX QUANTILE FOR CENSORED HAS BEEN TURNED OFF."
        #try:
            #float(options.maxQuantileforCensored)
            #print "MAX QUANTILE FOR CENSORED HAS BEEN TURNED OFF."

        script_writer.write("maxQuantileforCensored={0},".format(options.maxQuantileforCensored))
        #except:
        #    #print "MAX QUANTILE FOR CENSORED HAS BEEN TURNED OFF."
        #
        #    #script_writer.write("maxQuantileforCensored=NULL,")


    if "quantile" in options.normalization:
        script_writer.write("normalization=\"quantile\",")
    elif "equalizeMedians" in options.normalization:
        script_writer.write("normalization=\"equalizeMedians\",")
    elif "globalStandards" in options.normalization:
        script_writer.write("normalization=\"globalStandards\",")
        norm_peptide_set=set(combined_results[combined_results["StandardType"]=="globalStandard"]["PeptideSequence"])
        peptide_norm_string="c(\""+"\",\"".join(norm_peptide_set)+"\")"
        script_writer.write("nameStandards="+peptide_norm_string+",")
    elif "false" in options.normalization:
        script_writer.write("normalization=FALSE,")

    #script_writer.write("betweenRunInterferenceScore=TRUE,") #because why not have them?  # because it crashes the analysis on certain large projects...
    

    if "all" in options.feature_subset:
        script_writer.write("featureSubset=\"all\",")
    elif "topN" in options.feature_subset:
        if options.feature_subset_N==3:
            script_writer.write("featureSubset=\"top3\",")
        else:
            script_writer.write("featureSubset=\"topN\",")
            script_writer.write("featureSubset=\"n_top_feature="+str(options.feature_subset_N)+",")
    elif "highQuality" in options.feature_subset:
        script_writer.write("featureSubset=\"highQuality\",")

    

    script_writer.write("summaryMethod=\"TMP\",")

    if options.noimpute:
        script_writer.write("censoredInt=NULL,MBimpute=FALSE,")
    else:
        if int(options.censoredInt)==0:
            script_writer.write("censoredInt=0,")
        elif "NULL" in options.censoredInt:
            script_writer.write("censoredInt=\"NA\",")


        if options.MBimpute:
            script_writer.write("MBimpute=TRUE,")
        else:
            script_writer.write("MBimpute=FALSE,")
            

        if "minFeature" in options.cutoffCensored:
            script_writer.write("cutoffCensored=\"minFeature\",")
        elif "minRun" in options.cutoffCensored:
            script_writer.write("cutoffCensored=\"minRun\",")
        elif "minFeatureNRun" in options.cutoffCensored:
            script_writer.write("cutoffCensored=\"minFeatureNRun\",")

        if options.remove50missing:
            script_writer.write("remove50missing=TRUE,")
        else:
            script_writer.write("remove50missing=FALSE,")

    script_writer.write("logTrans=2)\n")
    #script_writer.write("logTrans=2,")
    #script_writer.write("skylineReport=TRUE)\n")
    ######################################################
    #Data is now stored in "TMP_result$ProcessedData"
    #TMP_result$RunlevelData
    #TMP_result$SummaryMethod
    #TMP_result$PredictBySurvival (if MBimpute with AFT)
    #In this section, we'll write out the QC type plots... and then make our comparison based on the experimental structure info... (which is stored in "group_information")

    script_writer.write("write.csv(TMP_result$ProcessedData,file=\"TMP_dataProcess_output.csv\")\n")
    script_writer.write("quant_matrix<-quantification(TMP_result,type=\"Sample\",format=\"matrix\")\n")
    script_writer.write("write.csv(quant_matrix,file=\"quantification.csv\")\n")
    script_writer.write("cond_quant_matrix<-quantification(TMP_result,type=\"Group\",format=\"matrix\")\n")
    script_writer.write("write.csv(cond_quant_matrix,file=\"condition_quantification.csv\")\n")

    script_writer.write("dataProcessPlots(data = TMP_result, type=\"QCplot\", ylimUp=35)\n")

    script_writer.write("dataProcessPlots(data = TMP_result, type=\"ProfilePlot\", ylimUp=35, featureName=\"NA\",width=7,height=7)\n")

    script_writer.write("dataProcessPlots(data = TMP_result, type=\"ConditionPlot\",save_condition_plot_result=TRUE)\n")

    #dataProcessPlots(data=ideal,type="QCPlot")
    #dataProcessPlots(data=ideal,type="ProfilePlot")
    #dataProcessPlots(data=ideal,type="ConditionPlot")

    ######################################################
    ### After we have those plots, we're going to ask for the comparisons themselves... (again, "group_information")
    
    #### Note that the conditions are arranged in alphabetical order <--------------- !!!!!!
    #
    #
    #Let's go ahead and make strings to describe our comparison matrix for MSstats!  Let's grab all of the biological conditions from group_information dataframe!

    #groups_abet=sorted(group_information['Biological Condition'].unique().tolist())  #This should be the same order that MSstats sees! -- it appears that it is not...
    #groups_abet=natsorted(group_information['Biological Condition'].unique().tolist(),alg=ns.IGNORECASE | ns.REAL)  #This should be the same order that MSstats sees! -- it appears that it is not... alphabetical ignore case
    #Natsort fails when numbers come into play, so we'll just put this to bed and use R to sort.
    with open("R_sorter.Rscript",'wb') as sorter_script_writer:
        sorter_script_writer.write("sort(c(\""+"\",\"".join(group_information['Biological Condition'].unique().tolist())+"\"))")
    groups_abet_cmdout=subprocess.check_output(["Rscript","R_sorter.Rscript"])
    print groups_abet_cmdout,"\nACTUALOUTPUT==================="
    groups_abet=[x.strip("\"").rstrip("\"") for x in groups_abet_cmdout.split()[1:]]
    new_abet=[]
    for each_group in groups_abet:
        if each_group.startswith('[') and each_group.endswith(']'):  #These groups are actually just line information from the R output.  We'll remove them.
            print "Group ",each_group," is not a real group... dropping..."
        else:
            new_abet.append(each_group)
    groups_abet=new_abet
    print groups_abet,"This should be alphabetical!"
    num_groups=len(groups_abet)

    group_to_MSstats_order={}
    MSstats_order_to_group={}
    control_index_list=[]
    test_index_list=[]
    for i in xrange(0,num_groups):
        if "C" in group_information[group_information['Biological Condition'] == groups_abet[i]]['Test or Control'].unique().tolist():
            control_index_list.append(i)
        else:
            print "was looking for C in...",group_information[group_information['Biological Condition'] == groups_abet[i]]['Test or Control'].unique().tolist()
            test_index_list.append(i)
        group_to_MSstats_order[groups_abet[i]]=i
        MSstats_order_to_group[i]=groups_abet[i]
    comparisons=[]
    comparison_names=[]
    control_fraction=1.0/float(len(control_index_list))
    print "MSstats o to g",MSstats_order_to_group

    for i in xrange(0,num_groups):  #A line for each comparison
        this_comparison_name=""
        if i in control_index_list:
            continue
        #Otherwise...
        new_comparison=[]
        for j in xrange(0,num_groups):   #A position in the line for each comparison
            if i==j:
                new_comparison.append(str(1))
                #this_comparison_name=MSstats_order_to_group[i]+" vs."+this_comparison_name
                this_comparison_name=MSstats_order_to_group[i]+" "#+" vs."+this_comparison_name
            elif j in control_index_list:
                new_comparison.append(str(-1.0*control_fraction))
                #this_comparison_name=this_comparison_name+" "+MSstats_order_to_group[j]+"&"
            else:
                new_comparison.append(str(0))
        this_comparison_name=this_comparison_name[:-1]
        comparison_names.append(this_comparison_name)
        comparisons.append(new_comparison)

    ######################### OKAY #######################
    ########## Now we have the comparisons built and we're
    ########## going to take those comparisons and build
    ########## the strings that MSstats needs to make its
    ########## analysis go!
    ########## comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
    ########## comparison2<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
    ########## comparison3<-matrix(c(-1,0,0,0,0,0,0,0,1,0),nrow=1)
    ########## comparison<-rbind(comparison1,comparison2, comparison3)
    ########## row.names(comparison)<-c("T3-T1","T7-T1","T9-T1")

    i=1
    temp_rbind_string=""
    for each_comparison_list in comparisons:
        script_writer.write("comparison"+str(i)+"<-matrix(c("+",".join(each_comparison_list)+"),nrow=1)\n")
        temp_rbind_string+="comparison"+str(i)+","
        i+=1
    temp_rbind_string=temp_rbind_string[:-1]
    #if len(temp_rbind_string.split(","))>1:
    script_writer.write("comparison<-rbind("+temp_rbind_string+")\n")
    #else:
    #    script_writer.write("comparison<-comparison1\n")
    temp_comparison_string=""
    for each_comp in comparison_names:
        temp_comparison_string+="\""+each_comp+"\""+","
    temp_comparison_string=temp_comparison_string[:-1]
    script_writer.write("row.names(comparison)<-c("+temp_comparison_string+")\n")
    # OKAY, so at this point we've now constructed our "comparison" object in the R script...
    # Now we will have to prepare the groupComparison!
    
    script_writer.write("comparisonResult<-groupComparison(contrast.matrix=comparison,data=TMP_result)\n")  #May need TMP_result$ProcessedData
    script_writer.write("print(levels(TMP_result$ProcessedData$GROUP_ORIGINAL))\n")
    # The heavy lifting is done!!! Let's write out a csv before we go crazy with plotting...

    script_writer.write("write.csv(comparisonResult$ComparisonResult,file=\"comparisonResult_output.csv\")\n")
    #script_writer.write("save(comparisonResult,file=\"image.RData\")\n")


    script_writer.write("quantification_csv<-read.csv(\"quantification.csv\",check.names=FALSE)\n")
    script_writer.write("comparison_csv<-read.csv(\"comparisonResult_output.csv\",check.names=FALSE)\n")
    script_writer.write("condition_plot_csv<-read.csv(\"ConditionPlot_value.csv\",check.names=FALSE)\n")
    #script_writer.write("experiment_design<-read.csv(\""+options.experiment_file+"\",sep=\"\\t\",check.names=FALSE)\n")


    script_writer.write("save.image(file=\"image.RData\")\n")
    #OKAY! So, now we're going to write out the plots... This may take a bit...
    #So, first, let's check if we can output a heatmap (number of comparisons >2)
    #if len(temp_comparison_string.split(","))>2:
    #    script_writer.write("groupComparisonPlots(data=comparisonResult$ComparisonResult,type="Heatmap", logBase.pvalue=2, sig="+str(options.significance)+", FCcutoff="+options.FCthreshold


    pass


#OKAY.... The R Script has been written!
#We're going to execute the R script now!
print "Copying RScript back to Galaxy..."
shutil.copy('MSstats_Script.R',options.RScriptOutput)
#subprocess.check_call(['Rscript', 'MSstats_Script.R'],shell=False,stderr=sys.stdout.fileno())



#A word to the wise... Don't run MSstats through Rscript... At least as of 2017/05/02, it causes some serious issues...  Much better to run straight into R itself, for whatever reason...

if os.name=="nt":
    #subprocess.check_call(["Rscript", 'MSstats_Script.R'],shell=True,stderr=sys.stdout.fileno())
    ##subprocess.check_call(["R",'<', 'MSstats_Script.R','--vanilla'],shell=True,stderr=sys.stdout.fileno())
    subprocess.check_call(["R < MSstats_Script.R --vanilla"],shell=True,stderr=sys.stdout.fileno())
else:
    #subprocess.check_call(['Rscript', 'MSstats_Script.R'],shell=True,stderr=sys.stdout.fileno())
    ##subprocess.check_call(['R','<', 'MSstats_Script.R','--vanilla'],shell=True,stderr=sys.stdout.fileno())
    subprocess.check_call(['R < MSstats_Script.R --vanilla'],shell=True,stderr=sys.stdout.fileno())

print "Moving files to final output locations...."
shutil.copy('TMP_dataProcess_output.csv',options.processedOutput)
shutil.copy('quantification.csv',options.quantificationOutput)
shutil.copy('condition_quantification.csv',options.quantificationConditionOutput)
shutil.copy('comparisonResult_output.csv',options.comparisonOutput)
shutil.copy('ConditionPlot_value.csv',options.conditionPlotCSVOutput)
shutil.copy('image.RData',options.RDataOutput)
shutil.copy('QCPlot.pdf',options.QCplotOutput)
shutil.copy('ProfilePlot.pdf',options.profilePlotOutput)
shutil.copy('ConditionPlot.pdf',options.conditionPlotOutput)


print "All done!"
