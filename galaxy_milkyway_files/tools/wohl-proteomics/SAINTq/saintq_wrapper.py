import os, sys, re
import optparse
import shutil
import pandas
import numpy
import gc
import subprocess
import uniprot as uni
from natsort import natsorted, ns

#####################################
#This is a script to combine the output reports from
#Skyline, in preparation for MSstats!  Let's get started.
#
#VERSION 0.1
version="0.1"
#DATE: 5/18/2017
date="5/18/2017"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the SAINTq wrapper for Galaxy, Wohlschlegel Lab UCLA"
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
parser.add_option("--mprophet_q",action="store",type="float",dest="mprophet_q") #We'll throw out things above this q-value threshold.
parser.add_option("--input_file",action="store",type="string",dest="input_file")
parser.add_option("--output_file",action="store",type="string",dest="output_file")
parser.add_option("--simple_output_file",action="store",type="string",dest="simple_output_file")
parser.add_option("--exp_file",action="store",type="string",dest="exp_file")
parser.add_option("--quant_level",action="store",type="string",dest="quant_level")#"peptide","protein",or "fragment"
parser.add_option("--best_prop_pep",action="store",type="float",dest="best_prop_pep")
parser.add_option("--min_n_pep",action="store",type="int",dest="min_n_pep")#Minimum number of peptides
parser.add_option("--best_prop_frag",action="store",type="float",dest="best_prop_frag")
parser.add_option("--min_n_frag",action="store",type="int",dest="min_n_frag")#Minimum number of fragments...
parser.add_option("--compress_n_ctrl",action="store",type="int",dest="compress_n_ctrl")
parser.add_option("--compress_n_rep",action="store",type="int",dest="compress_n_rep")
parser.add_option("--normalize_control",action="store",type="string",dest="normalize_control")#"true" or "false"
parser.add_option("--remove_repeated_peptides",action="store_true",dest="remove_repeated_peptides")
parser.add_option("--fill_missing_features",action="store_true",dest="fillMissingFeatures")

(options,args) = parser.parse_args()


group_information = pandas.read_csv(options.exp_file,sep='\t')


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


if options.quant_level!="protein":    
    input_data=pandas.read_csv(options.input_file,sep=',',index_col=False,low_memory=False)

    input_data.rename(columns={'Protein Name':'ProteinName'},inplace=True)
    combined_results=input_data[numpy.invert(input_data.ProteinName.str.contains("Decoys"))]
    combined_results.rename(columns={'ProteinName':'Protein Name'},inplace=True)

    combined_results.sort_values(by='Protein Name',inplace=True)

if options.mprophet_q is not None:
    combined_results.loc[combined_results['annotation_QValue']>options.mprophet_q,'Area']=0.0

if options.quant_level=="peptide":
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    combined_results['unique_name']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']+"_"+combined_results['File Name']+"_"+combined_results['Protein Name']
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(int)
    groups=combined_results.groupby(by=['unique_name'],as_index=False).agg({'Area':numpy.sum})
    combined_results.drop('Area',1,inplace=True)
    merged_results=pandas.merge(combined_results,groups,on=['unique_name'])
    merged_results.drop_duplicates(subset='unique_name',inplace=True)#change cols to subset for newer pandas...
    merged_results['Fragment Ion']="sum"
    combined_results=merged_results
    column_list=combined_results.columns.tolist()
    column_list.append(column_list.pop(column_list.index('Standard Type')))
    column_list.append(column_list.pop(column_list.index('Truncated')))
    combined_results.reindex(columns=column_list)
    combined_results.drop('unique_name',1,inplace=True)
    combined_results.fillna(value=0.0,inplace=True)

    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    
if options.remove_repeated_peptides:
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    combined_results['unique_name']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']+"_"+combined_results['Fragment Ion']+"_"+combined_results['File Name']
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(int)
    print "Removing duplicate peptides...",len(combined_results)
    combined_results.drop_duplicates(subset='unique_name',keep=False,inplace=True)
    print "Done!",len(combined_results)
    combined_results.drop('unique_name',1,inplace=True)

    
#if options.rename:
#    combined_results.rename(columns={'Protein Name':'ProteinName','Peptide Modified Sequence':'PeptideSequence','Precursor Charge':'PrecursorCharge','Fragment Ion':'FragmentIon','Product Charge':'ProductCharge','Isotope Label Type':'IsotopeLabelType','Standard Type':'StandardType','File Name':'Run','Area':'Intensity'},inplace=True)

if options.fillMissingFeatures:
    bioreplicate_dict={}
    condition_dict={}

    for each_run in combined_results['File Name'].unique():
        temp=combined_results[combined_results['File Name']==each_run]
        bioreplicate_dict[each_run]=temp['BioReplicate'].unique()[0]
        condition_dict[each_run]=temp['Condition'].unique()[0]

    grouped_df=combined_results.groupby(["Peptide Modified Sequence","Precursor Charge","Fragment Ion","Product Charge"])
    concat_list=[]
    correct_length=len(bioreplicate_dict.keys())
    for name,eachgroup in grouped_df:
        if len(eachgroup)!=correct_length:
            for each_name in bioreplicate_dict.keys():#name_list:
                if each_name not in eachgroup['File Name'].unique():
                    new_row=eachgroup.head(n=1).copy(deep=True)
                    new_row['File Name']=each_name
                    new_row['Intensity']=numpy.nan
                    new_row['Condition']=condition_dict[each_name]
                    new_row['BioReplicate']=bioreplicate_dict[each_name]
                    concat_list.append(new_row)
    concat_list.append(combined_results)
    combined_results=pandas.concat(concat_list)

    #combined_results=pandas.concat([combined_results,new_row])

if options.quant_level=="peptide":
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    combined_results['Product Charge']=combined_results['Product Charge'].astype(str)
    combined_results['unique_nofile_name']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']+"_"+combined_results['Protein Name']+"_"+combined_results["Fragment Ion"]+"_"+combined_results["Product Charge"]
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(int)
    combined_results['Product Charge']=combined_results['Product Charge'].astype(int)
    combined_results['file']=combined_results['Condition']+"_"+combined_results['File Name']
    combined_results.fillna(value=0.0,inplace=True)
    pivot_df=combined_results[['unique_nofile_name','file','Area']]
    pivot_df=pivot_df.pivot(index='unique_nofile_name',columns='file',values='Area')

    combined_results=combined_results[['Protein Name','Peptide Modified Sequence','Precursor Charge','unique_nofile_name']]
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    combined_results['peptide']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']
    combined_results=combined_results.drop_duplicates(subset=['unique_nofile_name'])
    combined_results=combined_results.merge(right=pivot_df, left_on='unique_nofile_name',right_index=True)
    combined_results=combined_results.drop('unique_nofile_name', 1)
    combined_results=combined_results.drop('Precursor Charge', 1)
    combined_results=combined_results.drop('Peptide Modified Sequence', 1)
    combined_results.sort_values(by=['peptide'],inplace=True)
    combined_results.rename(columns={'Protein Name':'Protein'},inplace=True)


elif options.quant_level=="fragment":
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    combined_results['Product Charge']=combined_results['Product Charge'].astype(str)
    combined_results['unique_nofile_name']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']+"_"+combined_results['Protein Name']+"_"+combined_results["Fragment Ion"]+"_"+combined_results["Product Charge"]
    combined_results['unique_name']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']+"_"+combined_results['Protein Name']+"_"+combined_results["Fragment Ion"]+"_"+combined_results["Product Charge"]+"_"+combined_results['File Name']
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(int)
    combined_results['Product Charge']=combined_results['Product Charge'].astype(int)
    combined_results['file']=combined_results['Condition']+"_"+combined_results['File Name']
    combined_results.fillna(value=0.0,inplace=True)
    pivot_df=combined_results[['unique_nofile_name','file','Area','unique_name']]
    pivot_df.drop_duplicates(subset='unique_name',inplace=True)
    pivot_df.drop('unique_name',1,inplace=True)
    pivot_df=pivot_df.pivot(index='unique_nofile_name',columns='file',values='Area')

    combined_results=combined_results[['Protein Name','Peptide Modified Sequence','Precursor Charge','Fragment Ion','Product Charge','unique_nofile_name']]
    combined_results['Precursor Charge']=combined_results['Precursor Charge'].astype(str)
    combined_results['Product Charge']=combined_results['Product Charge'].astype(str)
    combined_results['peptide']=combined_results['Peptide Modified Sequence']+"_"+combined_results['Precursor Charge']
    combined_results['fragment']=combined_results['Fragment Ion']+"_"+combined_results['Product Charge']
    combined_results=combined_results.drop_duplicates(subset=['unique_nofile_name'])
    combined_results=combined_results.merge(right=pivot_df, left_on='unique_nofile_name',right_index=True)
    combined_results=combined_results.drop('unique_nofile_name', 1)
    combined_results=combined_results.drop('Precursor Charge', 1)
    combined_results=combined_results.drop('Peptide Modified Sequence', 1)
    combined_results=combined_results.drop('Product Charge', 1)
    combined_results=combined_results.drop('Fragment Ion', 1)
    combined_results.sort_values(by=['peptide'],inplace=True)
    combined_results.rename(columns={'Protein Name':'Protein'},inplace=True)

elif options.quant_level=="protein":
    combined_results=pandas.read_csv(options.input_file,sep=',',index_col=0,low_memory=False)
    #combined_results=combined_results
    combined_results.update(combined_results[combined_results.select_dtypes(include=['number']).columns].applymap(numpy.exp2),overwrite=True)
    combined_results.fillna(value=0.0,inplace=True)
    #print combined_results
    #pass

combined_results.to_csv("SAINTq_combined_input.tsv",sep="\t",index=False)

print "We have now written the combined SAINTq intensity input file to SAINTq_combined_input.csv"

#Let's start by reading in the experiment structure.
#We're going to have to add the Control/Test headers for SAINTq.
header=True
with open("SAINTq_combined_input.tsv",'rb') as tab_reader:
    with open("final_saint_input_file.tsv",'wb') as final_writer:
        for each_row in tab_reader:
            if header==True:
                header=False
                columns=each_row.split("\t")
                
                if options.quant_level=="protein":
                    starting_index=1 #We skip the Protein column
                elif options.quant_level=="peptide":
                    starting_index=2 #Skip protein and Peptide columns
                elif options.quant_level=="fragment":
                    starting_index=3 #Skip protein, peptide and fragment columns

                #print "all columns...",columns
                #print "starting index is...",starting_index
                run_samples=columns[starting_index:]#We skip the Protein column
                run_T_or_C_vector=[]  #populate with tabs, C, or T's
                run_condition_vector=[]   #populate with tabs, and conditions


                for x in xrange(0,starting_index): #prepare the number of tabs needed...
                    run_T_or_C_vector.append("")
                    run_condition_vector.append("")

                for each_sample in run_samples:
                    if options.quant_level=="protein":
                        #print each_sample
                        run_condition=each_sample.rsplit("_",1)[0]
                        #print run_condition,"this is run cond"
                        #print group_information[group_information['Biological Condition']==run_condition]
                        if "T" in group_information[group_information['Biological Condition']==run_condition]['Test or Control'].tolist()[0]:
                            run_T_or_C_vector.append("T")
                        else:
                            run_T_or_C_vector.append("C")
                    else:
                        run_condition= group_information[group_information['Original File Name']==each_sample.rsplit(".",1)[0].split("_",1)[1]]['Biological Condition'].tolist()[0]
                        if "T" in group_information[group_information['Biological Condition']==run_condition]['Test or Control'].unique().tolist():
                            run_T_or_C_vector.append("T")
                        else:
                            run_T_or_C_vector.append("C")

                    run_condition_vector.append(run_condition)
                run_condition_vector[-1]+="\n"
                run_T_or_C_vector[-1]+="\n"
                final_writer.write("\t".join(run_T_or_C_vector))
                final_writer.write("\t".join(run_condition_vector))
            final_writer.write(each_row)
                    


with open("saintq_configuration",'wb') as config_writer:
    config_writer.write("normalize_control="+options.normalize_control+"\n")
    config_writer.write("input_filename=final_saint_input_file.tsv\n")
    config_writer.write("input_level="+options.quant_level+"\n")

    config_writer.write("protein_colname=Protein\n")
    if options.quant_level=="peptide" or options.quant_level=="fragment":
        config_writer.write("pep_colname=peptide\n")
    if options.quant_level=="fragment":
        config_writer.write("frag_colname=fragment\n")
    
    config_writer.write("compress_n_ctrl={0}\n".format(options.compress_n_ctrl))
    config_writer.write("compress_n_rep={0}\n".format(options.compress_n_rep))

    if options.quant_level=="peptide" or options.quant_level=="fragment":
        config_writer.write("min_n_pep={0}\n".format(options.min_n_pep))
        config_writer.write("best_prop_pep={0}\n".format(options.best_prop_pep))
    if options.quant_level=="fragment":
        config_writer.write("min_n_frag={0}\n".format(options.min_n_frag))
        config_writer.write("best_prop_frag={0}\n".format(options.best_prop_frag))


os.system("saintq saintq_configuration")

for each_file in os.listdir(os.getcwd()):
    if "scores_list_" in each_file:
        break
#copy raw output...
shutil.copy(each_file,options.output_file)

raw_df=pandas.read_csv(options.output_file,sep='\t')

pivot_output=raw_df.pivot(index='Prey',columns='Bait',values='AvgP')
pivot_output.fillna(value=0.0,inplace=True)
pivot_output.to_csv(options.simple_output_file,sep="\t")
#print pivot_output




#shutil.copy(simple_output,options.simple_output_file)

print "All done!"
