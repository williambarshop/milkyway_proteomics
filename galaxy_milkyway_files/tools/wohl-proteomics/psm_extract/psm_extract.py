import os, sys, re
import optparse
import shutil
import pandas
import numpy
import gc
import multiprocessing
from joblib import Parallel, delayed


parser = optparse.OptionParser()
#For psm extractor, we're going to filter by 1. PSM q-value and 2. FIDO q-values...
parser.add_option("--qthresh",action="store", type="float", dest="q_threshold")
parser.add_option("--fidoq",action="store", type="float", dest="fido_q_threshold")

(options,args) = parser.parse_args()

#### Check for FIDO q-filter
global gfido_q_threshold
gfido_q_threshold = options.fido_q_threshold
print "We're going to filter by Fido Q-Value of ",options.fido_q_threshold
q_threshold = options.q_threshold
print "We'll also filter at a PSM (percolator) q-value of ",options.q_threshold


experimental_data=pandas.read_csv("experimental_groups.tsv",sep='\t',index_col=False)

for index,eachrow in experimental_data.iterrows():
    eachrow['Original File Name']=eachrow['Original File Name']+".mzML" #Correcting the file names...

global idx_to_run_dict # Bad coding style, but makes for easy 'applies' when dealing with pandas dataframes...
idx_to_run_dict = {}

for index,row in experimental_data.iterrows():
    thisrun = row['Original File Name']
    thisidx = str(row['Crux File Integer'])
    idx_to_run_dict[thisidx]=thisrun

    
def makeRoundedModSequence(x):
    #group.loc[index,'unmodified sequence']=re.sub("[\(\[].*?[\)\]]", "", eachrow['sequence'])
    #mod_mass_list=re.findall(r"[-+]?\d*\.\d+|\d+",x['sequence'])
    mod_mass_list=re.findall(r'[\d\.\d]+',x['sequence'])
    rounded_mass_list=iter([str(int(round(float(y)))) for y in mod_mass_list])

    def replace_func(replace_list):
        def replace_inner(match):
            return next(replace_list)
        return replace_inner

    #print list(rounded_mass_list),"rounded list"
    #print x['sequence'],"sequence"

    #return re.sub(r"[-+]?\d*\.\d+|\d+",replace_func(rounded_mass_list),x['sequence'])
    return re.sub(r'[\d\.\d]+',replace_func(rounded_mass_list),x['sequence'])

    
def addFractionationGroup(x):
    return str(x['File Name'].split('.')[0].rsplit("-",1)[0])#[:-3])  # FIX THIS TO RETURN FRAC GROUP

def makeUniquePinName(x):
    #group.loc[index,'unmodified sequence']=re.sub("[\(\[].*?[\)\]]", "", eachrow['sequence'])
    return str(idx_to_run_dict[x['SpecId'].split("_")[1]]+".mzML."+str(x['ScanNr']))+"."+re.sub("[\(\[].*?[\)\]]", "", x['Peptide']).split(".")[1].strip()
    
def makeUniquePercoName(x):
    return str(x['file']+"."+str(x['scan']))+"."+x['unmodified sequence'].strip()

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
            unfiltered.loc[index,'fido_filter']=1 #This means we'll keep it when we do our df filter... We'll end $
        else:
            unfiltered.loc[index,'fido_filter']=0

        #print changed,type(changed),"this is changed..."
    return unfiltered



#combined_fname_scans=combined_groups.apply(combineFileNameAndScanPerco,axis=1)
#combined_groups['unique_name']=combined_fname_scans



infiles = []
for root, subFolders, files in os.walk(os.getcwd()):
    for eachfile in files:
        if str(eachfile).endswith(".target.psms.txt"):
            infiles.append(str(os.path.join(root,eachfile)))
            
pinfiles = []
for root, subFolders, files in os.walk(os.getcwd()):
    for eachfile in files:
        #print eachfile,"yo",str(eachfile)
        if str(eachfile).endswith(".pin"):
            pinfiles.append(str(os.path.join(root,eachfile)))

            
            


pin_columns=['SpecId','Label','ScanNr','ExpMass','CalcMass','RawScore','DeNovoScore','ScoreRatio','Energy','lnEValue','IsotopeError','lnExplainedIonCurrentRatio','lnNTermIonCurrentRatio','lnCTermIonCurrentRatio','lnMS2IonCurrent','Mass','PepLen','dM','absdM','MeanErrorTop7','sqMeanErrorTop7','StdevErrorTop7','Charge2','Charge3','Charge4','Charge5','Charge6','enzN','enzC','enzInt','ptm','A-Freq','C-Freq','D-Freq','E-Freq','F-Freq','G-Freq','H-Freq','I-Freq','K-Freq','L-Freq','M-Freq','N-Freq','P-Freq','Q-Freq','R-Freq','S-Freq','T-Freq','V-Freq','W-Freq','Y-Freq','B-Freq','Z-Freq','J-Freq','X-Freq','U-Freq','O-Freq','ppm','absppm','Peptide']

for eachfile in pinfiles:
    with open(eachfile,'rb') as header_reader:
        first_file_header=header_reader.next().split("\t")
        #print "these are in the header..."
        #print first_file_header
        break
first_file_header=[x.strip() for x in first_file_header]
#print "before filter",pin_columns
pin_columns = [col for col in pin_columns if col in first_file_header]
#print "after filter",pin_columns

pin_filtered=[]

for eachfile in pinfiles:
    
    #newdf=pandas.DataFrame.from_csv(eachfile,sep='\t',index_col=False,skiprows=[1])
    newdf=pandas.read_csv(eachfile,sep="\t",index_col=False,skiprows=[1],usecols=pin_columns,engine='python')
    #print newdf
    newdf_masked=newdf[newdf['Label'] == 1].copy(deep=True) # We want the targets only.
    pin_unique_names=newdf_masked.apply(makeUniquePinName,axis=1)
    #print pin_unique_names[1]
    newdf_masked['unique_name']=pin_unique_names
    pin_filtered.append(newdf_masked)

combined_pins=pandas.concat(pin_filtered)
#print combined_pins['unique_name']
combined_pins.drop('Peptide',axis=1,inplace=True) #Get rid of it, we already have sequence... jsut needed it to provide a unique key...
#print combined_pins.columns,"Combined Pins COLUMNS...."
#print "------------------------------------<<<<<<<" 
            
dataframe_vector=[]
#print infiles,"these are infiles"
for eachfile in infiles:
    newdf=pandas.read_csv(eachfile,sep='\t',index_col=False)
    perco_unique_names=newdf.apply(makeUniquePercoName,axis=1)
    newdf['unique_name']=perco_unique_names
    #print newdf
    dataframe_vector.append(newdf)
    del newdf
    gc.collect()

combined_results=pandas.concat(dataframe_vector)
#print combined_results.columns,"comb resus (perco) columns"
#combined_results.to_csv("/galaxy-central/tools/wohl-proteomics/psm_extract/perco_df.csv")
#combined_pins.to_csv("/galaxy-central/tools/wohl-proteomics/psm_extract/pins_df.csv")
merged_combined_results=combined_results.merge(right=combined_pins,how="inner",on="unique_name")#,indicator=True)        
merged_combined_results.drop('unique_name',1,inplace=True)
#print merged_combined_results,"after mergex"

####################################################
####################################################
####################################################
#
#This is the section where we have our dataframe with everything important--
#and before filtering.  So from this, we will dump out our QC plots.
#
####################################################
####################################################
####################################################






merged_combined_results=merged_combined_results[merged_combined_results['percolator q-value']<=q_threshold]
#merged_combined_results=merged_combined_results[merged_combined_results['percolator q-value']<=q_threshold]

merged_combined_results['fido_filter']=0
#print merged_combined_results,"pre group"
merged_combined_results['protein q-values']=merged_combined_results['protein q-values'].astype(str)########################################## ADDED TO HANDLE FLOAT SPLIT ISSUE
merged_combined_results_grp=merged_combined_results.groupby(numpy.arange(len(merged_combined_results))//multiprocessing.cpu_count())
merged_combined_results=applyParallelHalf(merged_combined_results_grp,fido_filter)
#print merged_combined_results,"pre fido filter"
merged_combined_results=merged_combined_results[merged_combined_results['fido_filter']==1]
merged_combined_results=merged_combined_results.drop('fido_filter',axis=1)

#merged_combined_results_grp=merged_combined_results.groupby(numpy.arange(len(merged_combined_results))//multiprocessing.cpu_count())
#merged_combined_results=applyParallelHalf(merged_combined_results_grp,makeRoundedModSequence)
skyline_mod_sequences=merged_combined_results.apply(makeRoundedModSequence,axis=1)
merged_combined_results['Skyline Modified Sequence']=skyline_mod_sequences
#makeRoundedModSequence
merged_combined_results.to_csv("MERGED_OUTPUT.tsv",sep='\t',index=False)


#if 'total matches/spectrum' in merged_combined_results.columns:
#    merged_combined_results.drop('total matches/spectrum',1,inplace=True)#Just cleaning up...
#print set(merged_combined_results['_merge']),"that was the set..."
#print merged_combined_results
#sys.exit(0)
