import os, sys, re
import optparse
import shutil
import pandas
import numpy
import subprocess
#import time

#####################################
#This is the wrapper for luciphor2.
#It will be responsible for taking command line arguments and building a
#configuration file for running of luciphor2 to help localize PTMs on
#identified modified peptides, and ultimately producing an estimate of
#false localization of the modifciations.
#
#VERSION 1.10A
version="1.20A"
#DATE: 12/22/2015
date="12/22/2015"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the lucphor2 wrapper for Galaxy, Wohlschlegel Lab UCLA"
print "Written by William Barshop"
print "Version: ",version
print "Date: ",date


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

#print sys.argv,"THESE ARE THE ARGS"

parser = optparse.OptionParser()
parser.add_option("--pout",action="store",type="string",dest="operation_folder")
parser.add_option("--hcd",action="store_true",dest="hcd")
parser.add_option("--cid",action="store_true",dest="cid")
parser.add_option("--ppm",action="store_true",dest="ppm")
parser.add_option("--ms2tol",action="store", type="float", dest="ms2tol")
parser.add_option("--minmz",action="store", type="float", dest="min_mz")
parser.add_option("--fixed",action="store", type="string", dest="fixed_mods")
parser.add_option("--target",action="store",type="string",dest="target_mods")
parser.add_option("--variable",action="store",type="string",dest="variable_mods")
parser.add_option("--nl",action="store",type="string",dest="neutral_loss")
parser.add_option("--dm",action="store",type="float",dest="decoy_mass")
parser.add_option("--dnl",action="store",type="string",dest="decoy_neutral_loss")
parser.add_option("--mc",action="store",type="int",dest="max_charge")
parser.add_option("--min_psms",action="store",type="int",dest="min_psms")
parser.add_option("--ml",action="store",type="int",dest="max_length")
parser.add_option("--mp",action="store",type="int",dest="max_num_permutation")
parser.add_option("--mst",action="store",type="float",dest="model_score_threshold")
parser.add_option("--st",action="store",type="float",dest="score_threshold")
parser.add_option("--nt",action="store",type="int",dest="num_threads")
parser.add_option("--mzml",action="store",type="string",dest="mzml_files")
parser.add_option("--expgroups",action="store",type="string",dest="exp_group_file")


(options,args) = parser.parse_args()




aa_masses = {'A' : 71.037114, 'R' : 156.101111, 'N' : 114.042927, 'D' : 115.026943, 'C' : 103.009185, 'E' : 129.042593, 'Q' : 128.058578, 'G' : 57.021464, 'H' : 137.058912, 'I' : 113.084064, 'L' : 113.084064, 'K' : 128.094963, 'M' : 131.040485, 'F' : 147.068414, 'P' : 97.052764, 'S' : 87.032028, 'T' : 101.047679, 'U' : 150.95363, 'W' : 186.079313, 'Y' : 163.06332, 'V' : 99.068414 }
######### I scraped these AA masses from Matrix, who proudly bring us mascot....


#########HERE I NEED TO
#########CRAWL THROUGH THE POUT FOLDER AND
#########GET A LIST OF ALL THE SPECTRAL FILES AND SUUUCH
#OKKKKKAYYYYY Let's get a list of all the MZML files we got!
mzml_dict={}
for eachfile in options.mzml_files.split(","):
    if eachfile in mzml_dict:
        print "Somehow you gave me the same file twice.... Just a warning, but you might want to take a look at what's going on!"
    else:
        mzml_dict[eachfile]={}

#### Okay, now we'll need to read in the percolator experimental groups so that we
#### can build the input files for Luciphor2 appropriately! (that means get the groups right...)
if options.exp_group_file is None:
    print "I need to know how to group the experiments to properly build my Luciphor input files... Gimmie the file!"
    sys.exit(2)

# I'm lazy, so let's just read the file into a dataframe.
#print options.exp_group_file,type(options.exp_group_file)

#with open(options.exp_group_file,'r') as expfile:
group_information = pandas.read_csv(options.exp_group_file,sep='\t')


#print group_information,"this was group info...."
#print type(group_information)
run_dict={}

#And then for easy access, we'll make a dict of it.
for index,row in group_information.iterrows():
    run_dict[row['Crux File Integer']]=row
#for each in run_dict.keys():
#    print type(each)

infiles = []
for root, subFolders, files in os.walk(options.operation_folder):
    for eachfile in files:
        if 'target.psms.txt' in eachfile:
            infiles.append(str(os.path.join(root,eachfile)))

dataframe_vector=[]
for eachfile in infiles:
    newdf=pandas.DataFrame.from_csv(eachfile,sep='\t',index_col=False)
    dataframe_vector.append(newdf)
    del newdf


combined_results=pandas.concat(dataframe_vector)
#del dataframe_vector

#######Okay, so now we have all our percolator PSMs loaded into a 
#######single pandas dataframe.  We'll take that concatetnated dataframe
#######and then iterate over all the results and place in the file names
#######while constructing the input files.

theruns=[]
#for each_run in group_information['Original File Name']:
#    thisrun=each_run+".mzML"
#    theruns.append(thisrun)

combined_results['file_idx']=combined_results['file_idx'].astype(str)
combined_results['file']=combined_results['file'].astype(str)

#for each_run in group_information['Original File Name']:
#    thisrun=str(each_run)+".mzML"


new_results=[]
    

for each_idx in group_information['Crux File Integer']:
    #print type(each_idx),type(combined_results['file_idx'])
    mask = combined_results[(combined_results.file_idx == str(each_idx))]
    mask['file']=run_dict[int(each_idx)]['Original File Name']+".mzML"
    new_results.append(mask)
    #print run_dict[int(each_idx)]['Original File Name']+".mzML"
    
combined_results=pandas.concat(new_results)
output_results=pandas.concat(new_results)

####### Okay, now we've taken care of file names, let's just go ahead and build us some luciphor inputs!

luci_input_writers={}
luci_input_files=[] #This is just a list of the files we'll have to run luciphor on!
for eachgroup in set(group_information['Fractionation Group ID String']):
    mywriter=open(options.operation_folder+eachgroup+".pin_out/crux-output/"+"luciphor_spectra_"+eachgroup+".tsv",'w')
    luci_input_files.append(options.operation_folder+eachgroup+".pin_out/crux-output/"+"luciphor_spectra_"+eachgroup+".tsv")
    mywriter.write("srcFile\tscanNum\tcharge\tPSMscore\tpeptide\tmodSites\n")
    luci_input_writers[eachgroup]=mywriter


run_to_group={}
for index,row in group_information.iterrows(): #['Original File Name']:
    thisrun = row['Original File Name']+".mzML"
    thisgroup = row['Fractionation Group ID String']
    run_to_group[thisrun]=thisgroup



target_masses = []
non_decimal = re.compile(r'[^\d.]+')
for eachmod in options.target_mods.split(","):
    if eachmod is not "":
        target_masses.append(non_decimal.sub('',eachmod))
target_masses=list(set(target_masses))
global target_mass_conversion #CRAPPY WORKAROUND....
target_mass_conversion={}


for each in target_masses:
    target_mass_conversion[float(each)]=str(each)
#We're also going to actually put the rest of our mod masses in there, too...

for each in options.fixed_mods.split(","):
    #print each
    if each is not "":
        target_mass_conversion[float(non_decimal.sub('',each))]=str(non_decimal.sub('',each))
for each in options.variable_mods.split(","):
    #print each
    if each is not "":
        target_mass_conversion[float(non_decimal.sub('',each))]=str(non_decimal.sub('',each))



rounded_target_masses=[str(float(each)) for each in target_masses]
#print "These are the rounded targets...",rounded_target_masses

combined_results['percolator q-value']=1.0-combined_results['percolator q-value']   #Luciphor expects PSM scores to increase...  So we'll go ahead and give the score as 1-q.

print target_mass_conversion,"this is target mass conv."
regex = re.compile('[^a-zA-Z]')
for each_idx in group_information['Crux File Integer']:
    mask = combined_results[(combined_results['file_idx'] == str(each_idx))]
    #print '|'.join(target_masses),"This is target masses...."
    mask2=mask[(mask['sequence'].str.contains('|'.join(rounded_target_masses)))] #numpy.any(mods in mask.sequence for mods in target_masses)]
    #mask=mask.mask(applyme)
    #print mask,"THIS IS MASK"
    #mask2 = mask[mask.sequence.contains('|'.join(target_masses))]
    #print mask2,"this is mask2<------------------"
    ### OKAY, NOW YOU'LL NEED TO WRITE THE OUTPUT FOR THIS FILE
    for index,eachrow in mask2.iterrows():
        thispep=eachrow['sequence']
        fullseqlen=len(thispep)
        unmodpeplen=len(regex.sub('',thispep))
        newpep=[]
        moddict={}
        adjustment=0
        actualindex=0
        nummodsadded=0
        #print "---------------------------"
        #print "Working on :",thispep
        i=0
        while i < len(thispep):
            #print "now...",thispep[i]
            #print "i is...",i
            #print "actual index is...",actualindex
            if thispep[i].isalpha():
                newpep.append(thispep[i])
                actualindex+=1
            else:
                buf=""
                while (i < len(thispep) and (not thispep[i].isalpha())):
                    buf+=thispep[i]
                    i+=1
                #print target_mass_conversion
                for eachmod in target_mass_conversion:
                    if str(eachmod) in buf:
                        buf=buf.replace(str(eachmod),target_mass_conversion[eachmod])
                        break
                #print "Adding ",buf,"at",actualindex-1+nummodsadded
                #print str(float(buf.replace('[','').replace(']',''))),"and",+aa_masses[regex.sub('',thispep)[actualindex-1+nummodsadded]]
                moddict[actualindex-1+nummodsadded]=str(float(buf.replace('[','').replace(']',''))+aa_masses[regex.sub('',thispep)[actualindex-1+nummodsadded]])  #,regex.sub('',thispep)[actualindex-1+nummodsadded])   #### I REPLACED realstart with ACTUALINDEX here
                nummodsadded+=1
            i+=1
        #print "---------------------------"
        modstr=""
        for eachmodposition in sorted(moddict.keys()):
            modseq=moddict[eachmodposition]
            modstr+=str(eachmodposition)+"="+modseq+","
        modstr=modstr[:-1]
        luci_input_writers[run_to_group[eachrow['file']]].write(str(eachrow['file'])+'\t'+str(eachrow['scan'])+'\t'+str(eachrow['charge'])+'\t'+str(eachrow['percolator q-value'])+'\t'+regex.sub('',eachrow['sequence'])+'\t'+modstr+"\n")
    
#mask = combined_results


#for eachgroup in set(group_information['Fractionation Group ID String']):
#    myreader=open(options.operation_folder+"luciphor_spectra_"+eachgroup+".tsv",'rb')
#    i=0
#    for eachline in myreader:
#        i+=1
#        #print eachline
#        if i>50:
#            break


##### THIS WILL WRITE THE COMMON SECTION OF THE CONFIGURATION FILE

with open("shared_config.txt",'w') as configfile:
    configfile.write("SPECTRUM_SUFFIX = mzML\n")
    configfile.write("INPUT_TYPE = 1\n")

    #### HCD OR CID ####
    if options.hcd:
        configfile.write("ALGORITHM = 1\n")
    elif options.cid:
        configfile.write("ALGORITHM = 0\n")
    else:
        print "You gotta tell me if it's HCD or CID data..."
        sys.exit(2)
    ####################

    configfile.write("TSV_HEADER = 1\n")
    

    #### MODEL PSMS REQUIREMENTS ####
    if options.min_psms:
        configfile.write("MIN_NUM_PSMS_MODEL = "+str(options.min_psms)+"\n")
    #################################


    #### MASS ACCURACY STUFFS ####
    configfile.write("MS2_TOL = "+str(options.ms2tol)+"\n")
    if options.ppm:
        configfile.write("MS2_TOL_UNITS = 1\n")
    else:
        configfile.write("MS2_TOL_UNITS = 0\n")
    ##############################

    #### MIN MZ ####
    configfile.write("MIN_MZ = "+str(options.min_mz)+"\n")
    #################################




    #configfile.write("WRITE_MATCHED_PEAKS_FILE = 1\n")


    #### FIXED MODS ####
    for eachmod in options.fixed_mods.split(","):
        if eachmod is not "":
            configfile.write("FIXED MOD = "+eachmod.strip()+"\n")
    ####################

    #### VARIABLE MODS, UNSCORED ####
    for eachmod in options.variable_mods.split(","):
        if eachmod is not "":
            configfile.write("VAR_MOD = "+eachmod.strip()+"\n")
    #################################

    #### TARGET MODS, SCORED ####
    notargets=True
    for eachmod in options.target_mods.split(","):
        if eachmod is not "":
            notargets=False
            configfile.write("TARGET_MOD = "+eachmod.strip()+"\n")
            
    if notargets:
        print "What? No target modifications!? Give me something to score!"
        sys.exit(2)
    #############################

    #### NEUTRAL LOSSES, SCORED ####
    for eachnl in options.neutral_loss.split(","):
        if eachnl is not "":
            configfile.write("NL = "+eachnl.strip()+"\n")
    ################################

    #### DECOY MASS SETTING ####
    configfile.write("DECOY_MASS = "+str(options.decoy_mass)+"\n")
    if options.decoy_mass is None:
        print "No decoy mass? How am I supposed to calculate an FLR without a decoy mass!?"
        sys.exit(2)
    ############################

    #### DECOY NEUTRAL LOSSES ####
    for eachnl in options.decoy_neutral_loss.split(","):
        if eachnl is not "":
            configfile.write("DECOY_NL = "+eachnl+"\n")
    ##############################

    configfile.write("MAX_CHARGE_STATE = "+str(options.max_charge)+"\n")

    configfile.write("MAX_PEP_LEN = "+str(options.max_length)+"\n")

    configfile.write("MAX_NUM_PERM ="+str(options.max_num_permutation)+"\n")
    
    configfile.write("SELECTION_METHOD = 0\n") #THIS SCRIPT IS HARD CODED FOR USE WITH PERCOLATOR SCORES...

    configfile.write("MODELING_SCORE_THRESHOLD = "+str(1.0-float(options.model_score_threshold))+"\n")

    configfile.write("SCORING_THRESHOLD = "+str(1.0-float(options.score_threshold))+"\n")

    configfile.write("MIN_NUM_PSMS_MODEL = 50\n")

    configfile.write("MOD_PEP_REP = 0\n")

    configfile.write("NUM_THREADS = "+str(options.num_threads)+"\n")

    configfile.write("RUN_MODE = 0\n")

##### ABOVE THIS, FOR EACH RUN, WE WILL HAVE TO APPEND THE
##### SPECTRUM_PATH AND INPUT_DATA SECTIONS AND OUTPUT_FILE
##### now: maybe just input_data and output_file?
## the plan:
#foreachgroup in thegroups
## 1chdir to the folder
## 2write full config
## 3run
## 4fix psms output file


basedir=os.getcwd()

#We'll copy the mzML files to the proper folders...
#files = [f for f in os.listdir('.')]
#for f in files:
#    print f

for eachrun in run_to_group:
    #thislink=os.readlink(eachrun)
    #shutil.copy(thislink,options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)
    #thislink=os.readlink(eachrun)
    shutil.copy(eachrun,options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)


#And we'll go ahead and start the ball rolling!

os.chdir(basedir)


for eachgroup in set(group_information['Fractionation Group ID String']):
    os.chdir(basedir)
    shutil.copy("shared_config.txt",options.operation_folder+eachgroup+".pin_out/crux-output/luciphor_cfg_"+eachgroup+".txt")

for eachwriter in luci_input_writers:
    luci_input_writers[eachwriter].close()


for eachgroup in set(group_information['Fractionation Group ID String']):
    os.chdir(basedir)
    shutil.copy("shared_config.txt",options.operation_folder+eachgroup+".pin_out/crux-output/luciphor_cfg_"+eachgroup+".txt")
    os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
    with open("luciphor_cfg_"+eachgroup+".txt",'a') as file_writer:
        filestr="luciphor_spectra_"+eachgroup+".tsv"
        outstr="luciphor_output_"+eachgroup+".tsv"
        file_writer.write("SPECTRUM_PATH="+str(os.getcwd())+"\n")
        file_writer.write("INPUT_DATA="+filestr+"\n")
        file_writer.write("OUTPUT_FILE="+outstr+"\n")



processes=[]
#for eachgroup in set(group_information['Fractionation Group ID String']):
#    os.chdir(basedir)
#    #shutil.copy("shared_config.txt",options.operation_folder+eachgroup+".pin_out/crux-output/luciphor_cfg_"+eachgroup+".txt")
#    os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
#    #with open("luciphor_cfg_"+eachgroup+".txt",'a') as file_writer:
#    #    filestr="luciphor_spectra_"+eachgroup+".tsv"
#    #    outstr="luciphor_output_"+eachgroup+".tsv"
#    #    file_writer.write("SPECTRUM_PATH="+str(os.getcwd())+"\n")
#    #    file_writer.write("INPUT_DATA="+filestr+"\n")
#    #    file_writer.write("OUTPUT_FILE="+outstr+"\n")
#    command = "java -jar /galaxy-central/tools/wohl-proteomics/luciphor2/lucXor.jar "+"luciphor_cfg_"+eachgroup+".txt"
#    print "RUNNING: ",command,"in folder",str(os.getcwd())
#    #os.system(command)
#    p=(subprocess.Popen(command,shell=True)
#    p.wait()
#    if p is not 0:
#        p=subprocess.Popen(command,shell=True)
#        p.wait()
#    if p is not 0:
#        p=subprocess.Popen(command,shell=True)
#        p.wait()


with open(os.devnull, "w") as fnull:
    for eachgroup in set(group_information['Fractionation Group ID String']):
        os.chdir(basedir)
        os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
        #processes.append(subprocess.Popen(command,shell=True)
        command = "java -jar /galaxy-central/tools/wohl-proteomics/luciphor2/lucXor.jar "+"luciphor_cfg_"+eachgroup+".txt"
        processes.append(subprocess.Popen(command.split()))
        #processes.append(subprocess.Popen(command.split(),stdout=fnull, stderr=fnull)) # These pipes will hide the output from luciphor... this is only for DEBUG DEBUG DEBUG purposes...
        #out, err = p.communicate()
        #p=os.system(command)
        #print "Luciphor exited with a status of ",str(os.system(command))
        #time.sleep(20)
    

    for each_proc in processes:
        each_proc.wait()


###### INPUT-FILE CLEANUP ######
os.chdir(basedir)
for eachrun in run_to_group:
##    #thislink = os.readlink(eachrun)
##    os.unlink(options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)
    os.remove(options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)


##### READ IN OUTPUT FILES #####

output_results['file_idx']=output_results['file_idx'].astype(str)
output_results['file']=output_results['file'].astype(str)

def extractScanNum(x):
    return str(x['specId'].split(".")[1])

def extractFileName(x):
    return str(x['specId'].split(".")[0])

            

group_to_luci_results={}
for eachgroup in set(group_information['Fractionation Group ID String']):
    os.chdir(basedir)
    os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
    luci_result=pandas.read_csv("luciphor_output_"+str(eachgroup)+".tsv",sep='\t')
    #scan_extracted=luci_result['specId'].split(".")[1]
    scan_extracted=luci_result.apply(extractScanNum,axis=1)
    fname_extracted=luci_result.apply(extractFileName,axis=1)
    luci_result['scan']=scan_extracted
    luci_result['file']=fname_extracted+".mzML"
    #luci_result.assign(scan = lambda x: x['specId'].split(".")[1])
    #for each,row in luci_result.iterrows():
    #    print each,row
    #sys.exit(2)###################################################
    group_to_luci_results[eachgroup]=luci_result




   
group_to_runids={}
for index,row in group_information.iterrows():
    if row['Fractionation Group ID String'] in group_to_runids:
        group_to_runids[row['Fractionation Group ID String']].append(str(row['Crux File Integer']))
    else:
        #print row,"THIS IS ROW!2"
        group_to_runids[row['Fractionation Group ID String']]=[str(row['Crux File Integer'])]


##### We'll read in the raw inputs, and then append extra luciphor information to them! #####
groups_output={}
for each_group in group_to_runids:
    groups_output[each_group]=output_results[(output_results['file_idx'].str.contains('|'.join(group_to_runids[each_group])))]
    
def replaceMods(input):
    modDict=target_mass_conversion
    #print "input",input
    output=input
    #print modDict,"This is modDict!"
    for eachmod in modDict:
        if str(eachmod) in output and str(modDict[eachmod]) not in output:
            if float(eachmod)>0.0:
                output=output.replace(str(eachmod),"+"+str(modDict[eachmod]))
            else:
                output=output.replace(str(eachmod),str(modDict[eachmod]))
    #print "output",output
    return output



##### We'll use this dict to convert fileidx to file names
fileidx_to_file={}
for index,row in group_information.iterrows():
    thisrun = row['Original File Name']
    thisidx = str(row['Crux File Integer'])
    fileidx_to_file[thisidx]=thisrun

modsite_to_modmass={}
modsite_to_modmass_fixed={}
modsite_to_modmass_variable={}
for eachmod in options.target_mods.split(","):
    modsite_to_modmass[eachmod[:1]]=non_decimal.sub('',eachmod)
for eachmod in options.fixed_mods.split(","):
    modsite_to_modmass_fixed[eachmod[:1]]=non_decimal.sub('',eachmod)

for eachmod in options.variable_mods.split(","):                        ###### NB: At the moment, I'm not sure how LuciPhor handles these...
    modsite_to_modmass_variable[eachmod[:1]]=non_decimal.sub('',eachmod)######     So I'll make this for later, but I don't know what to do with it... yet...


#print "This is the modsite dict!",modsite_to_modmass




##### We'll have to cram the luciphor results into the appropriate files... and modify the outputs!
modified_output={}

#### DEBUG LOGGING:
#logwriter=open(options.operation_folder+"luci_log.txt",'wb')





##########################################NEWNEWNEWNEWNEWNEWNEWNEW
def combineFileNameAndScanLuci(x):
    return str(x['specId'].split(".")[0]+".mzML."+x['specId'].split(".")[1])
def combineFileNameAndScanPerco(x):
    return str(x['file']+"."+str(x['scan']))

newlist=[]
for eachgroup in groups_output:
    newlist.append(groups_output[eachgroup]) 
combined_groups=pandas.concat(newlist)#This is the combined percolator results...
combined_fname_scans=combined_groups.apply(combineFileNameAndScanPerco,axis=1)
combined_groups['unique_name']=combined_fname_scans
del newlist

newlist=[]
for eachgroup in group_to_luci_results:
    newlist.append(group_to_luci_results[eachgroup])

combined_luci_groups=pandas.concat(newlist)#This is the combined luciphor results...
combined_fname_scans_luci=combined_luci_groups.apply(combineFileNameAndScanLuci,axis=1)
combined_luci_groups['unique_name']=combined_fname_scans_luci
del newlist


combined_luci_groups=combined_luci_groups.drop('scan',1)
combined_luci_groups=combined_luci_groups.drop('file',1)


merged_df=pandas.merge(combined_groups,combined_luci_groups, on="unique_name",sort=True,how="outer",indicator=False)
merged_df['numPPS'].fillna("0")
merged_df['numRPS'].fillna("0")
merged_df=merged_df.drop_duplicates(subset='unique_name')
merged_df['luci_sequence']=""

#merge_numPPS_mask=merged_df[merged_df['numPPS']>=1]

for index,eachrow in merged_df.iterrows():
    if eachrow['numRPS'] < 1:
        zz=0
        newPep=[]
        while zz<len(eachrow['sequence']):
            if not eachrow['sequence'][zz]=='[':
                newPep.append(eachrow['sequence'][zz])
            else:
                buffer=newPep.pop()
                while zz<len(eachrow['sequence']) and eachrow['sequence'][zz]!=']':
                    buffer+=eachrow['sequence'][zz]
                    zz+=1
                buffer+=eachrow['sequence'][zz]
                buffer=buffer[:2]+"+"+target_mass_conversion[buffer[2:-1]]+"]" #We'll cram in the full mass accuracy mod as given in the settings...
                newPep.append(buffer)
            zz+=1

        merged_df.loc[index,'luci_sequence']=''.join(newPep)
    elif eachrow['numRPS'] >= 1:
        zz=0
        newPep=[]
        while zz<len(eachrow['predictedPep1']):
            if not eachrow['predictedPep1'][zz].islower():
                site=eachrow['predictedPep1'][zz]
                if eachrow['predictedPep1'][zz] in modsite_to_modmass_fixed:
                    if float(modsite_to_modmass_fixed[site])>0.0:
                        newPep.append(site+"[+"+modsite_to_modmass_fixed[site]+"]")
                    else:
                        newPep.append(site+"["+modsite_to_modmass_fixed[site]+"]")
                else:
                    newPep.append(eachrow['predictedPep1'][zz])
            else:

                modsite=eachrow['predictedPep1'][zz].upper()
                if modsite in modsite_to_modmass:
                    if float(modsite_to_modmass[modsite])>0.0:
                        newPep.append(modsite+"[+"+modsite_to_modmass[modsite]+"]")
                    else:
                        newPep.append(modsite+"["+modsite_to_modmass[modsite]+"]")
                    if modsite in modsite_to_modmass_fixed:
                        if float(modsite_to_modmass_fixed[modsite])>0.0:
                            newPep.append("[+"+modsite_to_modmass_fixed[modsite]+"]")
                        else:
                            newPep.append("["+modsite_to_modmass_fixed[modsite]+"]")
                elif site in modsite_to_modmass_variable:
                    #if site in modsite_to_modmass_variable: ##### IS THIS HOW LUCIPHOR HANDLES VARIABLE MODS?
                    if float(modsite_to_modmass_variable[site])>0.0:
                        newPep.append("[+"+modsite_to_modmass_variable[site]+"]")
                    else:
                        newPep.append("["+modsite_to_modmass_variable[site]+"]")


            zz+=1
        merged_df.loc[index,'luci_sequence'] = ''.join(newPep)

merged_df.rename(columns={"specId":"luci_specId","peptide":"luci_peptide","pepProphet":"q-1","predictedPep1":"luci_predictedPep1","predictedPep2":"luci_predictedPep2","numPPS":"luci_numPPS","numRPS":"luci_numRPS","deltaScore":"luci_deltaScore","pep1score":"luci_pep1score","pep2score":"luci_pep2score","globalFLR":"luci_globalFLR","localFLR":"luci_localFLR"}, inplace=True)
merged_df=merged_df.drop('unique_name',1)


for each_idx in set(merged_df['file_idx']):
    modified_output[each_idx]=merged_df[merged_df['file_idx']==str(each_idx)]







##########################################NEWNEWNEWNEWNEWNEWNEWNEW
##########################################NEWNEWNEWNEWNEWNEWNEWNEW




'''
for each_group_df in groups_output:
    modify_me=groups_output[each_group_df].copy(deep=True)
    fixed_sequences=modify_me['sequence'].apply(replaceMods)
    modify_me['sequence']=fixed_sequences
    for each_idx in set(modify_me['file_idx']):
        one_run_modify=modify_me[modify_me['file_idx']==str(each_idx)].copy(deep=True)#deep=True)  # 10/14/2015 added copy ....
        one_run_modify['scan']=one_run_modify['scan'].astype(int)
        temp=group_to_luci_results[each_group_df].copy(deep=True)
        temp['scan']=temp['scan'].astype(int)
        print set(temp['file']),"runs in group",each_group_df
        one_run_luciphor=temp[temp['file'].str.contains(fileidx_to_file[str(each_idx)])]
	print one_run_luciphor,"this is after filtering for<------------"+fileidx_to_file[str(each_idx)]
        merged_df=pandas.merge(one_run_modify,one_run_luciphor, on="scan",sort=True,how="outer",indicator=False)
        merged_df['numPPS'].fillna("0")
        merged_df['numRPS'].fillna("0")
        os.chdir(basedir)

        #We'll iterate first over every scan in the merged output.
        for index,eachrow in merged_df.iterrows():
            if eachrow['numPPS'] >= 1:
                zz=0
                newPep=[]
                while zz<len(eachrow['predictedPep1']):
                    if not eachrow['predictedPep1'][zz].islower():
                        site=eachrow['predictedPep1'][zz]
                        if eachrow['predictedPep1'][zz] in modsite_to_modmass_fixed:
                            if float(modsite_to_modmass_fixed[site])>0.0:
                                newPep.append(site+"[+"+modsite_to_modmass_fixed[site]+"]")
                            else:
                                newPep.append(site+"["+modsite_to_modmass_fixed[site]+"]")
                        else:
                            newPep.append(eachrow['predictedPep1'][zz])
                    else:
                       
                        modsite=eachrow['predictedPep1'][zz].upper()
                        if modsite in modsite_to_modmass:
                            if float(modsite_to_modmass[modsite])>0.0:
                                newPep.append(modsite+"[+"+modsite_to_modmass[modsite]+"]")
                            else:
                                newPep.append(modsite+"["+modsite_to_modmass[modsite]+"]")
                            if modsite in modsite_to_modmass_fixed:
                                if float(modsite_to_modmass_fixed[modsite])>0.0:
                                    newPep.append("[+"+modsite_to_modmass_fixed[modsite]+"]")
                                else:
                                    newPep.append("["+modsite_to_modmass_fixed[modsite]+"]")
                        elif site in modsite_to_modmass_variable:
                            #if site in modsite_to_modmass_variable: ##### IS THIS HOW LUCIPHOR HANDLES VARIABLE MODS?
                            if float(modsite_to_modmass_variable[site])>0.0:
                                newPep.append("[+"+modsite_to_modmass_variable[site]+"]")
                            else:
                                newPep.append("["+modsite_to_modmass_variable[site]+"]")


                    zz+=1
                merged_df.loc[index,'sequence'] = ''.join(newPep)
        merged_df.rename(columns={"specId":"luci_specId","peptide":"luci_peptide","pepProphet":"q-1","predictedPep1":"luci_predictedPep1","predictedPep2":"luci_predictedPep2","numPPS":"luci_numPPS","numRPS":"luci_numRPS","deltaScore":"luci_deltaScore","pep1score":"luci_pep1score","pep2score":"luci_pep2score","globalFLR":"luci_globalFLR","localFLR":"luci_localFLR"}, inplace=True)
        modified_output[each_idx]=merged_df
'''
###########################################################################





idx_to_group={}
for index,row in group_information.iterrows(): #['Original File Name']:
    thisidx = str(row['Crux File Integer'])
    thisgroup = row['Fractionation Group ID String']
    idx_to_group[thisidx]=thisgroup




assembled_output={}

for eachkey in modified_output:
    group=idx_to_group[eachkey]
    if group in assembled_output:
        assembled_output[group].append(modified_output[eachkey])
    else:
        assembled_output[group]=[modified_output[eachkey]]
#print assembled_output.keys(),"these are the keys"
###### Now we'll go ahead and write out each group to a file in the approp. location!
for eachgroup in assembled_output:
    os.chdir(basedir)
    os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
    groupdf=pandas.concat(assembled_output[eachgroup])
    shutil.move(eachgroup+".percolator.target.psms.txt",eachgroup+".percolator.target.psms_uncorrected.txt")
    groupdf.to_csv(eachgroup+".percolator.target.psms.txt",sep="\t",index=False)
    
#for eachgroup in set(group_information['Fractionation Group ID String']):
#    os.chdir(basedir)
#    shutil.copy("shared_config.txt",options.operation_folder+eachgroup+".pin_out/crux-output/luciphor_cfg_"+eachgroup+".txt")
#    os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
#    with open("luciphor_cfg_"+eachgroup+".txt",'a') as file_writer:




#####for eachgroup in set(group_information['Fractionation Group ID String']):
#####    pass
#for eachgroup in set(group_information['Fractionation Group ID String']):
#    mywriter=open(options.operation_folder+eachgroup+".pin_out/crux-output/"+"luciphor_spectra_"+eachgroup+".tsv",'w')
#    luci_input_files.append(options.operation_folder+eachgroup+".pin_out/crux-output/"+"luciphor_spectra_"+eachgroup+".tsv")

#shutil.copy("shared_config.txt","/home/galaxy/shared_config.txt")
    



#parser.print_usage()

#print options.max_length
#for each in options:
#    print each
#logwriter.close()
print "-----------------------------------------------------------------------"

