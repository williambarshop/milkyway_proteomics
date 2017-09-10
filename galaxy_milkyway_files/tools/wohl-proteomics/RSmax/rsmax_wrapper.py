import os, sys, re
import optparse
import shutil
import pandas
import numpy
import subprocess
import math
import pyopenms
import gc
import operator
from pyopenms import *
from tqdm import *
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import cpu_count
#from pyopenms.sysinfo import free_mem

import time


#####################################
#This is the wrapper for RSMax... Adapted from the Luciphor2 wrapper.
#It will be responsible for taking command line arguments and building a
#input file for running of RSMax to help localize PTMs on
#identified modified peptides, and ultimately producing an estimate of
#confidence of PTM localization for the peptides
#Now automatically adds unscored mods!
#
#VERSION 0.95A
version="0.95A"
#DATE: 08/10/2016
date="08/10/2016"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the RSmax wrapper for Galaxy, Wohlschlegel Lab UCLA"
print "Written by William Barshop"
print "Version: ",version
print "Date: ",date

basedir=os.getcwd()

##############
#file = MzMLFile()
#experiment = MSExperiment()
#file.load(inputfile, experiment)
#That's temporary storage for pyopenms stuff...
##for spectrum in experiment:
##    if spectrum.getMSLevel() == 2:
##        precursor = spectrum.getPrecursors()
##        precursorMZ = precursor[0].getMZ()
##        precursorCharge = precursor[0].getCharge()
##        precursorUnchargedMass = precursor[0].getUnchargedMass()
##        #print precursorUnchargedMass
##        #pause = raw_input("Pause")
##        #precursorMZ=precursor[0][0]
##        spec_array = spectrum.get_peaks()
##        instrumentSettings = spectrum.getInstrumentSettings()
##        acquisitionInfo = spectrum.getAcquisitionInfo()
##        nativeID = spectrum.getNativeID()

def applyParallel(dfGrouped, func):
    retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group) for name, group in dfGrouped)
    return pandas.concat(retLst)


def fix_modification_parallel(group):
    for index,eachrow in group.iterrows():
        if ((eachrow['numModSites'] < 1) and (']' in eachrow['sequence'])):
            #pass
            #print eachrow
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
                    buffer=buffer[:2]+buffer[2:-1]+"]" #We'll cram in the full mass accuracy mod as given in the settings...
                    newPep.append(buffer)
                zz+=1
            #print group.loc[index,'pRS_sequence'],"BEFORE - 1"
            group.loc[index,'pRS_sequence']=''.join(newPep)
            #print group.loc[index,'pRS_sequence'],"AFTER - 1"

        elif eachrow['numModSites'] >= 1:
            #print group.loc[index,'pRS_sequence'],"BEFORE - 2"
            numPPS=0
            #print "----------------------------------"
            for each_aa in target_mod_AA:
                #print "looking for ",each_aa
                numPPS+=eachrow['sequence'].count(each_aa)
                #print "found ",str(eachrow['sequence'].count(each_aa))

            #print numPPS, "in peptide",eachrow['sequence']
            #print "----------------------------------"

            #original_peptide=eachrow['sequence']
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
                    #print buffer
                    if massStr_to_modID[buffer[3:-1]]=='1': # if this mass is the target mod, then we'll go ahead and drop it and only add back the AA.
                        buffer=buffer[:1]
                    else:
                        buffer=buffer[:-1]+"]" #We'll cram in the full mass accuracy mod as given in the settings...
                    newPep.append(buffer)
                zz+=1

            #now we will read the phosphoRS modification string and determine which positions to place the target modification!
            site_scores_split=eachrow['pRS_siteProbabilities'].split(';')
            site_dict={}
            for eachsite in site_scores_split:
                score=float(eachsite.split()[1])
                site_index = int(eachsite.split()[0].split('(')[1].split(')')[0])-1 #Corrected for index positions (now)!
                site_dict[site_index]=score
            sorted_sites = sorted(site_dict.items(), key=operator.itemgetter(1), reverse=True)
            #print sorted_sites,"THIS IS SORTED SITES!"
            #now we're going to slip target masses into each appropriate index!

            number_of_sites=eachrow['numModSites']
            ns = 0
            while ns<number_of_sites:
                newPep[int(sorted_sites[ns][0])]+="[+"+str(target_mod_mass)+"]"
                #print sorted_sites[ns],"this will be the ",ns,"th site to add a target mod!"
                ns+=1
            #print "replacing "+eachrow['sequence']+" with "+''.join(newPep)
            group.loc[index,'pRS_sequence']=''.join(newPep)
            group.loc[index,'pRS_numPPS']=numPPS
            #print group.loc[index,'pRS_sequence'],"AFTER - 2"
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

#print sys.argv,"THESE ARE THE ARGS"

parser = optparse.OptionParser()
parser.add_option("--pout",action="store",type="string",dest="operation_folder")
parser.add_option("--tool_dir",action="store",type="string",dest="tool_dir")
parser.add_option("--ms2tol",action="store", type="float", dest="ms2tol")
parser.add_option("--target",action="store",type="string",dest="target_mods")
#parser.add_option("--variable",action="store",type="string",dest="variable_mods")
parser.add_option("--nl",action="store",type="string",dest="neutral_loss")
parser.add_option("--nt",action="store",type="int",dest="num_threads")
parser.add_option("--mzml",action="store",type="string",dest="mzml_files")
parser.add_option("--expgroups",action="store",type="string",dest="exp_group_file")
parser.add_option("--mp",action="store",type="int",dest="max_permutations")
parser.add_option("--scorenl",action="store_true",dest="score_nl")
parser.add_option("--filter",action="store", type="float", dest="filter")
#"--scorenl"

(options,args) = parser.parse_args()



aa_masses = {'A' : 71.037114, 'R' : 156.101111, 'N' : 114.042927, 'D' : 115.026943, 'C' : 103.009185, 'E' : 129.042593, 'Q' : 128.058578, 'G' : 57.021464, 'H' : 137.058912, 'I' : 113.084064, 'L' : 113.084064, 'K' : 128.094963, 'M' : 131.040485, 'F' : 147.068414, 'P' : 97.052764, 'S' : 87.032028, 'T' : 101.047679, 'U' : 150.95363, 'W' : 186.079313, 'Y' : 163.06332, 'V' : 99.068414 }
######### I scraped these AA masses from Matrix, who proudly bring us Mascot!(ew)


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

filter = 1.0 #This allows everything as the default!
if options.filter is not None:
    filter=options.filter


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
    newdf=pandas.read_csv(eachfile,sep='\t',index_col=False)
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
combined_results['percolator q-value']=combined_results['percolator q-value'].astype(float)

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
    mywriter=open(options.operation_folder+eachgroup+".pin_out/crux-output/"+eachgroup+"_rsmax.txt",'w')
    luci_input_files.append(options.operation_folder+eachgroup+".pin_out/crux-output/"+eachgroup+"_rsmax.txt")
    #mywriter.write("specId\tpeptide\tmodPosition\tmodMap\tactivationType\tcharge\tpeak list\tprecursormz\n")
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
        target_masses.append(non_decimal.sub('',eachmod.split()[1]))
        #print "appended",non_decimal.sub('',eachmod),"from the input",eachmod.split()[1]
target_masses=list(set(target_masses))



#modId,ScoredmodMassDouble,NLmassDouble,ScoredModificationName,ScoredShortModName,ScoredModSiteAAs,ScoredElementalModificationStr
#We'll take in the info about the target mod...
target_mod_split = options.target_mods.split()
target_mod_AA = target_mod_split[0]
target_mod_mass = target_mod_split[1]
target_mod_elements = target_mod_split[2]

#print options.neutral_loss,type(options.neutral_loss),"yo"
if not options.neutral_loss is "":
    target_mod_NL_split = options.neutral_loss.split()
else:
    target_mod_NL_split = "A 0.00 H".split() #We're not using a NL...
target_mod_NL_AA = target_mod_NL_split[0]
target_mod_NL_mass = str(math.fabs(float(target_mod_NL_split[1])))
target_mod_NL_elements = target_mod_NL_split[2]


#We'll set up symlinks for each file...
for eachrun in run_to_group:
    thislink=os.readlink(eachrun)
    shutil.copy(thislink,options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)





#rounded_target_masses=str(target_mod_mass)
#print "These are the rounded targets...",rounded_target_masses


#print target_mass_conversion,"this is target mass conv."
regex = re.compile('[^a-zA-Z]')


#We'll set up the mass conversions before we get started on anything else...
#We're going to make a reference dictionary to convert the mass (as str) to the modid... which we'll go ahead and assign now!
#First, we'll make the modId 1 (scored modification)
massStr_to_modID = {str(target_mod_mass) : '1'}
#And we'll iterate that value up one for the next modification we introduce...
id_itr=2

#We're going to make an unscored_mod_dict, which will store key is mass(as str), value is AAs (as [str,str,str])
unscored_mod_dict={}



for each_idx in group_information['Crux File Integer']:

    file = MzMLFile()
    experiment = MSExperiment()
    print "loading file "
    print run_dict[int(each_idx)]['Original File Name']+".mzML"#run_dict[int(each_idx)]
    print ""
    file.load(os.path.join(basedir,run_dict[int(each_idx)]['Original File Name']+".mzML"), experiment)
    #That's temporary storage for pyopenms stuff...
    ##for spectrum in experiment:
    ##    if spectrum.getMSLevel() == 2:
    ##        precursor = spectrum.getPrecursors()
    ##        precursorMZ = precursor[0].getMZ()
    ##        precursorCharge = precursor[0].getCharge()
    ##        precursorUnchargedMass = precursor[0].getUnchargedMass()
    ##        #print precursorUnchargedMass
    ##        #pause = raw_input("Pause")
    ##        #precursorMZ=precursor[0][0]
    ##        spec_array = spectrum.get_peaks()
    ##        instrumentSettings = spectrum.getInstrumentSettings()
    ##        acquisitionInfo = spectrum.getAcquisitionInfo()
    ##        nativeID = spectrum.getNativeID()






    mask = combined_results[(combined_results['file_idx'] == str(each_idx))]
    #print '|'.join(target_masses),"This is target masses...."
    mask2=mask[(mask['sequence'].str.contains('|'.join(target_masses)))] #numpy.any(mods in mask.sequence for mods in target_masses)]
    #print mask2,"this is before the q-value filter..."
    mask3=mask2[(mask2['percolator q-value'] <= filter)]
    #print mask3,"THIS WAS THE FILTERED FINAL MASK FOR THIS RUN..."
    #sorted_mask2=mask2.sort_values('scan',ascending=True)
    #mask=mask.mask(applyme)
    #print mask,"THIS IS MASK"
    #mask2 = mask[mask.sequence.contains('|'.join(target_masses))]
    #print mask2,"this is mask2<------------------"
    ### OKAY, NOW YOU'LL NEED TO WRITE THE OUTPUT FOR THIS FILE
    print "Generating input files...."
    for index,eachrow in tqdm(mask3.iterrows()):
        thispep=eachrow['sequence']
        fullseqlen=len(thispep)
        unmodpeplen=len(regex.sub('',thispep))
        newpep=[]
        adjustment=0
        actualindex=0
        nummodsadded=0
        #print "---------------------------"
        #print "Working on :",thispep

        ################ So we're going to construct the ModMap for this PSM.

        #We'll iterate over the peptide sequence, reading the modifications and making a modmap of the PhosphoRS standard format.
        # 0.000100200100.0
        modmap="0."
        target_modsites=[]
        i=0
        while i < len(thispep):
            if thispep[i].isalpha():
                modmap+="0"
                #newpep.append(thispep[i])
                actualindex+=1
            elif thispep[i]=="[":
                #Now we know we're in a mod.
                mod_start_pos=i
                mod_end_pos=i #changed these from actualindex...
                actualindex-=1 #Because we're seeing a mod for the AA right before this.
                i-=1
                while thispep[mod_end_pos]!="]":
                    mod_end_pos+=1
                mod=thispep[mod_start_pos:mod_end_pos+1]
                #print "THIS IS MOD",mod
                diff=mod_end_pos+1-mod_start_pos
                #i=i-diff#correction for the peptide losing stuff
                #print thispep, "BEFORE"
                thispep=thispep[:mod_start_pos]+thispep[mod_end_pos+1:]
                #print thispep, "AFTER"
                modmap=modmap[:-1]#We'll remove the last 0 of modmap to place in the correct modificaiton info...
                #print massStr_to_modID,"this is the dictionary of interest..."
                try:
                    #print "Handling mod",mod[2:-1]
                    modmap+=massStr_to_modID[mod[2:-1]]
                    if massStr_to_modID[mod[2:-1]]!='1':
                        if thispep[i] in unscored_mod_dict[mod[2:-1]]:
                            pass #The AA is already known for this modification.
                        else:
                            unscored_mod_dict[mod[2:-1]].append(thispep[i])
                            print "Adding the AA "+thispep[i]+" to the modification of mass "+mod[2:-1],"from the peptide",thispep
                except: #Mod not yet in dict... let's add it!
                    unscored_mod_dict[mod[2:-1]]=[thispep[i]]
                    #And we'll add it to the modID dict, too, and iterate the id_itr for the next possible mod...
                    massStr_to_modID[str(mod[2:-1])]=str(id_itr)
                    id_itr+=1

                    modmap+=massStr_to_modID[mod[2:-1]]
                    print "We've added the unscored mod of mass "+mod[2:-1]+" for the AA "+thispep[i],"from the peptide",thispep
                    pass #Add mod here!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if massStr_to_modID[mod[2:-1]] == '1':
                    target_modsites.append(str(len(modmap)-3))
                    #target_modsites.append(str(actualindex))
            i+=1
            actualindex+=1 # THIS WAS ADDED 12/15/2015 TO CORRECT FOR INDEX SHIFTING
        modmap+=".0"
        #print modmap,"This is modmap!"
        mod_positions=",".join(target_modsites)
        specID=eachrow['file'].replace('.mzML','')+"."+str(eachrow['scan'])+'.'+str(eachrow['scan'])+'.'+str(eachrow['charge'])
        #print specID,"this is specID"
        spec_peaks=experiment[eachrow['scan']-1].get_peaks()#THIS WILL HAVE TO BE FORMATTED CORRECTLY....
        #print spec_peaks,"These are the peaks in a RAW format..."
        formatted_spec_peaks=""
        #for eachpeak in spec_peaks:
        #    formatted_spec_peaks+=str(round(eachpeak[0],6))+":"+str(round(eachpeak[1],6))+"," #6 decimal places of float accuracy...
        for i in xrange(len(spec_peaks[0])):
            formatted_spec_peaks+=str(round(spec_peaks[0][i],6))+":"+str(round(spec_peaks[1][i],6))+"," #6 decimal places of float accuracy...
        formatted_spec_peaks=formatted_spec_peaks[:-1]
        charge=eachrow['charge']
        activation_method="HCD"#by default
        activation_type=experiment[eachrow['scan']-1].getPrecursors()[0].getActivationMethods()
        precursorMZ=experiment[eachrow['scan']-1].getPrecursors()[0].getMZ()
        if 0 in activation_type:
            activation_method="CID"
        elif 5 in activation_type:
            activation_method="ECD"
        elif 8 in activation_type:
            activation_method="HCD"
        elif 11 in activation_type:
            activation_method="ETD"
        else:
            print "I DON\'T KNOW WHAT THE ACTIVATION METHOD IS.... SO WE\'LL CALL IT HCD...."


        #print "---------------------------"
        #FOR PHSOPHORS, thispep = peptide sequence, cleaned
        #               modmap = modmap
        #               mod_positions = positions of modifications "14,18"
        #mywriter.write("specId\tpeptide\tmodPosition\tmodMap\tactivationType\tcharge\tpeak list\tprecursormz\n")

        luci_input_writers[run_to_group[eachrow['file']]].write(specID+'\t'+thispep+'\t'+mod_positions+'\t'+modmap+'\t'+activation_method+'\t'+str(eachrow['charge'])+'\t'+formatted_spec_peaks+'\t'+str(precursorMZ)+'\n')################### ACTIVATION TYPE, PEAKS, PRECURSORMZ
        #luci_input_writers[run_to_group[eachrow['file']]].write(str(eachrow['file'])+'\t'+str(eachrow['scan'])+'\t'+str(eachrow['charge'])+'\t'+str(eachrow['percolator q-value'])+'\t'+regex.sub('',eachrow['sequence'])+'\t'+modstr+"\n")
    del file
    del experiment
    gc.collect()
#mask = combined_results


#And finally, we'll make the mod strings for running pRS via rsmax.exe
unscored_mod_str=""
print "This is the unscored mod dict...",unscored_mod_dict
for each_mod in unscored_mod_dict:
    unscored_mod_str+=massStr_to_modID[str(each_mod)]+","+str(each_mod)+","+"0.0"+","+"unscored"+massStr_to_modID[str(each_mod)]+","+"mod"+massStr_to_modID[str(each_mod)]+","+"".join(unscored_mod_dict[each_mod])+","+"XXXXX "





os.chdir(basedir)

for eachwriter in luci_input_writers:
    luci_input_writers[eachwriter].close()

if options.score_nl:
    score_neutrals="true"
else:
    score_neutrals="false"


processes=[]
#RSmax.exe /path/to/inputs/ scoreNLBoolean massToleranceDouble modId,ScoredmodMassDouble,NLmassDouble,ScoredModificationName,ScoredShortModName,ScoredModSiteAAs,ScoredElementalModificationStr
#modId,ScoredmodMassDouble,NLmassDouble,ScoredModificationName,ScoredShortModName,ScoredModSiteAAs,ScoredElementalModificationStr

target_modstr=massStr_to_modID[str(target_mod_mass)]+","+str(target_mod_mass)+","+target_mod_NL_mass+","+"targetMod"+","+"targ"+","+target_mod_AA+","+target_mod_elements
####COME BACK

with open(os.devnull, "w") as fnull:
    for eachgroup in tqdm(set(group_information['Fractionation Group ID String'])):
        os.chdir(basedir)
        os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
        #processes.append(subprocess.Popen(command,shell=True)
        ###command = "mono /nfs/galaxy/galaxy/tools/wohl-proteomics/RSmax/RSmax.exe "+os.getcwd()+" "+eachgroup+" "+score_neutrals+" "+str(options.ms2tol)+" "+target_modstr+" "+unscored_mod_str
        command = "mono "+options.tool_dir+"/RSmax.exe "+os.getcwd()+" "+eachgroup+" "+score_neutrals+" "+str(options.ms2tol)+" "+target_modstr+" "+unscored_mod_str
        print "running command... ",command
        #processes.append(subprocess.Popen(command.split()))
        #proc=subprocess.Popen(command.split())
        #proc.wait()
        ############subprocess.call(command.split(), stderr=sys.stdout.fileno())  #### SEEMED TO WORK ONCE UP ON A TIME....
        #processes.append(subprocess.Popen(command.split(),stdout=fnull, stderr=fnull)) # These pipes will hide the output from luciphor... this is only for DEBUG DEBUG DEBUG purposes...
        #out, err = p.communicate()
        p=os.system(command)
        #print "Luciphor exited with a status of ",str(os.system(command))
        time.sleep(5)
        print "Finished scoring this command...!"

    #for each_proc in processes:
    #    each_proc.wait()


###### INPUT-FILE CLEANUP ######
#This gets rid of the mzML files...
print "Beginning file cleanup..."
os.chdir(basedir)
for eachrun in run_to_group:
    os.remove(options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)
###    #thislink = os.readlink(eachrun)
###    os.unlink(options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)


##### READ IN OUTPUT FILES #####

output_results['file_idx']=output_results['file_idx'].astype(str)
output_results['file']=output_results['file'].astype(str)

#def extractScanNum(x):
#    return str(x['specId'].split(".")[1])

#def extractFileName(x):
#    return str(x['specId'].split(".")[0])

            

group_to_results={}
print "reading results in..."
for eachgroup in tqdm(set(group_information['Fractionation Group ID String'])):
    os.chdir(basedir)
    os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
    rsmax_result=pandas.read_csv("output_"+str(eachgroup)+"_rsmax.txt",sep='\t')
    #scan_extracted=rsmax_result['specId'].split(".")[1]
    #scan_extracted=rsmax_result.apply(extractScanNum,axis=1)
    #fname_extracted=rsmax_result.apply(extractFileName,axis=1)
    rsmax_result['scan']=rsmax_result['scan'].astype(int)
    #rsmax_result.assign(scan = lambda x: x['specId'].split(".")[1])
    #for each,row in rsmax_result.iterrows():
    #    print each,row
    #sys.exit(2)###################################################
    group_to_results[eachgroup]=rsmax_result


   
group_to_runids={}
group_to_files={}
for index,row in group_information.iterrows():
    if row['Fractionation Group ID String'] in group_to_runids:
        group_to_runids[row['Fractionation Group ID String']].append(str(row['Crux File Integer']))
    else:
        #print row,"THIS IS ROW!2"
        group_to_runids[row['Fractionation Group ID String']]=[str(row['Crux File Integer'])]
    if row['Fractionation Group ID String'] in group_to_files:
        group_to_files[row['Fractionation Group ID String']].append(str(row['Original File Name']))
    else:
        #print row,"THIS IS ROW!2"
        group_to_files[row['Fractionation Group ID String']]=[str(row['Original File Name'])]


##### We'll read in the raw inputs, and then append extra PhosphoRS information to them! #####
groups_output={}
for each_group in group_to_files:
    groups_output[each_group]=output_results[(output_results['file'].str.contains('|'.join(group_to_files[each_group])))]
    

##### We'll use this dict to convert fileidx to file names
fileidx_to_file={}
for index,row in group_information.iterrows():
    thisrun = row['Original File Name']
    thisidx = str(row['Crux File Integer'])
    fileidx_to_file[thisidx]=thisrun

modsite_to_modmass={}
modsite_to_modmass_fixed={}
modsite_to_modmass_variable={}
print options.target_mods,"target mods!"
print options.variable_mods,"unscored mods!"

#for eachmod in options.target_mods.split():
#    modsite_to_modmass[eachmod[:1]]=non_decimal.sub('',eachmod)

tmod_aa=options.target_mods.split()[0]
tmod_mass=options.target_mods.split()[1]
for each_aa in tmod_aa:
    modsite_to_modmass[each_aa]=tmod_mass


for eachmod in options.variable_mods.split(","):
    mod_aa=eachmod.split()[0]
    mod_mass=eachmod.split()[1]
    for each_aa in mod_aa:
        modsite_to_modmass[each_aa]=mod_mass
print "This is the modsite dict!",modsite_to_modmass


#for eachmod in options.fixed_mods.split(","):
#    modsite_to_modmass_fixed[eachmod[:1]]=non_decimal.sub('',eachmod)




##### unscored_mod_dict={}#key is mass(as str), value is AAs (as str)  #This exists here in the pRS implementation!


#for eachmod in options.variable_mods.split(","):                        ###### NB: At the moment, I'm not sure how LuciPhor handles these...
#    modsite_to_modmass_variable[eachmod[:1]]=non_decimal.sub('',eachmod)######     So I'll make this for later, but I don't know what to do with it... yet...




#####
#####
##### We'll have to cram the phosphoRS results into the appropriate files... and modify the outputs!
modified_output={}

#### DEBUG LOGGING:
#logwriter=open(options.operation_folder+"luci_log.txt",'wb')





##########################################NEWNEWNEWNEWNEWNEWNEWNEW
#def combineFileNameAndScanLuci(x):
#    return str(x['specId'].split(".")[0]+".mzML."+x['specId'].split(".")[1])
def combineFileNameAndScanPerco(x):
    return str(x['file']+"."+str(x['scan']))

newlist=[]
for eachgroup in groups_output:
    newlist.append(groups_output[eachgroup]) 
combined_groups=pandas.concat(newlist)#This is the combined percolator results...
#combined_fname_scans=combined_groups['mzML']+"."+str(combined_groups['scan'])
print "Generating unique names..."
combined_fname_scans=combined_groups.apply(combineFileNameAndScanPerco,axis=1)
combined_groups['unique_name']=combined_fname_scans
del newlist

newlist=[]
for eachgroup in group_to_results:
    newlist.append(group_to_results[eachgroup])

combined_luci_groups=pandas.concat(newlist)#This is the combined luciphor results...
#combined_fname_scans_luci=combined_luci_groups.apply(combineFileNameAndScanLuci,axis=1)
combined_fname_scans_luci=combined_luci_groups['mzML']+"."+combined_luci_groups['scan'].astype(str)
#print combined_fname_scans_luci
combined_luci_groups['unique_name']=combined_fname_scans_luci
combined_luci_groups.rename(columns={"sequence":"unmodified sequence","sequenceProbability":"pRS_peptideLocalizationProbability","phosphoRSscore":"pRS_score","siteProbabilities":"pRS_siteProbabilities"},inplace=True)

del newlist

combined_luci_groups=combined_luci_groups.drop('scan',1)
combined_luci_groups=combined_luci_groups.drop('mzML',1)


if 'unmodified sequence' in combined_luci_groups:
    combined_luci_groups=combined_luci_groups.drop('unmodified sequence',1)  ### IF WE ERROR ON UNMOD SEQUENCE OR UNMOD SEQUENCE IS BLANK, WE SHOULD DROP IT FROM COMBINED_GROUPS AND NOT THIS DF




merged_df=pandas.merge(combined_groups,combined_luci_groups, on="unique_name",sort=True,how="outer",indicator=False)
merged_df['numModSites'].fillna(0,inplace=True)
merged_df['numModSites']=merged_df['numModSites'].astype(int)
merged_df['pRS_sequence']=""
merged_df['pRS_numPPS']=0

#merged_df['numRPS'].fillna("0")
merged_df=merged_df.drop_duplicates(subset='unique_name')
#print merged_df,"THIS IS THE MERGED DATAFRAME."
#merge_numPPS_mask=merged_df[merged_df['numPPS']>=1]

merged_df['pRS_sequence']=merged_df['sequence']
#filtered_frame=merged_df[merged_df['sequence'].str.contains(']')]
#print filtered_frame,"filtered......!!!!!!!!!!"

print "Replacing mod positions from pRS outputs... In Parallel!"
#################### THIS WHOLE LOOP COULD PROBABLY BE PARALLELIZED WITH A GROUPBY....
#print "length before...",str(len(merged_df.index))
grouped_df=merged_df.groupby(np.arange(len(merged_df))//cpu_count())
merged_df=applyParallel(grouped_df,fix_modification_parallel)
#print "length after...",str(len(merged_df.index))


merged_df=merged_df.drop('unique_name',1)


for each_idx in set(merged_df['file_idx']):
    modified_output[each_idx]=merged_df[merged_df['file_idx']==str(each_idx)]



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
    
print "-----------------------------------------------------------------------"

