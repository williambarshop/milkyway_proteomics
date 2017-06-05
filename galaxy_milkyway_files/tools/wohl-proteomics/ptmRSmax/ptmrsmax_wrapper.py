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
import uniprot as uni
from pyteomics import mass
from pyopenms import *
from tqdm import *
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import cpu_count
import modification
#from pyopenms.sysinfo import free_mem

import time


#####################################
#This is the wrapper for ptmRSMax... Adapted from the Luciphor2 wrapper, and the phosphoRS wrapper (RSmax).
#It will be responsible for taking command line arguments and building a
#input file for running of ptmRSMax to help localize PTMs on
#identified modified peptides, and ultimately producing an estimate of
#confidence of PTM localization for the peptides
#
#VERSION 0.4A
version="0.4A"
#DATE: 04/07/2017
date="04/07/2017"
#ptmRSmax XML configuration file should be named inputFileBase+"_ptmRS_config.xml"

#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the ptmRSmax wrapper for Galaxy, Wohlschlegel Lab UCLA"
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
    #print dfGrouped,type(dfGrouped)
    #retLst = [x for x in [group.apply(fix_modification_parallel) for name, group in dfGrouped]]
    retLst = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(func)(group) for name, group in dfGrouped)
    #retLst = Parallel(n_jobs=8)(delayed(func)(group) for name, group in dfGrouped)
    #print retLst[0].columns,"before concat..."

    my_tmp=pandas.concat(retLst)
    return_df = my_tmp.reindex_axis(retLst[0].columns, axis=1)
    #print return_df.columns,"after concat..."
    return return_df


def fix_modification_parallel(group):

    #print group.columns,"in parallel at start..."
    #global numMod_list
    #global siteProbabilities_list
    filtered_modification_list=[]
    all_target_aa=[]
    mod_specific_aa={}
    first_time_setup=True
    for eachmod in modification_list:
        all_target_aa.extend(eachmod.getAAs())
        mod_specific_aa[eachmod.getName()]=eachmod.getAAs()


    for index,eachrow in group.iterrows():
        if ((eachrow['ptmRS_numMods'] < 1) and (']' in eachrow['sequence'])):
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
            group.loc[index,'ptmRS_sequence']=''.join(newPep)
            #print group.loc[index,'pRS_sequence'],"AFTER - 1"

        elif eachrow['ptmRS_numMods'] >= 1:
            #print group.loc[index,'pRS_sequence'],"BEFORE - 2"
            mod_specific_PPS={}
            total_numPPS=0
            #print "----------------------------------"
            for each_aa in all_target_aa:
                #print "looking for ",each_aa
                total_numPPS+=eachrow['sequence'].count(each_aa)
                #print "found ",str(eachrow['sequence'].count(each_aa))
            for each_mod in mod_specific_aa:
                for each_aa in mod_specific_aa[each_mod]:
                    if not each_mod in mod_specific_PPS:
                        mod_specific_PPS[each_mod]=0
                    mod_specific_PPS[each_mod]+=eachrow['sequence'].count(each_aa)

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
                    #print buffer,buffer[3:-1]
                    if buffer[3:-1] in modMass_targets: # if this mass is in our list of target mods, then we'll go ahead and drop it and only add back the AA.
                        #print "IT IS IN MOD MASS TO TARGETS..",buffer[3:-1]
                        buffer=buffer[:1]
                    else:
                        #print "IT IS NOT IN MOD MASS TO TARGETS..",buffer[3:-1],type(buffer[3:-1])
                        #for each in modMass_targets:
                        #    print each,type(each),"in modmass targets..."
                        buffer=buffer[:-1]+"]" #We'll cram in the full mass accuracy mod as given in the settings...
                    newPep.append(buffer)
                zz+=1
            #print newPep,"should have tmods removed..."
            #We have stripped out the target modifications.
            #now we will read the ptmRS modification string (FOR EACH TARGET MOD) and determine which positions to place the current target modification
            
            #print eachrow,eachrow.index 
            #print eachrow.index.tolist()
            if first_time_setup:
                #mod_set=[]
                index_localprob=eachrow.index.tolist().index('ptmRS_peptideLocalizationProbability')
                #print "index is ...",index_localprob
                remaining_list=eachrow.index.tolist()[index_localprob+1:-3]
                mod_set=set([x.split("_")[0] for x in remaining_list])
                #for each in remaining_list:
                #print "MOD SET IS",mod_set
                #print remaining_list
                for eachmod in modification_list:
                    if eachmod.getName() in mod_set:
                        filtered_modification_list.append(eachmod)
                        #######modmap_per_run[run_dict[int(each_idx)]['Original File Name']+".mzML"][eachrow['scan']]=modmap  # THIS IS NOW AVAILABLE...

                        #### Check ModificationManagement getScoredEquivalent() whatever it's called... maybe that's the key to
                        #### intelligently knowing which mods are grouped together and recounting nummods correctly...






                first_time_setup=False
            #print eachrow,"yo"
            #print filtered_modification_list,"this is the filtered list..."

            #for eachMod in modification_list:#this should be the target mods only...
            for eachMod in filtered_modification_list:#this should be the target mods only...
                #print "I'm in - 1"
                this_mod_name=eachMod.getName()
                this_mod_number_of_sites=int(eachrow[this_mod_name+"_numMods"])
                #if isinstance(eachrow[this_mod_name+"_siteProbabilities"],float):
                #    print ">>>>>>>>>>>>>>>>>>>>>"
                #    print eachrow
                #    print "<<<<<<<<<<<<<<<<<<<<<~~~~~~~~~~~~~~~~"
                this_mod_site_dict={}
                #print this_mod_number_of_sites,"this is the number of sites...."
                #print this_mod_name,"this is the name...."
                if this_mod_number_of_sites > 0:
                    this_mod_site_scores_split=eachrow[this_mod_name+"_siteProbabilities"].split(';')
                    #print this_mod_site_scores_split,"split scores..."

                    for eachsite in this_mod_site_scores_split:
                        #print eachsite,"eachsite"
                        score=float(eachsite.split()[1])
                        factor=int(eachsite.split()[0].rsplit('x',1)[1].replace(":",""))
                        site_index=str(int(eachsite.split()[0].split('(')[1].split(')')[0])-1)+"_"+str(factor)
                        this_mod_site_dict[site_index]=score
                    this_mod_sorted_sites = sorted(this_mod_site_dict.items(), key=operator.itemgetter(1), reverse=True)
                    this_mod_sorted_sites=this_mod_sorted_sites[:this_mod_number_of_sites]
                    #print this_mod_sorted_sites,"sorted sites, filtered... before"
                    this_mod_sorted_sites=[(x[0].split("_")[0],x[1]) for x in this_mod_sorted_sites]                    
                    #print this_mod_sorted_sites,"sorted sites, filtered... after"

                    collapse_dict={}
                    for eachmod in this_mod_sorted_sites:
                        #print eachmod,"mod sorted sites eachmod"
                        if eachmod[0] in collapse_dict:
                            collapse_dict[eachmod[0]]+=float(eachMod.getMass())
                        else:
                           collapse_dict[eachmod[0]]=float(eachMod.getMass())
                    #print collapse_dict,"collapse dict"
                    #I NEED TO COLLAPSE "POLYMERIC MODS" DOWN IF THEY ARE HIGH ENOUGH IN THE RANKING TO BE PLACED ON AMINOACIDS.

                    for each_mod in collapse_dict:
                    
                        #newPep[int(this_mod_sorted_sites[ns][0])]+="[+"+str(eachMod.getMass()*this_mod_factor_dict[int(this_mod_sorted_sites[ns][0])])+"]"
                        if collapse_dict[each_mod]>=0:
                            newPep[int(each_mod)]+="[+"+str(collapse_dict[each_mod])+"]"
                        else:
                            newPep[int(each_mod)]+="["+str(collapse_dict[each_mod])+"]"
                        #print "moving right along..",total_numPPS
                        #print sorted_sites[ns],"this will be the ",ns,"th site to add a target mod!"
                        #ns+=this_mod_factor_dict[int(this_mod_sorted_sites[ns][0])]
                        #ns+=1  #CHANGE THIS TO BE THE FACTOR
            #print "replacing "+eachrow['sequence']+" with "+''.join(newPep)
            group.loc[index,'ptmRS_totalNumPPS']=total_numPPS
            #print "placed totalnumpps",total_numPPS
            group.loc[index,'ptmRS_sequence']=''.join(newPep)
            #print "did the sequence, too"
            #print group.loc[index,'ptmRS_sequence'],"AFTER - 2"
   #print group.columns,"in parallel at end..."

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
parser.add_option("--ms2tol",action="store", type="float", dest="ms2tol")
parser.add_option("--scored",action="store",type="string",dest="scored")
parser.add_option("--unscored",action="store",type="string",dest="unscored")
#parser.add_option("--target",action="store",type="string",dest="target_mods")
#parser.add_option("--targetnames",action="store",type="string",dest="target_mod_names")
#parser.add_option("--variable",action="store",type="string",dest="variable_mods")
parser.add_option("--silac",action="store",type="string",dest="silac_mods")
#parser.add_option("--di",action="store",type="string",dest="diagnostic_ions")
#parser.add_option("--nl",action="store",type="string",dest="neutral_loss")
parser.add_option("--nt",action="store",type="int",dest="num_threads")
parser.add_option("--mzml",action="store",type="string",dest="mzml_files")
parser.add_option("--expgroups",action="store",type="string",dest="exp_group_file")
parser.add_option("--mp",action="store",type="int",dest="max_permutations")
parser.add_option("--scorenl",action="store_true",dest="score_nl")
parser.add_option("--filter",action="store", type="float", dest="filter")
parser.add_option("--xml",action="store",type="string",dest="xml_file")
#"--scorenl"

(options,args) = parser.parse_args()



aa_masses = {'A' : 71.037114, 'R' : 156.101111, 'N' : 114.042927, 'D' : 115.026943, 'C' : 103.009185, 'E' : 129.042593, 'Q' : 128.058578, 'G' : 57.021464, 'H' : 137.058912, 'I' : 113.084064, 'L' : 113.084064, 'K' : 128.094963, 'M' : 131.040485, 'F' : 147.068414, 'P' : 97.052764, 'S' : 87.032028, 'T' : 101.047679, 'U' : 150.95363, 'W' : 186.079313, 'Y' : 163.06332, 'V' : 99.068414 }
######### I scraped these AA masses from Matrix, who proudly bring us Mascot!
global modification_list#To copy into Parallel instances by joblib...
#global modification_list_scored#To copy into Parallel instances by joblib...
modification_list=[]
#modification_list_scored=[]

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
    print "I need to know how to group the experiments to properly build my ptmRS input files... Gimmie the file!"
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
ptmrs_config_writers={}

if options.xml_file:
    xml_locations=[]


for eachgroup in set(group_information['Fractionation Group ID String']):
    mywriter=open(options.operation_folder+eachgroup+".pin_out/crux-output/"+eachgroup+"_rsmax.txt",'w')
    if not options.xml_file:
        ptmrswriter=open(options.operation_folder+eachgroup+".pin_out/crux-output/"+eachgroup+"_ptmRS_config.xml",'w')
        ptmrs_config_writers[eachgroup]=ptmrswriter
    else:
        xml_locations.append(options.operation_folder+eachgroup+".pin_out/crux-output/"+eachgroup+"_ptmRS_config.xml")
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
for thismod in options.scored.split("___"):
    for eachmod in thismod.split(",")[0]:
        if eachmod is not "":
            try:
                mod_mass=float(eachmod)
            except:
                mod_mass=round(mass.calculate_mass(formula=eachmod),6)
            target_masses.append(mod_mass)
            #target_masses.append(non_decimal.sub('',eachmod))
target_masses=list(set(target_masses))



#modId,ScoredmodMassDouble,NLmassDouble,ScoredModificationName,ScoredShortModName,ScoredModSiteAAs,ScoredElementalModificationStr
#We'll take in the info about the target mod...
target_masses_list=[]
mods_by_mass={}
terminal_aa_re=re.compile(r'([\D]+$)')
for each_mod in options.scored.split("___"):
    this_mod_split=each_mod.split(",")
    if len(this_mod_split)<4:
        continue #We'll just skip it... Wrong formatting.
    newMod=modification.Modification()
    if len(this_mod_split)>4 and len(this_mod_split[5])>0:
        #print this_mod_split[5],"THIS IS THE NL STRING...<--------------"
        NL_mass=-1.0*float(non_decimal.sub('',this_mod_split[5]))
        NL_AAs=terminal_aa_re.search(this_mod_split[5]).groups(0)[0]
        #print NL_AAs,"These are the NL AAs... and all of them"
        #print "raw string...",this_mod_split[5]
        NL_string=str(NL_AAs)+" "+str(NL_mass)
    else:
        NL_string="-"
    try:
        mod_mass=float(this_mod_split[0])
    except:
        mod_mass=round(mass.calculate_mass(formula=this_mod_split[0]),6)
        pass #Handle conversion of elemental composition to mass... pyteomics?
    print this_mod_split
    if str(mod_mass) not in mods_by_mass:
        newMod.create(this_mod_split[4],this_mod_split[1]+" "+str(mod_mass),NL_string)#nameString,modString,NLstring
        mods_by_mass[str(mod_mass)]=newMod
        target_masses_list.append(str(newMod.getMass()))
        modification_list.append(newMod)
    else:
        update_me=mods_by_mass[str(mod_mass)]
        mass=update_me.getMass()
        aas=update_me.getAAs()
        NLs=update_me.getNL()
        #print NLs,"these are NLs..."
        for each_aa in this_mod_split[1]:
            if each_aa not in aas:
                aas.append(each_aa)
        update_me.setAAs(aas)
        if len(this_mod_split)>4 and len(this_mod_split[5])>0 and this_mod_split[5]!="-":
            NL_mass=-1.0*float(non_decimal.sub('',this_mod_split[5]))
            NL_AAs=terminal_aa_re.search(this_mod_split[5]).groups(0)[0]

            for each_AA in NL_AAs:
                if each_AA not in NLs:
                    NLs[each_AA]=NL_mass
        update_me.setNL(NLs)
        #print aas, mass,NLs,"after merge..."
        #pass#Update this with more AAs, etc



#target_mod_names_split = [x for x in options.target_mod_names.split(",")] # Name Strings
#target_mod_split = [x for x in options.target_mods.split(",")] # Taret Mod String
#target_mod_NL_split = [x for x in options.neutral_loss.split(",")] # NL Strings
##target_mod_DI_split = [x for x in options.diagnostic_ions.split(",")] # DIstring "diagnostic_ions"

##zipped_target_mods = zip(target_mod_names_split,target_mod_split,target_mod_NL_split)#,target_mod_DI_split)
print "We'll score ",str(len(modification_list)),"mods (ptmRS)"
for each_mod in modification_list:
    print each_mod.mass
    print each_mod.AAs
    print "------------"
#target_masses_list=[]
#for each_mod in zipped_target_mods:
#    newMod=modification.Modification()
#    newMod.create(each_mod[0],each_mod[1],each_mod[2])#,each_mod[3])
#    target_masses_list.append(str(newMod.getMass()))
#    modification_list.append(newMod)
#    #modification_list_scored.append(newMod)
##print modification_list
##sys.exit(2)


#target_mod_split = [x.split() for x in options.target_mods.split(",")] ###############################################<-~~~~~~~~~~~~~~~~~~~~~~~
#target_mod_AA = target_mod_split[0]
#target_mod_mass = target_mod_split[1]
#target_mod_elements = target_mod_split[2]

#print options.neutral_loss,type(options.neutral_loss),"yo"
#if not options.neutral_loss is "":
#    target_mod_NL_split = [x.split() for x in options.neutral_loss.split(",")]
#else:
#    target_mod_NL_split = ["A 0.00 H".split()] #We're not using a NL...

#target_mod_NL_AA = target_mod_NL_split[0]
#target_mod_NL_mass = str(math.fabs(float(target_mod_NL_split[1])))
#target_mod_NL_elements = target_mod_NL_split[2]





########### THIS IS WHERE YOU LEFT OFF ON 4/7/2017 --- CONTINUE ON WITH ADAPTING THE UNSCORED MODS TO THE NEW SYSTEM! YOU'RE ALMOST THERE.
unscored_mod_dict={}#key is mass(as str), value is AAs (as str)

if options.unscored!="None":
    unscored_mod_split = options.unscored.split("___")
    for each_mod in unscored_mod_split:
        this_mod_split=each_mod.split(",")
        print this_mod_split
        if len(this_mod_split)<4:
            continue
        try:
            mod_mass=float(this_mod_split[0])
        except:
            mod_mass=round(mass.calculate_mass(formula=this_mod_split[0]),6)
        mod_aas=this_mod_split[1]
        unscored_mod_dict[mod_mass]=mod_aas
    

#modId,ScoredmodMassDouble,NLmassDouble,ScoredModificationName,ScoredShortModName,ScoredModSiteAAs,ScoredElementalModificationStr


#We're going to make a reference dictionary to convert the mass (as str) to the modid... which we'll go ahead and assign now!
id_itr=1
modID_targets=[] # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ptmRS this will hold a list of the score target mods...
global modMass_targets
modMass_targets=[]
massStr_to_modID={}
for each_mod in modification_list:
    massStr_to_modID[str(each_mod.getMass())]=str(id_itr)
    modID_targets.append(str(id_itr))
    modMass_targets.append(str(each_mod.getMass()))
    each_mod.setID(id_itr)
    id_itr+=1
#print modMass_targets
#print massStr_to_modID
#for each_mod in modification_list:
#    print each_mod,each_mod.getID()
#sys.exit(2)


target_mods_less_than_this=id_itr

#Now, we'll go ahead and iterate over all the other mods...
#id_itr=2
for each_mod in unscored_mod_dict:
    massStr_to_modID[str(each_mod)]=str(id_itr)
    id_itr+=1
unscored_mod_str=""
for each_mod in unscored_mod_dict:
    unscored_mod_str+=massStr_to_modID[str(each_mod)]+","+str(each_mod)+","+"0.0"+","+"unscored"+massStr_to_modID[str(each_mod)]+","+"mod"+massStr_to_modID[str(each_mod)]+","+unscored_mod_dict[each_mod]


for eachrun in run_to_group:
    #print eachrun,"we're dealing with this run..."
    ##thislink=os.readlink(eachrun)
    ##shutil.copy(thislink,options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)
    shutil.copy(eachrun,options.operation_folder+run_to_group[eachrun]+".pin_out/crux-output/"+eachrun)




#global target_mass_conversion
#target_mass_conversion={}


#target_mass_conversion[str(round(float(target_mod_mass),2))]=str(target_mod_mass)
#We're also going to actually put the rest of our mod masses in there, too...

#for eachmod in options.variable_mods.split(","):
#    mod_aa=eachmod.split()[0]
#    mod_mass=eachmod.split()[1]
#    target_mass_conversion[str(mod_mass)]=str(mod_mass)



#if options.variable_mods is not None:
#    for each in unscored_mod_dict.keys():
#        if each is not "":
#            target_mass_conversion[str(non_decimal.sub('',each))]=str(non_decimal.sub('',each))



#rounded_target_masses=str(target_mod_mass)
#print "These are the rounded targets...",rounded_target_masses
global modmap_per_run
modmap_per_run={}

#print target_mass_conversion,"this is target mass conv."
regex = re.compile('[^a-zA-Z]')
for each_idx in group_information['Crux File Integer']:

    file = MzMLFile()
    experiment = MSExperiment()
    print "loading file "
    print run_dict[int(each_idx)]['Original File Name']+".mzML"#run_dict[int(each_idx)]
    print ""
    file.load(os.path.join(basedir,run_dict[int(each_idx)]['Original File Name']+".mzML"), experiment)
    if not run_dict[int(each_idx)]['Original File Name']+".mzML" in modmap_per_run:
        modmap_per_run[run_dict[int(each_idx)]['Original File Name']+".mzML"]={}
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
    mask2=mask[(mask['sequence'].str.contains('|'.join(target_masses_list)))] #numpy.any(mods in mask.sequence for mods in target_masses)]
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
                mod_end_pos=i#changed these from actualindex...
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
                modmap+=massStr_to_modID[mod[2:-1]]
                if massStr_to_modID[mod[2:-1]] in modID_targets:
                    target_modsites.append(str(massStr_to_modID[mod[2:-1]])+":"+str(len(modmap)-3))
                    #target_modsites.append(str(actualindex))
            i+=1
            actualindex+=1 # THIS WAS ADDED 12/15/2015 TO CORRECT FOR INDEX SHIFTING
        modmap+=".0"
        modmap_per_run[run_dict[int(each_idx)]['Original File Name']+".mzML"][eachrow['scan']]=modmap
        #print modmap,"This is modmap!"
        mod_positions=",".join(target_modsites)
        specID=eachrow['file'].replace('.mzML','')+"."+str(eachrow['scan'])+'.'+str(eachrow['scan'])+'.'+str(eachrow['charge'])
        #print specID,"this is specID"
        spec_peaks=experiment[eachrow['scan']-1].get_peaks()#THIS WILL HAVE TO BE FORMATTED CORRECTLY....
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
        elif 13 in activation_type: ################################################################ WHAT IS THE ETHCD ACTIVATION TYPE INT?
            activation_method="EThcD"
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


os.chdir(basedir)

for eachwriter in luci_input_writers:
    luci_input_writers[eachwriter].close()

if options.score_nl:
    score_neutrals="true"
else:
    score_neutrals="false"


processes=[]
#RSmax.exe /path/to/inputs/ scoreNLBoolean massToleranceDouble modId,ScoredmodMassDouble,NLmassDouble,ScoredModificationName,ScoredShortModName,ScoredModSiteAAs,ScoredElementalModificationStr
#modId,ScoredmodMassDouble,NLmassDouble,NL_aa,ScoredModificationName,ScoredShortModName,ScoredModSiteAAs
target_modstr=""
for eachMod in modification_list:
    NL_aas=""
    NL_mass=""
    AAstr=""
    for eachNL in eachMod.getNL():
        NL_aas+=eachNL
        NL_mass=str(eachMod.getNL()[eachNL])
    for eachAA in eachMod.getAAs():
        AAstr+=eachAA
    target_modstr+=eachMod.getID()+","+str(eachMod.getMass())+","+NL_mass+","+NL_aas+","+eachMod.getName()+","+eachMod.getName()+","+AAstr+" "
#print target_modstr
#sys.exit(2)
#target_modstr=massStr_to_modID[str(target_mod_mass)]+","+str(target_mod_mass)+","+target_mod_NL_mass+","+"targetMod"+","+"targ"+","+target_mod_AA+","+target_mod_elements
####COME BACK
if not options.xml_file:
    for eachwriter in ptmrs_config_writers:
        thiswriter=ptmrs_config_writers[eachwriter]

        thiswriter.write("<AnyPTM>\n")
        for eachMod in modification_list:
            thiswriter.write("<modification name=\"{0}\" abbreviation=\"{1}\" searchdefined=\"TRUE\" mass=\"{2}\" unimodId=\"{3}\" >\n".format(eachMod.getName(),eachMod.getName(),str(eachMod.getMass()),eachMod.getID()))
            for eachAA in eachMod.getAAs():
                thiswriter.write("<target aminoacid=\"{0}\" />\n".format(eachAA))

            #for eachNL in eachMod.getNL(): #float is value, key is AA
            nl_dict=eachMod.getNL()
            #print nl_dict

            if isinstance(nl_dict,dict):
                #print "i'm in!"
                if len(nl_dict.keys())>0:
                    theAAs=nl_dict.keys()
                    for eachAA in theAAs:
                        nl_mass=nl_dict[eachAA]
                        if nl_mass > 0.0:
                            thiswriter.write("<neutralloss abbreviation=\"{0}\" mass=\"{1}\">\n".format(eachMod.getName()+"_nl_"+eachAA,nl_mass))
                            #for eachAA in theAAs:
                            thiswriter.write("<target aminoacid=\"{0}\" />\n".format(eachAA))
                            #thiswriter.write("<target aminoacid=\"{0}\" factor=\"1\"/>\n".format(eachAA))
                            #thiswriter.write("<target aminoacid=\"{0}\" factor=\"2\"/>\n".format(eachAA)) 
                            #thiswriter.write("<target aminoacid=\"{0}\" factor=\"3\"/>\n".format(eachAA))
                            thiswriter.write("</neutralloss>\n")
                
            
                #for eachAA in nl_dict.keys():
                 
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            thiswriter.write("<equivalentmodification name=\"{0}\" factor=\"1\" new=\"FALSE\"/>\n".format(eachMod.getName()))
            thiswriter.write("</modification>\n")

        thiswriter.write("<FragmentIonCompositionPreference>\n<FragmentIonComposition ActivationType=\"CID\" FragmentIonComposition=\"b,y\" NeutralLossFragmentIonComposition=\"\"/>\n<FragmentIonComposition ActivationType=\"HCD\" FragmentIonComposition=\"b,y\" NeutralLossFragmentIonComposition=\"b,y\"/>\n<FragmentIonComposition ActivationType=\"ETHcD\" FragmentIonComposition=\"b,y,c,zPrime,zRadical\" NeutralLossFragmentIonComposition=\"b,y\"/>\n<FragmentIonComposition ActivationType=\"ETD\" FragmentIonComposition=\"c,zPrime,zRadical\" NeutralLossFragmentIonComposition=\"\"/>\n</FragmentIonCompositionPreference>\n")



        thiswriter.write("</AnyPTM>\n")

    for eachwriter in ptmrs_config_writers:
        ptmrs_config_writers[eachwriter].close()


else:
    for each in xml_locations:
        shutil.copy(options.xml_file,each)






with open(os.devnull, "w") as fnull:
    for eachgroup in tqdm(set(group_information['Fractionation Group ID String'])):
        os.chdir(basedir)
        os.chdir(options.operation_folder+eachgroup+".pin_out/crux-output/")
        #processes.append(subprocess.Popen(command,shell=True)
        command = "mono /galaxy-central/tools/wohl-proteomics/ptmRSmax/ptmRSmax.exe "+os.getcwd()+" "+eachgroup+" "+str(options.num_threads)+" "+str(options.max_permutations)+" "+score_neutrals+" "+str(options.ms2tol)+" "+str(target_mods_less_than_this-1)+" "+target_modstr+" "+unscored_mod_str
        #command = "mono /home/galaxy/wohl-proteomics-backups/2016-07-25/wohl-proteomics/ptmRSmax/ptmRSmax.exe "+os.getcwd()+" "+eachgroup+" "+str(options.num_threads)+" "+str(options.max_permutations)+" "+score_neutrals+" "+str(options.ms2tol)+" "+str(target_mods_less_than_this-1)+" "+target_modstr+" "+unscored_mod_str
        print "running command... ",command
        #processes.append(subprocess.Popen(command.split()))
        #proc=subprocess.Popen(command.split())
        #proc.wait()
        ############subprocess.call(command.split(), stderr=sys.stdout.fileno())  #### SEEMED TO WORK ONCE UP ON A TIME....
        #processes.append(subprocess.Popen(command.split(),stdout=fnull, stderr=fnull)) # These pipes will hide the output from luciphor... this is only for DEBUG DEBUG DEBUG purposes...
        #out, err = p.communicate()
        p=os.system(command)
        #print "Luciphor exited with a status of ",str(os.system(command))
        time.sleep(1)
        print "Finished scoring this command...!"
        time.sleep(1)

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
    try:
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
    except:
        print "WARNING ::::: Couldn't find results for ",str(eachgroup)

   
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

#modsite_to_modmass={}
#modsite_to_modmass_fixed={}
#modsite_to_modmass_variable={}
#print options.target_mods,"target mods!"
#print options.variable_mods,"unscored mods!"

#for eachmod in options.target_mods.split():
#    modsite_to_modmass[eachmod[:1]]=non_decimal.sub('',eachmod)

#tmod_aa=options.target_mods.split()[0]
#tmod_mass=options.target_mods.split()[1]
#for each_aa in tmod_aa:
#    modsite_to_modmass[each_aa]=tmod_mass


#for eachmod in options.variable_mods.split(","):
#    mod_aa=eachmod.split()[0]
#    mod_mass=eachmod.split()[1]
#    for each_aa in mod_aa:
#        modsite_to_modmass[each_aa]=mod_mass
#print "This is the modsite dict!",modsite_to_modmass


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
    try:
        newlist.append(group_to_results[eachgroup])
    except:
        print "We didn't find a result file from ",eachgroup

combined_luci_groups=pandas.concat(newlist)#This is the combined results
#combined_fname_scans_luci=combined_luci_groups.apply(combineFileNameAndScanLuci,axis=1)
combined_fname_scans_luci=combined_luci_groups['mzML']+"."+combined_luci_groups['scan'].astype(str)
#print combined_fname_scans_luci
combined_luci_groups['unique_name']=combined_fname_scans_luci
print "column names are"
combined_luci_groups.rename(columns={"sequence":"unmodified sequence","sequenceProbability":"ptmRS_peptideLocalizationProbability","ptmRSscore":"ptmRS_score","all_numMods":"ptmRS_numMods"},inplace=True)
print combined_luci_groups.columns.values.tolist()

del newlist

combined_luci_groups=combined_luci_groups.drop('scan',1)
combined_luci_groups=combined_luci_groups.drop('mzML',1)


#These are global so that they will make their way into the joblib Parallel instances...
global numMod_list
global siteProbabilities_list

numMod_list=[]
siteProbabilities_list=[]
for each in combined_luci_groups.columns.values.tolist():
    if "_numMods" in each:
        numMod_list.append(each)
    elif "_siteProbabilities" in each:
        siteProbabilities_list.append(each)



print combined_luci_groups.columns,"These are the ptmRS group columns..."
merged_df=pandas.merge(combined_groups,combined_luci_groups, on="unique_name",sort=True,how="outer",indicator=False)  #CHANGED BACK TO OUTER BECAUSE WE WERE THROWING AWAY UNMOD PEPTIDES! (3/2/16)
#merged_df=combined_groups.merge(combined_luci_groups, on="unique_name",sort=True,how="inner",indicator=False)   #CHANGED TO INNER 2/22/2016
#print merged_df.columns.values.tolist(),"fwd merged cols"
#print merged_df.columns,"columns, too... jsut after..."
#firestick
#reverse_merged_df=pandas.merge(combined_luci_groups,combined_groups, on="unique_name",sort=True,how="outer",indicator=False)
#print reverse_merged_df.columns.values.tolist(),"rev merged cols"

merged_df['ptmRS_numMods'].fillna(0,inplace=True)
merged_df['ptmRS_numMods']=merged_df['ptmRS_numMods'].astype(int)
merged_df['ptmRS_sequence']=""
merged_df['ptmRS_numPPS']=0

#merged_df['numRPS'].fillna("0")
merged_df=merged_df.drop_duplicates(subset='unique_name')
#merge_numPPS_mask=merged_df[merged_df['numPPS']>=1]
#print target_mass_conversion,"this is the mass conv dict"

merged_df['ptmRS_sequence']=merged_df['sequence']
#filtered_frame=merged_df[merged_df['sequence'].str.contains(']')]
#print filtered_frame,"filtered......!!!!!!!!!!"

print "Replacing mod positions from ptmRS outputs... In Parallel!"
#################### THIS WHOLE LOOP COULD PROBABLY BE PARALLELIZED WITH A GROUPBY....
#print "length before...",str(len(merged_df.index))

grouped_df=merged_df.groupby(np.arange(len(merged_df))//cpu_count())
####print merged_df,type(merged_df)
####merged_df = merged_df.apply(fix_modification_parallel,axis=0)
#print merged_df.columns,"columns, too... BEFORE PARALLEL"
merged_df['ptmRS_totalNumPPS']=np.nan
merged_df=applyParallel(grouped_df,fix_modification_parallel)
print merged_df.columns,"columns, too... AFTER PARALLEL"
#print "length after...",str(len(merged_df.index))

#print merged_df
merged_df=merged_df.drop('unique_name',1)



def stripMods(x):
    return re.sub("[\(\[].*?[\)\]]", "", x['sequence'])
#grouped_df2=merged_df.groupby(np.arange(len(merged_df))//cpu_count())
merged_df['unmodified sequence']=merged_df.apply(stripMods,axis=1)
#merged_df=applyParallel(grouped_df2,stripMods)

if "unmodified sequence_x" in merged_df.columns.tolist():
    #merged_df['unmodified sequence']=merged_df['unmodified sequence_x']+merged_df['unmodified sequence_y']
    merged_df=merged_df.drop('unmodified sequence_x',1)
    merged_df=merged_df.drop('unmodified sequence_y',1)
    #merged_df.rename(columns={"unmodified sequence_y":"unmodified sequence"},inplace=True)
    #merged_df['unmodified sequence']=merged_df['sequence']



#print merged_df,"THIS IS THE MERGED DATAFRAME AFTER PROCESSING....."
print merged_df.columns,"columns, too... END"

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

