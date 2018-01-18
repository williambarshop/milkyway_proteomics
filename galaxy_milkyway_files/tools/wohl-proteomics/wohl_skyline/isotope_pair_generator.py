import os, sys, re, time
import collections
import optparse
from copy import deepcopy
from xml.etree import ElementTree
import xmltodict


#####################################
#This is the wrapper for Paired Isotope Addition from two variable mods on a skyline file.
#
#VERSION 0.50A
version="0.50A"
#DATE: 1/17/2018
date="1/17/2018"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the Skyline Isotope Pair Corrector, Wohlschlegel Lab UCLA"
print "Written by William Barshop"
print "Version: ",version
print "Date: ",date




parser = optparse.OptionParser()

parser.add_option("--skylineFile",action="store",type="string",dest="skylineFile")
parser.add_option("--outputFile",action="store",type="string",dest="outputFile")
parser.add_option("--lightMass",action="store", type="float", dest="lightMass")
parser.add_option("--heavyMass",action="store", type="float", dest="heavyMass")
parser.add_option("--targetAA",action="store", type="string", dest="target_AA")
parser.add_option("--removeUnmod",action="store_true", dest="removeUnmod")


(options,args) = parser.parse_args(sys.argv)


global proton_mass
proton_mass=1.007825035
global C13_mass
C13_mass=13.00335483-12.00#it's the isotope delta



def correct_Isotope_Partners(inputProtein,lightMass,heavyMass,target_AA,removeUnmod):
    #This method takes in a protein and modification information-- copies the protein, finds unpaired heavy and light mods, and adds the pairs.
    #
    #NOTE: Please make sure to handle situations where one modification state has different charge states identified than the other!
    #
    this_new_protein=deepcopy(inputProtein)
    
    sequences_with_heavy={}
    sequences_with_light={}

    sequences_with_heavy_bibliospec={}
    sequences_with_light_bibliospec={}
    
    if 'dict' in str(type(this_new_protein['peptide'])) or 'Dict' in str(type(this_new_protein['peptide'])):
        this_new_protein_peptides=[this_new_protein['peptide']]
       
    else:
        this_new_protein_peptides=this_new_protein['peptide']
    if True: #shitty code workaround from old format of if else, embarrassing but lazy
        #for each_peptide in this_new_protein['peptide']:
        for each_peptide in this_new_protein_peptides:
            #print each_peptide[0],"this is each peptide[0]\n==============================="
            if str(round(lightMass,4)) in each_peptide['@modified_sequence']:
                if 'dict' in str(type(each_peptide['precursor'])) or "Dict" in str(type(each_peptide['precursor'])):
                    the_precursors=[each_peptide['precursor']]
                else:
                    the_precursors=each_peptide['precursor']
                if each_peptide['@modified_sequence'] in sequences_with_light:
                    pass
                else:
                    sequences_with_light[each_peptide['@modified_sequence']]=[]
                    sequences_with_light_bibliospec[each_peptide['@modified_sequence']]={}
                for each_precursor in the_precursors:
                    this_charge=each_precursor['@charge']
                    sequences_with_light[each_peptide['@modified_sequence']].append(this_charge)
                    if 'bibliospec_spectrum_info' in each_precursor:
                        sequences_with_light_bibliospec[each_peptide['@modified_sequence']][this_charge]=each_precursor['bibliospec_spectrum_info']
                
            if str(round(heavyMass,4)) in each_peptide['@modified_sequence']:
                if 'dict' in str(type(each_peptide['precursor'])) or "Dict" in str(type(each_peptide['precursor'])):
                    the_precursors=[each_peptide['precursor']]
                else:
                    the_precursors=each_peptide['precursor']
                if each_peptide['@modified_sequence'] in sequences_with_heavy:
                    pass
                else:
                    sequences_with_heavy[each_peptide['@modified_sequence']]=[]
                    sequences_with_heavy_bibliospec[each_peptide['@modified_sequence']]={}

                for each_precursor in the_precursors:
                    this_charge=each_precursor['@charge']
                    sequences_with_heavy[each_peptide['@modified_sequence']].append(this_charge)
                    if 'bibliospec_spectrum_info' in each_precursor:
                        sequences_with_heavy_bibliospec[each_peptide['@modified_sequence']][this_charge]=each_precursor['bibliospec_spectrum_info']

    
    
    
    peptides_with_mods=sequences_with_heavy.keys()
    peptides_with_mods.extend(sequences_with_light.keys())
    peptides_with_mods_charges={}
    for each_dict in [sequences_with_heavy,sequences_with_light]: ##### CONSIDER ALL THIS CODE FOR IF YOU ADD MULTI-CONCURRENT LABEL SUPPORT
        for each_peptide_key in each_dict:
            for each_charge in each_dict[each_peptide_key]:
                if each_peptide_key not in peptides_with_mods_charges:
                    peptides_with_mods_charges[each_peptide_key]=[]
                if each_peptide_key.replace(str(round(heavyMass,4)),str(round(lightMass,4))) not in peptides_with_mods_charges:
                    peptides_with_mods_charges[each_peptide_key.replace(str(round(heavyMass,4)),str(round(lightMass,4)))]=[]
                if each_peptide_key.replace(str(round(lightMass,4)),str(round(heavyMass,4))) not in peptides_with_mods_charges:
                    peptides_with_mods_charges[each_peptide_key.replace(str(round(lightMass,4)),str(round(heavyMass,4)))]=[]
                
                if each_charge not in peptides_with_mods_charges[each_peptide_key]:
                    peptides_with_mods_charges[each_peptide_key].append(each_charge)
                if each_charge not in peptides_with_mods_charges[each_peptide_key.replace(str(round(heavyMass,4)),str(round(lightMass,4)))]:
                    peptides_with_mods_charges[each_peptide_key.replace(str(round(heavyMass,4)),str(round(lightMass,4)))].append(each_charge)
                if each_charge not in peptides_with_mods_charges[each_peptide_key.replace(str(round(lightMass,4)),str(round(heavyMass,4)))]:
                    peptides_with_mods_charges[each_peptide_key.replace(str(round(lightMass,4)),str(round(heavyMass,4)))].append(each_charge)
                
    
    #print peptides_with_mods_charges, "These are the peptides which we will need to pair with their proper isotope partner!"
    completed_peptides=[] #This will be a list of peptides_with_mods which have already been completed.
    if len(peptides_with_mods)==0:
        if not removeUnmod:
            return inputProtein #There are no peptides to fix! We can call it a day on this protein.
        else:
            return None
    
    #So, we'll iterate over these peptides again, and stop whenever we see a peptide in our list of peptides_with_mods.
    #If we find that peptide, we'll set it up as a heavy/light pair, and iterate back over all of our peptides to grab the other state if it exists.
    #If it exists, then it's easy to set up the pair by copying it into place with the appropriate changes to make it a heavy mod.
    #We should make sure that the same charge states are available for both isotope labels!
    #Once that's done, the peptide should be removed from the list "peptides_with_mods" and added to the "completed_peptides" list to make sure that we
    #don't process it again.
    if 'dict' in str(type(this_new_protein['peptide'])) or 'Dict' in str(type(this_new_protein['peptide'])):
        #print this_new_protein['peptide'],"impending crash"
        #for each_key in this_new_protein['peptide']:
        #    each_peptide=this_new_protein['peptide'][each_key]
        #    if each_peptide['@modified_sequence'] in peptides_with_mods:
        #        print "it was a dict and a modified sequence! Do something here!!!<====================="
        #        pass #do something here because the peptide is modified
        iterate_over_me=[this_new_protein['peptide']]
    else:
        iterate_over_me=this_new_protein['peptide']
    if True:#shitty, again
        ctr=0
        new_peptide_list=[]
        #iterate_over_me=this_new_protein['peptide']
        for each_peptide in iterate_over_me:
            light=False
            heavy=False
            opposite_in_mods=False
            #print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
            #print each_peptide
            if each_peptide['@modified_sequence'] in completed_peptides or each_peptide['@modified_sequence'].replace(str(round(heavyMass,4)),str(round(lightMass,4))) in completed_peptides or each_peptide['@modified_sequence'].replace(str(round(lightMass,4)),str(round(heavyMass,4))) in completed_peptides:
                ##consider this section for multi-concurrent label support
                #print each_peptide['@modified_sequence'],each_peptide,"is supposedly already in",completed_peptides
                #print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                continue
            #print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
            modcount=0
            modcount+=each_peptide['@modified_sequence'].count(str((round(lightMass,4))))
            modcount+=each_peptide['@modified_sequence'].count(str((round(heavyMass,4))))
            if modcount==0:
                if not removeUnmod:
                    new_peptide_list.append(each_peptide)
                else:
                    print "This peptide has no modsite and will be removed(0): {0}".format(each_peptide['@modified_sequence'])
                    #print each_peptide['@modified_sequence'].count(str((round(lightMass,4))))
                    #print each_peptide['@modified_sequence'].count(str((round(heavyMass,4))))
                    #print str((round(heavyMass,4))),str((round(lightMass,4)))
            if each_peptide['@modified_sequence'] in peptides_with_mods:
                if modcount >1:
                    print "This peptide has more than one modsite and will be removed: {0}".format(each_peptide['@modified_sequence'])
                    continue
                elif modcount==0:
                    if not removeUnmod:
                        new_peptide_list.append(each_peptide)
                    else:
                        print "This peptide has no modsite and will be removed(1): {0}".format(each_peptide['@modified_sequence'])
                    continue
                #is what we're looking at heavy or light?
                if str(round(lightMass,4)) in each_peptide['@modified_sequence']:
                    light=True
                    if each_peptide['@modified_sequence'].replace(str(round(lightMass,4)),str(round(heavyMass,4))) in peptides_with_mods:
                        opposite_in_mods=True
                       
                else:
                    heavy=True #if not, must be heavy.
                    if each_peptide['@modified_sequence'].replace(str(round(heavyMass,4)),str(round(lightMass,4))) in peptides_with_mods:
                        opposite_in_mods=True

                if not heavy and not light:#it's not modified! just add it and move on
                    if not removeUnmod:
                        new_peptide_list.append(each_peptide)
                    else:
                        print "This peptide has no modsite and will be removed(2): {0}".format(each_peptide['@modified_sequence'])
                    continue
                #now that we know if we're looking at heavy or light, and if the opposite label for this mod state  of peptide is in the list, we can
                #go ahead and add the necessary info and get things combined.
                
                #If we start with heavy, we're going to replace that in a couple of places with the light version first.
                if heavy:
                    each_peptide['@calc_neutral_pep_mass']=str(float(each_peptide['@calc_neutral_pep_mass'])-(heavyMass-lightMass)) 
                    each_peptide['@modified_sequence']=each_peptide['@modified_sequence'].replace(str(round(heavyMass,4)),str(round(lightMass,4)))
                    mod_aa_indicies=[]
                    #print str(type(each_peptide['variable_modifications'])),"this is the heavy mod type"
                    #print each_peptide['variable_modifications']
                    if 'dict' in str(type(each_peptide['variable_modifications'])) or "Dict" in str(type(each_peptide['variable_modifications'])):
                        for each_modification_key in each_peptide['variable_modifications']:
                            each_modification=each_peptide['variable_modifications'][each_modification_key]
                            #print "looking for a heavy mod...."
                            if each_modification['@mass_diff']=="+"+str(round(heavyMass,1)): #We're assuming that it's positive...
                                mod_aa_indicies.append(each_modification['@index_aa'])
                                each_modification['@mass_diff']="+"+str(round(lightMass,1))
                                each_modification['@modification_name']=str(round(lightMass,1))
                                #print "found the mod... Let's change it to light!",each_modification
                            
                    heavy_mod_dict=collections.OrderedDict()
                    internal_list=[]
                    for each_aa_index in mod_aa_indicies:
                        mod_dict=collections.OrderedDict()
                        mod_dict['@index_aa']=each_aa_index
                        mod_dict['@modification_name']=str(round(heavyMass-lightMass,1))
                        mod_dict['@mass_diff']="+"+str(int(round(heavyMass-lightMass,0)))
                        internal_list.append(mod_dict)
                    heavy_mod_dict['explicit_heavy_modifications']={}
                    heavy_mod_dict['explicit_heavy_modifications']['explicit_modification']=internal_list
                    #internal_list
                    each_peptide_new=each_peptide.__class__()
                    for key,value in each_peptide.items():
                        each_peptide_new[key]=value
                        #print "the key is now",key
                        if key=='variable_modifications' and len(internal_list)>0:
                            #print "found it!"
                            each_peptide_new['explicit_modifications']=heavy_mod_dict
                            #print each_peptide_new['explicit_modifications'],"added"
                            #each_peptide_new['explicit_heavy_modifications']=internal_list
                    #each_peptide_new['explicit_modifications']=heavy_mod_dict_new
                    this_new_protein['peptide'][ctr]=each_peptide_new
                    each_peptide=this_new_protein['peptide'][ctr]
                else:
                    #it's light
                    mod_aa_indicies=[]
                    for each_modification_key in each_peptide['variable_modifications']:
                        each_modification=each_peptide['variable_modifications'][each_modification_key]
                        #print each_modification
                        if each_modification['@mass_diff']=="+"+str(round(lightMass,1)): #We're assuming that it's positive...
                            mod_aa_indicies.append(each_modification['@index_aa'])
                            each_modification['@mass_diff']="+"+str(round(lightMass,1))
                            each_modification['@modification_name']=str(round(lightMass,1))
                            #print "found the mod... Let's change it to light!"
                            
                    light_mod_dict=collections.OrderedDict()
                    internal_list=[]
                    for each_aa_index in mod_aa_indicies:
                        mod_dict=collections.OrderedDict()
                        mod_dict['@index_aa']=each_aa_index
                        mod_dict['@modification_name']=str(round(heavyMass-lightMass,1))
                        mod_dict['@mass_diff']="+"+str(int(round(heavyMass-lightMass,0)))
                        internal_list.append(mod_dict)
                    light_mod_dict['explicit_heavy_modifications']={}
                    light_mod_dict['explicit_heavy_modifications']['explicit_modification']=internal_list
                    #light_mod_dict['explicit_heavy_modifications']=internal_list
                    each_peptide_new=each_peptide.__class__()
                    for key,value in each_peptide.items():
                        each_peptide_new[key]=value
                        #print "the key is now",key
                        if key=='variable_modifications' and len(internal_list)>0:
                            #print "found it!"
                            each_peptide_new['explicit_modifications']=light_mod_dict
                            #print each_peptide_new['explicit_modifications'],"added"
                    #light_mod_dict['explicit_heavy_modifications']=internal_list
                    #each_peptide['explicit_modifications']=light_mod_dict
                    this_new_protein['peptide'][ctr]=each_peptide_new
                    each_peptide=this_new_protein['peptide'][ctr]

                    
                    
                desired_charges=peptides_with_mods_charges[each_peptide['@modified_sequence']]
                if 'dict' in str(type(each_peptide['precursor'])) or "Dict" in str(type(each_peptide['precursor'])):    
                    iterable_precursors=[each_peptide['precursor']]
                    #pass #<---------------------------------------------------------=========================================<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CHECK THIS LATER!!!
                else: #should be a dict. This seems to be the most common case.             
                    iterable_precursors=each_peptide['precursor']
                    
                example_precursor=iterable_precursors[0]
                isotope_delta_dict={}
                isotope_mass_dict={}
                transition_ctr=0
                charge=int(example_precursor['@charge'])
                for each_transition in example_precursor['transition']:
                    if "@mass_index" in each_transition:
                        isotope=int(each_transition['@mass_index'])
                    else:
                        isotope=0
                        isotope_zero_mass=example_precursor['transition'][transition_ctr]['product_mz']
                    isotope_mass_dict[str(float(isotope))]=example_precursor['transition'][transition_ctr]['product_mz']
                    transition_ctr+=1
                    
                for each_mass_index in isotope_mass_dict:
                    isotope_mass=isotope_mass_dict[each_mass_index]
                    isotope_delta_dict[each_mass_index]=(float(isotope_mass_dict[each_mass_index])-float(isotope_mass_dict["0.0"]))*float(charge)
                
                
                
                final_precursors=[]
                for each_desired_charge in desired_charges:
                    #print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                    #print "We're now working on light-charge",each_desired_charge
                    original_charge=example_precursor['@charge']
                    new_precursor=deepcopy(example_precursor)
                    if light:
                        try:
                            if each_desired_charge in sequences_with_light_bibliospec[each_peptide['@modified_sequence']]:
                                new_precursor['bibliospec_spectrum_info']=sequences_with_light_bibliospec[each_peptide['@modified_sequence']][each_desired_charge]
                        except: #It wasn't in the dict.
                            pass
                    #for each_key in new_precursor:
                    #    print each_key,new_precursor[each_key]
                    
                    if light and original_charge==each_desired_charge:
                        #print "were gonna make a ship on",new_precursor
                        final_precursors.append(new_precursor)
                        continue
                    new_precursor['@charge']=each_desired_charge
                    if light:
                        new_precursor['@precursor_mz']=str(round((float(new_precursor['@calc_neutral_mass'])+(float(each_desired_charge)*proton_mass))/float(each_desired_charge),6))
                        
                    if heavy:
                        new_precursor['@precursor_mz']=str(round(((float(new_precursor['@calc_neutral_mass'])-(heavyMass-lightMass))+(float(each_desired_charge)*proton_mass))/float(each_desired_charge),6))
                        new_precursor['@modified_sequence']=new_precursor['@modified_sequence'].replace(str(round(heavyMass,4)),str(round(lightMass,4)))
                        new_precursor['@calc_neutral_mass']=str(float(new_precursor['@calc_neutral_mass'])-(heavyMass-lightMass))
                    
                    
                    #Now we need to correct mz values of the transitions and transition precursors.
                    #print new_precursor['transition']
                    transition_ctr=0
                    charge=int(new_precursor['@charge'])
                    for each_transition in new_precursor['transition']:
                        if "@mass_index" in each_transition:
                            isotope=int(each_transition['@mass_index'])
                        else:
                            isotope=0
                        #isotope_delta_dict[str(isotope)]
                        #isotope=int(each_transition['@isotope_dist_rank'])-1
                        #print "changing the value of product",new_precursor['transition'][transition_ctr]['product_mz']
                        ##new_precursor['transition'][transition_ctr]['product_mz']=str(float(new_precursor['@precursor_mz'])+((float(isotope)*C13_mass)/charge))
                        new_precursor['transition'][transition_ctr]['product_mz']=str(round(float(new_precursor['@precursor_mz'])+((float(isotope_delta_dict[str(float(isotope))])/charge)),6))
                        #print "after:",new_precursor['transition'][transition_ctr]['product_mz']
                        new_precursor['transition'][transition_ctr]['precursor_mz']=str(float(new_precursor['@precursor_mz']))

                        transition_ctr+=1
                    
                    final_precursors.append(new_precursor)
                    
                    
                for each_desired_charge in desired_charges:
                    #print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                    #print "We're now working on heavy-charge",each_desired_charge
                    original_charge=example_precursor['@charge']
                    new_precursor=deepcopy(example_precursor)
                    new_precursor['@isotope_label']="heavy"
                    if heavy:
                        try:
                            if each_desired_charge in sequences_with_heavy_bibliospec[each_peptide['@modified_sequence']]:
                                new_precursor['bibliospec_spectrum_info']=sequences_with_heavy_bibliospec[each_peptide['@modified_sequence']][each_desired_charge]
                        except:
                            pass#also was not in the dict
                    #for each_key in new_precursor:
                    #    print each_key,new_precursor[each_key]
                    
                    if heavy and original_charge==each_desired_charge:
                        final_precursors.append(new_precursor)
                        continue
                    new_precursor['@charge']=each_desired_charge
                    if light:
                        new_precursor['@precursor_mz']=str(round(((float(new_precursor['@calc_neutral_mass'])+(heavyMass-lightMass))+(float(each_desired_charge)*proton_mass))/float(each_desired_charge),6))
                        new_precursor['@modified_sequence']=new_precursor['@modified_sequence'].replace(str(round(lightMass,4)),str(round(heavyMass,4)))
                        new_precursor['@calc_neutral_mass']=str(float(new_precursor['@calc_neutral_mass'])+(heavyMass-lightMass))
                    if heavy:
                        new_precursor['@precursor_mz']=str(round(((float(new_precursor['@calc_neutral_mass']))+(float(each_desired_charge)*proton_mass))/float(each_desired_charge),6))

                    transition_ctr=0
                    charge=int(new_precursor['@charge'])
                    #print new_precursor
                    for each_transition in new_precursor['transition']:
                        if "@mass_index" in each_transition:
                            isotope=int(each_transition['@mass_index'])
                        else:
                            isotope=0
                        #isotope=int(each_transition['@isotope_dist_rank'])-1
                        #print "changing the value of product",new_precursor['transition'][transition_ctr]['product_mz']
                        ##new_precursor['transition'][transition_ctr]['product_mz']=str(float(new_precursor['@precursor_mz'])+((float(isotope)*C13_mass)/charge))
                        #print "after:",new_precursor['transition'][transition_ctr]['product_mz'],"from a shift of",((float(isotope)*C13_mass)/charge),isotope,C13_mass,charge
                        ##new_precursor['transition'][transition_ctr]['precursor_mz']=str(float(new_precursor['@precursor_mz']))#+((heavyMass-lightMass)/charge))
                        new_precursor['transition'][transition_ctr]['product_mz']=str(round(float(new_precursor['@precursor_mz'])+((float(isotope_delta_dict[str(float(isotope))])/charge)),6))
                        #print "after:",new_precursor['transition'][transition_ctr]['product_mz']
                        new_precursor['transition'][transition_ctr]['precursor_mz']=str(float(new_precursor['@precursor_mz']))
                        transition_ctr+=1
                        
                    
                    final_precursors.append(new_precursor)
                    #charges_in_current=[]
                    #print each_peptide['precursor']
                    
                    #for each_precursor_key in each_peptide['precursor']: #I actually shouldn't be iterating over every key here. I should just be grabbing the key that I need to change the precursors.
                    #    
                    #    each_precursor=each_peptide['precursor'][each_precursor_key]
                    #    #print "here is the precursor key!: ",each_precursor_key,each_precursor['@charge']
                    #    pass
                    #    #
                    #    #NOTE, I AM SKIPPING ANY CHANGES TO THE PRECURSOR SECTION HERE IN HOPES THAT SKYLINE CAN FIGURE IT OUT ON ITS OWN!
                    #    #
                    #    #
                    #print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                each_peptide['precursor']=final_precursors
                pass #do something here because the peptide is modified
            if removeUnmod and modcount==0:
                pass
            else:
                new_peptide_list.append(each_peptide)
            completed_peptides.append(each_peptide['@modified_sequence'])
            ctr+=1
            pass
        this_new_protein['peptide']=new_peptide_list
    
    
    return this_new_protein
    
    
######################################## END OF PROTEIN MODIFICATION SECTION ########################################
    
    
def skyline_isotope_mod_generator(skylineFile,outputFile,lightMass,heavyMass,target_AA,removeUnmod):
    isotope_delta=heavyMass-lightMass

    # Open the Skyline File
    with open(skylineFile,'rb') as fd:
        doc = xmltodict.parse(fd.read())

    # Modification Initial setup - one time
    # We're going to start with two variable mods, and switch to one variable mod, with one heavy isotope modificaiton
    #print doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']
    #doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']
    
    #If there's no heavy mods already, let's set one up.
    doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications'] = collections.OrderedDict()
    new_dict=collections.OrderedDict()
    #doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications']['5']='five'
    #print doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']
    #print doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications']
    #print type(doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications'])
    new_dict['static_modification']=collections.OrderedDict()
    doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications'] = new_dict
    #print doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications']
    #Now it exists. We can populate it!
    
    doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications']['static_modification']['@name']=str(round(heavyMass-lightMass,1))#"heavy {0} {1}".format(str(isotope_delta),target_AA)
    doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications']['static_modification']['@aminoacid']=target_AA
    doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications']['static_modification']['@massdiff_monoisotopic']=isotope_delta
    doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['heavy_modifications']['static_modification']['@massdiff_average']=isotope_delta
    
    #Now that we have the heavy isotope mod set up, let's make sure that we don't have a mod for the full heavy mass still.  If so, we'll remove it.
    
    new_mods=[]
    for each_modification in doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['static_modifications']['static_modification']:
        #print each_modification_key,"key"
        #print "whole thing",doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['static_modifications'][each_modification_key]
        #each_modification=doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['static_modifications']['static_modification'][each_modification_key]
        #print each_modification,"this is a modificaiton!"
        #for index, item in each_modification[0].items():
        #    print index, item
        if '@massdiff_monoisotopic' in each_modification:
            #print each_modification["@massdiff_monoisotopic"], type(each_modification["@massdiff_monoisotopic"])
            if float(each_modification["@massdiff_monoisotopic"])==heavyMass and target_AA in each_modification["@aminoacid"]:
                #print each_modification
                
                if target_AA==each_modification['@aminoacid']:
                    #del doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['static_modifications']
                    continue #in this case, we just won't add it to "new_mods"
                    
                else:
                    each_modification['@aminoacid']=each_modification['@aminoacid'].replace(target_AA,"")
                    new_mods.append(each_modification)
            else:
                each_modification['@name']=str(round(lightMass,1))
                new_mods.append(each_modification)
        else:
            new_mods.append(each_modification)
    #print new_mods,"this is new mods!"
    doc['srm_settings']['settings_summary']['peptide_settings']['peptide_modifications']['static_modifications']['static_modification']=new_mods
    #del new_mods
    
    #Now we'll iterate over each of the peptides and find those with modifications matching either the heavy or light mod masses!
    new_proteins=[]
    #print doc['srm_settings']['settings_summary']
    i=0
    for each_protein in doc['srm_settings']['protein']:
        this_new_protein=correct_Isotope_Partners(each_protein,lightMass,heavyMass,target_AA,removeUnmod)
        #if this_new_protein
        #check to make sure that everything has a mod pep again?
        if this_new_protein is not None:
            new_proteins.append(this_new_protein)
        #print new_proteins
        #if i==2:
        #    print new_proteins
        #    sys.exit(2)
        i+=1
    
    doc['srm_settings']['protein']=new_proteins
    
    
    outFile = open(skylineFile[:-4]+'-temp.sky', 'w')
    outFile.write(xmltodict.unparse(doc, pretty=True))
    outFile.close()

    # Simply to prettify the skyline file so that the xml syntax is readable by Skyline
    #os.remove(skylineFile)
    with open(skylineFile[:-4]+'-temp.sky', 'rt') as f:
        tree = ElementTree.parse(f)
    tree.write(outputFile)
    #shutil.copy(skylineFile,"C:\\Users\\JamesLab\\LUMOSSKYLINETEST_WITHDECOY.sky")
    os.remove(skylineFile[:-4]+'-temp.sky')

#    sys.exit(2)

    
print "Now running our script to change heavy and light variable modifications into paired heavy/light isotope combos for skyline..."
print "We'll be looking at modifications on the amino acid \"{0}\".".format(options.target_AA)
print "The light mass to consider will be {0}, and the heavy {1}.".format(options.lightMass,options.heavyMass)
print "We will thusly infer that the isotope delta mass is {0}".format(options.heavyMass-options.lightMass)
print "Lastly, we will be operating on the file {0} and producing the output file {1}".format(options.skylineFile,options.outputFile)
if options.removeUnmod:
    print "We will be removing all unmodified peptides."

skyline_isotope_mod_generator(options.skylineFile,options.outputFile,options.lightMass,options.heavyMass,options.target_AA,options.removeUnmod)
