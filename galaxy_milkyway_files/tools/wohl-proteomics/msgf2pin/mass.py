import re
from orangecontrib.bio.ontology import OBOParser
import orangecontrib
from pyteomics import mass
import sys
mod_bracket_re='\[|\]'
mod_dict={}

pepseq="VS[UNIMOD:21]KLKNWEY[UNIMOD:21]R"
mod_obo="/galaxy-central/tools/wohl-proteomics/MSGFcrux/unimod.obo"
with open(mod_obo,"r") as unimod_obo:
    obo_parser=orangecontrib.bio.ontology.OBOParser(unimod_obo)
    if "UNIMOD" in pepseq:
        print pepseq,"working on this..."
        split_peptide=re.split(mod_bracket_re,pepseq)
        #print split_peptide
        new_peptide=[]
        mod_mass=0.0
        unmod_pepmass=0.0
        for each_split in split_peptide:
            if "UNIMOD" in each_split:
                if each_split in mod_dict:
                    mod_mass+=mod_dict[each_split]
                else:
                    thismod_mass=""
                    print each_split," was not in the dictionary!"
                    trigger=False
                    unimod_obo.seek(0)
                    for event,value in obo_parser:
                        if "TAG_VALUE" in event and not trigger:
                            print event,value
                            if each_split in value[1]:
                                trigger=True
                        elif trigger:
                            if "delta_mono_mass" in value[1]:
                                print value[1],"val1 should be the mass"
                                thismod_mass=value[1].split("\"")[1]
                                trigger=False
                                break
                            else:
                                continue
                    if thismod_mass=="":
                        print "ERROR: ERROR: ERROR: THE MASS FOR THIS MOD WAS NOT FOUND IN THE UNIMOD OBO FILE..."
                        sys.exit(2)
                    mod_dict[each_split]=float(thismod_mass)
                    mod_mass+=mod_dict[each_split]


                #print each_split,"with unimod..."
            else:
                new_peptide.extend(each_split)
        unmod_pepmass=mass.fast_mass(new_peptide)
        print mod_mass,"mod mass"
        print unmod_pepmass+mod_mass,"This is the combined mass..."
        print new_peptide




