### This file will contain an obo file parser
### The code from this parser is based on the work of Uli Kohler, as noted below.
### The OBO file parser will replace MilkyWay's dependence upon the Orange-Bioinformatics package.
### The class construction and _minor_ changes were added.

from __future__ import with_statement
from collections import defaultdict

__author__    = "Uli Kohler"
__copyright__ = "Copyright 2013 Uli Kohler"
__license__   = "Apache v2.0"
__version__   = "1.1"

class obo_parser:
    #unimod_to_mass_dict={}
    #unimod_to_modstr_dict={}

    def __init__(self):
        self.unimod_to_mass_dict={}
        self.unimod_to_modstr_dict={}
        self.unimod_int_to_modstr_dict={}
        self.unimod_int_to_mass_dict={}

    def processGOTerm(self,goTerm):
        """
        In an object representing a GO term, replace single-element lists with
        their only member.
        Returns the modified object as a dictionary.
        """
        ret = dict(goTerm) #Input is a defaultdict, might express unexpected behaviour
        for key, value in ret.items():
            if len(value) == 1:
                ret[key] = value[0]
        return ret

    def parseGOOBO(self,filename):
        """
        Parses a Gene Ontology dump in OBO v1.2 format.
        Yields each 
        Keyword arguments:
            filename: The filename to read
        """
        with open(filename, "r") as infile:
            currentGOTerm = None
            for line in infile:
                line = line.strip()
                if not line: continue #Skip empty
                if line == "[Term]":
                    if currentGOTerm: yield self.processGOTerm(currentGOTerm)
                    currentGOTerm = defaultdict(list)
                elif line == "[Typedef]":
                    #Skip [Typedef sections]
                    currentGOTerm = None
                else: #Not [Term]
                    #Only process if we're inside a [Term] environment
                    if currentGOTerm is None: continue
                    key, sep, val = line.partition(":")
                    currentGOTerm[key].append(val.strip())
            #Add last term
            if currentGOTerm is not None:
                yield self.processGOTerm(currentGOTerm)


    def consumeGO(self, go_file):
        for each in self.parseGOOBO(go_file):
            if not 'xref' in each.keys():
                continue
            mass=[x.split(" ",1)[1].replace("\"","") for x in each['xref'] if "delta_mono_mass" in x]
            if len(mass) < 1:
                continue #Nothing there...
            else:
                mass_float=float(mass[0])
            if float(mass_float)<0:
                mass_str="[{0}]".format(mass_float)
            else:
                mass_str="[+{0}]".format(mass_float)
            #print mass_str
            #print each['id']
            #print mass_float
            self.unimod_to_mass_dict[each['id']]=mass_float
            self.unimod_to_modstr_dict[each['id']]=mass_str
            self.unimod_int_to_mass_dict[each['id'].rsplit(":",1)[1]]=mass_float
            self.unimod_int_to_modstr_dict[each['id'].rsplit(":",1)[1]]=mass_str
        pass #Here's the place to add them to the unimod_to_mass_dict and modstr_dict
