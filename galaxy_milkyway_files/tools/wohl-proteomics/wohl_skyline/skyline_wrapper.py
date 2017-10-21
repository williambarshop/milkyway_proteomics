import os, sys, re, time
import optparse
import shutil
import pandas
import numpy
import subprocess,tempfile,glob
import datetime
import sqlite3
import string
from pyteomics import parser as pyparser
from pyteomics import fasta
from operator import itemgetter
import collections
#from tqdm import *
import uniprot as uni
import os.path
import xmltodict
import collections
import warnings
from random import randint
from copy import deepcopy
from xml.etree import ElementTree


###### Redirection of Pandas warning messages to STDOUT instead of STDERR for Galaxy's sake...  ######
def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))

warnings.showwarning = customwarn
#####                CODE IMPLEMENTED FROM STACK OVERFLOW ANSWER #858916                        ######



#warnings.filterwarnings('ignore')

#import subprocesswin32 as subprocess
#import time

#####################################
#This is the wrapper for SkylineRunner.
#
#VERSION 0.70A
version="0.70A"
#DATE: 7/21/2017
date="7/21/2017"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the SkylineRunner wrapper for Galaxy, Wohlschlegel Lab UCLA"
print "Written by William Barshop and Hee Jong Kim"
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

print sys.argv,"THESE ARE THE ARGS"
# some versions of galaxy append junk onto the end of the commands that we have to get rid of before argument parsing.
# so let's take care of that now...
arguments=sys.argv
final_args=[]
for each_arg in sys.argv:
    if each_arg==";":
        break
    else:
        final_args.append(each_arg)
print "final argument set:",final_args


parser = optparse.OptionParser()
parser.add_option("--sslgz",action="store",type="string",dest="ssl_gz")
parser.add_option("--mzml",action="store",type="string",dest="mzml")
parser.add_option("--fractions",action="store_true",dest="fractions")
parser.add_option("--makeblib",action="store_true",dest="makeblib")
parser.add_option("--db",action="store",type="string",dest="db") #PROTEINS TO QUANTIATE (ACCEPTED)
parser.add_option("--add_mods",action="store_true",dest="add_mods")
parser.add_option("--mode",action="store",type="string",dest="mode")
parser.add_option("--ms1_res",action="store", type="int", dest="ms1_res")
parser.add_option("--ms1_res_mz",action="store", type="float", dest="ms1_res_mz")
parser.add_option("--ms2_res",action="store", type="int", dest="ms2_res")
parser.add_option("--ms2_res_mz",action="store", type="float", dest="ms2_res_mz")
parser.add_option("--skyline_file",action="store",type="string",dest="skyline_file")
parser.add_option("--rt_tolerance",action="store", type="float", dest="rt_tolerance")
parser.add_option("--ms1_template",action="store", type="string", dest="ms1_template")
parser.add_option("--ms2_template",action="store", type="string", dest="ms2_template")
parser.add_option("--dia_template",action="store", type="string", dest="dia_template")
parser.add_option("--msstats",action="store", type="string", dest="msstats")
parser.add_option("--protdb",action="store", type="string", dest="protdb")
parser.add_option("--mprophet",action="store",type="string",dest="mprophet") #Should be "secondbest" or "decoy" or None
parser.add_option("--output",action="store", type="string", dest="output")
parser.add_option("--background_db",action="store", type="string", dest="background_db") #if no protdb provided, the backgrounddb FASTA should be provided insteads
parser.add_option("--enzyme",action="store", type="string", dest="enzyme") #if no protdb provided, the enzyme
parser.add_option("--missedcleave",action="store", type="int", dest="missedcleave") #if no protdb provided, the max missed cleaves
parser.add_option("--exp_file",action="store", type="string", dest="exp_file") # Wohl-Experimental Design File
parser.add_option("--chromatograms_file",action="store", type="string", dest="chromatograms") # Wohl-Experimental Design File
parser.add_option("--quantiva_transitions",action="store", type="string", dest="quantiva_transitions") # Export quantiva transitions?
parser.add_option("--peak_boundaries",action="store", type="string", dest="peak_boundaries")
parser.add_option("--num_procs",action="store", type="int", dest="num_procs")

(options,args) = parser.parse_args(final_args)

###First things first, let's unzip the tar.gz file...
working_directory = os.getcwd()


#print "Testing SKR NOSHELL"
#test=subprocess.Popen(args="SkylineCmd")
#returncode = test.wait()
#print "done testing SKR, now batch"
#with open("test.bat",'wb') as tempwriter:
#    tempwriter.write("SkylineCmd")

#test2=subprocess.Popen("test.bat")
#returncode2=test2.wait()
#print "done testing SKR-batch, now moving on..."




global tmp_stderr_name
global tmp_stdout_name
tmp_stderr_name = tempfile.NamedTemporaryFile(dir = working_directory, suffix = '.stderr').name
tmp_stdout_name = tempfile.NamedTemporaryFile(dir = working_directory, suffix = '.stdout').name


def correctBackgroundProteome(skylineFile):
    with open(skylineFile,'rb') as fd:
        doc = xmltodict.parse(fd.read())
        time_string=str(datetime.datetime.now().time()).replace(":",'').replace('.','')
	if 'background_proteome' in doc['srm_settings']['settings_summary']['peptide_settings'].keys():
            #doc['srm_settings']['settings_summary']['background_proteome']
	    doc['srm_settings']['settings_summary']['peptide_settings']['background_proteome']['@name']="background_proteome"+time_string
            doc['srm_settings']['settings_summary']['peptide_settings']['background_proteome']['@database_path']="background_proteome"+time_string+".protdb"
	else:
            #print doc['srm_settings']['settings_summary']
	    doc['srm_settings']['settings_summary']['peptide_settings']['background_proteome']={}
            doc['srm_settings']['settings_summary']['peptide_settings']['background_proteome']['@name']="background_proteome"+time_string
            doc['srm_settings']['settings_summary']['peptide_settings']['background_proteome']['@database_path']="background_proteome"+time_string+".protdb"

    # Initial temporary file to store the modified skyline file
    outFile = open(skylineFile[:-4]+'-temp.sky', 'w')
    outFile.write(xmltodict.unparse(doc, pretty=True))
    outFile.close()

    # Simply to prettify the skyline file so that the xml syntax is readable by Skyline
    with open(skylineFile[:-4]+'-temp.sky', 'rt') as f:
        tree = ElementTree.parse(f)
    tree.write(skylineFile)
    os.remove(skylineFile[:-4]+'-temp.sky')
    #move the background proteome to the correct name locations...
    with open(skylineFile,'rt') as f:
        line_ctr=0
        digest_line=None
        background_proteome_line=None
        background_proteome_line_content=None
        for each_line in f:
            if 'digest_settings' in each_line:
                digest_line=line_ctr
            elif '<background_proteome' in each_line:
                background_proteome_line_content=each_line
                background_proteome_line=line_ctr
            line_ctr+=1            

    with open(skylineFile,'rt') as f:
        line_ctr=0
        new_file=[]
        for each_line in f:
            if line_ctr<=digest_line:
                new_file.append(each_line)
            elif line_ctr==digest_line+1:
                new_file.append(background_proteome_line_content)
                if '<background_proteome' not in each_line:
                    new_file.append(each_line)
            elif line_ctr==background_proteome_line:
                pass
            elif '<background_proteome':
                new_file.append(each_line)
            line_ctr+=1
    os.remove(skylineFile)
    with open(skylineFile,'w') as fwriter:
        for each_line in new_file:
            fwriter.write(each_line)

    return time_string #We will need this later to rename the background proteome...

#def renameBackgroundProteome(time_string):
#    os.rename("background_proteome.protdb","background_proteome"+time_string+".protdb")

########################### BEGIN METHODS FOR DECOY PEPTIDE GENERATION ###########################
# Skyline Decoy Peptide Generator with mass shifts
def skyline_decoy_generator(skylineFile):

    #shutil.copy(skylineFile,"C:\\Users\\JamesLab\\LUMOSSKYLINETEST_NODECOY.sky")


    # Open the Skyline File
    attempt=0
    #while attempt<10:
    #    time.sleep(30)
    #    try:
    with open(skylineFile,'rb') as fd:
        doc = xmltodict.parse(fd.read())
    #    except Exception as e:
    #        shutil.copyfile(skylineFile,"temporary")
    #        os.remove(skylineFile)
    #        shutil.move("temporary",skylineFile)
    #        print "moved the file back and forth to temporary"
    #        print "FAILED",e.message,e.args
    #        attempt+=1

    #doc=convert(doc)
    max_mz=float(doc['srm_settings']['settings_summary']['transition_settings']['transition_instrument']['@max_mz'])
    min_mz=float(doc['srm_settings']['settings_summary']['transition_settings']['transition_instrument']['@min_mz'])

    # Decoy Entry Initial setup - one time
    doc['srm_settings']['peptide_list'] = collections.OrderedDict()
    doc['srm_settings']['peptide_list']['@label_name'] = 'Decoys'
    doc['srm_settings']['peptide_list']['@websearch_status'] = 'X'
    doc['srm_settings']['peptide_list']['@auto_manage_children'] = 'false'
    doc['srm_settings']['peptide_list']['@decoy'] = 'true'

    # Prepare the decoy peptide placeholder
    doc['srm_settings']['peptide_list']['peptide'] = []

    # Duplicate the peptide level info for generating the decoy dataset
    for each_protein in doc['srm_settings']['protein']:
        #print each_protein['peptide']
        if 'dict' in str(type(each_protein['peptide'])) or 'Dict' in str(type(each_protein['peptide'])):
            doc['srm_settings']['peptide_list']['peptide'].append(deepcopy(each_protein['peptide']))
            continue
        else:
            for each_peptide in each_protein['peptide']:
                #print type(each_peptide),each_peptide,"\n\n"
                #sys.exit(2)
                doc['srm_settings']['peptide_list']['peptide'].append(deepcopy(each_peptide))
    #sys.exit(2)
    # Generate the Decoy Peptides with shifted m/z
    for each_decoy_peptide in doc['srm_settings']['peptide_list']['peptide']:
        if 'list' in str(type(each_decoy_peptide)):
            print each_decoy_peptide,"this is a list and we expect a dict of somekind -- SKIPPING"
            continue
            #sys.exit(2)
        if 'str' in str(type(each_decoy_peptide)) or 'unicode' in str(type(each_decoy_peptide)):
            print each_decoy_peptide,"this is a str/unicode and we expect a dict of somekind -- SKIPPING"
            continue
            #sys.exit(2)
            
        outofrange=True
        while outofrange:
            outofrange=False
            random_decoy_mass_shift = randint(1,5) * (randint(0,1)*2-1)
            min_result=0
            max_result=0

            if 'list' not in str(type(each_decoy_peptide['precursor'])):
                each_precursor=each_decoy_peptide['precursor']
                tmp = float(each_precursor['@precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                if tmp < min_result or min_result == 0:
                    min_result=tmp
                if tmp > max_result or max_result == 0:
                    max_result=tmp
                if min_result != 0 and min_result < min_mz:
                    outofrange=True
                    break
                if max_result != 0 and max_result > max_mz:
                    outofrange=True
                    break
                    
                # Iterate transition
                if 'list' not in str(type(each_precursor['transition'])):
                    each_transition=each_precursor['transition']
                    tmp = float(each_transition['precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                    if tmp < min_result or min_result == 0:
                        min_result=tmp
                    if tmp > max_result or max_result == 0:
                        max_result=tmp
                    if min_result != 0 and min_result < min_mz:
                        outofrange=True
                        break
                    if max_result != 0 and max_result > max_mz:
                        outofrange=True
                        break
                    tmp = float(each_transition['product_mz']) + (random_decoy_mass_shift * 1.000454)
                    if tmp < min_result or min_result == 0:
                        min_result=tmp
                    if tmp > max_result or max_result == 0:
                        max_result=tmp
                    if min_result != 0 and min_result < min_mz:
                        outofrange=True
                        break
                    if max_result != 0 and max_result > max_mz:
                        outofrange=True
                        break
                else:
                    for each_transition in each_precursor['transition']:#.items():
                        tmp = float(each_transition['precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                        if tmp < min_result or min_result == 0:
                            min_result=tmp
                        if tmp > max_result or max_result == 0:
                            max_result=tmp
                        if min_result != 0 and min_result < min_mz:
                            outofrange=True
                            break
                        if max_result != 0 and max_result > max_mz:
                            outofrange=True
                            break
                        tmp = float(each_transition['product_mz']) + (random_decoy_mass_shift * 1.000454)
                        if tmp < min_result or min_result == 0:
                            min_result=tmp
                        if tmp > max_result or max_result == 0:
                            max_result=tmp
                        if min_result != 0 and min_result < min_mz:
                            outofrange=True
                            break
                        if max_result != 0 and max_result > max_mz:
                            outofrange=True
                            break
            else:
                for each_precursor in each_decoy_peptide['precursor']:#.items():
                    tmp = float(each_precursor['@precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                    if tmp < min_result or min_result == 0:
                        min_result=tmp
                    if tmp > max_result or max_result == 0:
                        max_result=tmp
                    if min_result != 0 and min_result < min_mz:
                        outofrange=True
                        break
                    if max_result != 0 and max_result > max_mz:
                        outofrange=True
                        break
                    # Iterate transition
                    for each_transition in each_precursor['transition']:
                        tmp = float(each_transition['precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                        if tmp < min_result or min_result == 0:
                            min_result=tmp
                        if tmp > max_result or max_result == 0:
                            max_result=tmp
                        if min_result != 0 and min_result < min_mz:
                            outofrange=True
                            break
                        if max_result != 0 and max_result > max_mz:
                            outofrange=True
                            break
                        tmp = float(each_transition['product_mz']) + (random_decoy_mass_shift * 1.000454)
                        if tmp < min_result or min_result == 0:
                            min_result=tmp
                        if tmp > max_result or max_result == 0:
                            max_result=tmp
                        if min_result != 0 and min_result < min_mz:
                            outofrange=True
                            break
                        if max_result != 0 and max_result > max_mz:
                            outofrange=True
                            break
        
        # Remove a few tags...
        del each_decoy_peptide['@start']
        del each_decoy_peptide['@end']
        del each_decoy_peptide['@prev_aa']
        del each_decoy_peptide['@next_aa']

        # Add Decoy specific tags
        each_decoy_peptide['@auto_manage_children'] = 'false'
        each_decoy_peptide['@decoy'] = 'true'

        #If we're this far in, that means the mass shift was okay.  Let's go ahead and move forward!
        if 'list' not in str(type(each_decoy_peptide['precursor'])):
            #print "THIS IS NOT A LIST OF PRECURSOR",each_decoy_peptide['precursor']
            each_precursor=each_decoy_peptide['precursor']
            each_precursor['@auto_manage_children'] = "false"
            #print "FAILED ON",each_precursor['@auto_manage_children']

            #sys.exit(2)
            each_precursor['@decoy_mass_shift'] = random_decoy_mass_shift
            each_precursor['@precursor_mz'] = float(each_precursor['@precursor_mz']) + (random_decoy_mass_shift * 1.000454)
            # Iterate transition
            if 'list' not in str(type(each_precursor['transition'])):
                #print each_precursor['transition'], "NOT LIST PREC, AND NOT LIST TRANSITION"
                #sys.exit(2)
                each_transition=each_precursor['transition']
                # Add decoy_mass_shift tag
                each_transition['@decoy_mass_shift'] = random_decoy_mass_shift
                # Calculate & replace the shifted m/z
                each_transition['precursor_mz'] = float(each_transition['precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                each_transition['product_mz'] = float(each_transition['product_mz']) + (random_decoy_mass_shift * 1.000454)
            else:
                #print each_precursor['transition'],"WORKING..............."
                #sys.exit(2)
                for each_transition in each_precursor['transition']:#.items():
                    # Add decoy_mass_shift tag
                    each_transition['@decoy_mass_shift'] = random_decoy_mass_shift
                    # Calculate & replace the shifted m/z
                    each_transition['precursor_mz'] = float(each_transition['precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                    each_transition['product_mz'] = float(each_transition['product_mz']) + (random_decoy_mass_shift * 1.000454)
                    # each_transition['collision_energy'] = PUT THE CALCULATED VALUE HERE IF SHIFTED COLLISION ENERGY IS NEEDED
        else:
            #print "THIS IS A LIST OF PRECURSOR!", each_decoy_peptide['precursor']
            for each_precursor in each_decoy_peptide['precursor']:#.items():

                # Add decoy_mass_shift tags

                each_precursor['@auto_manage_children'] = "false"
                #print "FAILED ON",each_precursor['@auto_manage_children']

                #sys.exit(2)
                each_precursor['@decoy_mass_shift'] = random_decoy_mass_shift
                each_precursor['@precursor_mz'] = float(each_precursor['@precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                # Iterate transition
                for each_transition in each_precursor['transition']:#.items():
                    # Add decoy_mass_shift tag
                    each_transition['@decoy_mass_shift'] = random_decoy_mass_shift
                    # Calculate & replace the shifted m/z
                    each_transition['precursor_mz'] = float(each_transition['precursor_mz']) + (random_decoy_mass_shift * 1.000454)
                    each_transition['product_mz'] = float(each_transition['product_mz']) + (random_decoy_mass_shift * 1.000454)
                    # each_transition['collision_energy'] = PUT THE CALCULATED VALUE HERE IF SHIFTED COLLISION ENERGY IS NEEDED

    # Initial temporary file to store the decoy-included skyline file
    outFile = open(skylineFile[:-4]+'-temp.sky', 'w')
    outFile.write(xmltodict.unparse(doc, pretty=True))
    outFile.close()

    # Simply to prettify the skyline file so that the xml syntax is readable by Skyline
    os.remove(skylineFile)
    with open(skylineFile[:-4]+'-temp.sky', 'rt') as f:
        tree = ElementTree.parse(f)
    tree.write(skylineFile)
    #shutil.copy(skylineFile,"C:\\Users\\JamesLab\\LUMOSSKYLINETEST_WITHDECOY.sky")
    os.remove(skylineFile[:-4]+'-temp.sky')

###########################  END METHODS FOR DECOY PEPTIDE GENERATION  ###########################




########################### BEGIN METHODS FOR PROTDB GENERATION ###########################

def create_protdb(db_filename):
    db = sqlite3.connect(db_filename)
    cursor = db.cursor()
    cursor.execute('''
        CREATE TABLE ProteomeDbDigestedPeptide (
            Id  integer PRIMARY KEY AUTOINCREMENT,
            Digestion BIGINT  NOT NULL,
            Sequence  TEXT  NOT NULL,
            UNIQUE (
            Digestion,
            Sequence
            ),
            CONSTRAINT FK33FE12D9FB0ADC88 FOREIGN KEY (
            Digestion
            )
            REFERENCES ProteomeDbDigestion
        )
    ''')
    cursor.execute('''
        CREATE TABLE ProteomeDbDigestedPeptideProtein (
            Id  integer PRIMARY KEY AUTOINCREMENT,
            Peptide BIGINT  NOT NULL,
            Protein BIGINT  NOT NULL,
            UNIQUE (
            Peptide,
            Protein
            ),
            CONSTRAINT FK6C498DEE7B29BF71 FOREIGN KEY (
            Peptide
            )
            REFERENCES ProteomeDbDigestedPeptide,
            CONSTRAINT FK6C498DEE6EF722BA FOREIGN KEY (
            Protein
            )
            REFERENCES ProteomeDbProtein
        )
    ''')
    cursor.execute('''
        CREATE TABLE ProteomeDbDigestion (
            Id  integer PRIMARY KEY AUTOINCREMENT,
            Version INT  NOT NULL,
            Name      TEXT  UNIQUE,
            MinSequenceLength INT,
            MaxSequenceLength INT,
            Description TEXT
        )
    ''')
    cursor.execute('''
        CREATE TABLE ProteomeDbProtein (
            Id   integer PRIMARY KEY AUTOINCREMENT,
            Version  INT     NOT NULL,
            Sequence TEXT   UNIQUE
        )
    ''')
    cursor.execute('''
        CREATE TABLE ProteomeDbProteinAnnotation (
            Id  integer PRIMARY KEY AUTOINCREMENT,
            Version INT  NOT NULL,
            Protein BIGINT,
            Name    TEXT,
            Value   TEXT,
            UNIQUE (
            Protein,
            Name
            ),
            CONSTRAINT FKB67691996EF722BA FOREIGN KEY (
            Protein
            )
            REFERENCES ProteomeDbProtein
        )
    ''')
    cursor.execute('''
        CREATE TABLE ProteomeDbProteinName (
            Id  integer PRIMARY KEY AUTOINCREMENT,
            Version  INT  NOT NULL,
            Protein  BIGINT,
            IsPrimary     BOOL,
            Name    TEXT    UNIQUE,
            Description  TEXT,
            PreferredName   TEXT,
            Accession     TEXT,
            Gene    TEXT,
            Species  TEXT,
            WebSearchStatus TEXT,
            CONSTRAINT FK50469A296EF722BA FOREIGN KEY (
            Protein
            )
            REFERENCES ProteomeDbProtein
        )
    ''')
    cursor.execute('''
        CREATE TABLE ProteomeDbSchemaVersion (
            Id   integer PRIMARY KEY AUTOINCREMENT,
            SchemaVersionMajor INT,
            SchemaVersionMinor INT
        )
    ''')
    db.close()

def cleave_trunc(sequence, enzyme, missed, enzyme_n_terminal_re, enzyme_c_terminal_exception,enzyme_c_terminal_re):
    # Actual digestion
    cleaved_peptide_set = pyparser.cleave(sequence, enzyme, missed)

    # Clean up the list for only skyline needs
    peptide_sequence_holder =""
    processed_peptide_list = []
    for digested_peptide in sorted(cleaved_peptide_set):
        truncated_peptide_sequence = digested_peptide[:7]
        if truncated_peptide_sequence != peptide_sequence_holder and len(digested_peptide) >= 4:
            processed_peptide_list.append(truncated_peptide_sequence)
            peptide_sequence_holder = truncated_peptide_sequence

    # Reorder the digested peptides corresponding to the protein sequence
    each_peptide_info=[]
    for each_peptide in processed_peptide_list:
        # locate each digested peptide's location and add index into the list
        peptide_position = [m.start() for m in re.finditer(enzyme_n_terminal_re[:-1] + each_peptide + ')', sequence)]
        peptide_position = [x + 1 for x in peptide_position]
        # Beginning of the protein's index 
        if sequence.startswith(each_peptide):
            peptide_position = [0]
        # Handle multiple occurrence of same peptide sequence 
        for each_peptide_position in peptide_position:
            # Take care of ending peptide mismatch
            if len(each_peptide) < 7 and re.search(enzyme_c_terminal_re, each_peptide) and each_peptide_position != len(sequence) - len(each_peptide):
                continue
            # Take care of c-terminal cleavage exception: e.g. Trypsin KR|'P'
            if (each_peptide_position+len(each_peptide)) < len(sequence):
                if len(each_peptide) < 7 and sequence[each_peptide_position+len(each_peptide)] == enzyme_c_terminal_exception:
                    continue
            each_peptide_info.append([each_peptide_position,each_peptide])
    reorganized_peptide_list = sorted(each_peptide_info, key=itemgetter(0))

    # Remove index
    final_peptide_list = []
    for each_peptide in reorganized_peptide_list:
        final_peptide_list.append(each_peptide[1])

    return final_peptide_list


def fasta_processing(fasta_filename, enzyme, missed, enzyme_n_terminal_re, enzyme_c_terminal_exception,enzyme_c_terminal_re):
    # fasta processing for digested peptide, protein sequence, and their relationship
    # Define Dictionaries
    DbProteinName = collections.OrderedDict()
    DbProtein = collections.OrderedDict()
    DbDigestedPeptide = collections.OrderedDict()
    DbDigestedPeptideProtein = collections.OrderedDict()
    DbDigestedPeptideProtein_duplicate_check = []

    protein_index = 1
    peptide_index = 1
    peptideProtein_index = 1

    print "* Protein Digestion and Mapping..."
    fasta_tuples = fasta.read(fasta_filename)
    for descr, protein_seq in fasta_tuples:
        # * containing sequence conversion with X amino acid notation
        protein_seq = string.replace(protein_seq, '*', 'X')

        # Validate the protein sequence uniqueness
        if not protein_seq in DbProtein.values():
            # ID - Accesion & description
            DbProteinName[protein_index] = [descr.split(' ')[0], ' '.join(descr.split(' ')[1:])]
            # Full Sequence storage
            DbProtein[protein_index] = protein_seq
            digested_peptides = cleave_trunc(protein_seq, enzyme, missed, enzyme_n_terminal_re, enzyme_c_terminal_exception, enzyme_c_terminal_re)
            for each_peptide in digested_peptides:
                if not each_peptide in DbDigestedPeptide:
                    DbDigestedPeptide[each_peptide] = peptide_index
                    DbDigestedPeptideProtein[peptideProtein_index] = [peptide_index, protein_index]
                    DbDigestedPeptideProtein_duplicate_check.append(hash(str(peptide_index)+','+str(protein_index)))
                    peptideProtein_index += 1
                    peptide_index += 1
                elif hash(str(DbDigestedPeptide[each_peptide])+','+str(protein_index)) in DbDigestedPeptideProtein_duplicate_check:
                    pass
                else:
                    DbDigestedPeptideProtein[peptideProtein_index] = [DbDigestedPeptide[each_peptide], protein_index]
                    DbDigestedPeptideProtein_duplicate_check.append(hash(str(DbDigestedPeptide[each_peptide])+','+str(protein_index)))
                    peptideProtein_index += 1

            protein_index +=1
            DbDigestedPeptideProtein_duplicate_check = []
    print "* Protein Digestion and Mapping has been completed."
    return (DbProteinName, DbProtein, DbDigestedPeptide, DbDigestedPeptideProtein)



def insert_digestion_info(db_filename, enzyme_name):
    db = sqlite3.connect(db_filename)
    cursor = db.cursor()
    # Insert ProteomeDbDigestion Info
    digestionInfo = (1, enzyme_name, 4, 7)
    cursor.execute('''INSERT INTO ProteomeDbDigestion(Version, Name, MinSequenceLength, MaxSequenceLength) VALUES (?,?,?,?)''', digestionInfo)
    db.commit()
    db.close()

def insert_schema_version(db_filename):
    db = sqlite3.connect(db_filename)
    cursor = db.cursor()
    # Insert ProteomeDbSchemaVersion Info
    cursor.execute('''INSERT INTO ProteomeDbSchemaVersion(SchemaVersionMajor, SchemaVersionMinor) VALUES (0,1)''')
    db.commit()
    db.close()

def insert_protein_sequence(db_filename, DbProtein):
    db = sqlite3.connect(db_filename)
    cursor = db.cursor()

    # Protein List Prep for multiple insertion
    protein_insertion_list = []
    for protein_id, protein_seq in DbProtein.items():
        protein_insertion_list.append((protein_id, "1", protein_seq))
    # INSERT PROTEIN SEQUENCE
    print "* Inserting protein sequences..."
    cursor.executemany('''
        INSERT INTO ProteomeDbProtein(Id, Version, Sequence)
        VALUES(?,?,?)''', protein_insertion_list)
    print "* protein sequences have been inserted."
    db.commit()
    db.close()

def insert_peptide_sequence(db_filename, DbDigestedPeptide):
    db = sqlite3.connect(db_filename)
    cursor = db.cursor()
    
    # Digested Peptide List Prep for multiple insertion
    digested_peptide_list=[]
    for pep, pep_id in DbDigestedPeptide.items():
        digested_peptide_list.append((pep_id, "1", pep))
    
    # INSERT DIGESTED PEPTIDE SEQUENCE
    print "* Inserting digested peptide sequences..."
    cursor.executemany('''
        INSERT INTO ProteomeDbDigestedPeptide(Id, Digestion, Sequence)
        VALUES(?,?,?)''', digested_peptide_list)
    print "* digested peptide sequences have been inserted."
    db.commit()
    db.close()

def insert_relationship(db_filename, DbDigestedPeptideProtein):
    db = sqlite3.connect(db_filename)
    cursor = db.cursor()
    
    # Relationship List Prep for multiple insertion
    relationship_insertion_list=[]
    for relationship_id, mapping_list in DbDigestedPeptideProtein.items():
        relationship_insertion_list.append((relationship_id, mapping_list[0], mapping_list[1]))

    # INSERT PEPTIDE PROTEIN RELATIONSHIP
    print "* Inserting peptide-protein mapping..."
    cursor.executemany('''
        INSERT INTO ProteomeDbDigestedPeptideProtein(Id, Peptide, Protein)
        VALUES(?,?,?)''', relationship_insertion_list)
    print "* Peptide-protein mapping has been completed."
    db.commit()
    db.close()

def uniprot_insert_protein_info(db_filename, DbProteinName):
    db = sqlite3.connect(db_filename)
    cursor = db.cursor()
    # INSERT PEPTIDE PROTEIN RELATIONSHIP
    print "Inserting Uniprot Information..."

    for protein_id, description in DbProteinName.items():
        protein_info = uni.retrieve(description[0])
        WebSearchStatus = 'X#U' + description[0]
        for line in protein_info.splitlines():
            if line.split('   ')[0] == 'ID':
                PreferredName = line.split('   ')[1]
            if line.split('   ')[0] == 'GN':
                Gene = line.split('   ')[1].replace('Name=', '').replace('; Synonyms=', ' ').replace(',', '').replace(';', '')
            if line.split('   ')[0] == 'OS':
                Species = line.split('   ')[1].replace('.', '')
        cursor.execute('''
            INSERT INTO ProteomeDbProteinName(Version, Protein, IsPrimary, Name, Description, PreferredName, Accession, Gene, Species, WebSearchStatus)
            VALUES(?,?,?,?,?,?,?,?,?,?)''', ('2',protein_id,'1',description[0],description[1],PreferredName,description[0],Gene,Species,WebSearchStatus,))
        db.commit()
    db.close()


###########################  END METHODS FOR PROTDB GENERATION  ###########################









def read_stderr(): #THIS CODE COMES FROM THE MSCONVERT PACKAGE FROM PROTK/GALAXY-P/BGRUENING(probably)/J.CHILTON
    stderr = ''
    if(os.path.exists(tmp_stderr_name)):
        with open(tmp_stderr_name, 'rb') as tmp_stderr:
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read(buffsize)
                    if not stderr or len(stderr) % buffsize != 0:
                        break
            except OverflowError:
                pass
    return stderr





#cmd="tar xzvf "+options.ssl_gz
#cmd="7z e "+options.ssl_gz+" && 7z x "+options.ssl_gz.rsplit(".",1)[0]

#shutil.copy(options.ssl_gz,"archive.tar.gz")
cmd="7z x "+options.ssl_gz+" -so -tgzip | 7z x -aoa -si -ttar"
print cmd,"This is the cmd..."

with open(tmp_stderr_name, 'wb') as tmp_stderr:
    with open(tmp_stdout_name, 'wb') as tmp_stdout:
        proc = subprocess.Popen(args=cmd, shell=True, stderr=tmp_stderr.fileno(), stdout=tmp_stdout.fileno(), env=os.environ)
        returncode = proc.wait()
        if returncode != 0:
            raise Exception, "Program returned with non-zero exit code %d. stderr: %s" % (returncode, read_stderr())

os.chdir("ssl_files")
os.mkdir("skyline_output/")
output_dir=os.path.join(os.getcwd(),"skyline_output/")

basedir=os.getcwd()

list_of_folders_with_ssl=set(folder for folder, subfolders, files in os.walk(basedir) for file_ in files if os.path.splitext(file_)[1] == '.ssl')

def getTemplateFile(mode,ms1_template,ms2_template,dia_template,provided):
    if "dda_ms1" in mode:
        return ms1_template
    elif "ms2_prm" in mode:
        return ms2_template
    elif "dia" in mode:
        return dia_template
    elif "provided_file" in mode:
        return provided
    else:
        print "I don't know what this mode is!"
        sys.exit(0)

def prependLine(filename,line):# From StackOverflow Q#5914627
    with open(filename,'r+') as myfile:
        content= myfile.read()
        myfile.seek(0,0)
        myfile.write(line.rstrip('\r\n')+'\n'+content)




if options.fractions:
    #We have fractionated data.  This means that we'll make a Skyline file for each folder (read: fraction) of the dataset.
    #Let's go ahead and build the shared portion of the command...
    skyline_runner_cmd = "SkylineCmd --batch-commands="
    #skyline_runner_cmd = "SkylineRunner --batch-commands="
    os.chdir(basedir)
    cfgwriter=open("skyline_batch",'wb')
    cfgwriter.write("--timestamp\n")
    if options.add_mods:
        cfgwriter.write("--import-search-add-mods\n")
    if "dda_ms1" in options.mode:
        cfgwriter.write("--full-scan-precursor-res=\n")
        cfgwriter.write("--full-scan-precursor-res-mz=\n")
        cfgwriter.write("--full-scan-rt-filter-tolerance="+str(options.rt_tolerance)+"\n")
    elif "ms2_prm" in options.mode:
        cfgwriter.write("--full-scan-product-res=\n")
        cfgwriter.write("--full-scan-product-res-mz=\n")
        cfgwriter.write("--full-scan-rt-filter-tolerance="+str(options.rt_tolerance)+"\n")
        pass
    elif "dia" in options.mode:
        cfgwriter.write("--full-scan-precursor-res=\n")
        cfgwriter.write("--full-scan-precursor-res-mz=\n")
        cfgwriter.write("--full-scan-product-res=\n")
        cfgwriter.write("--full-scan-product-res-mz=\n")
        cfgwriter.write("--full-scan-rt-filter-tolerance="+str(options.rt_tolerance)+"\n")
        pass
    elif "provided_file" in options.mode:
        cfgwriter.write("--full-scan-rt-filter-tolerance="+str(options.rt_tolerance)+"\n")

        pass
    else:
        print "I can\'t interpret the mode I\'m operating in!"
        sys.exit(1)
    
    cfgwriter.write("--import-fasta "+options.db+"\n")
    cfgwriter.close()

    for each_fraction_fldr in list_of_folders_with_ssl:
        os.chdir(basedir)
        shutil.copy("skyline_batch",each_fraction_fldr+"/skyline_batch")
        if options.protdb is not None:
            shutil.copy(options.protdb,each_fraction_fldr+"/background_proteome.protdb")
        
        os.chdir(each_fraction_fldr)
        skyline_filename=each_fraction_fldr+".sky"
        protdb_timestring=correctBackgroundProteome(skyline_filename)

        prependLine("skyline_batch","--in="+skyline_filename+"\n")
        #Let's go ahead and build a Blib library for the search results in this folder! EDIT: I THINK THIS IS HANDLED BY SKYLINE NOW.... COMMENTED OUT.
        
        list_of_ssl_files=[]

        for file in glob.glob("*.ssl"):
            list_of_ssl_files.append(file)
        #file_list="+".join(list_of_ssl_files)
        #blib_cmd="BlibBuild "+file_list+" -i "+each_fraction_fldr+" -a milkyway-galaxy"
        #split_blib_cmd=blib_cmd.split()
        #proc = subprocess.Popen(args=split_blib_cmd, shell=True, env=os.environ)
        #returncode = proc.wait()
        #if returncode != 0:
        #    raise Exception, "Program returned with non-zero exit code %d." % (returncode)
        #filter_cmd="BlibFilter "+each_fraction_fldr+".blib "+each_fraction_fldr+".filtered.blib"
        #split_filter_cmd=filter_cmd.split()
        #proc = subprocess.Popen(args=split_filter_cmd, shell=True, env=os.environ)
        #returncode = proc.wait()
        #if returncode != 0:
        #    raise Exception, "Program returned with non-zero exit code %d." % (returncode)


        single_cfgwriter=open("skyline_batch",'ab')

        for each_ssl_file in list_of_ssl_files:
            single_cfgwriter.write("--import-search-file="+file+"\n")


        if not "provided_file" in options.mode:
            shutil.copy(os.path.join(working_directory,getTemplateFile(options.mode,ms1_template,ms2_template,dia_template,options.skyline_file)),each_fraction_fldr+".sky")
        else:
            shutil.copy(os.path.join(working_directory,options.skyline_file),each_fraction_fldr+".sky")





        #final_cmd=skyline_runner_cmd+each_fraction_fldr+".sky "+skyline_runner_cmd_2

        single_cfgwriter.write("--reintegrate-create-model\n")
        single_cfgwriter.write("--reintegrate-model-name="+each_fraction_fldr+"\n")
        #single_cfgwriter.write("--reintegrate-annotate-scoring\n")
        single_cfgwriter.write("--reintegrate-overwrite-peaks\n")
        if options.msstats is not None:
            single_cfgwriter.write("--report-name=Wohl_MSstats_Input\n")
            single_cfgwriter.write("--report-file=MSstats_"+each_fraction_fldr+".csv\n")
        single_cfgwriter.write("--out=docker_protection_temporary.sky\n")
        #single_cfgwriter.write("--save\n")
        single_cfgwriter.close()
        final_cmd=skyline_runner_cmd+"skyline_batch"
        #final_cmd=final_cmd.split()
        proc = subprocess.Popen(args=final_cmd, shell=True, env=os.environ)
        returncode = proc.wait()
        if returncode != 0:
            raise Exception, "Program returned with non-zero exit code %d." % (returncode)
        os.remove(skyline_filename)
        shutil.copy("docker_protection_temporary.sky",skyline_filename)
        os.remove("docker_protection_temporary.sky")



    pass
else:
    #This is an unfractionated dataset.  We will import everything together into a single Skyline file.
    # REWRITE THIS SECTION FOR NON-FRAC DATASETS AND THEN TEST
    skyline_filename="combined_analysis.sky"

    if not "provided_file" in options.mode:
        shutil.copy(os.path.join(working_directory,getTemplateFile(options.mode,ms1_template,ms2_template,dia_template,options.skyline_file)),skyline_filename)
    else:
        if options.skyline_file is None:
            print "You need to give me a skyline file!"
            sys.exit(2)
        print "Copying file",str(os.path.join(working_directory,options.skyline_file))
        print "now just",skyline_filename
        print "Currently in folder",str(os.getcwd())
        shutil.copy(options.skyline_file,skyline_filename)
        pass
    #we'll make sure the background proteome is correct...
    protdb_timestring=correctBackgroundProteome(skyline_filename)


    if options.protdb is not None:
        print "Copying background protdb"
        shutil.copy(options.protdb,"background_proteome.protdb")
        print "Completed..."

    #if options.protdb is None:
    else:
        if options.background_db is not None:
            fasta_filename=options.background_db
            pass # GENERATE THE PROTDB FILE USING THE BACKGROUND PROVIDED BY options.background_db
        elif options.db is not None:
            fasta_filename=options.db
            print "USING THE FULL BACKGROUND AS ACCEPTED PROTEINS... USE WITH CAUTION!"
            pass #GENERATE PROTDB from options.db AS FALLBACK
        else:
            print "I need either a protDB background file, or a background database / database to work with..."
            sys.exit(2)
        if options.missedcleave is not None:
            missed=options.missedcleave
        else:
            print "I need you to give me a value for --missedcleave"
            sys.exit(2)
        ####### MAKE THE PROTDB NOW #######
        protdb_name="background_proteome"+protdb_timestring+".protdb"
        if options.enzyme == 'TrypsinP':
            enzyme_name = 'Trypsin/P'
            enzyme = '(?<=[RK])'
            enzyme_n_terminal_re = '(?=[RK])'
            enzyme_c_terminal_exception = None
            enzyme_c_terminal_re = '[ACDEFGHILMNPQSTVWY]$'
        elif options.enzyme == 'Trypsin':
            enzyme_name = 'Trypsin'
            enzyme = '(?<=[RK])(?=[^P])'
            enzyme_n_terminal_re = '(?=[RK])'
            enzyme_c_terminal_exception = 'P'
            enzyme_c_terminal_re = '[ACDEFGHILMNPQSTVWY]$'
        elif options.enzyme == 'Chymotrypsin':
            enzyme_name = 'Chymotrypsin'                    # Should follow the same name convention as Skyline does
            enzyme = '(?<=[FWYL])(?=[^P])'                     # Regular Expression of digestion rule
            enzyme_n_terminal_re = '(?=[FWYL])'                # N-terminal recognizing amino acids
            enzyme_c_terminal_exception = 'P'                # C-terminal amino acids that prevent the digestion
            enzyme_c_terminal_re = '[ACDEGHIKMNPQRSTV]$'    # All amino acids except N-terminal recognizing amino acids
        else:
            print "You have to choose a predefined enzyme."
            sys.exit(2)

        create_protdb(protdb_name)
        print "* protDB file " + protdb_name + " has been generated."
        insert_digestion_info(protdb_name, enzyme_name)
        print "* Digestion information has been inserted."
        insert_schema_version(protdb_name)
        print "* protDB schema version has been inserted."
        DbProteinName, DbProtein, DbDigestedPeptide, DbDigestedPeptideProtein = fasta_processing(fasta_filename, enzyme, missed, enzyme_n_terminal_re, enzyme_c_terminal_exception,enzyme_c_terminal_re)
        # Insert protein sequence, digested peptide sequence, and their relationship
        insert_protein_sequence(protdb_name, DbProtein)
        insert_peptide_sequence(protdb_name, DbDigestedPeptide)
        insert_relationship(protdb_name, DbDigestedPeptideProtein)

        #renameBackgroundProteome(protdb_timestring)

    #skyline_runner_cmd = "SkylineCmdSlow --batch-commands="
    #skyline_runner_cmd = "SkylineCmdSlow "
    #os.chdir(output_dir)
    cfgwriter=open("skyline_batch",'wb')
    cfgwriter.write("--timestamp\n")
    cfgwriter.write("--in="+skyline_filename+"\n")




    if options.add_mods:
        cfgwriter.write("--import-search-add-mods\n")
    if "dda_ms1" in options.mode:
        cfgwriter.write("--full-scan-precursor-res=\n")
        cfgwriter.write("--full-scan-precursor-res-mz=\n")
        cfgwriter.write("--full-scan-rt-filter-tolerance="+str(options.rt_tolerance)+"\n")
    elif "ms2_prm" in options.mode:
        cfgwriter.write("--full-scan-product-res=\n")
        cfgwriter.write("--full-scan-product-res-mz=\n")
        cfgwriter.write("--full-scan-rt-filter-tolerance="+str(options.rt_tolerance)+"\n")
        pass
    elif "dia" in options.mode:
        cfgwriter.write("--full-scan-precursor-res=\n")
        cfgwriter.write("--full-scan-precursor-res-mz=\n")
        cfgwriter.write("--full-scan-product-res=\n")
        cfgwriter.write("--full-scan-product-res-mz=\n")
        cfgwriter.write("--full-scan-rt-filter-tolerance="+str(options.rt_tolerance)+"\n")
        pass
    elif "provided_file" in options.mode:
        cfgwriter.write("--full-scan-rt-filter-tolerance="+str(options.rt_tolerance)+"\n")

        pass
    else:
        print "I can\'t interpret the mode I\'m operating in!"
        sys.exit(1)

    #list_of_ssl_files=[]

    if not options.makeblib:
        if os.path.isfile("filtered.blib"):
            print "USING FILTERED.BLIB"
            cfgwriter.write("--add-library-path=filtered.blib\n")
        elif os.path.isfile("combined_analysis.blib"):
            print "COULDN\'T FIND FILTERED. USING COMBINED ANALYSIS"
            cfgwriter.write("--add-library-path=combined_analysis.blib\n")
        elif os.path.isfile("combined_spectral_lib.blib"):
            cfgwriter.write("--add-library-path=combined_spectral_lib.blib\n")
        else:
            print "I can\'t find a BLIB spectral library to use! It\'s not named \'combined_spectral_lib.blib\' nor \'combined_analysis.blib\' ... rename it to one of these, or throw the --makeblib flag!"
            sys.exit(2)
    else:
        print "YOU ASKED TO MAKE BLIB FILES. KNOW THAT THIS IS NOT WHAT YOU WANT IF YOU\'RE DOING PTM ANALYSES, AS IT WILL INCLUDE TYPE 2\'S"
        for file in glob.glob("*.ssl"):
            cfgwriter.write("--import-search-file="+file+"\n")
    


    cfgwriter.write("--import-fasta="+options.db+"\n")

    if options.mprophet is not None and "decoy" in options.mprophet:

        cfgwriter.write("--out=docker_protection_temporary.sky\n")
        cfgwriter.close()
        
        with open('skyline_batch','rb') as batchreader:
            whole=batchreader.read().replace("\r"," ").replace("\n"," ")

        final_cmd="SkylineCmd "+whole
        #final_cmd="SkylineRunner "+whole

        final_cmd=final_cmd.split()
        print "execution 1:"
        print "About to execute command ",' '.join(final_cmd)

        notConnectedSkyline=True
        attempt=1
        while notConnectedSkyline:
            if attempt == 5:
                print "I can't get a connection to Skyline... something wrong!\nQuitting!~"
                sys.exit(2)
            time.sleep(60)
            print "MAKING ATTEMPT NUMBER ",str(attempt)
            proc = subprocess.Popen(args=final_cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)#, env=os.environ)
            #proc = subprocess.Popen("execute.bat")#, env=os.environ)
            output_communication=proc.communicate()[0]
            print output_communication
            returncode = proc.wait()
            if "Error: Could not connect to Skyline." in output_communication:
                attempt+=1
                continue
            if returncode != 0:
                raise Exception, "Program returned with non-zero exit code %d." % (returncode)
            else:
                notConnectedSkyline=False


        ######################### NOW I CALL TO BUILD THE DECOY PEPTIDES ############################
        time.sleep(30)
        #skyline_decoy_generator("combined_analysis.sky")
        print "Building Skyline Mass-Shifted Decoys..."
        skyline_decoy_generator("docker_protection_temporary.sky")
        print "finished adding decoy peptides!"
        shutil.move(skyline_filename,"ORIGINAL.sky")
        shutil.copy("docker_protection_temporary.sky",skyline_filename)
        os.remove("docker_protection_temporary.sky")
        time.sleep(30)

        ################################# END DECOYS #########################################
        #Now we're going to strip out all the results so that we can reimport with the decoys!
        strip_cmd="SkylineCmd --in="+skyline_filename+" --remove-all=. --out=docker_protection_temporary.sky"
        notConnectedSkyline=True
        attempt=1
        while notConnectedSkyline:
            if attempt == 5:
                print "I can't get a connection to Skyline... something wrong!\nQuitting!~"
                sys.exit(2)
            time.sleep(60)
            print "MAKING ATTEMPT NUMBER ",str(attempt)
            proc = subprocess.Popen(args=strip_cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)#, env=os.environ)
            #proc = subprocess.Popen("execute.bat")#, env=os.environ)
            output_communication=proc.communicate()[0]
            print output_communication
            returncode = proc.wait()
            if "Error: Could not connect to Skyline." in output_communication:
                attempt+=1
                continue
            if returncode != 0:
                raise Exception, "Program returned with non-zero exit code %d." % (returncode)
            else:
                notConnectedSkyline=False

        time.sleep(20)

        os.remove(skyline_filename)
        shutil.copy("docker_protection_temporary.sky",skyline_filename)
        os.remove("docker_protection_temporary.sky")

        shutil.move("skyline_batch","before_decoys")

        cfgwriter=open("skyline_batch",'wb')
        cfgwriter.write("--timestamp\n")
        cfgwriter.write("--in="+skyline_filename+"\n")
        #cfgwriter.write("--remove-all=.\n")

    #cfgwriter.write("--out=docker_protection_temporary.sky\n")


    #We're going to add the decoy peptides, save out the document, and after that we'll import the data in parallel and finally handle mProphet modeling
    if options.mprophet is not None and "decoy" not in options.mprophet:
        if "reverse" in options.mprophet:
            cfgwriter.write("--decoys-add=reverse\n")
        elif "shuffle" in options.mprophet:
            cfgwriter.write("--decoys-add=shuffle\n")
        else:
            pass


    #Now we'll have to save the document and then we can reopen it to import in parallel
    ###########################~~~~
        cfgwriter.write("--out=docker_protection_temporary.sky\n")
        cfgwriter.close()
        print "Setting up the file so that we can use parallel import!"
        with open('skyline_batch','rb') as batchreader:
            whole=batchreader.read().replace("\r"," ").replace("\n"," ")

        setup_cmd="SkylineCmd "+whole

        setup_cmd=setup_cmd.split()
        print "execution 2:"
        print "About to execute command ",' '.join(setup_cmd)

        notConnectedSkyline=True
        attempt=1
        while notConnectedSkyline:
            if attempt == 5:
                print "I can't get a connection to Skyline... something wrong!\nQuitting!~"
                sys.exit(2)
            time.sleep(60)
            print "MAKING ATTEMPT NUMBER ",str(attempt)
            proc = subprocess.Popen(args=setup_cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)#, env=os.environ)
            output_communication=proc.communicate()[0]
            print output_communication
            returncode = proc.wait()
            if "Error: Could not connect to Skyline." in output_communication:
                attempt+=1
                continue
            if returncode != 0:
                raise Exception, "Program returned with non-zero exit code %d." % (returncode)
            else:
                notConnectedSkyline=False
        
        os.remove(skyline_filename)
        shutil.copy("docker_protection_temporary.sky",skyline_filename)
        os.remove("docker_protection_temporary.sky")


        os.rename("skyline_batch","skyline_batch_setup")
        cfgwriter=open("skyline_batch",'wb')
        cfgwriter.write("--timestamp\n")
        cfgwriter.write("--in="+skyline_filename+"\n")


    #Now that the file has had all the targets set up, and resaved, the file should now be ready for import with parallelism.
    if options.num_procs is not None:
        cfgwriter.write("--import-process-count={0}\n".format(options.num_procs))

    
    if options.mprophet is not None:
        if "secondbest" in options.mprophet:
            cfgwriter.write("--import-all=.\n")
            cfgwriter.write("--reintegrate-model-name=2NDBESTcombinedanalysis"+str(datetime.datetime.now().time()).replace(":",'').replace('.','')+"\n")
            cfgwriter.write("--reintegrate-create-model\n")
            cfgwriter.write("--reintegrate-model-second-best\n")
            #cfgwriter.write("--reintegrate-annotate-scoring\n")
            cfgwriter.write("--reintegrate-overwrite-peaks\n")

        elif "reverse" in options.mprophet:
            #cfgwriter.write("--decoys-add=reverse\n")
            cfgwriter.write("--import-all=.\n")
            cfgwriter.write("--reintegrate-model-name=REVERSEcombinedanalysis"+str(datetime.datetime.now().time()).replace(":",'').replace('.','')+"\n")
            cfgwriter.write("--reintegrate-create-model\n")
            #cfgwriter.write("--reintegrate-annotate-scoring\n")
            cfgwriter.write("--reintegrate-overwrite-peaks\n")

        elif "shuffle" in options.mprophet:
            #cfgwriter.write("--decoys-add=shuffle\n")
            cfgwriter.write("--import-all=.\n")
            cfgwriter.write("--reintegrate-model-name=SHUFFLEcombinedanalysis"+str(datetime.datetime.now().time()).replace(":",'').replace('.','')+"\n")
            cfgwriter.write("--reintegrate-create-model\n")
            #cfgwriter.write("--reintegrate-annotate-scoring\n")
            cfgwriter.write("--reintegrate-overwrite-peaks\n")

        elif "decoy" in options.mprophet:
            cfgwriter.write("--import-all=.\n")
            cfgwriter.write("--reintegrate-model-name=DECOYcombinedanalysis"+str(datetime.datetime.now().time()).replace(":",'').replace('.','')+"\n")
            cfgwriter.write("--reintegrate-create-model\n")
            #cfgwriter.write("--reintegrate-annotate-scoring\n")
            cfgwriter.write("--reintegrate-overwrite-peaks\n")
            #print "DECOYS DO NOT WORK YET!"
            #cfgwriter.close()
            #sys.exit(2)
        else:
            pass
            cfgwriter.write("--import-all=.\n")


    else:
        pass
        cfgwriter.write("--import-all=.\n")

    

    #cfgwriter.write("--out=docker_protection_temporary.sky\n")
    cfgwriter.close()

    single_cfgwriter=open("skyline_batch",'ab')



    #final_cmd=skyline_runner_cmd+each_fraction_fldr+".sky "+skyline_runner_cmd_2

    if options.msstats is not None:
        shutil.copy("C:\\skyline\\WOHL_MSSTATS_REPORT.skyr","WOHL_MSSTATS_REPORT.skyr")
        single_cfgwriter.write("--report-add=WOHL_MSSTATS_REPORT.skyr\n")
        #single_cfgwriter.write("--report-add=\"/galaxy-central/tools/wohl-proteomics/wohl_skyline/WOHL_MSSTATS_REPORT.skyr\"\n")
        single_cfgwriter.write("--report-conflict-resolution=overwrite\n")
        single_cfgwriter.write("--report-name=Wohl_MSstats_Input\n")
        single_cfgwriter.write("--report-file=MSstats_input.csv\n")


    single_cfgwriter.write("--chromatogram-file={0}\n".format(options.chromatograms))
    single_cfgwriter.write("--chromatogram-precursors\n")
    single_cfgwriter.write("--chromatogram-products\n")
    single_cfgwriter.write("--chromatogram-base-peaks\n")
    single_cfgwriter.write("--chromatogram-tics\n")
    #shutil.copy("C:\\Users\\JamesLab\\Skyline\\peak_boundaries.skyr","peak_boundaries.skyr")
    #single_cfgwriter.write("--report-add=peak_boundaries.skyr\n")
    #single_cfgwriter.write("--report-name=Peak_Boundaries\n")
    #single_cfgwriter.write("--report-file={0}\n".format("peak_boundaries.csv"))#options.peak_boundaries))
    single_cfgwriter.write("--out=docker_protection_temporary.sky\n")




    single_cfgwriter.close()
    #final_cmd=skyline_runner_cmd+"skyline_batch"

    with open('skyline_batch','rb') as batchreader:
        whole=batchreader.read().replace("\r"," ").replace("\n"," ")

    final_cmd="SkylineCmd "+whole
    #final_cmd="SkylineRunner "+whole

    final_cmd=final_cmd.split()
    print "execution 3:"
    print "About to execute command ",' '.join(final_cmd)

    notConnectedSkyline=True
    attempt=1
    while notConnectedSkyline:
        if attempt == 5:
            print "I can't get a connection to Skyline... something wrong!\nQuitting!~"
            sys.exit(2)
        time.sleep(60)
        print "MAKING ATTEMPT NUMBER ",str(attempt)
        proc = subprocess.Popen(args=final_cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)#, env=os.environ)
        #proc = subprocess.Popen("execute.bat")#, env=os.environ)
        output_communication=proc.communicate()[0]
        print output_communication
        returncode = proc.wait()
        if "Error: Could not connect to Skyline." in output_communication:
            attempt+=1
            continue
        if returncode != 0:
            raise Exception, "Program returned with non-zero exit code %d." % (returncode)
        else:
            notConnectedSkyline=False
        
        #stdout, stderr = subprocess.Popen("cmd.exe", stdout=subprocess.PIPE,startupinfo=sysuser).communicate()

    os.remove(skyline_filename)
    shutil.copy("docker_protection_temporary.sky",skyline_filename)
    os.remove("docker_protection_temporary.sky")


    shutil.copy("C:\\skyline\\peak_boundaries.skyr","peak_boundaries.skyr")

    peak_cfgwriter=open("peak_skyline_batch",'ab')
    peak_cfgwriter.write("--timestamp\n")
    peak_cfgwriter.write("--in="+skyline_filename+"\n")
    peak_cfgwriter.write("--report-conflict-resolution=overwrite\n")
    peak_cfgwriter.write("--report-add=peak_boundaries.skyr\n")
    peak_cfgwriter.write("--report-name=Peak_Boundaries\n")
    peak_cfgwriter.write("--report-file={0}\n".format(options.peak_boundaries))
    peak_cfgwriter.close()

    with open('peak_skyline_batch','rb') as batchreader:
        whole=batchreader.read().replace("\r"," ").replace("\n"," ")

    final_cmd="SkylineCmd "+whole

    final_cmd=final_cmd.split()
    print "execution 4:"
    print "About to execute command ",' '.join(final_cmd)


    notConnectedSkyline=True
    attempt=1
    while notConnectedSkyline:
        if attempt == 5:
            print "I can't get a connection to Skyline... something wrong!\nQuitting!~"
            sys.exit(2)
        time.sleep(60)
        print "MAKING ATTEMPT NUMBER ",str(attempt)
        proc = subprocess.Popen(args=final_cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)#, env=os.environ)
        #proc = subprocess.Popen("execute.bat")#, env=os.environ)
        output_communication=proc.communicate()[0]
        print output_communication
        returncode = proc.wait()
        if "Error: Could not connect to Skyline." in output_communication:
            attempt+=1
            continue
        if returncode != 0:
            raise Exception, "Program returned with non-zero exit code %d." % (returncode)
        else:
            notConnectedSkyline=False

    pass



if options.msstats is not None:
    if options.fractions:
        pass#put the MSstats gathering stuff here for fractionated data
    else:
        os.chdir(basedir)
        if options.exp_file is not None:
            msstats_df=pandas.read_csv("MSstats_input.csv")
            exp_structure=pandas.read_csv(options.exp_file,sep="\t",low_memory=False)
            exp_structure['Original File Name']=exp_structure['Original File Name']+".mzML"
            exp_structure.set_index(keys='Original File Name',inplace=True)
            name_to_biocond={}
            group_to_biorep={}
            #k=1
            for index,each in exp_structure.iterrows():
                #print each,"each"
                #print index,"index"
                name_to_biocond[index]=each['Biological Condition']
                group_to_biorep[each['Fractionation Group ID String']]=each['BioReplicate']
                #k+=1
            print "name to biocond",name_to_biocond
            #for eachkey in name_to_biocond:
            #    msstats_df.loc[msstats_df['File Name'] == eachkey, 'Condition'] = name_to_biocond[eachkey]
            msstats_df['Condition']=msstats_df['File Name'].map(name_to_biocond)
            #    #filter=msstats_df[msstats_df['File Name']==eachkey]
            #    #filter['Condition']=name_to_biocond[eachkey]
            #msstats_df['BioReplicate']=msstats_df['File Name'].map(group_to_biorep)
            print "group to biorep",group_to_biorep
            for eachgroup in group_to_biorep:
                msstats_df.loc[msstats_df['File Name'].str.contains(eachgroup), 'BioReplicate'] = group_to_biorep[eachgroup]
            #    #filter=msstats_df[msstats_df['File Name'].str.contains(eachgroup)]
            #    #filter['BioReplicate']=group_to_biorep[eachgroup]
            os.remove("MSstats_input.csv")
            msstats_df.to_csv("MSstats_input.csv",index=False)
        
        shutil.move("MSstats_input.csv",options.msstats)
        


        #Done with the *.raw and *.mzML files...
        delete_list = [ f for f in os.listdir(".") if f.endswith(".raw") or f.endswith(".mzML") ]
        for f in delete_list:
            os.remove(f)
            print "removing ",f
        #shutil.copy("peak_boundaries.csv",options.peak_boundaries)

#MOVE ThE OUTPUT TARGZ
os.chdir(basedir)



#zip_cmd="7z a -ttar -so skyline_output"+str(datetime.datetime.now()).replace(" ","_").replace(":","_")+".tar ssl_files/ | 7z a -si OUTPUTARCHIVE.tar.gz"
#zip_cmd="7z a -ttar -so skyline_output"+str(datetime.datetime.now()).replace(" ","_").replace(":","_")+".tar * | 7z a -si OUTPUTARCHIVE.tar.gz"
os.chdir("..")
zip_cmd="7z a -ttar -so skyline_output"+str(datetime.datetime.now()).replace(" ","_").replace(":","_")+".tar ssl_files/ | pigz.exe > OUTPUTARCHIVE.tar.gz"
subprocess.call(args=zip_cmd,shell=True)

shutil.move("OUTPUTARCHIVE.tar.gz",options.output)

print "-------------------------------------------------------------------------"
print "We're all done here."
print "-------------------------------------------------------------------------"
sys.exit()

