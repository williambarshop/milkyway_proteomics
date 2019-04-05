import warnings
warnings.filterwarnings('ignore') #because of plotly in pymzml
import os, sys, pandas, pyteomics, glob, numpy, subprocess, shutil, re, optparse
from itertools import izip_longest
from subprocess import Popen, STDOUT
from array import array
import numpy
import sqlite3, zlib
from pymzml import run

#####################################
#This is the script to run Specter on MilkyWay SSL converter outputs.
#
#
#VERSION 0.05
version="0.05"
#DATE: 4/4/2019
date="4/4/2019"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the Specter wrapper for Galaxy, Wohlschlegel Lab UCLA"
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


#commands=["sh SpecterStandalone.sh {0} {1} 0 start {2} ".format(each_run.replace(".mzML",""),spectral_library,num_processors,instrument_type,search_fragment_mass_tolerance_in_ppm) for each_run in glob.glob("*.mzML")]

#os.path.join(jar_path,"MSPLIT-DIAv"+search_version+".jar")

parser = optparse.OptionParser()
#parser.add_option("--sslgz",action="store",type="string",dest="sslgz")
parser.add_option("--mzmlfiles",action="store",type="string",dest="dda_mzml_files")
parser.add_option("--diafiles",action="store",type="string",dest="diafiles")
parser.add_option("--psmtsv",action="store",type="string",dest="psm_tsv")
parser.add_option("--ms2searchppm",action="store", type="float", dest="ms2searchppm")
parser.add_option("--diafolder",action="store",type="string",dest="diamzmlfolder")
parser.add_option("--fdrtofilter",action="store", type="float", dest="fdrtofilter")
parser.add_option("--num_processors",action="store", type="int", dest="num_processors")
parser.add_option("--tool_path",action="store",type="string",dest="tool_path")
parser.add_option("--instrument_type",action="store",type="string",dest="instrument_type")

(options,args) = parser.parse_args()

spectral_library="filtered.blib" #This, is by definition, the spectral library to use in MilkyWay from the SSL Converter.
specter_path=os.path.join(options.tool_path,"Specter/")

print "The Specter path is {0}".format(specter_path)
num_processors=options.num_processors
instrument_type=options.instrument_type
fdr_to_filter=options.fdrtofilter

search_fragment_mass_tolerance_in_ppm=options.ms2searchppm
dia_folder=options.diamzmlfolder


basedir=os.getcwd()


#Alright, we'll start by unzipping the ssl file tar.gz
#shutil.copy(options.sslgz,"ssl.tar.gz")
#subprocess.check_output("tar xvf {0}".format("ssl.tar.gz"))

os.chdir(os.path.join(basedir,"ssl_files"))
dia_files=[]
os.mkdir(os.path.join(basedir,"dia_temp"))
dia_temp_folder=os.path.join(basedir,"dia_temp")
for each_dia_file in glob.glob("*.mzML"):
    os.rename(each_dia_file,os.path.join(dia_temp_folder,each_dia_file))
    dia_files.append(os.path.join(dia_temp_folder,each_dia_file))
for each_file in options.dda_mzml_files.split(","):
    os.rename(os.path.join(basedir,each_file),os.path.join(os.getcwd(),each_file))

ssl_files=glob.glob("*.ssl")
ssl_dfs=[]
for each_ssl in ssl_files:
    ssl_dfs.append(pandas.read_csv(each_ssl,sep="\t"))
ssl_df=pandas.concat(ssl_dfs)
del ssl_dfs


#Read the psm table, but only the columns we need...
psm_df=pandas.read_csv(options.psm_tsv,sep="\t",usecols=['file','scan','sequence','spectrum precursor m/z','charge','protein id'])

#Let's load all the mzML files we need.

#THIS IS WHEN WE BRING THE DIA MZML FILES IN!
for each_run in ssl_df['file'].unique():
    os.unlink(each_run)

for each_dia_mzml in dia_files:
    os.rename(each_dia_mzml,os.path.join(os.getcwd(),each_dia_mzml.rsplit("/",1)[1]))

parallel_limit=num_processors
commands=["python {0}Specter_Standalone.py {1} {2} 0 start {3} {4} {5}".format(specter_path,each_run.replace(".mzML",""),spectral_library.replace(".blib",""),num_processors,instrument_type,search_fragment_mass_tolerance_in_ppm) for each_run in glob.glob("*.mzML")]
groups = [(Popen(cmd, stderr=subprocess.STDOUT,shell=True)
          for cmd in commands)] * parallel_limit

print "commands are:",commands
for processes in izip_longest(*groups):
    for p in filter(None, processes):
        print "waiting..."
        p.wait()

error_codes={}


for each_msplit_output in glob.glob("*_msplit_search.out"):
    print "working on",each_msplit_output
    filtered_output=each_msplit_output.replace("_search.out","_filtered.out")
    print "running","java -Xmx{0}G -cp {1} UI.SWATHFilter -r {2} -o {3} -fdr {4} -rt {5}".format(max_gigs_mem_per_run,os.path.join(jar_path,"MSPLIT-DIAv"+search_version+".jar"),each_msplit_output,filtered_output,fdr_to_filter,"1" if filter_by_rt else "0")
    
    error_codes[each_msplit_output.replace("_msplit_search.out",".mzML")]=os.system("java -Xmx{0}G -cp {1} UI.SWATHFilter -r {2} -o {3} -fdr {4} -rt {5}".format(max_gigs_mem_per_run,os.path.join(jar_path,"MSPLIT-DIAv"+search_version+".jar"),each_msplit_output,filtered_output,fdr_to_filter,"1" if filter_by_rt else "0"))
#for each_key in error_codes:
#    with open(each_key.replace(".mzML","_msplit_filter_log.txt"),'wb') as logwriter:
#        logwriter.write(error_codes[each_key])

#final_rt_filtered=glob.glob("*rtFiltered.txt")




specter_final_outputs=[]
for each_file in final_rt_filtered:
    specter_final_outputs.append(pandas.read_csv(each_file,sep='\t',index_col=False))

specter_df=pandas.concat(specter_final_outputs)
specter_df['Peptide']=specter_df['Peptide'].apply(lambda x: re.sub(r'([+-]\d+\.\d+)',r'[\1]',x))

decimal_extractor=re.compile(r"(\d+.\d+)")

def rounder(match): #Skyline rounds to a single decimal for mods...
    return"{:.1f}".format(float(match.group()))

specter_df['Peptide']=specter_df['Peptide'].apply(lambda x: re.sub(decimal_extractor,rounder,x))
#specter_df['Peptide']=specter_df['Peptide'].apply(lambda x: re.sub(r'([+-]\d+\.\d+)',r'[\1]',x))

unique_peptides=specter_df['Peptide'].unique()

spectral_lib="filtered.blib"

for each_run in specter_df['#File'].unique():
    this_run_df=specter_df[specter_df['#File']==each_run]
    z=0
    #with mzml.MzML(os.path.join(basedir,"ssl_files",each_run),use_index=True,iterative=False) as this_mzml: #build_id_cache=True ,retrieve_refs=True
    #with mzml.read(os.path.join(basedir,"ssl_files",each_run),use_index=False,iterative=False) as this_mzml: #build_id_cache=True ,retrieve_refs=True
    #with mzml.PreIndexedMzML(os.path.join(basedir,"ssl_files",use_index=True,each_run),iterative=False) as this_mzml: #build_id_cache=True ,retrieve_refs=True
    mzml_reader=run.Reader(os.path.join(basedir,"ssl_files",each_run))
    #print mzml_reader[5].peaks
    print "Finished loading the mzML file:",each_run
    #tqdm.write("Finished loading the mzML file: "+each_run)
    with open(each_run.rsplit(".",1)[0]+"_specter.ssl",'wb') as ssl_writer:
        ssl_writer.write("file\tscan\tcharge\tsequence\n")
        with open(each_run.rsplit(".",1)[0]+".ms2",'wb') as ms2_writer:
            ms2_writer.write("H\tSource file\t{0}\n".format(each_run))
            for each_spec_tuple_key in decompressed_spectra:
                each_spec_tuple=decompressed_spectra[each_spec_tuple_key]
                this_pep=spec_mapping[each_spec_tuple[1][0]]
                #print "this pep is",this_pep
                #break
                this_pep_run_df=this_run_df[this_run_df['seq_charge']==this_pep]
                if this_pep not in run_pepcharges[each_run]:
                    continue
                spectral_ms2_string=""
                zipped_data=each_spec_tuple[0]
                for each in zipped_data:
                    spectral_ms2_string+="{0}\t{1}\n".format(each[0],each[1])

                #Finish this selection for the seq_charge
                for index, each_ID in this_pep_run_df.iterrows():
                    z+=1

                    ms2_writer.write("S\t{0}\t{0}\t{1}\n".format(z,spec_mapping_precursor_mass[this_pep]))
                    #Insert RT here, by taking scan# and reading out the mzML file absolute RT...
                    #ms2_writer.write("I\tRTime\t{0}\n".format(this_mzml["controllerType=0 controllerNumber=1 scan={0}".format(each_ID['Scan#'])]["scanList"]["scan"][0]['scan start time']))
                    #ms2_writer.write("I\tRTime\t{0}\n".format(this_mzml.get_by_id("controllerType=0 controllerNumber=1 scan={0}".format(each_ID['Scan#']))["scanList"]["scan"][0]['scan start time']))
                    ms2_writer.write("I\tRTime\t{0}\n".format(mzml_reader[each_ID['Scan#']]['scan start time']))
                    ms2_writer.write("Z\t{0}\t{1}\n".format(this_pep.split("_")[1],(float(spec_mapping_precursor_mass[this_pep])*float(this_pep.split("_")[1]))-((float(this_pep.split("_")[1])-1)*1.007276)   ))
                    ms2_writer.write(spectral_ms2_string)

                    ssl_writer.write("{0}\t{1}\t{2}\t{3}\n".format(each_run.replace(".mzML",".ms2"),z,this_pep.split("_")[1],this_pep.split("_")[0]))


#print "{0} -a milkyway.auth *_specter.ssl specter_filtered.blib".format(os.path.join(jar_path,"blibbuild\\BlibBuild.exe"))
os.system("{0} -a milkyway-galaxy *_specter.ssl combined_spectral_lib.blib 2>&1".format("wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibBuild.exe"))
os.system("{0} combined_spectral_lib.blib filtered.blib 2>&1".format("wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibFilter.exe"))

print "All done!"
sys.exit(0)
