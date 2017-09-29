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
#This is the script to run MSPLIT-DIA on MilkyWay SSL converter outputs.
#
#
#VERSION 0.05
version="0.05"
#DATE: 9/28/2017
date="9/28/2017"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the MSPLIT-DIA wrapper for Galaxy, Wohlschlegel Lab UCLA"
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

parser = optparse.OptionParser()
#parser.add_option("--sslgz",action="store",type="string",dest="sslgz")
parser.add_option("--mzmlfiles",action="store",type="string",dest="dda_mzml_files")
parser.add_option("--diafiles",action="store",type="string",dest="diafiles")
parser.add_option("--psmtsv",action="store",type="string",dest="psm_tsv")
parser.add_option("--ms2daltons",action="store", type="float", dest="ms2daltons")
parser.add_option("--cyclescans",action="store", type="int", dest="cyclescans")
parser.add_option("--maxgigsmemory",action="store", type="float", dest="maxgigsmemory")
parser.add_option("--ms2searchppm",action="store", type="float", dest="ms2searchppm")
parser.add_option("--daltonparenttolerance",action="store", type="float", dest="search_parent_mass_tolerance_da")
parser.add_option("--diafolder",action="store",type="string",dest="diamzmlfolder")
parser.add_option("--variablewindowfile",action="store",type="string",dest="variablewindowfile")
parser.add_option("--ms1start",action="store", type="float", dest="ms1start")
parser.add_option("--ms1end",action="store", type="float", dest="ms1end")
parser.add_option("--fdrtofilter",action="store", type="float", dest="fdrtofilter")
parser.add_option("--filterbyrt",action="store_true", dest="filterbyrt")
parser.add_option("--tool_path",action="store",type="string",dest="tool_path")

search_version="07192015"
version="02102015"
(options,args) = parser.parse_args()


jar_path="{0}".format(options.tool_path)
variable_windows_file=options.variablewindowfile
ms1_start=options.ms1start
ms1_end=options.ms1end
fdr_to_filter=options.fdrtofilter
filter_by_rt=options.filterbyrt


num_scans_per_cycle=options.cyclescans
max_gigs_mem_per_run=options.maxgigsmemory
search_fragment_mass_tolerance_in_ppm=options.ms2searchppm
search_patent_mass_tolerance_da=options.search_parent_mass_tolerance_da
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

#Let's load all the mzML files we need.  We'll do it one at a time so that we 
output_libraries=[]
for each_mzml in ssl_df['file'].unique():
    #print each_mzml
    write_spectra=[]
    subset_psm=psm_df[psm_df['file']==each_mzml]
    subset_psm=subset_psm[subset_psm['scan'].isin(ssl_df[ssl_df['file']==each_mzml]['scan'].tolist())]
    with open(each_mzml.rsplit(".",1)[0]+"_msplit.mgf",'wb') as mgf_writer:
        output_libraries.append(each_mzml.rsplit(".",1)[0]+"_msplit.mgf")
        this_mzml=run.Reader(os.path.join(basedir,"ssl_files",each_mzml))
        #with mzml.MzML(os.path.join(basedir,"ssl_files",each_mzml),build_id_cache=True,use_index=True,iterative=False,retrieve_refs=True) as this_mzml:
        for index,each_scan in subset_psm.iterrows():
            this_spectrum = this_mzml[each_scan['scan']]
            psm_row=subset_psm[subset_psm['scan']==each_scan['scan']].iloc[0]
            psm_charge=psm_row['charge']
            psm_mz=psm_row['spectrum precursor m/z']
            #print psm_charge
            #print psm_mz
            #print psm_row
            #break
            #print this_spectrum['precursorList']

            #print psm_df.loc[numpy.logical_and(psm_df['file']==each_scan['file'],psm_df['scan']==each_scan['scan'])]['charge'].iloc[0]

            mgf_writer.write("BEGIN IONS\n")
            mgf_writer.write("TITLE=Scan Number: {0} Retention Time: PT{1}S PROTEIN: {2}\n".format(each_scan['scan'],str(this_spectrum['scan start time']*60.0),each_scan['protein id'].split(",")[0]))
            mgf_writer.write("CHARGE={0}{1}\n".format(psm_charge,'+' if psm_charge>0 else "-"  ))
            mgf_writer.write("PEPMASS={0}\n".format(psm_mz))
            mgf_writer.write("SEQ={0}\n".format(each_scan['sequence'].replace("[","").replace("]","")))
            mgf_writer.write("RTINSECONDS={0}\n".format(str(this_spectrum['scan start time']*60.0)))
            mgf_writer.write("SCAN={0}\n".format(each_scan['scan']))
            for each_tuple in this_spectrum.peaks:#zip(this_spectrum['m/z array'],this_spectrum.i):
                mass=each_tuple[0]
                intensity=each_tuple[1]
                mgf_writer.write("{0} {1}\n".format(mass,intensity))
            mgf_writer.write("END IONS\n")


fragment_mass_tolerance_in_da=options.ms2daltons
decoy_libraries=[]
for each_library in output_libraries:
    print "Processing ",each_library," to add decoys..."
    print subprocess.check_output("java -cp "+os.path.join(jar_path,"MSPLIT-DIAv"+version+".jar")+" org.Spectrums.DecoySpectrumGenerator {0} {1} {2}".format(each_library,each_library.rsplit(".",1)[0]+"_wdecoys.mgf",str(fragment_mass_tolerance_in_da)))
    decoy_libraries.append(each_library.rsplit(".",1)[0]+"_wdecoys.mgf")


#We have to join each decoy lib together into a megalibrary.
os.system("cat *msplit_wdecoys.mgf > final_msplit_lib.mgf")


#THIS IS WHEN WE BRING THE DIA MZML FILES IN!
for each_run in ssl_df['file'].unique():
    os.unlink(each_run)
#
#
#BRING IN THE DIA FILES HERE
#
#
#
#for each_dia_mzml in options.diafiles.split(","):#glob.glob("{0}/*.mzML".format(dia_folder)):
#    os.rename(os.path.join(basedir,each_dia_mzml),os.path.join(os.getcwd(),each_dia_mzml))
for each_dia_mzml in dia_files:
    os.rename(each_dia_mzml,os.path.join(os.getcwd(),each_dia_mzml.rsplit("/",1)[1]))
with open(variable_windows_file,'rb') as freader:
    flines=freader.read().replace("\r","").split("\n")
    flines_fixed=[]
    for each_line in flines:
        if each_line=="":
            continue
        this_line=['MS2']
        this_line.extend(each_line.replace(",",'\t').split('\t'))
        flines_fixed.append('\t'.join(this_line)+"\n")
    final_flines_fixed=['#Scan\twindowBegin\twindowEnd\n','MS1\t{0}\t{1}\n'.format(ms1_start,ms1_end)]
    final_flines_fixed.extend(flines_fixed)

variable_windows_final=variable_windows_file.rsplit('.',1)[0]+"_msplit_windows.tsv"
with open(variable_windows_file.rsplit('.',1)[0]+"_msplit_windows.tsv",'wb') as fwriter:
    for each_line in final_flines_fixed:
        fwriter.write(each_line)

parallel_limit=8
commands = ["java -Xmx{0}G -jar {1} {2} {3} {4} {5} {6} {7} {8}".format(max_gigs_mem_per_run,os.path.join(jar_path,"MSPLIT-DIAv"+search_version+".jar"),search_patent_mass_tolerance_da,search_fragment_mass_tolerance_in_ppm,num_scans_per_cycle,each_run,"final_msplit_lib.mgf",each_run.rsplit(".",1)[0]+"_msplit_search.out",str(variable_windows_final) if variable_windows_file!=None else "") for each_run in glob.glob("*.mzML")]
groups = [(Popen(cmd, stderr=subprocess.STDOUT)
          for cmd in commands)] * parallel_limit

print "commands are:",commands
for processes in izip_longest(*groups):
    for p in filter(None, processes):
        print "waiting..."
        p.wait()

output_logs={}
for each_msplit_output in glob.glob("*_msplit_search.out"):
    print "working on",each_msplit_output
    filtered_output=each_msplit_output.replace("_search.out","_filtered.out")
    print "running","java -Xmx{0}G -cp {1} UI.SWATHFilter -r {2} -o {3} -fdr {4} -rt {5}".format(max_gigs_mem_per_run,os.path.join(jar_path,"MSPLIT-DIAv"+search_version+".jar"),each_msplit_output,filtered_output,fdr_to_filter,"1" if filter_by_rt else "0")
    
    output_logs[each_msplit_output.replace("_msplit_search.out",".mzML")]=subprocess.check_output("java -Xmx{0}G -cp {1} UI.SWATHFilter -r {2} -o {3} -fdr {4} -rt {5}".format(max_gigs_mem_per_run,os.path.join(jar_path,"MSPLIT-DIAv"+search_version+".jar"),each_msplit_output,filtered_output,fdr_to_filter,"1" if filter_by_rt else "0"))
for each_key in output_logs:
    with open(each_key.replace(".mzML","_msplit_filter_log.txt"),'wb') as logwriter:
        logwriter.write(output_logs[each_key])

final_rt_filtered=glob.glob("*rtFiltered.txt")

msplit_final_outputs=[]
for each_file in final_rt_filtered:
    msplit_final_outputs.append(pandas.read_csv(each_file,sep='\t',index_col=False))

msplit_df=pandas.concat(msplit_final_outputs)
msplit_df['Peptide']=msplit_df['Peptide'].apply(lambda x: re.sub(r'([+-]\d+\.\d+)',r'[\1]',x))

decimal_extractor=re.compile(r"(\d+.\d+)")

def rounder(match): #Skyline rounds to a single decimal for mods...
    return"{:.1f}".format(float(match.group()))

msplit_df['Peptide']=msplit_df['Peptide'].apply(lambda x: re.sub(decimal_extractor,rounder,x))
#msplit_df['Peptide']=msplit_df['Peptide'].apply(lambda x: re.sub(r'([+-]\d+\.\d+)',r'[\1]',x))

unique_peptides=msplit_df['Peptide'].unique()

spectral_lib="filtered.blib"

blib_connection=sqlite3.connect(spectral_lib)


blib_cursor=blib_connection.cursor()

blib_cursor.execute("SELECT * FROM RefSpectra")
refSpectra=blib_cursor.fetchall()

#Before we get the refSpectraPeaks we need, we're going to gather a list of the peptide/charges we need for
#each run.
run_pepcharges={}
msplit_df['seq_charge']=msplit_df['Peptide']+"_"+msplit_df['z'].astype(str)
for each_run in msplit_df['#File'].unique():
    run_pepcharges[each_run]=[]
    run_view=msplit_df[msplit_df["#File"]==each_run]
    no_duplicates=run_view['seq_charge'].drop_duplicates().tolist()
    run_pepcharges[each_run].extend(no_duplicates)
    #print run_pepcharges[each_run]
    #break

spec_mapping={}
spec_mapping_precursor_mass={}
for each_spectrum in refSpectra:
    spec_mapping[each_spectrum[0]]=each_spectrum[4]+"_"+str(each_spectrum[3])
    spec_mapping_precursor_mass[each_spectrum[4]+"_"+str(each_spectrum[3])]=each_spectrum[2]
    #print each_spectrum
    #break

#NOTE: BlibBuild will store the raw byte array in zlib compression didn't do anything... So we'll add in
#some excepts to handle those situations...
blib_cursor.execute("SELECT * FROM RefSpectraPeaks")
refSpectraPeaks=blib_cursor.fetchall()
#for each_spectra_set in blib_cursor.fet(50000):
    ##print type(each_spectra_set)
    
    
    
decompressed_spectra={}
for i in refSpectraPeaks:
    this_pep=spec_mapping[i[0]]
    if this_pep not in run_pepcharges[each_run]:
        continue

    try:
        mzarray=numpy.fromstring(zlib.decompress(str(i[1])))
    except:
        mzarray=numpy.array(array('d',str(i[1])))
    try:
        intensityarray=numpy.fromstring(str(i[2]),dtype=numpy.float32)
        if numpy.isnan(intensityarray).any() or not (intensityarray>=0).all():
            intensityarray=numpy.fromstring(zlib.decompress(str(i[2])),dtype=numpy.float32)
            #break
    except:
        intensityarray=numpy.fromstring(zlib.decompress(str(i[2])),dtype=numpy.float32)
        if numpy.isnan(intensityarray).any() or not (intensityarray>=0).all():
            intensityarray=numpy.fromstring(str(i[2]),dtype=numpy.float32)

    zipped_data=zip(mzarray,intensityarray)
    decompressed_spectra[i[0]]=(zipped_data,i)
    
    
    
    
    
    

for each_run in msplit_df['#File'].unique():
    this_run_df=msplit_df[msplit_df['#File']==each_run]
    z=0
    #with mzml.MzML(os.path.join(basedir,"ssl_files",each_run),use_index=True,iterative=False) as this_mzml: #build_id_cache=True ,retrieve_refs=True
    #with mzml.read(os.path.join(basedir,"ssl_files",each_run),use_index=False,iterative=False) as this_mzml: #build_id_cache=True ,retrieve_refs=True
    #with mzml.PreIndexedMzML(os.path.join(basedir,"ssl_files",use_index=True,each_run),iterative=False) as this_mzml: #build_id_cache=True ,retrieve_refs=True
    mzml_reader=run.Reader(os.path.join(basedir,"ssl_files",each_run))
    #print mzml_reader[5].peaks
    print "Finished loading the mzML file:",each_run
    #tqdm.write("Finished loading the mzML file: "+each_run)
    with open(each_run.rsplit(".",1)[0]+"_msplit.ssl",'wb') as ssl_writer:
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


#print "{0} -a milkyway.auth *_msplit.ssl msplit_filtered.blib".format(os.path.join(jar_path,"blibbuild\\BlibBuild.exe"))
subprocess.check_output("{0} -a milkyway-galaxy *_msplit.ssl combined_spectral_lib.blib".format("wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibBuild.exe"))
subprocess.check_output("{0} combined_spectral_lib.blib filtered.blib".format("wine /galaxy-central/tools/wohl-proteomics/ssl_converter/skylineblib/BlibFilter.exe"))

print "All done!"
sys.exit(0)
