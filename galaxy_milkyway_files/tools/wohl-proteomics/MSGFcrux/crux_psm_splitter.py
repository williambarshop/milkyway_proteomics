import sys
import getopt
import os
import pandas as pd

###############################################
#CRUX PSM SPLITTER
#William Barshop - 2015/07/13
date="2015/07/13"
#Laboratory of James A. Wohlschlegel, UCLA
#
#VERSION:--------------------------------------
#0.51 ALPHA
version=0.51
#
#DESCRIPTION:----------------------------------
#crux_psm_splitter.py is a script designed to split crux percolator output psm files by the input analysis file.  It can also assemble NSAF and raw spectral counts through crux spectral-counts
#
#The script will require as INPUT a "mapping" file with two columns, the data analysis file name, and the corresponding
#crux-internal integer label.  This .cruxmap file can be obtained from the CRUX PIN FIXER script.
#
#You must give an output FOLDER in which all the separated psms files will be generated.
#The script will also 
#
#SYNTAX AND USAGE:-----------------------------
#python crux_psm_splitter.py -i <INPUT-TAB.target.psms.txt> -o <OUTPUT-FOLDER> -m <MAP-FILE.cruxmap> -s <TARGET.PROTEINS.TXT> -q <PROTEIN Q-VALUE CUTOFF> -f <PROTEIN-DB.fasta> -r
#
###############################################

def help():
    print "Syntax and usage as follows:"
    print "python crux_psm_splitter.py -i <INPUT-TAB.target.psms.txt> -o <OUTPUT-FOLDER> -m <MAP-FILE> -s -q <PROTEIN Q-VALUE FILTER> -f <PROTEIN_FASTA> -r"
    print "-r is for RAW spec counts vs NSAF"
    print "Version:",version
    print "Dated:",date

input = None
output = None
mapfile = None
speccount = False
proteinfile = None
protdb = None
raw = False
dnsaf = False
qvalue = .01 # DEFAULT Q-VALUE IS .01


try:
    opts, args = getopt.getopt(sys.argv[1:], "ho:i:m:s:q:f:rd", ["help", "output=","input=","map=","s=","q=","f="])
except getopt.GetoptError as err:
    # print help information and exit:
    print str(err)
    print "-------------------------------------------------"
    print "Invalid syntax or option!"
    help()
    sys.exit(2)

for option, value in opts:
    if option in ("-h","--help"):
        help()
        sys.exit()
    elif option in ("-o","--output"):
        output = value
        if output[-1] != "/":
            output+="/"
    elif option in ("-i","--input"):
        input = value
    elif option in ("-m","--map"):
        mapfile = value
    elif option in ("-s","--speccount"):
        speccount=True
        proteinfile = value
    elif option in ("-q","--qvalue"):
        qvalue = float(value)
    elif option in ("-f","--fasta"):
        protdb = value
    elif option in ("-r","--raw"):
        raw=True
    elif option in ("-d","--dnsaf"):
        dnsaf=True

if input is None:
    print "Come on... at least give me an input...."
    help()
    sys.exit(2)

if mapfile is None:
    mapfile = input.rsplit(".",1)[0]
    mapfile +=".cruxmap"

if output is None:
    if not os.path.exists('psm-split-output'):
        os.makedirs('psm-split-output')
    output='psm-split-output/'
    #output = input.rsplit(".",1)[0]
    #output +=".crux_ready.pin"

if speccount:
    print "Q-value filtering will be at",qvalue
    print "Protein information will be pulled from",proteinfile
    if protdb is None and not raw:
        print "If you're planning on calculating NSAF values, you'll need a protein FASTA."
        sys.exit(2)
    if protdb is not None:
        print "Protein FASTA to be loaded:",protdb
    elif raw:
        print "RAW Spectral counts will be generated!"
        

print "Input file will be read from",input
print "Mapfile will be read from",mapfile
print "Output file will be written to",output


print "--------------------------------------"
print "Reading mappings from mapfile..."



mapdict={}
with open(mapfile,'r') as openfile:
    for each in openfile:
        line=each.split("\t") #assuming tab delim
        if "Original File Name" in line[1]:
            continue
        else:
            mapdict[line[0]]=str(line[1])[:-1]


run_writers={}
for each_run in mapdict.keys():
    run_writers[each_run]=open(output+str(mapdict[each_run])+".psms.txt","wb")
    

with open(input,'r') as inputfile:
    for each_line in inputfile:
        line=each_line.split("\t")
        if "file" in line[0] and "file_idx" in line[1]:
            for each_writer in run_writers:
                run_writers[each_writer].write(each_line)
            continue
        #print line
        run_writers[line[1]].write(each_line)
        

print "Cleaning up and closing files..."
for each_run in run_writers.keys():
    run_writers[each_run].close()

run_readers={}
if speccount: #We'll do the spectral counting stuff allllllll inside here
    #Let's go ahead and run the crux spectral-counting




    for each_run in mapdict.keys():
        print mapdict[each_run]
        run_command = "crux spectral-counts"
        if raw:
            run_command += " --measure RAW"
        elif dnsaf:
            run_command += " --measure dNSAF"
            run_command += " --protein-database "+protdb
        else:
            run_command += " --measure NSAF"
            run_command += " --protein-database "+protdb
        run_command += " --overwrite T"
        run_command += " --fileroot "+str(mapdict[each_run])
        run_command += " "+output+str(mapdict[each_run])+".psms.txt"


        print "Running command...",run_command
        print "--------------------------------------"
        os.system(run_command)



    #Now we'll go ahead and read in those outputs from the crux-output folder.
    #OK!

    print "-----------------------------------------"

    print "Let's grab the list of proteins passing our q-value filter!"

    passing_proteins={}
    with open(proteinfile,'r') as protein_file:
        for each_line in protein_file:
            line = each_line.split("\t")
            if "percolator score" in line[0]:
                continue
            elif float(line[2]) <= qvalue:
                passing_proteins[line[3]]=float(line[2])

    #print passing_proteins

    #raw_input("yo")
    print "ALLRIGHT -- now we're going to go ahead and read all the spectral counts into a dict of dicts, and maybe make that into a dataframe..."

    for each_run in mapdict.keys():
        run_readers[mapdict[each_run]]=open("crux-output/"+str(mapdict[each_run])+".spectral-counts.target.txt","r")
        print "Added file","crux-output/"+str(mapdict[each_run])+".spectral-counts.target.txt"


    run_dict={}

    for each_file in run_readers:
        for each_line in run_readers[each_file]:
            #print run_readers[each_file],"yo"
            line=each_line.split("\t")
            line[1]=line[1][:-1]
            #print line
            if "protein id" in line[0]:
                continue
            if each_file not in run_dict:
                protein_dict={}
                run_dict[each_file]=protein_dict
            if line[0] in passing_proteins:
                if not raw:
                    run_dict[each_file][line[0]]=float(line[1])
                else:
                    run_dict[each_file][line[0]]=int(line[1])
            else:
                continue
            

            #if "file" in line[0] and "file_idx" in line[1]:
            #    continue

    #print run_dict.keys()
    if raw:
        dataframe=pd.DataFrame(run_dict,dtype='int32')
        dataframe=dataframe.fillna(0).astype(int)
    else:
        dataframe=pd.DataFrame(run_dict,dtype='float').fillna(0.0)
    #dataframe=dataframe.convert_objects(convert_numeric=True)
    #print dataframe.columns.values
    if raw:
        speccount_output=output+"filtered-spectral-counts-RAW.csv"
    else:
        speccount_output=output+"filtered-spectral-counts-NSAF.csv"

    with open(speccount_output,'wb') as output_writer:
        dataframe.to_csv(path_or_buf=output_writer,sep="\t")
    #print dataframe
    #print dataframe.dtypes
    print "Spectral Counts dataframe written to disk at ",speccount_output
    

print "Finished and closed!"
