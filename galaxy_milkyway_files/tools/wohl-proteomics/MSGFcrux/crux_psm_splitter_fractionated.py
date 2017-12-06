import sys
import getopt
import os
import pandas as pd
import numpy as np
import shutil
import Bio
from Bio import SeqIO

###############################################
#CRUX PSM SPLITTER
#William Barshop - 2015/07/13
#Updated for fractionated data processing...
date="2015/10/22"
#Laboratory of James A. Wohlschlegel, UCLA
#
#VERSION:--------------------------------------
#0.96c ALPHA
version=0.96
#
#DESCRIPTION:----------------------------------
#crux_psm_splitter_fractionated.py is a script designed to split crux percolator output psm files by the input analysis file.  It can also assemble NSAF and raw spectral counts through crux spectral-counts
#
#The script will require as INPUT a "mapping" file with two columns, the data analysis file name, and the corresponding
#crux-internal integer label.  This .cruxmap file can be obtained from the CRUX PIN FIXER script.
#
#You must give an output FOLDER in which all the separated psms files will be generated.
#The script will also 
#
#SYNTAX AND USAGE:-----------------------------
#python crux_psm_splitter_fractionated.py -i <INPUT-FOLDER> -o <OUTPUT-FOLDER> -m <MAP-FILE.cruxmap> -s <FIDO_OUTPUT> -q <PROTEIN Q-VALUE CUTOFF> -f <PROTEIN-DB.fasta> -r
#
###############################################

def help():
    print "Syntax and usage as follows:"
    print "python crux_psm_splitter_fractionated.py -i <INPUT-FOLDER> -o <OUTPUT-FOLDER> -m <MAP-FILE> -s <FIDO_OUTPUT> -q <PROTEIN Q-VALUE FILTER> -f <PROTEIN_FASTA> -r --saint"
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
saint = False
threshold = False
groupquant = False
quant_level="protein"
decoy_str = "Reverse_"
qvalue = .01 # DEFAULT PROTEIN LEVEL Q-VALUE IS .01
ptm_threshold_type=None
ptm_threshold=None

try:
    opts, args = getopt.getopt(sys.argv[1:], "ho:i:m:s:q:f:rdg:t:l:", ["help", "output=","input=","map=","s=","q=","f=","saint=","decoystr=","groupquant=",'threshold=','quantlevel=','ptmThresholdType=','ptmThreshold='])
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
    elif option in ("--saint"):
        saint=True
        saint_interact,saint_bait,saint_prey=value.split(",")
    elif option in ("--decoystr"):
        decoy_str=value
    elif option in ("--groupquant"):
        groupquant=value
    elif option in ("--threshold","-t"):
        threshold=True
        psm_threshold=value
    elif option in ("--quantlevel","-l"):
        quant_level=value
    elif option in ("--ptmThresholdType"):
        ptm_threshold_type=value
    elif option in ("--ptmThreshold"):
        ptm_threshold=value

if saint and protdb is None:
    print "You must give me a protein database to produce SAINT file inputs..."
    sys.exit(2)

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
        print "RAW Spectral counts will be generated (instead of NSAF)!"
        

print "Input folder will be read from",input

infiles=[]
for root, subFolders, files in os.walk(input):
    for eachfile in files:
        if 'target.psms.txt' in eachfile:
            infiles.append(str(os.path.join(root,eachfile)))
#print infiles,"THIS IS INFILES"


print "Mapfile will be read from",mapfile
print "Output file will be written to",output


print "--------------------------------------"
print "Reading mappings from mapfile..."


###################################################################################################
#
#   THIS BLOCK OF CODE IS REALLY POORLY DESIGNED. IT ASSUMES A CONSTANT COLUMN POSITIONING, WHICH IS
#   PRONE TO CAUSE PROBLEMS.  BELOW IS A MUCH BETTER WAY OF EXTRACTING DATA FROM A DATAFRAME BY COLUMN
#   IDENTIFICATION.  YOU SHOULD SWITCH TO THE LATTER FORM.
#   -- it's me from the future. dataframes are the way to go.

mapdict={} # keys are strings of file_idx, value is run name (without mzML)
groupIDdict={} # keys are strings of file_idx, value is GROUP NAME
with open(mapfile,'r') as openfile:
    for each in openfile:
        line=each.split("\t") #assuming tab delim
        if "Original File Name" in line[1]:
            continue
        else:
            mapdict[line[0]]=str(line[1]) #[:-1] ########### I TOOK THIS OFF... WHY IS IT HERE? !!!!!!!!!!!!!!!!!!!!!!!!
            #mapdict[line[0]]=str(line[1]) #[:-1] ########### I TOOK THIS OFF... WHY IS IT HERE? !!!!!!!!!!!!!!!!!!!!!!!!
            groupIDdict[line[0]] = line[2]            
print mapdict,"THIS IS MAPDICT!!!"
print groupIDdict,"THIS IS IDDICT!!!"
#
#
#
#
###################################################################################################



group_information = pd.read_csv(mapfile,sep='\t')

if saint:
    saint_run_dict={} # Key is file_idx, value is file_name.mzML
    saint_group_to_run_dict={} # Key is group, value is [1, 2, 3, 4] list of file_idx belonging to runs in the group...
    saint_group_to_file_name={} # key is group, value is ["xxxx.mzML", "xxxx.mzML"]
    saint_group_to_biol_cond={} # key is  group, value is "biological condition"
    saint_group_to_test={} # key is  group, value is "T" or "C"

    for index,row in group_information.iterrows():
        saint_run_dict[str(row['Crux File Integer'])]=row['Original File Name']+".mzML"
        if row['Fractionation Group ID String'] in saint_group_to_run_dict:
            saint_group_to_run_dict[row['Fractionation Group ID String']].append(str(row['Crux File Integer']))
        else:
            saint_group_to_run_dict[row['Fractionation Group ID String']] = [str(row['Crux File Integer'])]
        if row['Fractionation Group ID String'] in saint_group_to_file_name:
            saint_group_to_file_name[row['Fractionation Group ID String']].append(str(row['Original File Name'])+".mzML")
        else:
            saint_group_to_file_name[row['Fractionation Group ID String']] = [str(row['Original File Name'])+".mzML"]
        if row['Fractionation Group ID String'] in saint_group_to_biol_cond:
            #saint_group_to_file_name[row['Fractionation Group ID String']].append(str(row['Original File Name'])+".mzML")
            if saint_group_to_biol_cond[row['Fractionation Group ID String']] != str(row['Biological Condition']):
                print "ERROR: SAINT FRACTIONATION GROUPS ARE INCONSISTENT... (biological conditions)"
                print saint_group_to_biol_cond[row['Fractionation Group ID String']],"and",str(row['Biological Condition'])
                sys.exit(2)
        else:
            saint_group_to_biol_cond[row['Fractionation Group ID String']] = str(row['Biological Condition'])
        if row['Fractionation Group ID String'] in saint_group_to_test:
            #saint_group_to_file_name[row['Fractionation Group ID String']].append(str(row['Original File Name'])+".mzML")
            if saint_group_to_test[row['Fractionation Group ID String']] != str(row['Test or Control']):
                print "ERROR: SAINT FRACTIONATION GROUPS ARE INCONSISTENT... (test conditions)"
                print saint_group_to_test[row['Fractionation Group ID String']],"and",str(row['Test or Control'])
                sys.exit(2)
        else:
            saint_group_to_test[row['Fractionation Group ID String']] = str(row['Test or Control'])

    ######## Now I'm going to make a list of dictionaties for each fractionation group id string with biol_cond, and test
    dict_list=[]
    for eachgroup in saint_group_to_biol_cond:
        dict_list.append({"Fractionation Group ID String":eachgroup,"Biological Condition":saint_group_to_biol_cond[eachgroup],"Test or Control":saint_group_to_test[eachgroup]})

    saint_key_df=pd.DataFrame(dict_list)
    #print saint_key_df
    #print "that was it...."
    #sys.exit(2)




#groupdict={}
#for eachfile in infiles:
#    with open(eachfile,'r') as openfile:
#         pass####################################################



run_writers={}
###for each_run in mapdict.keys():
###    run_writers[each_run]=open(output+str(mapdict[each_run])+".psms.txt","wb")

already_made={}
for each_group in groupIDdict.keys():
    if groupIDdict[each_group] not in already_made:
        already_made[groupIDdict[each_group]]=open(output+str(groupIDdict[each_group])+".psms.txt","wb")
        run_writers[each_group]=already_made[groupIDdict[each_group]]
        ###run_writers[each_group]=open(output+str(groupIDdict[each_group])+".psms.txt","wb")
    else:
        run_writers[each_group]=already_made[groupIDdict[each_group]]
print run_writers,"run writers"
header=None
columns=0
col_name_to_header_index={}
for eachfile in infiles:
    with open(eachfile,'rb') as inputfile:
        print "Opening ",eachfile
        for each_line in inputfile:
            #print each_line,"eachline"
            line=each_line.split("\t")
            #print line, "THESE ARE THE HEADERS"
            if "file" in line and "file_idx" in line:
                if header is None:
                    header=line
                    columns=len(header)
                    for eachoutput in already_made:
                        already_made[eachoutput].write("\t".join(header))
                i=0
                for each_column in header:  
                    col_name_to_header_index[header[i]]=i
                    i+=1
                #header=line
                continue

            #print col_name_to_header_index,"CONVERSION DICT"
            #### THIS IS WHERE I WILL FILTER THE FILES FOR PTM LOCALIZATION SCORES IF DESIRED!
            #if filter
            #if this-type-of-filter
            #if not good for threshold
            #then continue, else pass
            if ptm_threshold_type is not None:
                if ptm_threshold_type=="ptmRS":
                    if float(line[col_name_to_header_index['ptmRS_peptideLocalizationProbability']])<ptm_threshold:
                        continue
                elif ptm_threshold_type=="phoshpoRS":
                    if float(line[col_name_to_header_index['pRS_peptideLocalizationProbability']])<ptm_threshold:
                        continue
                elif ptm_threshold_type=="LuciPHOr":
                    if float(line[col_name_to_header_index['luci_globalFLR']])>ptm_threshold:
                        continue


            if header is None:
                #run_writers[line[1]].write("\t".join(header))
                #header=None
                print "THIS SHOULDN'T HAPPEN! Header issue"
            #print mapdict,line
            #line[col_name_to_header_index['file']]=mapdict[line[col_name_to_header_index['file']]]+".mzML"   #Is this necessary anymore?
            ###line[0]=mapdict[line[1]]+".mzML"
            newline="\t".join(line)
            if len(line)!=columns:
                print line,"PROBLEM LINE!"
            else:
                ###run_writers[line[1]].write(newline) # OTHERWISE, INSTEAD OF NEWLINE, GIVE JUST EACH LINE
                #print run_writers, "these are the writers..."
                run_writers[line[col_name_to_header_index['file_idx']]].write(newline) # OTHERWISE, INSTEAD OF NEWLINE, GIVE JUST EACH LINE
        

print "Cleaning up and closing files..."
for each_run in run_writers.keys():
    run_writers[each_run].close()


    








run_readers={}
if speccount: #We'll do the spectral counting stuff allllllll inside here
    #Let's go ahead and run the crux spectral-counting


    print mapdict,"THIS IS MAPDICT"
    print groupIDdict,"THIS IS GROUPIDDICT"

    ###for each_run in mapdict.keys():
    for each_run in set(groupIDdict.values()):
        ###print mapdict[each_run]
        print "Working on spectral counts for ",each_run,"group"
        run_command = "crux spectral-counts"
        if threshold:
            run_command += " --threshold "+str(psm_threshold)
        if raw:
            run_command += " --measure RAW"
        elif dnsaf:
            run_command += " --measure dNSAF"
            run_command += " --protein-database "+protdb
        else:
            run_command += " --measure NSAF"
            run_command += " --protein-database "+protdb
        run_command += " --quant-level "+quant_level
        run_command += " --overwrite T"
        ###run_command += " --fileroot "+str(mapdict[each_run])
        ###run_command += " "+output+str(mapdict[each_run])+".psms.txt"
        run_command += " --fileroot "+str(each_run)
        run_command += " "+output+str(each_run)+".psms.txt"


        print "Running command...",run_command
        os.system(run_command)
        print "--------------------------------------"
        print "--------------------------------------"



    #Now we'll go ahead and read in those outputs from the crux-output folder.
    #OK!

    print "-----------------------------------------"

    print "Let's grab the list of proteins passing our q-value filter!"

    #NBNBNBNBNBNBNBNB:    FOR THE FUTURE, THE BELOW SECTION SHOULD DEFINITELY BE CHANGED TO PANDAS READING IN THE PROTEIN FILE
    #                     AND THEN FILTERING AFTER THAT....

    passing_proteins={} # These are actually protein groups...
    with open(proteinfile,'r') as protein_file:
        for each_line in protein_file:
            line = each_line.split("\t")
            if "q-value" in line[0]:
                continue
            elif float(line[0]) <= qvalue: #and decoy_str not in line[1]:
                if "," in line[3]:
                    temp=line[3].replace("{","").replace("}","").split(",")
                    temp2=[]
                    for eachone in temp:
                        if not decoy_str in eachone:
                            temp2.append(eachone)
                    if len(temp2)>0:
                        passing_proteins[','.join(temp2)]=float(line[0])
                else:
                    if not decoy_str in line[3]:
                        passing_proteins[line[3]]=float(line[0])

    #print passing_proteins,"passing proteins"

    #raw_input("yo")
    print "ALLRIGHT -- now we're going to go ahead and read all the spectral counts into a dict of dicts, and maybe make that into a dataframe..."


    #print mapdict,"THIS IS MAPDICT!!"

    ###for each_run in mapdict.keys():
    for each_run in groupIDdict.keys():
        ###run_readers[mapdict[each_run]]=open("crux-output/"+str(mapdict[each_run])+".spectral-counts.target.txt","r")
        run_readers[groupIDdict[each_run]]=open("crux-output/"+str(groupIDdict[each_run])+".spectral-counts.target.txt","r")
        ###print "Added file","crux-output/"+str(mapdict[each_run])+".spectral-counts.target.txt"
        print "Added file","crux-output/"+str(groupIDdict[each_run])+".spectral-counts.target.txt"

    #shutil.copytree(output,"/home/galaxy/TESTINGOUTPUT/")

    run_dict={}

    for each_file in run_readers:
        protein_dict={}
        run_dict[each_file]=protein_dict


    for each_file in run_readers:
        for each_line in run_readers[each_file]:
            line=each_line.split("\t")
            line[1]=line[1][:-1]
            #print line,"MYLINE"
            if "protein id" in line[0]:
                continue
            for eachproteingroup in passing_proteins:
                #print line[0],eachproteingroup
                if line[0] in eachproteingroup:
                    if "," in eachproteingroup: #It's a protein group...
                        #print "I FOUND A PROTEIN GROUP!",line[0],eachproteingroup,each_file
                        #if eachproteingroup.replace("{","").replace("}","").strip() in run_dict[each_file]:
                        for each_file2 in run_readers:
                            if raw:
                                run_dict[each_file2][eachproteingroup.replace("{","").replace("}","").strip()]=['0']*(len(eachproteingroup.split(",")))
                                #run_dict[each_file2][eachproteingroup.replace("{","").replace("}","").strip()]=['0']*(len(eachproteingroup.split(","))-eachproteingroup.count(decoy_str)
                            else:
                                run_dict[each_file2][eachproteingroup.replace("{","").replace("}","").strip()]=['0.0']*(len(eachproteingroup.split(",")))
                                #run_dict[each_file2][eachproteingroup.replace("{","").replace("}","").strip()]=['0.0']*(len(eachproteingroup.split(","))-eachproteingroup.count(decoy_str))

    

    #for each_file in run_readers:
    #    #protein_dict={}
    #    print run_dict[each_file],"run dict for",each_file


    for each_file in run_readers:
        run_readers[each_file].seek(0)
        #print each_file,"in run readers..."


    for each_file in run_readers:
        for each_line in run_readers[each_file]:
            line=each_line.split("\t")
            line[1]=line[1][:-1]
            #print line,"MYLINE"
            if "protein id" in line[0] or "sequence" in line[0]:
                continue
            if each_file not in run_dict:
                protein_dict={}
                run_dict[each_file]=protein_dict
            if "peptide" in quant_level:
                if raw:
                    run_dict[each_file][line[0]]=int(line[1])
                else:
                    run_dict[each_file][line[0]]=float(line[1])
                continue

            for eachproteingroup in passing_proteins:
                #if line[0] in eachproteingroup:
                if line[0] in eachproteingroup:
                    
                    #if eachproteingroup.replace("{","").replace("}","").strip() in run_dict[each_file]:
                    #    print "old",run_dict[each_file][eachproteingroup.replace("{","").replace("}","").strip()],"new",line[1]

                    if "," in eachproteingroup: #It's a protein group...
                        if eachproteingroup.replace("{","").replace("}","").strip() in run_dict[each_file]:
                            j=0
                            splitgrp=eachproteingroup.split(",")
                            protgrp_list=run_dict[each_file][eachproteingroup.replace("{","").replace("}","").strip()]
                            while j < len(splitgrp)-1:
                                #if line[0] in splitgrp[j]:
                                if line[0] == splitgrp[j]:
                                    stop=j
                                    break
                                j+=1
                            #print protgrp_list,"yo"
                            if protgrp_list[j]!='0.0' and protgrp_list[j]!="0":
                                print "WARNING:::: REPLACING A VALUE UNEXPECTEDLY... old",protgrp_list[j],"new",line[1]
                            protgrp_list[j]=line[1]

                            run_dict[each_file][eachproteingroup.replace("{","").replace("}","").strip()]=protgrp_list                            
                        else:
                            j=0
                            splitgrp=eachproteingroup.split(",")
                            protgrp_list=['0.0']*(len(splitgrp))
                            while j < len(splitgrp):
                                #if line[0] in splitgrp[j]:
                                if line[0] == splitgrp[j]:
                                    stop=j
                                    break
                                j+=1
                            protgrp_list[j]=line[1]

                            run_dict[each_file][eachproteingroup.replace("{","").replace("}","").strip()]=protgrp_list
                        

                    else:
                        if not raw:
                            run_dict[each_file][eachproteingroup.replace("{","").replace("}","").strip()]=float(line[1])
                        else:
                            #print int(line[1]),"yo"
                            run_dict[each_file][line[0]]=int(line[1])
                    break

            #if line[0] in passing_proteins:
            #    if not raw:
            #        run_dict[each_file][line[0]]=float(line[1])
            #    else:
            #        run_dict[each_file][line[0]]=int(line[1])
            #else:
            #    continue
            

            #if "file" in line[0] and "file_idx" in line[1]:
            #    continue
    #print run_dict
    #print run_dict.keys()

    for eachrun in run_dict:
        thisrun=run_dict[eachrun]
        for each_prot in thisrun:
            thisprot=thisrun[each_prot]
            #if "," in each_prot:
            #    print "it's a list!",each_prot
            #print thisprot,"thisprot"
            #print each_prot,"each_prot"
            #print "----------------------"
            #if type(thisprot)!=list:
            #    print "THIS IS NOT A LIST",thisprot
            #    #pass
            #else:
            #    pass
            #    #print "THIS IS A LIST",thisprot,"<-----------------------------"
            if (type(thisprot)!=list) and ("," in each_prot):
                #print "I'M IN!<---------------------------------"
                if not raw:
                    thisrun[each_prot]=['0.0']*len(thisprot.split(","))
                else:
                    thisrun[each_prot]=['0']*len(thisprot.split(","))
            elif type(thisprot) is list:
                pass
                #for each in thisprot:
                #    if raw:
                #        if type(each) is not int:
                #            print "MORE CLEANUP!"


    if raw:
        dataframe=pd.DataFrame(run_dict,dtype='object')
        dataframe=dataframe.fillna("0")
    else:
        dataframe=pd.DataFrame(run_dict,dtype='object').fillna("0.0")
    #dataframe=dataframe.convert_objects(convert_numeric=True)
    #print dataframe.columns.values
    if raw:
        speccount_output=output+"filtered-spectral-counts-RAW.csv"
    else:
        speccount_output=output+"filtered-spectral-counts-NSAF.csv"



    if groupquant is not False:
        gqframe=dataframe.copy(deep=True)
        if "average" in groupquant:
            #We will collapse ambig. protein group values into the average of their counts/NSAFs
            #groups_output[each_group]=output_results[(output_results['file_idx'].str.contains('|'.join(group_to_runids[each_group])))]
            speccount_mask=dataframe[(dataframe.index.to_series().str.contains(','))]
            for eachindex,eachrow in speccount_mask.iterrows():
                group_length=eachindex.count(",")+1
                for eachcolumn in dict(eachrow):
                    #this_avg=np.mean(eachrow[eachcolumn])
                    if raw:
                        this_avg=int(round(np.mean(map(int,eachrow[eachcolumn]))))
                    else:
                        this_avg=np.mean(map(float,eachrow[eachcolumn]))
                    gqframe.ix[eachindex,eachcolumn]=this_avg



        with open(speccount_output,'wb') as output_writer:
            gqframe.to_csv(path_or_buf=output_writer,sep="\t",index_label="ProteinID")
    else:
        with open(speccount_output,'wb') as output_writer:
            dataframe.to_csv(path_or_buf=output_writer,sep="\t",index_label="ProteinID")
    print dataframe
    #print dataframe.dtypes
    print "Spectral Counts dataframe written to disk at ",speccount_output
    #shutil.copytree(output,"/home/galaxy/TESTINGOUTPUT/")
    #shutil.copytree(os.path.join(os.getcwd(),"psm-split-output"),"/home/galaxy/PSMSPLITOUTPUT/")

#saint_run_dict={} # Key is file_idx, value is file_name.mzML
#saint_group_to_run_dict={} # Key is group, value is [1, 2, 3, 4] list of file_idx belonging to runs in the group...
#saint_group_to_file_name={} # key is group, value is ["xxxx.mzML", "xxxx.mzML"]
#saint_group_to_biol_cond={} # key is  group, value is "biological condition"
#saint_group_to_test={} # key is  group, value is "T" or "C"

#saint_key_df <------ this is the dataframe with saint fraction group, biological condition and test information
#        dict_list.append({"Fractionation Group ID String":eachgroup,"Biological Condition":saint_group_to_biol_cond[eachgroup],"Test or Control":saint_group_to_test[eachgroup]})








if saint:
    #print "writing files to",saint_interact,saint_bait,saint_prey,"<--------------------------------"
    #Let's go ahead and open up our FASTA db...
    #We'll construct a dictionary of {Handle:Length} information
    protlen_dict={}
    with open(protdb,'rb') as openfasta:
        for record in SeqIO.parse(openfasta,"fasta"):
            protlen_dict[record.id]=len(record.seq)

    #print protlen_dict

    interact_writer=open(saint_interact,'wb')
    bait_writer=open(saint_bait,'wb')
    prey_writer=open(saint_prey,'wb')
    print "->>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<-"
    print "-------------------------------------------------"
    print "We'll begin, now, to process SAINT input files!"
    print "First, we'll make our BAIT file!"
    

    for index,each_row in saint_key_df.iterrows():
        #print each_row
        if each_row["Test or Control"].strip() is "C":
            bait_writer.write(each_row["Fractionation Group ID String"]+"\t"+each_row["Fractionation Group ID String"]+"_ctrl"+"\t"+each_row["Test or Control"]+"\n")
        else:    
            bait_writer.write(each_row["Fractionation Group ID String"]+"\t"+each_row["Biological Condition"]+"\t"+each_row["Test or Control"]+"\n")

    print "Now, we'll make out PREY file!"
    for each_protein in protlen_dict:
        prey_writer.write(each_protein+"\t"+str(protlen_dict[each_protein])+"\n")
    if groupquant is not False:
        prot_groups=gqframe[(gqframe.index.to_series().str.contains(','))]
        for thisindex,thisgroup in prot_groups.iterrows():
            prey_writer.write(thisindex+"\t"+str(int(round(np.mean([protlen_dict[prot] for prot in thisindex.split(',')]))))+"\n")


    print "Finally, we'll do the hard part... Building the INTERACT file!"

    if groupquant is False:
        iterateoverme=dataframe
    else:
        iterateoverme=gqframe

    for index,each_prot_group in iterateoverme.iterrows():
        if groupquant is False:
            if len(index.split(","))>1:
                i=0#This holds our list position.
                for each_protein in index.split(","):
                    if not decoy_str in each_protein:
                        for eachgroup in dataframe.columns.values.tolist():
                            if raw:
                                interact_writer.write(eachgroup+"\t"+saint_group_to_biol_cond[eachgroup]+"\t"+each_protein+"\t"+str(each_prot_group[eachgroup][i])+"\n")
                            else:
                                interact_writer.write(eachgroup+"\t"+saint_group_to_biol_cond[eachgroup]+"\t"+each_protein+"\t"+str(int(round(float(each_prot_group[eachgroup][i])*1000000.0)))+"\n")

                    i+=1
            #print index.split(",")
            #print each_prot_group
            else:
                for eachgroup in iterateoverme.columns.values.tolist():
                    #print eachgroup
                    if not decoy_str in index:
                        #print "THIS IS INDEX",index
                        if raw:
                            interact_writer.write(eachgroup+"\t"+saint_group_to_biol_cond[eachgroup]+"\t"+index+"\t"+str(each_prot_group[eachgroup])+"\n")
                        else:
                            interact_writer.write(eachgroup+"\t"+saint_group_to_biol_cond[eachgroup]+"\t"+index+"\t"+str(int(round(float(each_prot_group[eachgroup])*1000000.0)))+"\n")

        else:
            for eachgroup in iterateoverme.columns.values.tolist():
                if not decoy_str in index:
                    if raw:
                        interact_writer.write(eachgroup+"\t"+saint_group_to_biol_cond[eachgroup]+"\t"+index+"\t"+str(each_prot_group[eachgroup])+"\n")
                    else:
                        interact_writer.write(eachgroup+"\t"+saint_group_to_biol_cond[eachgroup]+"\t"+index+"\t"+str(int(round(float(each_prot_group[eachgroup])*1000000.0)))+"\n")


    interact_writer.close()
    bait_writer.close()
    prey_writer.close()






print "Finished and closed!"
