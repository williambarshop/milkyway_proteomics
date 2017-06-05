import orangecontrib
from orangecontrib.bio.ontology import OBOParser
import sys
import csv
import getopt
import re
import os
from os import listdir
from os.path import isfile, join



###############################################
#CRUX PIN FIXER FOR FRACTIONATED DATA
#William Barshop - 2015/08/24
date="2015/08/24"
#Laboratory of James A. Wohlschlegel, UCLA
#
#VERSION:--------------------------------------
#0.7.5 ALPHA
version="0.7.5"
#
#DESCRIPTION:----------------------------------
#crux_pin_fixer.py is a script designed to correct pin files generated from msgf2pin
#and reformat the PSMIDs to be crux friendly.
#This requires reformatting of the PSMids from msgf2pin (formatted as : 2015-05-26-wb-HEK293-Std-rtid2_SII_192_1_192_2_1 )
#to something more crux friendly -- ie replacing the filename with an integer fileID and placing in either
#"target" or "decoy" at the start of the PSM file ID... like ( target_0_35_1_5 )
#
#The script will OUTPUT a "mapping" file with two columns, the data analysis file name, and the corresponding
#crux-internal integer label.  If no filename is given for the output, the basename of the pin will be used
#and appended with the extension .cruxmap
#Additionally, if no output location is given, the reprocessed pin file will be output to the basename of the input file
#and appended with .pin
#
#SYNTAX AND USAGE:-----------------------------
#python crux_pin_fixer.py -i input_folder -o output_folder -m <MAP-FILE> -u <unimod.obo>
#
###############################################

def help():
    print "Syntax and usage as follows:"
    print "python crux_pin_fixer.py -i input_folder -o output_folder -m <MAP-FILE> -u <unimod.obo> -e <experiment_file>"
    print "Version:",version
    print "Dated:",date

input = None
output = None
mapfile = None
mod_obo = None
experiments = None
seleno=False

try:
    opts, args = getopt.getopt(sys.argv[1:], "ho:i:m:u:e:s", ["help", "output=","input=","map=","mod_obo=","exp="])
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
    elif option in ("-i","--input"):
        input_folder = value
        onlyfiles = [ f for f in listdir(input_folder) if isfile(join(input_folder,f)) and f.endswith(".tsv") ]
        out=""
        for each in onlyfiles:
            out+=each+"," #"output/"+
        input=out[:-1]
        #print input,"THIS IS INPUT"
    elif option in ("-m","--map"):
        mapfile = value
    elif option in ("-u","--unimod"):
        mod_obo = value
    elif option in ("-e","--exp"):
        experiments = value
    elif option in ("-s","--seleno"):
        seleno = True

if input_folder is None:
    print "Come on... at least give me an input...."
    help()
    sys.exit(2)

if mapfile is None:
    if len(input.split(","))>1:
        mapfile = "fractionated_analysis.cruxmap"
    else:
        mapfile = input.rsplit(".",1)[0]
        mapfile += ".cruxmap"

if output is None:
    if len(input.split(","))>1:
        output_folder="newoutput"
        os.mkdir(output_folder)
        output = []
        output_list=True
        for each in input.split(","):
            output.append(each.rsplit(".",1)[0]+".pin")
    else:
        output_list=False
        output = input.rsplit(".",1)[0]
        output +=".pin"

if mod_obo is None:
    print "mod_obo is needed! Please run with \'-u unimod_obo\' in the arguments "
    sys.exit(2)

output_list=True
output_folder=output
output=[]
for each in input.split(","):
    output.append(output_folder+each.rsplit(".",1)[0]+".pin")

print "Input file will be read from",input
if len(input.split(","))>1:
    print "Output files will be written to folder output"
else:
    print "Output file will be written to",output
print "Mapfile output file will be written to",mapfile

fileNameDictionary={}
file_int=0


mod_bracket_re='\[|\]'
mod_dict={}



zipped=zip(input.split(","),output)
print zipped,"this is zipped..."
print output,"this is output"



for eachfiles in zipped:
    with open(mod_obo,"r") as unimod_obo:
        peptide_index=0
        obo_parser=orangecontrib.bio.ontology.OBOParser(unimod_obo)
        with open("output/"+eachfiles[0],"r") as openfile:
            with open(eachfiles[1],"wb") as outputwriter:
                for eachline in openfile:
                    line=eachline.split("\t")
                    if peptide_index==0:
                        for each_column in line:
                            if each_column=="Peptide":
                                break
                            else:
                                #print each_column
                                peptide_index+=1
                      


                    if line[0]=="SpecId":
                        outputwriter.write(eachline)
                        continue
                    elif line[0]=="DefaultDirection":
                        outputwriter.write(eachline)
                        continue
                    file_name=line[0].rsplit("_",7)[0]
                    #file_name=line[0].split("_")[0]
                    #print file_name,"this is file name..."
                    #print "coming from...",line[0]
                    if file_name not in fileNameDictionary:
                        fileNameDictionary[file_name]=file_int
                        file_int+=1

                    if line[1] == "1":
                        target="target"
                    else:
                        target="decoy"

                    #print line, "THIS IS LINE"
                    #print line[0], "THIS IS LINE 0"






                    #print len(line),"this is line len"
                    #print line
                    #print "peptide index",peptide_index,len(line)
                    #print line[peptide_index]


                    ####new_peptide=line[peptide_index].replace("UNIMOD:","")  # WISHFUL THINKING


                    new_peptide=None
                    #print "LIN53",line[peptide_index]
                    #print line, peptide_index
                    #print "line len",len(line)
                    if seleno:
                        this_pep=line[peptide_index].replace("[UNIMOD:162]","")
                        line[peptide_index]=this_pep


                    if "UNIMOD" in line[peptide_index]:
                        split_peptide=re.split(mod_bracket_re,line[peptide_index])
                        #print split_peptide
                        new_peptide=[]
                        for each_split in split_peptide:
                            if "UNIMOD" in each_split:
                                if each_split in mod_dict:
                                    new_peptide.append(mod_dict[each_split])
                                else:
                                    mass=""
                                    print each_split," was not in the dictionary!"
                                    #print obo_parser,type(obo_parser),"this is obo parser"
                                    trigger=False
                                    unimod_obo.seek(0)
                                    for event,value in obo_parser:
                                        if "TAG_VALUE" in event and not trigger:
                                            #print event,value,"this is event"
                                            if each_split in value[1]:
                                                #print "YO",value[1],"and each split is",each_split
                                                trigger=True
                                                #pass
                                        elif trigger:
                                            #print event,value,"these are event and value..."
                                            if "delta_mono_mass" in value[1]:
                                                print value[1],"val1 should be the mass"

                                                mass="[+"
                                                mass+=value[1].split("\"")[1]
                                                #print mass
                                                if "-" in mass:
                                                    mass=mass.replace("+","")
                                                    #print "The mass was below zero..."
                                                mass+="]"
                                                #print mass,"THIS IS MAH MASS"
                                                trigger=False
                                                break
                                        else:
                                            continue
                                    if mass=="":
                                        print "ERROR: ERROR: ERROR: THE MASS FOR THIS MOD WAS NOT FOUND IN THE UNIMOD OBO FILE..."
                                        sys.exit(2)    
                                    mod_dict[each_split]=mass
                                    new_peptide.append(mass)
                                        
                            else:
                                new_peptide.append(each_split)
                        new_peptide=''.join(new_peptide)
                        

        
                    #print mod_dict,"this is the mod dict..."
        

                                    



                    newSpecID = target + "_" + str(fileNameDictionary[file_name]) + "_" + str(line[2]) + "_" + str(line[0].split("_")[5]) + "_" + str(line[0].split("_")[6])
                    #print newSpecID
                    #print line
                    #print eachline
                    line[0] = newSpecID
                    #sys.exit(2)
                    if not new_peptide is None:
                        line[peptide_index]=new_peptide #####CHECK IF THIS SHOULD BE 53
                    #print line,"AFTER"
                    writeLine = ""
                    for eachItem in line:
                        writeLine += eachItem
                        writeLine += "\t"
                    writeLine=writeLine[:-2]+"\n"
                    #print writeLine
                    outputwriter.write(writeLine)
        
                    #pause=raw_input("continue?")


#            if len(myline)>2:
#                if myline[2] in bio_cond_groups:
#                    bio_cond_groups[myline[1]].append(myline[2])
#                else:
#                    bio_cond_groups[myline[1]]=[myline[2]]
#print bio_cond_groups,"These are the biological conditions."

bio_cond_groups={}
experimentDict={}
groupDict={}
test_control_groups={}
biol_rep_groups={}
if experiments is not None:
    with open(experiments,'r') as openfile:
        for eachline in openfile:
            if len(eachline.split("___")) > 1:
                group_IDstring=eachline.strip().split("___")[0]
                group_name=eachline.strip().split("___")[1]
                biol_group=eachline.strip().split("___")[2]
                test_control=eachline.strip().split("___")[3]
                biol_replicate=eachline.strip().split("___")[4]

                groupDict[group_IDstring]=group_name
                bio_cond_groups[group_IDstring]=biol_group
                test_control_groups[group_IDstring]=test_control
                biol_rep_groups[group_IDstring]=biol_replicate
    for each in fileNameDictionary:
        for eachgroup in groupDict:
            #print "working on",eachgroup,"looking to match",groupDict[eachgroup],"to",each
            if eachgroup in each:
            #if groupDict[eachgroup] in each:
                #experimentDict[each]=groupDict[eachgroup]
                experimentDict[each]=eachgroup
                break
                
print "We're going to be looking through",fileNameDictionary
print "-----------------"
with open(mapfile,"wb") as mapfilewriter:
    if experiments is None:
        mapfilewriter.write("Crux File Integer\tOriginal File Name\n")
    else:
        mapfilewriter.write("Crux File Integer\tOriginal File Name\tFractionation Group ID String\tFractionation Group Name\tBiological Condition\tTest or Control\tBioReplicate\n")
    for each in fileNameDictionary:
        if experiments is None:
            mapfilewriter.write(str(fileNameDictionary[each])+"\t"+each+"\n")
        else:
            #print each,"each!<<<<<<<<<<<<<<"
            #print fileNameDictionary[each]
            #print experimentDict,"exp dict!<-===============-"
            #print experimentDict[each]
            #print groupDict
            #print groupDict[experimentDict[each]]
            #print bio_cond_groups[experimentDict[each]]
            #print test_control_groups[experimentDict[each]]

            #mapfilewriter.write(str(fileNameDictionary[each])+"\t"+each+"\t"+experimentDict[each]+"\t"+groupDict[experimentDict[each]]+"\t"+bio_cond_groups[experimentDict[each]]+"\t"+test_control_groups[experimentDict[each]]+"\t"+biol_rep_groups[experimentDict[each]]+"\n")
            mapfilewriter.write(str(fileNameDictionary[each])+"\t"+each+"\t"+groupDict[experimentDict[each]]+"\t"+groupDict[experimentDict[each]]+"\t"+bio_cond_groups[experimentDict[each]]+"\t"+test_control_groups[experimentDict[each]]+"\t"+biol_rep_groups[experimentDict[each]]+"\n")
print "Finished!"

