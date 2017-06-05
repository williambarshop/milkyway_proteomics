import sys, os
from glob import glob

#Mixed Label Stripping Tool
#This script will read through SSL files and strip
#PSMs with mixed labeling by SILAC searches...

arguments=sys.argv[1:]
labeled_aas=arguments[0].split(",")
input_folder=arguments[1]
filter_no_SILAC_aas=False
#print "argument length is",str(len(arguments))
if (len(arguments) < 3):
    filter_no_SILAC_aas=False
elif "True" in arguments[2]:
    print "We will filter out peptides with no labelable AAs..."
    filter_no_SILAC_aas=True

print arguments[2],type(arguments[2]),"This was the 2nd arg"

if (len(arguments) < 4):
    in_place=True
    output_file="temporary_ssl_conversion.ssl"
else:
    output_file=arguments[3]


modified_labeled_aas=[]

for each in labeled_aas:
    modified_labeled_aas.append(each+"[")


ssl_input_files = [y for x in os.walk(input_folder) for y in glob(os.path.join(x[0], '*.ssl'))]
#### Below, add code to handle each input file instead of one static input... loops!

for input_file in ssl_input_files:
    wrote_header=False
    exclude_count=0
    include_count=0
    nolabel_count=0
    writer=open(output_file,'wb')
    with open(input_file,'rb') as openfile:
        for eachline in openfile:
            if not wrote_header:
                writer.write(eachline)
                wrote_header=True
                continue
            thisline=eachline.strip().split("\t")
            exclude=False
            if any(x in thisline[3] for x in modified_labeled_aas): #This peptide has heavy labeled AAs...
                modlist={}
                alllist={}
                for each_mod_aa in modified_labeled_aas:
                    modlist[each_mod_aa]=thisline[3].count(each_mod_aa)
                    alllist[each_mod_aa]=thisline[3].count(each_mod_aa[:-1])
                    if modlist[each_mod_aa] != alllist[each_mod_aa]:
                        exclude=True
                    else:
                        pass
            if filter_no_SILAC_aas:
                if not any(x in thisline[3] for x in labeled_aas):
                    nolabel_count+=1
                    #print eachline,"THIS WAS EXCLUDED BECAUSE IT CONTAINED NO LABELED AAS"
                    exclude=True
            if not exclude:
                writer.write(eachline)
                include_count+=1
            else:
                #print "EXCLUDED:",eachline
                exclude_count+=1
    writer.close() #Close it up!            

    if in_place:
        os.remove(input_file)
        os.rename(output_file,input_file)

    print "All done with "+input_file+"! We removed "+str(exclude_count-nolabel_count)+" PSMs with mixed labels and kept "+str(include_count)+" PSMs"
    if filter_no_SILAC_aas:
        print "We also removed "+str(nolabel_count)+" peptides with no labelable amino acids."
