version=".92B"
date="7/31/2018"

import sys, os, subprocess, shutil
import pandas
from os import listdir
from os.path import isfile, join

def find_charge_range(input_pin_files): #This function is called to return information on the max an min charges considered among a collection of percolator 'pin' input files
    all_max_charge=2
    all_min_charge=2
    for each in onlyfiles:
        line_ctr=0
        with open(folder+"/"+each,'rb') as pin_reader:
            peptide_index=0
            index_ctr=0
            expmass_index=0
            specID_index=0
            scan_index=0
            label_index=0
            charge_indicies={}
            charges=[]
            pin_lines=[]
            for eachline in pin_reader:
                if line_ctr==0:
                    header=eachline.split("\t")
                    for each_item in header:
                        if each_item == "Peptide":
                            peptide_index=index_ctr
                        elif each_item == "SpecId":
                            specID_index=index_ctr
                        elif each_item == "ScanNr":
                            scan_index=index_ctr
                        elif each_item == "ExpMass":
                            expmass_index = index_ctr
                        elif each_item == "Label":
                            label_index=index_ctr
                        elif "Charge" in each_item:
                            charge_indicies[each_item]=index_ctr
                            charges.append(int(each_item.replace("Charge","")))
                        index_ctr+=1
                else:
                    break
            line_ctr=0
            charge_index_max=max(charge_indicies.values())
            charge_index_min=min(charge_indicies.values())
            charge_max=max(charges)
            charge_min=min(charges)
            if charge_max>all_max_charge:
                all_max_charge=charge_max
            if charge_min<all_min_charge:
                all_min_charge=charge_min
    return all_max_charge,all_min_charge



arguments=sys.argv
print arguments
arguments=arguments[1:]
independent_perco=arguments[0]
arguments=arguments[1:]
folder=arguments[0]
arguments=arguments[1:]
exp_group_file=arguments[0]
arguments=arguments[1:]

proton_mass=1.007276

onlyfiles = [ f for f in listdir(folder) if isfile(join(folder,f)) ]

startingdir=os.getcwd()

#cmd = " ".join(arguments)
#print " ".join(arguments)
processes=[]

os.mkdir("output")
#print onlyfiles,"onlyfiles"

group_information = pandas.read_csv(exp_group_file,sep='\t')

run_dict={} # Key is file_idx, value is file_name.mzML
run_dict_reverse={}
group_to_run_dict={} # Key is group, value is [1, 2, 3, 4] list of file_idx belonging to runs in the group...
group_to_file_name={} # key is group, value is ["xxxx.mzML", "xxxx.mzML"]
run_to_group_dict={} #Key is run, value is group for pin file
for index,row in group_information.iterrows():
    run_dict[str(row['Crux File Integer'])]=row['Original File Name']+".mzML"
    run_dict_reverse[row['Original File Name']+".mzML"]=row['Crux File Integer']
    run_to_group_dict[row['Original File Name']+".mzML"]=row['Fractionation Group ID String']
    if row['Fractionation Group ID String'] in group_to_run_dict:
        group_to_run_dict[row['Fractionation Group ID String']].append(str(row['Crux File Integer']))
    else:
        group_to_run_dict[row['Fractionation Group ID String']] = [str(row['Crux File Integer'])]
    if row['Fractionation Group ID String'] in group_to_file_name:
        group_to_file_name[row['Fractionation Group ID String']].append(str(row['Original File Name'])+".mzML")
    else:
        group_to_file_name[row['Fractionation Group ID String']] = [str(row['Original File Name'])+".mzML"]





#Here we actually run percolator
if independent_perco=="T": #If =="T" then we'll run percolator separately for each file or group
    for each in onlyfiles:
        #print each
        cmd=arguments[:]
        cmd.append("--fileroot")
        cmd.append(each.split(".")[0])
        cmd.append(each)
        os.mkdir("output/"+each+"_out")
        os.rename(folder+"/"+each,"output/"+each+"_out/"+each)
        os.chdir("output/"+each+"_out")
        print "---------------------------"
        print cmd,"this will run"
        print "---------------------------"
        processes.append(subprocess.Popen(cmd))
        os.chdir(startingdir)
    for each in processes:
        each.wait()

    for each in onlyfiles:  #here i need to take the group and filter it into the runs that belong TO THIS GROUP so use the group_to_run_dict to get the fileidx set
        clean_name=each.rsplit(".",1)[0]
        files_to_filter=[clean_name+'.percolator.decoy.peptides.txt',clean_name+'.percolator.target.peptides.txt',clean_name+'.percolator.target.psms.txt',clean_name+'.percolator.decoy.psms.txt']
        os.chdir(startingdir+"/output/"+each.replace(".mzML",".pin")+"_out/crux-output")##

        for each_newfile in group_to_file_name[each.replace(".pin","")]:#run_dict_reverse:
            print "handling {0} and {1}".format(each,each_newfile)
            if each_newfile.replace(".mzML","") != each.replace(".pin",""):
                if not os.path.isdir(startingdir+"/output/"+each_newfile.replace(".mzML",".pin")+"_out/"):
                    os.mkdir(startingdir+"/output/"+each_newfile.replace(".mzML",".pin")+"_out/")
                shutil.copytree(startingdir+"/output/"+each.replace(".mzML",".pin")+"_out/crux-output",startingdir+"/output/"+each_newfile.replace(".mzML",".pin")+"_out/crux-output")
                os.chdir(startingdir+"/output/"+each_newfile.replace(".mzML",".pin")+"_out/crux-output")
                #Now we'll load each of the following files and filter them out, and rename them....
                for each_filter_file in files_to_filter:
                    this_pin_df=pandas.read_csv(each_filter_file,sep='\t')
                    this_pin_df=this_pin_df[this_pin_df['file_idx']==run_dict_reverse[each_newfile]]
                    this_pin_df.to_csv(each.rsplit(".",1)[0]+'.'+'.'.join(each_filter_file.rsplit(".",4)[1:]),sep='\t',index=False)
                os.chdir(startingdir)#

        #cleanup
        #for each_filter_file in files_to_filter:
        #    os.remove(each_filter_file)
    #os.chdir(startingdir)
    #cleanup
    #os.chdir("output/")
    #os.system("tar -cvf - combined_out/ 2>/dev/null | pigz -9 -p 24 > combined_perco.tar.gz")
    #shutil.rmtree("combined_out/")
    os.chdir(startingdir)


else:  #This assumes =="F", and we'll run percolator ONCE on the AGGREGATE DATA.
    #2. Run crux percolator for the aggregated data
    #3. Read in the outputs [psms:[targets, decoys],peptides:[targets,decoys] as a single dataframes, and then split them by file and output them to proper folders]
    #4. Copy the XML files and other outputs to each folder for completeness.
    all_max_charge,all_min_charge=find_charge_range(onlyfiles)
    print "Max charge is {0} and min charge is {1}".format(all_max_charge,all_min_charge)
    os.mkdir("output/combined_out")
    first_file=True
    with open("output/combined_out/combined_input.pin",'wb') as pin_writer:
        #print run_dict_reverse,"new"
        #print onlyfiles,"old"
        for each in onlyfiles:
            line_ctr=0
            with open(folder+"/"+each,'rb') as pin_reader:
                peptide_index=0
                index_ctr=0
                expmass_index=0
                specID_index=0
                scan_index=0
                label_index=0
                charge_indicies={}
                charges=[]
                pin_lines=[]
                for eachline in pin_reader:
                    if line_ctr==0:
                        header=eachline.split("\t")
                        for each_item in header:
                            if each_item == "Peptide":
                                peptide_index=index_ctr
                            elif each_item == "SpecId":
                                specID_index=index_ctr
                            elif each_item == "ScanNr":
                                scan_index=index_ctr
                            elif each_item == "ExpMass":
                                expmass_index = index_ctr
                            elif each_item == "Label":
                                label_index=index_ctr
                            elif "Charge" in each_item:
                                charge_indicies[each_item]=index_ctr
                                charges.append(int(each_item.replace("Charge","")))
                            index_ctr+=1
                    else:
                        break
                line_ctr=0
                charge_index_max=max(charge_indicies.values())
                charge_index_min=min(charge_indicies.values())
                charge_max=max(charges)
                charge_min=min(charges)
                add_before=[]
                add_after=[]
                add_header_before=[]
                add_header_after=[]
                add_direction_before=[]
                add_direction_after=[]
                if charge_min>all_min_charge:
                    for i in xrange(all_min_charge,charge_min):
                        add_before.append("0")
                        add_header_before.append("Charge{0}".format(i))
                        add_direction_before.append("0")
                if charge_max<all_max_charge:
                    for i in xrange(charge_max+1,all_charge_max+1):
                        add_after.append("0")
                        add_header_after.append("Charge{0}".format(i))
                        add_direction_after.append("0")
                pin_reader.seek(0)
                for each_line in pin_reader:
                    if line_ctr <2 and first_file:
                        if line_ctr==0:
                            front_half=each_line.split("\t")[:charge_index_min]
                            charges=each_line.split("\t")[charge_index_min:charge_index_max+1]
                            back_half=each_line.split("\t")[charge_index_max+1:]
                            #After this segment, ensure that line is split properly, add charges missing into the charges middle section, and then re-concat them back together and replace each_line
                            front_half.extend(add_header_before)
                            front_half.extend(charges)
                            front_half.extend(add_header_after)
                            front_half.extend(back_half)
                            each_line="\t".join(front_half)
                            
                        else:
                            front_half=each_line.split("\t")[:charge_index_min]
                            charges=each_line.split("\t")[charge_index_min:charge_index_max+1]
                            back_half=each_line.split("\t")[charge_index_max+1:]
                            #After this segment, ensure that line is split properly, add charges missing into the charges middle section, and then re-concat them back together and replace each_line
                            front_half.extend(add_direction_before)
                            front_half.extend(charges)
                            front_half.extend(add_direction_after)
                            front_half.extend(back_half)
                            each_line="\t".join(front_half)
                        print each_line#################################################
                        pin_writer.write(each_line)
                    elif line_ctr >=2:
                        split_line=each_line.split("\t")

                        front_half=each_line.split("\t")[:charge_index_min]
                        charges=each_line.split("\t")[charge_index_min:charge_index_max+1]
                        back_half=each_line.split("\t")[charge_index_max+1:]
                        #After this segment, ensure that line is split properly, add charges missing into the charges middle section, and then re-concat them back together and replace each_line
                        front_half.extend(add_before)
                        front_half.extend(charges)
                        front_half.extend(add_after)
                        front_half.extend(back_half)
                        each_line="\t".join(front_half)

                        #each_line="\t".join(split_line)
                        pin_writer.write(each_line)

                    if line_ctr == 1 and first_file:
                        first_file=False #no more headers!
                    line_ctr+=1

                #We only want the header from the first file, top two lines.  Everything else we should just write in...
                #Then we'll execute percolator on the aggregate data
                #Then we'll split it and divy it up between all the different folders per file input



            os.mkdir("output/"+each+"_out") # We'll make the folders into which we'll end up sticking the final split outputs
            os.rename(folder+"/"+each,"output/"+each+"_out/"+each) #Move each separate pin file into its own folder, too
    cmd=arguments[:]
    cmd.append("--fileroot")
    cmd.append("combined_analysis")
    cmd.append("combined_input.pin")
    os.chdir("output/combined_out")
    print "Running....",cmd
    percolator_process=subprocess.Popen(cmd)
    percolator_process.wait()
    combined_out_folder=os.getcwd()
    files_to_filter=['combined_analysis.percolator.decoy.peptides.txt','combined_analysis.percolator.target.peptides.txt','combined_analysis.percolator.target.psms.txt','combined_analysis.percolator.decoy.psms.txt']
    for each in run_dict_reverse:
        os.chdir(combined_out_folder)
        shutil.copytree("crux-output",startingdir+"/output/"+each.replace(".mzML",".pin")+"_out/crux-output")
        os.chdir(startingdir+"/output/"+each.replace(".mzML",".pin")+"_out/crux-output")
        #Now we'll load each of the following files and filter them out, and rename them....
        for each_filter_file in files_to_filter:
            print "Reading file {0}".format(each_filter_file)
            this_pin_df=pandas.read_csv(each_filter_file,sep='\t')
            this_pin_df=this_pin_df[this_pin_df['file_idx']==run_dict_reverse[each.rsplit(".",1)[0]+".mzML"]]
            this_pin_df.to_csv(each.rsplit(".",1)[0]+'.'+'.'.join(each_filter_file.rsplit(".",4)[1:]),sep='\t',index=False)
            os.remove(each_filter_file)
    os.chdir(startingdir)
    #cleanup
    os.chdir("output/")
    os.system("tar -cvf - combined_out/ 2>/dev/null | pigz -9 -p 24 > combined_perco.tar.gz")
    shutil.rmtree("combined_out/")
    os.chdir(startingdir)
    #sys.exit(2)

    #First, we'll combine all the pin files together.
    #for each in onlyfiles:


for neweach in run_dict_reverse:
    each=neweach.replace(".mzML",".pin")
    if os.path.exists("output/"+each+"_out/"+each):
        pass #it already exists!
    else:
        #Let's go find the pin data...
        #first, what run group is this in?
        this_group=run_to_group_dict[neweach]
        if not os.path.isdir("output/"+each+"_out/"):
            os.mkdir("output/"+each+"_out/")
        with open("output/"+each+"_out/"+each,'wb') as filewriter:
            with open("output/"+this_group+".pin_out/"+this_group+".pin",'rb') as filereader:
                linectr=0
                peptide_index=0
                index_ctr=0
                expmass_index=0
                specID_index=0
                scan_index=0
                label_index=0
                charge_indicies={}
                pin_lines=[]
                for eachline in filereader:
                    if linectr==0:
                        header=eachline.split("\t")
                        for each_item in header:
                            if each_item == "Peptide":
                                peptide_index=index_ctr
                            elif each_item == "SpecId":
                                specID_index=index_ctr
                            elif each_item == "ScanNr":
                                scan_index=index_ctr
                            elif each_item == "ExpMass":
                                expmass_index = index_ctr
                            elif each_item == "Label":
                                label_index=index_ctr
                            elif "Charge" in each_item:
                                charge_indicies[each_item]=index_ctr
                            index_ctr+=1
                        #filewriter.write("SpecId\tlabel\tScanNr\tPeptide\tcharge\n")
                        filewriter.write(eachline)
                    elif linectr==1:
                        filewriter.write(eachline)
                    else:
                        #print eachline[specID_index]
                        if eachline.split("\t")[specID_index].split("_")[1]==str(run_dict_reverse[neweach]):
                            filewriter.write(eachline)
                    linectr+=1



#    os.system("tar -cvf - combined_out/ 2>/dev/null | pigz -9 -p 24 > combined_perco.tar.gz")
for each in onlyfiles:
    #if each not in [x.replace(".mzML",".pin") for x in run_dict_reverse.keys()]:
    if each.replace(".pin",".mzML") not in run_dict_reverse:
        #If it's not there, then we're going to mask the pin file into a tar.gz
        os.system("tar -cvf - output/{0}_out/{0} 2>/dev/null | pigz -9 -p 24 > output/{1}.tar.gz".format(each,each.replace(".pin","")))
        shutil.rmtree("output/{0}".format(each+"_out"))

        #os.mkdir("output/"+each+"_out")
        #os.rename(folder+"/"+each,"output/"+each+"_out/"+each)
        #os.chdir("output/"+each+"_out")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


################ To prepare for percolator file output, we're going to read in the percolator pin files
################ Once they're read in, we'll have to generate a unique ID of file_idx+"_"+scan to merge the pin sequence on
################ We're only going to bring in the sequence from the percolator input file.
################ ===== As of crux version 3.1, we will also begin to fix the charge column...  We also correct the spectral m/z column
os.chdir(startingdir)
pin_list=[]

spectral_expmz={}
#for each in onlyfiles:
for neweach in run_dict_reverse:
    each=neweach.replace(".mzML",".pin")
    os.chdir("output/"+each+"_out")
    with open("temp.tsv",'wb') as filewriter:
        with open(each,'rb') as filereader:
            linectr=0
            peptide_index=0
            index_ctr=0
            expmass_index=0
            specID_index=0
            scan_index=0
            label_index=0
            charge_indicies={}
            pin_lines=[]
            for eachline in filereader:
                if linectr==0:
                    header=eachline.split("\t")
                    for each_item in header:
                        if each_item == "Peptide":
                            peptide_index=index_ctr
                        elif each_item == "SpecId":
                            specID_index=index_ctr
                        elif each_item == "ScanNr":
                            scan_index=index_ctr
                        elif each_item == "ExpMass":
                            expmass_index = index_ctr
                        elif each_item == "Label":
                            label_index=index_ctr
                        elif "Charge" in each_item:
                            charge_indicies[each_item]=index_ctr
                        index_ctr+=1
                    filewriter.write("SpecId\tlabel\tScanNr\tPeptide\tcharge\n")
                elif linectr==1:
                    pass
                else:
                    thisline=eachline.split("\t")
                    peptide = thisline[peptide_index].split(".",1)[1]
                    peptide = peptide.rsplit(".",1)[0]
                    this_run=run_dict[thisline[specID_index].rsplit("_",4)[1]].rsplit(".",1)[0]
                    charge=0
                    for each_charge in charge_indicies:
                        if thisline[charge_indicies[each_charge]]=="1":
                            charge=each_charge.replace("Charge","")
                    if charge == 0 or charge =="0":
                        print "We couldnt find a charge for "+str(thisline[specID_index])+" so its defaulting to one! This should raise suspicions!"
                        charge="1"
                    
                    if this_run not in spectral_expmz:
                        spectral_expmz[this_run]={}
                    #print thisline,"charge is",charge
                    spectral_expmz[this_run][int(thisline[scan_index])]=str((float(thisline[expmass_index])+((float(charge)-1)*proton_mass))/float(charge))

                    filewriter.write(thisline[specID_index]+"\t"+thisline[label_index]+"\t"+thisline[scan_index]+"\t"+peptide+"\t"+str(charge)+"\n")
                linectr+=1

    new_in=pandas.read_csv("temp.tsv",sep='\t')
    os.remove("temp.tsv")
    #print new_in,"new csv..."
    #new_in['label']=new_in['label'].astype(int)
    new_in=new_in[new_in['label']==1]
    pin_list.append(new_in)
    os.chdir(startingdir)

print "Done fixing pin files..."

#target_1_58_58_2
def makeUniqueID(x):
    return x['SpecId'].split("_")[1]+"_"+str(x['ScanNr'])

def makeUniqueIDpsms(x):
    return str(x['file_idx'])+"_"+str(x['scan'])

pin_megaframe=pandas.concat(pin_list)
unique_names=pin_megaframe.apply(makeUniqueID,axis=1)
pin_megaframe['unique_name']=unique_names
pin_megaframe.drop('SpecId',axis=1,inplace=True)
pin_megaframe.drop('ScanNr',axis=1,inplace=True)
pin_megaframe.drop('label',axis=1,inplace=True)


os.chdir(startingdir)
#for each in onlyfiles:
for neweach in run_dict_reverse:
    each=neweach.replace(".mzML",".pin")

    os.chdir("output/"+each+"_out")
    target_psms=[ f for f in listdir("crux-output/") if (isfile(join("crux-output/",f)) and "target.psms.txt" in f)]
    decoy_psms=[ f for f in listdir("crux-output/") if (isfile(join("crux-output/",f)) and "decoy.psms.txt" in f)]
    #target_peptides=[ f for f in listdir("crux-output/") if (isfile(join("crux-output/",f)) and "target.peptides.txt" in f)]
    #decoy_peptides=[ f for f in listdir("crux-output/") if (isfile(join("crux-output/",f)) and "decoy.peptides.txt" in f)]
    #os.chdir("crux-output")
    corrected=[]
    for eachfile in target_psms:
        print "About to read in ",eachfile
        df=pandas.read_csv("crux-output/"+eachfile,sep='\t')
        df['file_idx']=df['file_idx'].astype(int)
        #print df,"this is the dataframe! <--------------------------"
        fileset=set(df['file_idx'])
        for eachone in fileset:
            eachint=int(eachone)
            mask=df[df['file_idx'] == eachone]
            mask['file']=run_dict[str(eachone)]
            #print run_dict[str(eachone)]
            corrected.append(mask)
            #print mask,"this is mask for",str(eachone)
    #print "I was handling ",each
    #print corrected,"the frame itself"
    corrected_df=pandas.concat(corrected)
    corrected_df.drop('charge',axis=1,inplace=True)
    unique_names_psms=corrected_df.apply(makeUniqueIDpsms,axis=1)
    corrected_df['unique_name']=unique_names_psms
    merged_df=corrected_df.merge(pin_megaframe,how='inner',on='unique_name')#
    corrected_df=merged_df
    corrected_df.drop('unique_name',axis=1,inplace=True)
    corrected_df.drop('sequence',axis=1,inplace=True)
    corrected_df.rename(columns={'Peptide':'sequence'},inplace=True)


    #Fixing precursor spectral m/z for targets
    corrected_df['scan']=corrected_df['scan'].astype(int)
    for index,each_group in corrected_df.groupby('file'):
        corrected_df.loc[corrected_df['file']==index,'spectrum precursor m/z']=corrected_df.loc[corrected_df['file']==index,'scan'].map(spectral_expmz[index.rsplit(".",1)[0]])



    #HERE WE CAN MAKE UNIQUE IDS
    #AND THEN DO A LEFT JOIN TO TAKE IN THE OLD (good) SEQUENCES
    #AND THEN GET RID OF THE BAD SEQUENCES AND REPLACE THEM WITH THE GOOD
    #SEQUENCES.
    corrected_df.to_csv("crux-output/"+eachfile,sep='\t',index=False)
    corrected_decoys=[]
    for eachfile in decoy_psms:
        df=pandas.read_csv("crux-output/"+eachfile,sep='\t')
        df['file_idx']=df['file_idx'].astype(int)
        #print df,"this is the dataframe!"
        fileset=set(df['file_idx'])
        for eachone in fileset:
            eachint=int(eachone)
            mask=df[df['file_idx'] == eachone]
            mask['file']=run_dict[str(eachone)]
            corrected_decoys.append(mask)
            #print mask,"this is mask for",str(eachone)

    #Fixing precursor spectral m/z for decoys
    corrected_df['scan']=corrected_df['scan'].astype(int)
    for index,each_group in corrected_df.groupby('file'):
        corrected_df.loc[corrected_df['file']==index,'spectrum precursor m/z']=corrected_df.loc[corrected_df['file']==index,'scan'].map(spectral_expmz[index.rsplit(".",1)[0]])
    corrected_df=pandas.concat(corrected_decoys)






    corrected_df.to_csv("crux-output/"+eachfile,sep='\t',index=False)
    #this_run=



    #onlyfiles = [ f for f in listdir(folder) if (isfile(join(folder,f)) and ".txt" in f)]
    os.chdir(startingdir)
