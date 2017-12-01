'''
2017-12-01
Author: WILLIAM BARSHOP
LABORATORY OF JAMES A WOHLSCHLEGEL

FIDO-Q SCRIPT.
This script will take the Fido raw output and produce Baysen "q-values" (pFDR) values for each protein!	


Usage: python fido_q.py $fido_output $decoy_string $output <"true" or "false" for FidoCT>
Caution... This script has some old code... Sorry.
'''

import sys, os,time

arguments = sys.argv[1:]
input_file=arguments[0]
decoy_string=arguments[1]
output_file=arguments[2]

useEmpirical=False
if "--useEmpirical" in arguments[3]:
   useEmpirical=True
fidoCT=False
clip_zero=False
if arguments[4]=="true":
    fidoCT=True
if len(arguments)>4:
    if arguments[5]=="true":
        clip_zero=True



ptoq_file="p_to_q_dict.csv"

def recalculateAverage(avg,new,count): #Look up this as "average on the go/fly" type calculations
    return ((avg * float(count) + new) / float(count + 1))





print "Decoy string is: ",decoy_string

pi0=1.0
#Let's calculate pi0 (ratio of targets to decoys)
target_set=set()
decoy_set=set()
header=True
with open(input_file,'r') as input_reader:
    for eachline in input_reader:
        if header:
            header=False
            continue
        #print eachline
        protein_name_list=eachline.split("{")[1].split("}")[0].strip().split(",")
        for protein_name in protein_name_list:
            #print protein_name
            if decoy_string in protein_name:
                decoy_set.add(protein_name)
            else:
                target_set.add(protein_name)
if len(decoy_set) <1:
    print "No decoys in the dataset... Can't estimate pi0!"
    sys.exit(1)
pi0=float(len(target_set))/float(len(decoy_set))
print "pi0 is set at: ",pi0," from ",str(len(target_set))," targets and ",str(len(decoy_set))," decoys..."




probability_dict={}
with open(input_file,'r') as openfile:
    for eachline in openfile:
        if not fidoCT:
            cap=eachline.split(" ",1)
        else:
            cap=eachline.split("\t",1)
        if "Posterior" in cap[0]: #This ensures compatibility with FidoCT output...
            continue
        if probability_dict.has_key(float(cap[0])):
            #print probability_dict[float(cap[0])]
            extend_with=cap[1].strip().replace("{","").replace("}","").strip().split(",")
            if not type(extend_with)==list:
                extend_with=[extend_with]
            for k in xrange(0,len(extend_with)):
                extend_with[k]=extend_with[k].strip()
            probability_dict[float(cap[0])].append(extend_with)
            #print "NOW IT'S ",probability_dict[float(cap[0])]
            pass
        else:            
            fixed=cap[1].strip().replace("{","").replace("}","").strip().split(",")
            if type(fixed)==list:
                #print type(fixed)
                pass
            else:
                fixed=[fixed]
                #print fixed
            for each in xrange(0,len(fixed)):
                fixed[each]=fixed[each].strip()
            probability_dict[float(cap[0])]=[fixed]






for each in probability_dict:
    a=probability_dict[each]
    #print a
    for each2 in a:
        count_decoy=0
        for each3 in each2:
            #print each3
            #print type(each3)
            if "Reverse_" in each3:
                #print "FOUND A DECOY!",each3
                count_decoy+=1
    #    if count_decoy<len(each2) and (count_decoy > 0):
    #        print(each2)
    #        print len(each2)," and ",count_decoy

    #if count_decoy<len(a) and (count_decoy > 0):
    #    print a
    #    print len(a)," and ",count_decoy
      
    #if len(a)>1:
    #    print len(a)," length"
    #    print a
    #    print "-------------------"



#################CALCULATION AS BAYSEAN pFDR q-value, NOT EMPIRICAL ##################
ptoq={}
prot_to_q={}
j=0
current=None
for each_probability in sorted(probability_dict.keys(),reverse=True):
    error_prob=1-each_probability
    prot_groups_at_prob=probability_dict[each_probability]

    if current is None:
        current = error_prob
    else:
        current = recalculateAverage(current,error_prob,j)

    ptoq[each_probability]=current

    for each_group in prot_groups_at_prob:
        if "list" in str(type(each_group)):
            for each_prot in each_group:
                prot_to_q[each_prot]=current
        else:
            prot_to_q[each_group]=current
    j+=1




empirical_ptoq={}
ptogFDR={}
qtogFDR={}
targetcount = 0
decoycount = 0
for each_probability in sorted(probability_dict.keys(),reverse=True):
    prob=each_probability
    prot_groups_at_prob=probability_dict[prob]
    for eachgroup in prot_groups_at_prob:
        decoy=True
        for eachprot in eachgroup:
            if decoy_string not in eachprot:
                decoy=False
                break
        if decoy:
            decoycount+=1
        else:
            targetcount+=1
    
    #current_global_FDR = (2.0*float(decoycount))/float(targetcount+decoycount) # FOR COMBINED TARGET-DECOY SEARCHES, UNADJUSTED
    #current_global_FDR = (float(decoycount))/float(targetcount+decoycount) # FOR SEPARATED TARGET-DECOY SEARCHES, UNADJUSTED
    current_global_FDR = (pi0*float(decoycount)/float(targetcount+decoycount)) # FOR SEPARATED TARGET-DECOY SEARCHES, ADJUSTED BY PI0 FOR UNEQUAL TARGET/DECOY POPULATION SIZES
    if targetcount>0:
        current_q = (pi0*float(decoycount)/float(targetcount+decoycount))#right now... this is gFDR...
        if current_q > 1.0:
            current_q=1.0
    else:
        current_q=0.0
    empirical_ptoq[prob]=current_q
    ptogFDR[prob]=current_global_FDR
    qtogFDR[current_q]=current_global_FDR
    #print "at ",prob," gFDR is ",current_global_FDR," q is ",current_q," numprot=",str(targetcount+decoycount)," numtarget=",targetcount," numdecoy=",decoycount

current_q=None
for each_probability in sorted(empirical_ptoq.keys(),reverse=False):
    #print "current prob:",prob
    prob=each_probability
    if current_q is None:
        current_q=empirical_ptoq[prob]
    elif current_q > empirical_ptoq[prob]:
        current_q=empirical_ptoq[prob]
    elif current_q < empirical_ptoq[prob]:
        #print "swapping",empirical_ptoq[prob],"with",current_q
        empirical_ptoq[prob]=current_q
    else:
        continue



reconstructed_probdict={}
ptoq_dict_file=open(ptoq_file,'wb')
ptoq_dict_file.write("FIDO Probability,q-value,statistical q-value,empirical q-value,protein group\n")
with open(output_file,'w') as writer:
    writer.write("q-value\tstatistical q-value\tempirical q-value\tprotein group\n")
    for each_probability in sorted(probability_dict.keys(),reverse=True):
        prot_groups_at_prob=probability_dict[each_probability]
        reconstructed_probdict[each_probability]=[]
        for eachgroup in prot_groups_at_prob:
            decoy=True
            for eachprot in eachgroup:
                if decoy_string not in eachprot:
                    decoy=False
                    break
            newlist=[]
            for eachprot in eachgroup:
                if not decoy:
                    if decoy_string in eachprot:
                        pass
                    else:
                        newlist.append(eachprot)   #Mixed ambiguous target/decoy groups go to the target.
            if decoy:
                newlist=eachgroup
            reconstructed_probdict[each_probability].append(newlist)
    for each_probability in sorted(reconstructed_probdict.keys(),reverse=True):
        prot_groups_at_prob=probability_dict[each_probability]
        if clip_zero and each_probability<=0.0:
            continue

        for eachgroup in prot_groups_at_prob:
            if useEmpirical:
                writer.write(str(empirical_ptoq[each_probability])+"\t"+str(ptoq[each_probability])+"\t"+str(empirical_ptoq[each_probability])+"\t"+",".join(eachgroup)+"\n")
                ptoq_dict_file.write(','.join([str(each_probability),str(empirical_ptoq[each_probability]),str(ptoq[each_probability]),str(empirical_ptoq[each_probability]),'\t'.join(eachgroup)])+'\n')
            else:
                writer.write(str(ptoq[each_probability])+"\t"+str(ptoq[each_probability])+"\t"+str(empirical_ptoq[each_probability])+"\t"+",".join(eachgroup)+"\n")
                ptoq_dict_file.write(','.join([str(each_probability),str(ptoq[each_probability]),str(ptoq[each_probability]),str(empirical_ptoq[each_probability]),'\t'.join(eachgroup)])+'\n')

            #print str(empirical_ptoq[each_probability])
            #print ','.join([str(each_probability),str(empirical_ptoq[each_probability]),'\t'.join(eachgroup)])+'\n'
    #for each_probability in sorted(probability_dict.keys(),reverse=True):
    #    prot_groups_at_prob=reconstructed_probdict[each_probability]
    #    
    #    #print reconstructed_probdict[each_probability]," at ",each_probability
#    #Maybe iterate over all the proteins and do a protein-Q file?  Maybe just make a data frame in the for loop above with all the data and dump it to disk

ptoq_dict_file.close()
