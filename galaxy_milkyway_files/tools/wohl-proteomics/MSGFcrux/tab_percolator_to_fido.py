import sys,os
import pandas
from Bio import SeqIO

"""
Converter script for Percolator output to create Fido-readable files.

By Jose Fernandez Navarro,
SLIGHTLY MODIFIED BY WILLIAM BARSHOP (LABORATORY OF JAMES A WOHLSCHLEGEL)
These modifications allow input of multiple tab delimited peptide-level percolator outputs.  These outputs will
be merged into a pandas dataframe.  Repeated peptides will have the best (1-PEP) used for Fido input.
-- Also have added in support for All PSMs mode...
2016-09-12

usage (from galaxy):
tab_percolator_to_fido.py percolator_output_parent_folder/ $target_fasta $decoy_fasta $graph_output $target_decoy_output <true/false for all_scores> <true/false for psm_level>
"""

import sys

from lxml import etree

class Hit(object):
    """Wrapper class"""
    def __init__(self, isdecoy=False,score=0.0,pep=0.0,q=0.0,p=0.0,name="",proteins=[]):
        self.isdecoy = isdecoy
        self.score = score
        self.pep = pep
        self.q = q
        self.p = p
        self.name = name
        self.proteins = proteins

class Perco2Fido(object):
    def __init__(self, args):
        in_folder = args[0]
        infiles=[]
        outfile = args[3]
        targetfasta = args[1]
        decoyfasta = args[2]
        tdoutput = args[4]
        all_scores=False
        if "true" in args[5] or "True" in args[5] or "TRUE" in args[5]:
            all_scores=True
        psm_level=False
        if "psms" in args[6] or "PSMS" in args[6] or "PSM" in args[6] or "psm" in args[6] or "true" in args[6] or "TRUE" in args[6]:
            psm_level=True



        targetfasta_seq = SeqIO.parse(open(targetfasta,'rb'),'fasta')        
        decoyfasta_seq = SeqIO.parse(open(decoyfasta,'rb'),'fasta')        

        
        targets=['{',' ']
        decoys=['{',' ']
        for each in targetfasta_seq:
            name, seq=each.id,str(each.seq)
            targets.append(each.id)
            targets.append(" , ")
        targets=targets[:-1]
        targets.append(" ")
        targets.append("}")

        for each in decoyfasta_seq:
            name, seq=each.id,str(each.seq)
            decoys.append(each.id)
            decoys.append(" , ")
        decoys=decoys[:-1]
        decoys.append(" ")
        decoys.append("}")

        with open(tdoutput,'w') as output:
            output.write(''.join(targets)+"\n")
            output.write(''.join(decoys)+"\n")
        





        for root, subFolders, files in os.walk(in_folder):
            for eachfile in files:
                if 'peptides.txt' in eachfile and not psm_level:
                    infiles.append(str(os.path.join(root,eachfile)))
                elif 'psms.txt' in eachfile and psm_level:
                    infiles.append(str(os.path.join(root,eachfile)))

        dataframe_vector=[]
        for eachfile in infiles:
            #newdf=pandas.DataFrame.from_csv(eachfile,sep='\t',index_col=False)
            newdf=pandas.read_csv(eachfile,sep='\t',index_col=False)
            if 'decoy' in eachfile:
                decoy=1
            else:
                decoy=0
            newdf.insert(0,'decoy',decoy)
            dataframe_vector.append(newdf)
    
        percolator_results=pandas.concat(dataframe_vector)
        #print dataframe_vector[0].columns.values
        del dataframe_vector
        peptides, pepunique = self.readPercolatorPeptides(percolator_results)
        if all_scores:
            self.writeFidoInputAllScores(peptides, outfile)
        else:
            self.writeFidoInput(pepunique, outfile)
        
        
    def writeFidoInput(self, peptides,fidoOutputFile):#,fidoOutputFile2): # current released fido only outputs one file
    #   graph input looks like this:
    #    e EEEMPEPK
    #    r SW:TRP6_HUMAN
    #    r GP:AJ271067_1
    #    r GP:AJ271068_1
    #    p 0.9849
    #    e LLEIIQVR
    #    r SW:TRP6_HUMAN
    #    r GP:AJ271067_1
    #    r GP:AJ271068_1
    #    p 0.0
    #    e FAFNNKPNLEWNWK
    #    r gi|1574458|gb|AAC23247.1|
    #    p 0.9750

    #{ target1 , target2 , target3 , ... }
    #{ decoy1 , decoy2 , decoy3 , ... }

        f = open(fidoOutputFile, "w")
        for peptide_key in peptides:
            peptide=peptides[peptide_key]
            #if(not peptide.isdecoy):
            prots = peptide.proteins
            pepname = peptide.name
            pepprob = 1 - peptide.pep
            f.write("e " + pepname + "\n")
            for prot in prots:
                f.write("r " + prot + "\n")
            f.write("p " + str(pepprob) + "\n")    
        f.close()

        for peptide_key in peptides:
            peptide=peptides[peptide_key]
            #if(not peptide.isdecoy):
            prots = peptide.proteins


        '''
        for peptide in peptides:
            #if(not peptide.isdecoy):
            prots = peptide.proteins
            pepname = peptide.name
            pepprob = 1 - peptide.pep
            f.write("e " + pepname + "\n")
            for prot in prots:
                f.write("r " + prot + "\n")
            f.write("p " + str(pepprob) + "\n")    
        f.close()
        '''
    def writeFidoInputAllScores(self, peptide_list,fidoOutputFile):#,fidoOutputFile2): # current released fido only outputs one file
    #   graph input looks like this:
    #    e EEEMPEPK
    #    r SW:TRP6_HUMAN
    #    r GP:AJ271067_1
    #    r GP:AJ271068_1
    #    p 0.9849
    #    e LLEIIQVR
    #    r SW:TRP6_HUMAN
    #    r GP:AJ271067_1
    #    r GP:AJ271068_1
    #    p 0.0
    #    e FAFNNKPNLEWNWK
    #    r gi|1574458|gb|AAC23247.1|
    #    p 0.9750

    #{ target1 , target2 , target3 , ... }
    #{ decoy1 , decoy2 , decoy3 , ... }
        

        f = open(fidoOutputFile, "w")
        for peptide in peptide_list:
            #peptide=peptides[peptide_key]
            #if(not peptide.isdecoy):
            prots = peptide.proteins
            pepname = peptide.name
            pepprob = 1 - peptide.pep
            f.write("e " + pepname + "\n")
            for prot in prots:
                f.write("r " + prot + "\n")
            f.write("p " + str(pepprob) + "\n")    
        f.close()

        for peptide in peptide_list:
            prots = peptide.proteins


        '''
        for peptide in peptides:
            #if(not peptide.isdecoy):
            prots = peptide.proteins
            pepname = peptide.name
            pepprob = 1 - peptide.pep
            f.write("e " + pepname + "\n")
            for prot in prots:
                f.write("r " + prot + "\n")
            f.write("p " + str(pepprob) + "\n")    
        f.close()
        '''
######## implement this possibly with a newer fido version:
#        f = open(fidoOutputFile2, "w")
#        f.write("{ ")
#        for protein in [x for x in proteins if not x.isdecoy]:
#            f.write(printable(protein.name) + " , ")
#        f.write("}\n")
#        f.write("{ ")
#        for protein in [x for x in proteins if x.isdecoy]:
#            f.write(printable(protein.name) + " , ")
#        f.write("}\n")
#        f.close()



    def readPercolatorPSM(self, elems):
        ##process percolator Peptides,ptms and proteins
        percolatorPSMs = []
        percolatorPSMdict = dict()
        for elem in elems.iter("{http://per-colator.com/percolator_out/13}psm"):
            decoy = True   
            if (elem.get("{http://per-colator.com/percolator_out/13}decoy") == "false"):
                decoy = False
            score = elem.findtext("{http://per-colator.com/percolator_out/13}svm_score")
            pep = elem.findtext("{http://per-colator.com/percolator_out/13}pep")
            q = elem.findtext("{http://per-colator.com/percolator_out/13}q_value")
            p = elem.findtext("{http://per-colator.com/percolator_out/13}q_value")
            name = elem.findall("{http://per-colator.com/percolator_out/13}peptide_seq")[0].get("seq")
            percolatorPSMs.append(Hit(decoy, float(score), float(pep), float(q), float(p), str(name), []))
            #if(not decoy): #REMOVED TO ALLOW REVERSE HITS FOR FIDO!
            if(not percolatorPSMdict.has_key(name)):
                percolatorPSMdict[name] = Hit(decoy, float(score), float(pep), float(q), float(p), str(name), [])
            elif(percolatorPSMdict[name].pep > float(pep)):
                percolatorPSMdict[name].pep = float(pep)
                percolatorPSMdict[name].score = float(score)
                percolatorPSMdict[name].q = float(q)
                percolatorPSMdict[name].p = float(p)
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
           
            return percolatorPSMs, percolatorPSMdict


    def readPercolatorPeptides(self, perco_dataframe):
        percolatorPeptides = []
        peptideUnique = dict()
        for index,eachitem in perco_dataframe.iterrows():
            #print eachitem,type(eachitem),"eachitem"
            decoy = True
            if (eachitem['decoy']==0):
                decoy = False
            score = eachitem['percolator score']
            pep = eachitem['percolator PEP']
            q = eachitem['percolator q-value']
            p = eachitem['percolator q-value']
            name = eachitem['sequence']
            proteins = eachitem['protein id'].split(",")
            percolatorPeptides.append(Hit(decoy,float(score),float(pep),float(q),float(p),str(name),proteins))

            #if(not decoy): #CHANGED TO ALLOW DECOYS THROUGH THE ANALYSIS
            if(not peptideUnique.has_key(str(name)) ):
                peptideUnique[str(name)] = Hit(decoy,float(score),float(pep),float(q),float(p),str(name),proteins)
            else:
                #print "Repeated peptide : " + name + " Found in Percolator files"
                #print "Checking PEP score... and decoy status"
                if (peptideUnique[str(name)].isdecoy and not decoy) or (decoy and not peptideUnique[str(name)].isdecoy):
                    prot_combined=peptideUnique[str(name)].proteins
                    prot_combined.extend(proteins)
                    #print prot_combined," after extension"
                    peptideUnique[str(name)].proteins=prot_combined
                    peptideUnique[str(name)].isdecoy=False
                    proteins=prot_combined

                if peptideUnique[str(name)].pep < pep:
                    #print "Replacing "+str(name)+"pep score with",pep
                    peptideUnique[str(name)].pep=float(pep)
                    peptideUnique[str(name)].q=float(q)
                    peptideUnique[str(name)].p=float(p)
                    peptideUnique[str(name)].score=float(score)
                    #c=set(peptideUnique[str(name)].proteins).union(set(proteins))
                    #d=set(peptideUnique[str(name)].proteins).intersection(set(proteins))
                    ##if len(list(set(peptideUnique[str(name)].proteins) - set(proteins)))>0:
                    #if len(list(c-d))>0:
                    #    print "THIS IS STRANGE-- PROTEIN LISTS WERE DIFFERENT IN LENGTH...."
                    #    print list(c-d),"differences"
                    #    print peptideUnique[str(name)].proteins,"CURRENT"
                    #    print proteins,"TO BE"
                    #peptideUnique[str(name)] = Hit(decoy,float(score),float(pep),float(q),float(p),str(name),proteins)
                else:
                    pass                    

            
        return percolatorPeptides,peptideUnique

    def readPercolatorProteins(self, elems):
        percolatorProteins = []
        for elem in elems.iter("{http://per-colator.com/percolator_out/13}protein"):
            decoy = True   
            if (elem.get("{http://per-colator.com/percolator_out/13}decoy") == "false"):
                decoy = False
            pep = elem.findtext("{http://per-colator.com/percolator_out/13}pep")
            q = elem.findtext("{http://per-colator.com/percolator_out/13}q_value")
            p = elem.findtext("{http://per-colator.com/percolator_out/13}p_value")
            name = elem.get("{http://per-colator.com/percolator_out/13}protein_id")
            qmp = elem.findtext("{http://per-colator.com/percolator_out/13}q_value_emp")
            if (not qmp):
                qmp = 0
            percolatorProteins.append(Hit(decoy,float(qmp),float(pep),float(q),float(p),str(name),[]))
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
       
        return percolatorProteins

def main():
    convert = Perco2Fido(sys.argv[1:])
        
if __name__ == '__main__':
    main()

