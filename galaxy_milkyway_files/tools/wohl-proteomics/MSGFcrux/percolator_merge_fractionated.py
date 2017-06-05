import sys, os, subprocess
from os import listdir
from os.path import isfile, join

arguments=sys.argv
print arguments
arguments=arguments[1:]
folder=arguments[0]
arguments=arguments[1:]

onlyfiles = [ f for f in listdir(folder) if isfile(join(folder,f)) ]
d=os.path.join(os.getcwd(),folder)
onlyfolders = [os.path.join(d,o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
everythingbase = [ f for f in listdir(os.getcwd())]
print everythingbase,"These should be everything here..."
print os.getcwd()
startingdir=os.getcwd()

os.chdir(join(startingdir,folder))
basedir=os.getcwd()

xml_files=[]
for eachfolder in onlyfolders:
    if eachfolder != "":
        print "EACH FOLDER IS:",eachfolder
        os.chdir(os.path.join(eachfolder,"crux-output/"))
        onlyxmlfiles=[ f for f in listdir(os.getcwd()) if (isfile(os.path.join(os.getcwd(),f)) and ".pout.xml" in f) ]
        for each in onlyxmlfiles:
            xml_files.append(os.getcwd()+"/"+each)
    os.chdir(basedir)


os.chdir(startingdir)

length_xml=len(xml_files)
i=0
with open('script_file','wb') as writer:
    for each in xml_files:
        if i<length_xml-1:
            writer.write(each+"\n")
            print each,"withnewline"
        else:
            writer.write(each)
            print each,"withOUTnewline"
        i+=1


