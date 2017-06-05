import subprocess, os
#import shutil
from subprocess import PIPE

import sys

class Fido(object):
    def __init__(self, args):
        #print args,"args!!!"
        infile = args[1]
        outfile = args[0]
        targetdecoy=args[2]
        
        
        options = self.parseOptions(args[3:])
        #print "so my options are... ",options
        self.runFido(options, infile, outfile, targetdecoy)

    def runFido(self, options, infile, outfile, targetdecoy):
        command = ['/galaxy-central/tools/wohl-proteomics/fido/FidoChooseParameters']
        command.append(options)
        #print command," just after append..."
        #shutil.copyfile(infile,"fidographchoose")
        command.append(infile)
        command.append(targetdecoy)
        #print command,"this is the final thing..."
        #command.append(" > ")
        fixed_command=[]
        for each in command:
            for each_one in each.split(" "):
                fixed_command.append(each_one)
        print "About to run",fixed_command
        print "instead of",command
        
        #error = subprocess.call(command, stdout=output2, stderr=sys.stdout)
        program = subprocess.Popen(fixed_command, stdout=PIPE, stderr=PIPE)
        output, err = program.communicate()
        returncode = program.returncode
        print output,"this was output..."
        print err,"this was err..."
        #print error,"is this the output?"
        splits = err[-30:].split()
        #print err[-30:].split(),"this is the split of the last 30 chars"
        #print a[-3:],"the last three... will this do?"
        #lastline=err[-13:-1]
        
        print splits,"this is my splits!"
        #splits=["0.5","0.1","0.01"]
        gamma,alpha,beta=splits[-3],splits[-2],splits[-1]
        print "gamma=",gamma," alpha=",alpha," beta=",beta,
        command = ['/galaxy-central/tools/wohl-proteomics/fido/Fido']
        #shutil.copyfile(infile,"fidograph")
        command.append(infile)
        command.append(gamma)
        command.append(alpha)       
        command.append(beta)
        fidooutput = open(outfile, 'w')
        error = subprocess.call(command, stdout=fidooutput, stderr=sys.stdout)       
        fidooutput.close()
        
        
    def parseOptions(self, options):
        #options = [x.split('=') for x in options]
        #options = dict([x.split('=') for x in options])
        
        # take out the options that have not been filled in
            
        # put dict back in a list
        #print options, " I'm starting with this..."
        #print type(options)
        options_out = ' '.join(options)
        
        return options_out

def main():
    fido = Fido(sys.argv[1:])

if __name__ == '__main__':
    main()
