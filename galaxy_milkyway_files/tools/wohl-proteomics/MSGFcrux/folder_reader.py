import os
import sys
from os import listdir
from os.path import isfile, join




onlyfiles = [ f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1],f)) ]

out=""
#first=True
for each in onlyfiles:
#    if first:
#        first=False
#        continue
    out+="output/"+each+","

print out[:-1]
