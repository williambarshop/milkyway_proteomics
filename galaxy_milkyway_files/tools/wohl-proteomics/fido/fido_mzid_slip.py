import pyteomics
from pyteomics import mzid
import sys, os


arguments=sys.argv

arguments=arguments[1:]

fido_groups=arguments[0]
input_mzid=arguments[1]

mzid_file=mzid.MzIdentML(input_mzid,build_id_cache=True)

etree=mzid_file.build_tree

root=etree.getroot()
for each in root:
    print each.tag,each.attrib

#for each in etree:
#    print each
