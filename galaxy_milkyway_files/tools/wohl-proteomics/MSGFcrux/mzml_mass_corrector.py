import sys
import pandas as pd
import numpy as np
from pyopenms import *

#argument structure will be python mzml_mass_corrector.py mass_correction_file filename_in minmslevel maxmslevel
arguments=sys.argv[1:]
corrections=pd.read_csv(arguments[0],index_col=0,names=['ppm_correction']) #this is the mass corrections file.
filename_in = arguments[1]
minMSlevel = int(arguments[2])
maxMSlevel = int(arguments[3])
filename_out = "output/tempout.mzML"



file = MzMLFile()
exp = MSExperiment()
options = file.getOptions()
options.setCompression(True)
#options.setMSLevels([minMSlevel,maxMSlevel])
file.setOptions(options)


ppm_correction=corrections.loc[filename_in,'ppm_correction']
actual_correct=float(1/(1.0+(ppm_correction/(10**6))))
file.load(filename_in, exp)

#print exp.getExperimentalSettings()
outexp = pyopenms.MSExperiment()
#print "Made outexp..."
#print hold,"hold"
#print "iterating over hold...."
#for each in hold:
#    print each

#exp.setExperimentalSettings(hold)

print "We'll correct",filename_in,"by",str(ppm_correction),"ppm"




#file.store(filename_out+"NO_CHANGE", exp)


for spectrum_num in xrange(0,exp.size()):
   spectrum=exp[spectrum_num]
   if spectrum.getMSLevel()<minMSlevel or spectrum.getMSLevel()>maxMSlevel:
      pass
   else:
      spec_array = spectrum.get_peaks()
      #print spec_array[0]
      #print type(spec_array)
      spec_tpose=spec_array[0]
      spec_tpose=spec_tpose.transpose()
      #spec_tpose[0]=spec_tpose[0]*float(1/(1.0+(ppm_correction/(10**6)))) #vectorized_massCorrect(spec_tpose[0],filePPMshift) # CHANGE THIS TO DOT...
      spec_tpose=np.dot(spec_tpose,actual_correct) #vectorized_massCorrect(spec_tpose[0],filePPMshift) # CHANGE THIS TO DOT...
      spectrum.set_peaks((spec_tpose.transpose(),spec_array[1]))
      spectrum.updateRanges()
   outexp.addSpectrum(spectrum)
   #exp[spectrum_num]=spectrum

outexp.updateRanges()
#exp.swap(outexp)
#exp.updateRanges()
file.store(filename_out, outexp)
