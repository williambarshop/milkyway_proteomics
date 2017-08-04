#from tqdm import tqdm
import os, sys
import shutil
import subprocess
import argparse

def DIA_Umpire_parameter_generator(paraName = "DIA_Umpire_params.txt", 
                                   Thread = 32, 
                                   RPmax = 25, 
                                   RFmax = 300, 
                                   CorrThreshold = 0.2, 
                                   DeltaApex = 0.6, 
                                   RTOverlap = 0.3, 
                                   AdjustFragIntensity = "true", 
                                   BoostComplementaryIon = "false",
                                   ExportPrecursorPeak = "true",
                                   SE_MS1PPM = 10,
                                   SE_MS2PPM = 15,
                                   SE_Resolution = 15000,
                                   SE_SN = 1.1,
                                   SE_MS2SN = 1.1,
                                   SE_EstimateBG = "false",
                                   SE_MinMSIntensity = 1,
                                   SE_MinMSMSIntensity = 1,
                                   SE_NoMissedScan = 2,
                                   SE_MaxCurveRTRange = 2,
                                   SE_RemoveGroupedPeaks = "true",
                                   SE_RemoveGroupedPeaksRTOverlap = 0.3,
                                   SE_RemoveGroupedPeaksCorr = 0.3,
                                   SE_MinNoPeakCluster = 2,
                                   SE_MaxNoPeakCluster = 4,
                                   SE_IsoPattern = 0.3,
                                   SE_MassDefectFilter = "true",
                                   SE_MassDefectOffset = 0.1,
                                   SE_StartCharge = 1,
                                   SE_EndCharge = 6,
                                   SE_MS2StartCharge = 1,
                                   SE_MS2EndCharge = 6,
                                   SE_MinFrag = 10,
                                   SE_StartRT = 0,
                                   SE_EndRT = 9999,
                                   SE_MinMZ = 200,
                                   SE_MinPrecursorMass = 600,
                                   SE_MaxPrecursorMass = 6000,
                                   WindowType = "V_SWATH",
                                   WindowSize = 10.0,
                                   WindowSetting = None):
    paraFile = open(paraName, 'w')
    paraFile.write("Thread = " + str(Thread) + "\n")
    paraFile.write("RPmax = " + str(RPmax) + "\n")
    paraFile.write("RFmax = " + str(RFmax) + "\n")
    paraFile.write("CorrThreshold = " + str(CorrThreshold) + "\n")
    paraFile.write("DeltaApex = " + str(DeltaApex) + "\n")
    paraFile.write("RTOverlap = " + str(RTOverlap) + "\n")
    paraFile.write("AdjustFragIntensity = " + AdjustFragIntensity + "\n")
    paraFile.write("BoostComplementaryIon = " + BoostComplementaryIon + "\n")
    paraFile.write("ExportPrecursorPeak = " + ExportPrecursorPeak + "\n")
    paraFile.write("SE.MS1PPM = " + str(SE_MS1PPM) + "\n")
    paraFile.write("SE.MS2PPM = " + str(SE_MS2PPM) + "\n")
    paraFile.write("SE.Resolution = " + str(SE_Resolution) + "\n")
    paraFile.write("SE.SN = " + str(SE_SN) + "\n")
    paraFile.write("SE.MS2SN = " + str(SE_MS2SN) + "\n")
    paraFile.write("SE.EstimateBG = " + SE_EstimateBG + "\n")
    paraFile.write("SE.MinMSIntensity = " + str(SE_MinMSIntensity) + "\n")
    paraFile.write("SE.MinMSMSIntensity = " + str(SE_MinMSMSIntensity) + "\n")
    paraFile.write("SE.NoMissedScan = " + str(SE_NoMissedScan) + "\n")
    paraFile.write("SE.MaxCurveRTRange = " + str(SE_MaxCurveRTRange) + "\n")
    paraFile.write("SE.RemoveGroupedPeaks = " + SE_RemoveGroupedPeaks + "\n")
    paraFile.write("SE.RemoveGroupedPeaksRTOverlap = " + str(SE_RemoveGroupedPeaksRTOverlap) + "\n")
    paraFile.write("SE.RemoveGroupedPeaksCorr = " + str(SE_RemoveGroupedPeaksCorr) + "\n")
    paraFile.write("SE.MinNoPeakCluster = " + str(SE_MinNoPeakCluster) + "\n")
    paraFile.write("SE.MaxNoPeakCluster = " + str(SE_MaxNoPeakCluster) + "\n")
    paraFile.write("SE.IsoPattern = " + str(SE_IsoPattern) + "\n")
    paraFile.write("SE.MassDefectFilter = " + SE_MassDefectFilter + "\n")
    paraFile.write("SE.MassDefectOffset = " + str(SE_MassDefectOffset) + "\n")
    paraFile.write("SE.StartCharge = " + str(SE_StartCharge) + "\n")
    paraFile.write("SE.EndCharge = " + str(SE_EndCharge) + "\n")
    paraFile.write("SE.MS2StartCharge = " + str(SE_MS2StartCharge) + "\n")
    paraFile.write("SE.MS2EndCharge = " + str(SE_MS2EndCharge) + "\n")
    paraFile.write("SE.MinFrag = " + str(SE_MinFrag) + "\n")
    paraFile.write("SE.StartRT = " + str(SE_StartRT) + "\n")
    paraFile.write("SE.EndRT = " + str(SE_EndRT) + "\n")
    paraFile.write("SE.MinMZ = " + str(SE_MinMZ) + "\n")
    paraFile.write("SE.MinPrecursorMass = " + str(SE_MinPrecursorMass) + "\n")
    paraFile.write("SE.MaxPrecursorMass = " + str(SE_MaxPrecursorMass) + "\n")
    paraFile.write("WindowType = " + WindowType + "\n")
    if WindowSetting == "SWATH":
        paraFile.write("WindowSize = " + str(WindowSize) + "\n")
        #pass
    #if WindowSetting == None and WindowType=="V_SWATH":
    #    #DEFAULT 55WINDOW SETTING
    #    #paraFile.write("399.75	410.55\n410.05	420.45\n419.95	429.05\n428.55	436.85\n436.35	443.95\n443.45	450.35\n449.85	456.75\n456.25	463.15\n462.65	469.55\n469.05	475.85\n475.35	482.25\n481.75	488.65\n488.15	495.05\n494.55	502.15\n501.65	508.55\n508.05	515.65\n515.15	522.05\n521.55	528.45\n527.95	535.55\n535.05	541.95\n541.45	549.05\n548.55	556.15\n555.65	563.25\n562.75	570.35\n569.85	578.25\n577.75	586.05\n585.55	593.15\n592.65	601.65\n601.15	608.75\n608.25	616.55\n616.05	625.15\n624.65	633.65\n633.15	642.15\n641.65	651.35\n650.85	662.05\n661.55	671.95\n671.45	682.65\n682.15	694.05\n693.55	705.35\n704.85	718.15\n717.65	733.05\n732.55	748.05\n747.55	764.35\n763.85	782.15\n781.65	801.35\n800.85	823.35\n822.85	846.05\n845.55	869.55\n869.05	896.55\n896.05	926.35\n925.85	961.85\n961.35	1012.35\n1011.85	1091.95\n1091.45	1217.65\n1217.15	1543.75\n")
    #    pass
    if WindowSetting == None:
        print "Leaving window settings blank..."
        paraFile.write("\n\n")
    else:
        print "Writing the window settings from the input provided!"
        paraFile.write("==window setting begin\n")
        with open(WindowSetting,'rb') as windowReader:
            for each_line in windowReader:
                each_line=each_line.replace(",","\t")
                paraFile.write(each_line)
        paraFile.write("==window setting end")

    	#print WindowSetting
    paraFile.close()
    print("DIA-Umpire parameter has been generated!")

########################################
VERSION = "1.0.0b"
UPDATED = "2017/04/15"
########################################
print "================================================="
print "== DIA-Umpire Wrapper"
print "== Written by Hee Jong Kim, William Barshop"
print "== Version: ", VERSION
print "== Last Updated: ", UPDATED
print "================================================="

parser = argparse.ArgumentParser(description="Complete DIA-Umpire Runner!")
# Set the location of parameter text
parser.add_argument('-paraName', action='store', dest='paraName', default="DIA_Umpire_params.txt", help="Set up the name and location of parameter text")
parser.add_argument('-paraOutput', action='store', dest='paraOutput', help="Set up the output location of parameter text")
# Galaxy tool folder, where we can start to find DIA-Umpire.
parser.add_argument('-tooldir',action='store',dest='tooldir',help="Galaxy tool folder")
# mzML Output
parser.add_argument('-mzMLoutput',action='store',dest='mzMLoutput',help="This is a required string to place the output mzML")
# Thread setting
parser.add_argument('-Thread', action='store', dest='Thread', default=32, type=int, help="Set up # of available cores")
# Precursor-fragments grouping parameters
parser.add_argument('-RPmax', action='store', dest='RPmax', default=25, type=int, help="Precursor-fragments grouping parameters")
parser.add_argument('-RFmax', action='store', dest='RFmax', default=300, type=int, help="Precursor-fragments grouping parameters")
parser.add_argument('-CorrThreshold', action='store', dest='CorrThreshold', default=0.2, type=float, help="Precursor-fragments grouping parameters")
parser.add_argument('-DeltaApex', action='store', dest='DeltaApex', default=0.6, type=float, help="Precursor-fragments grouping parameters")
parser.add_argument('-RTOverlap', action='store', dest='RTOverlap', default=0.3, type=float, help="Precursor-fragments grouping parameters")
#Fragment intensity adjustments
parser.add_argument('-AdjustFragIntensity', action='store', dest='AdjustFragIntensity', default="true", help="Fragment intensity adjustments")
parser.add_argument('-BoostComplementaryIon', action='store', dest='BoostComplementaryIon', default="false", help="Fragment intensity adjustments")
#Export detected MS1 features (output feature file can be loaded and mapped to RAW data in BatMass)
parser.add_argument('-ExportPrecursorPeak', action='store', dest='ExportPrecursorPeak', default="true", help="Export detected MS1 features")
#Signal extraction: mass accuracy and resolution
parser.add_argument('-MS1PPM', action='store', dest='MS1PPM', default=10, type=int, help="Signal extraction: mass accuracy and resolution")
parser.add_argument('-MS2PPM', action='store', dest='MS2PPM', default=15, type=int, help="Signal extraction: mass accuracy and resolution")
parser.add_argument('-Resolution', action='store', dest='Resolution', default=15000, type=int, help="Signal extraction: mass accuracy and resolution")
#Signal extraction: signal to noise filter
parser.add_argument('-SN', action='store', dest='SN', default=1.1, type=float, help="Signal extraction: signal to noise filter")
parser.add_argument('-MS2SN', action='store', dest='MS2SN', default=1.1, type=float, help="Signal extraction: signal to noise filter")
#Signal extraction: minimum signal intensity filter
parser.add_argument('-EstimateBG', action='store', dest='EstimateBG', default="false", help="Signal extraction: minimum signal intensity filter")
parser.add_argument('-MinMSIntensity', action='store', dest='MinMSIntensity', default=1, type=int, help="Signal extraction: minimum signal intensity filter")
parser.add_argument('-MinMSMSIntensity', action='store', dest='MinMSMSIntensity', default=1, type=int, help="Signal extraction: minimum signal intensity filter")
#Signal extraction: peak curve detection and isotope grouping
parser.add_argument('-NoMissedScan', action='store', dest='NoMissedScan', default=2, type=int, help="Signal extraction: peak curve detection and isotope grouping")
parser.add_argument('-MaxCurveRTRange', action='store', dest='MaxCurveRTRange', default=2, type=int, help="Signal extraction: peak curve detection and isotope grouping")
parser.add_argument('-RemoveGroupedPeaks', action='store', dest='RemoveGroupedPeaks', default="true", help="Signal extraction: peak curve detection and isotope grouping")
parser.add_argument('-RemoveGroupedPeaksRTOverlap', action='store', dest='RemoveGroupedPeaksRTOverlap', default=0.3, type=float, help="Signal extraction: peak curve detection and isotope grouping")
parser.add_argument('-RemoveGroupedPeaksCorr', action='store', dest='RemoveGroupedPeaksCorr', default=0.3, type=float, help="Signal extraction: peak curve detection and isotope grouping")
parser.add_argument('-MinNoPeakCluster', action='store', dest='MinNoPeakCluster', default=2, type=int, help="Signal extraction: peak curve detection and isotope grouping")
parser.add_argument('-MaxNoPeakCluster', action='store', dest='MaxNoPeakCluster', default=4, type=int, help="Signal extraction: peak curve detection and isotope grouping")
#Signal extraction: filtering of MS1 features
parser.add_argument('-IsoPattern', action='store', dest='IsoPattern', default=0.3, type=float, help="Signal extraction: filtering of MS1 features")
parser.add_argument('-MassDefectFilter', action='store', dest='MassDefectFilter', default="true", help="Signal extraction: filtering of MS1 features")
parser.add_argument('-MassDefectOffset', action='store', dest='MassDefectOffset', default=0.1, type=float, help="Signal extraction: filtering of MS1 features")
#Signal extraction: other 
parser.add_argument('-StartCharge', action='store', dest='StartCharge', default=1, type=int, help="Signal extraction: other")
parser.add_argument('-EndCharge', action='store', dest='EndCharge', default=6, type=int, help="Signal extraction: other")
parser.add_argument('-MS2StartCharge', action='store', dest='MS2StartCharge', default=1, type=int, help="Signal extraction: other")
parser.add_argument('-MS2EndCharge', action='store', dest='MS2EndCharge', default=6, type=int, help="Signal extraction: other")
parser.add_argument('-MinFrag', action='store', dest='MinFrag', default=10, type=int, help="Signal extraction: other")
parser.add_argument('-StartRT', action='store', dest='StartRT', default=0, type=int, help="Signal extraction: other")
parser.add_argument('-EndRT', action='store', dest='EndRT', default=9999, type=int, help="Signal extraction: other")
parser.add_argument('-MinMZ', action='store', dest='MinMZ', default=200, type=int, help="Signal extraction: other")
parser.add_argument('-MinPrecursorMass', action='store', dest='MinPrecursorMass', default=600, type=int, help="Signal extraction: other")
parser.add_argument('-MaxPrecursorMass', action='store', dest='MaxPrecursorMass', default=60000, type=int, help="Signal extraction: other")
#Isolation window setting
parser.add_argument('-WindowType', action='store', dest='WindowType', default="V_SWATH", help="Isolation window setting")
parser.add_argument('-WindowSize', action='store', dest='WindowSize', default=10.0, type=float, help="Isolation window setting")
parser.add_argument('-WindowSetting', action='store', dest='WindowSetting', default=None, help="Window setting")
# DIA-Umpire running associated setting
parser.add_argument('-inputFolder', action='store', dest='inputFolder', default=None, help="Set up the path of mzML file(s) location")


# DIA-Umpire parameter generator execution based on given parameters from commandline
parsed_args = parser.parse_args()
DIA_Umpire_parameter_generator(paraName = parsed_args.paraName, 
                                   Thread = parsed_args.Thread, 
                                   RPmax = parsed_args.RPmax, 
                                   RFmax = parsed_args.RFmax, 
                                   CorrThreshold = parsed_args.CorrThreshold, 
                                   DeltaApex = parsed_args.DeltaApex, 
                                   RTOverlap = parsed_args.RTOverlap, 
                                   AdjustFragIntensity = parsed_args.AdjustFragIntensity, 
                                   BoostComplementaryIon = parsed_args.BoostComplementaryIon,
                                   ExportPrecursorPeak = parsed_args.ExportPrecursorPeak,
                                   SE_MS1PPM = parsed_args.MS1PPM,
                                   SE_MS2PPM = parsed_args.MS2PPM,
                                   SE_Resolution = parsed_args.Resolution,
                                   SE_SN = parsed_args.SN,
                                   SE_MS2SN = parsed_args.MS2SN,
                                   SE_EstimateBG = parsed_args.EstimateBG,
                                   SE_MinMSIntensity = parsed_args.MinMSIntensity,
                                   SE_MinMSMSIntensity = parsed_args.MinMSMSIntensity,
                                   SE_NoMissedScan = parsed_args.NoMissedScan,
                                   SE_MaxCurveRTRange = parsed_args.MaxCurveRTRange,
                                   SE_RemoveGroupedPeaks = parsed_args.RemoveGroupedPeaks,
                                   SE_RemoveGroupedPeaksRTOverlap = parsed_args.RemoveGroupedPeaksRTOverlap,
                                   SE_RemoveGroupedPeaksCorr = parsed_args.RemoveGroupedPeaksCorr,
                                   SE_MinNoPeakCluster = parsed_args.MinNoPeakCluster,
                                   SE_MaxNoPeakCluster = parsed_args.MaxNoPeakCluster,
                                   SE_IsoPattern = parsed_args.IsoPattern,
                                   SE_MassDefectFilter = parsed_args.MassDefectFilter,
                                   SE_MassDefectOffset = parsed_args.MassDefectOffset,
                                   SE_StartCharge = parsed_args.StartCharge,
                                   SE_EndCharge = parsed_args.EndCharge,
                                   SE_MS2StartCharge = parsed_args.MS2StartCharge,
                                   SE_MS2EndCharge = parsed_args.MS2EndCharge,
                                   SE_MinFrag = parsed_args.MinFrag,
                                   SE_StartRT = parsed_args.StartRT,
                                   SE_EndRT = parsed_args.EndRT,
                                   SE_MinMZ = parsed_args.MinMZ,
                                   SE_MinPrecursorMass = parsed_args.MinPrecursorMass,
                                   SE_MaxPrecursorMass = parsed_args.MaxPrecursorMass,
                                   WindowType = parsed_args.WindowType,
                                   WindowSize = parsed_args.WindowSize,
                                   WindowSetting = parsed_args.WindowSetting)

if parsed_args.paraOutput is not None:
    shutil.copy(parsed_args.paraName,parsed_args.paraOutput)

# Actual execution part starts here
overwrite=True
DIAUmpire = parsed_args.tooldir+"/v2.1.2/DIA_Umpire_SE.jar"
parent_folder = parsed_args.inputFolder
print [x for x in os.walk(parent_folder)]
input_folder_list = [parent_folder]#list(set(folder for folder, subfolders, files in os.walk(parent_folder) for file_ in files if os.path.splitext(file_)[1] in ['.mzML']))#['.raw','.mzXML','.mzML']))
output_folder_name = "DIA_Umpire_outputs"
output_path = os.path.join(os.getcwd(),output_folder_name)
if not os.path.exists(output_path):
    os.mkdir(output_path)
print "We'll work in these folders:",input_folder_list
current_dir=os.getcwd()
# Iteration of each file: file conversion and DIA-Umpire execution
print "Working on folders: ",input_folder_list
for input_folder in input_folder_list: #tqdm(input_folder_list):
    raw_input=[current_dir+"/"+input_folder+f for f in os.listdir(input_folder) if f.endswith('.mzML')]
    #filtered_list = [x.rsplit(".",1)[0] for x in raw_input]
    #print filtered_list
    for each_file in raw_input:
        os.chdir(output_path)
        shutil.copy(each_file,output_path)
        #os.chdir(each_file.rsplit("/",1)[0])
        print os.getcwd(),"working dir of file..."
        if len([output_path+f for f in os.listdir(output_path) if f == '{0}.mzXML'.format(each_file.rsplit("/",1)[1].rsplit(".",1)[0])])<1:
            #if each_file.endswith('.raw'):
            #    subprocess.call("msconvert --mzXML --filter \"peakPicking true 1-\" {0}.raw".format(each_file.rsplit("/",1)[1].rsplit(".",1)[0]),shell=True)
            #else:
            subprocess.call("msconvert --mzXML --filter \"peakPicking true 1-\" {0}.mzML".format(output_path+"/"+each_file.rsplit("/",1)[1].rsplit(".",1)[0]),shell=True,stderr=sys.stdout.fileno())
            #shutil.move(each_file,each_file.split('.')[0]+"_original.mzML")
        print "About to run...\n","java -jar -Xmx32G {0} {1}.mzXML {2}".format(DIAUmpire,output_path+"/"+each_file.rsplit(".",1)[0].rsplit("/",1)[1],parsed_args.paraName)
        shutil.copy(current_dir+"/"+parsed_args.paraName,os.getcwd()+"/")
        if (overwrite and len([output_path+f for f in os.listdir(output_path) if '{0}_Q'.format(each_file.rsplit(".",1)[0]) and f.endswith('.mgf') and f!=f.rsplit("_Q",1)[0]+".mgf" and each_file.rsplit(".",1)[0].rsplit("/",1)[1] in f]) !=0) or len([output_path+"/"+f for f in os.listdir(output_path) if '{0}_Q'.format(each_file.rsplit(".",1)[0]) and f.endswith('.mgf') and f!=f.rsplit("_Q",1)[0]+".mgf" and each_file.rsplit(".",1)[0].rsplit("/",1)[1] in f]) == 0:
            print "current directory",os.getcwd()
            subprocess.call("java -jar -Xmx32G \"{0}\" \"{1}.mzXML\" {2}".format(DIAUmpire,output_path+"/"+each_file.rsplit(".",1)[0].rsplit("/",1)[1],parsed_args.paraName),shell=True)
        
        os.chdir(current_dir)

        mgf_files = [output_path+"/"+f for f in os.listdir(output_path) if '{0}_Q'.format(each_file.rsplit(".",1)[0]) and f.endswith('.mgf') and f!=f.rsplit("_Q",1)[0]+".mgf" and each_file.rsplit(".",1)[0].rsplit("/",1)[1] in f]
        print "Taking care of ...",mgf_files
        if len(mgf_files)>0:
            output=mgf_files[0].rsplit("_Q",1)[0]+".mgf"
        else:
            print "... There are no files here for {0}... I can't do anything.".format(each_file)
            continue
        writer=open(output,'wb')
        scan_ctr=1

        output_without_extension=output.split(".")[0]

        for each_file in mgf_files:
            with open(each_file,'r') as filereader:
                contents=filereader.readlines()
            for eachline in contents:
                if "TITLE=" in eachline:
                    writer.write("SCANS={0}\n".format(str(scan_ctr)))
                    charge=eachline.rsplit(".",1)[1]
                    writer.write("TITLE="+output_without_extension+"."+str(scan_ctr)+"."+str(scan_ctr)+"."+charge+"\n")
                    scan_ctr+=1
                else:
                    writer.write(eachline)
            print "Finished processing "+each_file+" ...."
        writer.close()    
        #cleanup Q1/Q2/Q3 files...
        for each_file in mgf_files:
            if "_Q1.mgf" in each_file or "_Q2.mgf" in each_file or "_Q3.mgf" in each_file:
                os.remove(each_file)
                print "removed ",each_file
        #current_dir=os.getcwd()
        #os.chdir(input_folder)
        #os.mkdir("mzML_output")
        #os.chdir("mzML_output")
        subprocess.call("msconvert --mzML --filter \"peakPicking true 1-\" {0} --outdir {1} --outfile {2}".format(output_path +"/"+output.rsplit("/",1)[1].rsplit(".",1)[0]+".mgf",parsed_args.mzMLoutput.rsplit("/",1)[0],os.path.basename(parsed_args.mzMLoutput)),shell=True,stderr=sys.stdout.fileno())
        os.rename(parsed_args.mzMLoutput+".mzML",parsed_args.mzMLoutput)
        #os.chdir(current_dir)
        print "Finsihed making mzML output"
        #shutil.copy(output_path +"/"+output.rsplit("/",1)[1].rsplit(".",1)[0]+".mzML",parsed_args.mzMLoutput)

