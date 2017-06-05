import os, sys, re
import optparse
import shutil
import pandas
import numpy
import gc
import subprocess

#####################################
#This is a script to combine the output reports from
#Skyline, in preparation for MSstats!  Let's get started.
#
#VERSION 0.70A
version="0.70A"
#DATE: 10/11/2016
date="10/11/2016"
#####################################
print "-----------------------------------------------------------------------"
print "Welcome to the MSstats wrapper for Galaxy, Wohlschlegel Lab UCLA"
print "Written by William Barshop"
print "Version: ",version
print "Date: ",date

basedir=os.getcwd()


####################################
#Argument parsing! So much fun!
#We'll use OptParse even though some
#people really rave about argparse...
#
#
# NB: With Optparse, if an option is
#     not specified, it will take a
#     value of None
####################################

parser = optparse.OptionParser()
parser.add_option("--experiment_file",action="store",type="string",dest="experiment_file")
parser.add_option("--folder",action="store",type="string",dest="operation_folder",default=".")
parser.add_option("--msstats-image-RData",action="store",type="string",dest="image_RData")
parser.add_option("--msstats-comparison-csv",action="store",type="string",dest="comparison_csv")

################# OUTPUTS ################################
parser.add_option("--comparisonPlotOutput",action="store",type="string",dest="comparisonPlotOutput")
parser.add_option("--heatmapOutput",action="store",type="string",dest="heatmapOutput")
parser.add_option("--volcanoPlotOutput",action="store",type="string",dest="volcanoPlotOutput")
parser.add_option("--RScriptOutput",action="store",type="string",dest="RScriptOutput")



################## BELOW THIS ARE PLOTTING OPTIONS ##############################   These are actually all going to be moved into a separate tool
#general options
parser.add_option("--significance",action="store",type="float",dest="significance") # For the volcano plots...
parser.add_option("--FCthreshold",action="store",type="float",dest="FCthreshold") # FC threshold For the volcano plots...
parser.add_option("--ylimUp",action="store",type="float",dest="ylimUp") # ylimUp threshold for the plots
parser.add_option("--ylimDown",action="store",type="float",dest="ylimDown") # ylimDown threshold for plots
parser.add_option("--xlimUp",action="store",type="float",dest="xlimUp") # xlimUp threshold for Volcano plots
parser.add_option("--autoAxes",action="store_true",dest="autoAxes")
parser.add_option("--xAxisSize",action="store",type="int",dest="xAxisSize")
parser.add_option("--yAxisSize",action="store",type="int",dest="yAxisSize")
parser.add_option("--width",action="store",type="int",dest="width",default=10)
parser.add_option("--height",action="store",type="int",dest="height",default=10)

#HeatMap
parser.add_option("--numProtein",action="store",type="int",dest="numProtein",default=180) # Number of proteins per heatmap... Max is 180
parser.add_option("--clustering",action="store",type="string",dest="clustering",default="protein") # clustering type for heatmap... Can be "protein", "comparison", "both"

#VolcanoPlot
parser.add_option("--dotSize",action="store",type="int",dest="dotSize",default=3)#volcanoplot
parser.add_option("--textSize",action="store",type="int",dest="textSize",default=4)#volcanoplot
parser.add_option("--proteinName",action="store_true",dest="proteinName") # On volcano plot, draw protein names?
parser.add_option("--legendSize",action="store",type="int",dest="legendSize",default=7)




(options,args) = parser.parse_args()
if options.autoAxes:
    xlimUp="FALSE"
    ylimUp="FALSE"
    ylimDown="FALSE"
else:
    xlimUp=options.xlimUp
    ylimUp=options.ylimUp
    ylimDown=options.ylimDown

if options.proteinName:
    proteinName="TRUE"
else:
    proteinName="FALSE"

print "Now we're going to prepare the R script for MSstats graphing..."

#Let's start by reading in the experiment structure.
group_information = pandas.read_csv(options.experiment_file,sep='\t')
comparison_df = pandas.read_csv(options.comparison_csv)
with open("MSstats_Script.R",'wb') as script_writer:
    script_writer.write("library(MSstats)\n")
    script_writer.write("setwd(\""+str(basedir)+"\")\n")   #We're going to set the current directory...
    script_writer.write("load(\""+str(options.image_RData)+"\")\n")
    #script_writer.write("comparisonResult<-read.csv(\""+str(options.comparison_csv)+"\")\n")   #We will load in the input CSV file! (In this case by absolute path, though that's not necessary...)

    #script_writer.write("write.csv(comparisonResult$ComparisonResult,file=\"comparisonResult_output.csv\")\n")

    #OKAY! So, now we're going to write out the plots... This may take a bit...
    #So, first, let's check if we can output a heatmap (number of comparisons >2)

    if len(comparison_df['Label'].unique().tolist())>=2:
        #script_writer.write("groupComparisonPlots(data=comparisonResult$ComparisonResult,type=\"Heatmap\", logBase.pvalue=2, sig="+str(options.significance)+", FCcutoff="+str(options.FCthreshold)+",ylimUp="+str(ylimUp)+",ylimDown="+str(ylimDown)+",xlimUp="+str(xlimUp)+",x.axis.size="+str(options.xAxisSize)+",y.axis.size="+str(options.yAxisSize)+",numProtein="+str(options.numProtein)+",clustering=\""+options.clustering+"\",width="+str(options.width)+",height="+str(options.height)+")\n") #add width, height, address
        script_writer.write("groupComparisonPlots(data=comparisonResult$ComparisonResult,type=\"Heatmap\", logBase.pvalue=2,x.axis.size="+str(options.xAxisSize)+",y.axis.size="+str(options.yAxisSize)+",numProtein="+str(options.numProtein)+",clustering=\""+options.clustering+"\",width="+str(options.width)+",height="+str(options.height)+")\n") #add width, height, address
        #pass

    script_writer.write("groupComparisonPlots(data=comparisonResult$ComparisonResult,ProteinName=\""+proteinName+"\",type=\"VolcanoPlot\", logBase.pvalue=2, sig="+str(options.significance)+", FCcutoff="+str(options.FCthreshold)+",ylimUp="+str(ylimUp)+",ylimDown="+str(ylimDown)+",xlimUp="+str(xlimUp)+",x.axis.size="+str(options.xAxisSize)+",dot.size="+str(options.dotSize)+",text.size="+str(options.textSize)+",legend.size="+str(options.legendSize)+",width="+str(options.width)+",height="+str(options.height)+",which.Comparison=\"all\")\n")
    script_writer.write("groupComparisonPlots(data=comparisonResult$ComparisonResult,type=\"ComparisonPlot\", sig="+str(options.significance)+",x.axis.size="+str(options.xAxisSize)+",dot.size="+str(options.dotSize)+",legend.size="+str(options.legendSize)+",width="+str(options.width)+",height="+str(options.height)+",which.Comparison=\"all\")\n")



#OKAY.... The R Script has been written!
#We're going to execute the R script now!
print "Copying RScript back to Galaxy..."
shutil.copy('MSstats_Script.R',options.RScriptOutput)
subprocess.check_call(['Rscript', 'MSstats_Script.R'],shell=False,stderr=sys.stdout.fileno())

print "Moving files to final output locations...."
#print os.listdir(os.getcwd())
#shutil.copy('TMP_dataProcess_output.csv',options.processedOutput)
#shutil.copy('comparisonResult_output.csv',options.comparisonOutput)

shutil.copy('VolcanoPlot.pdf',options.volcanoPlotOutput)
if len(comparison_df['Label'].unique().tolist())>2:
    shutil.copy('Heatmap.pdf',options.heatmapOutput)
shutil.copy('ComparisonPlot.pdf',options.comparisonPlotOutput)


print "All done!"
