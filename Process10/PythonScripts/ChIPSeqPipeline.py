#!/bin/env /home/mib-cri/software/PythonInstall/bin/python2.7
import argparse
import textwrap
import os
import ConfigParser
import sys
import subprocess

parser = argparse.ArgumentParser(
	prog='ChIPSeqPipeline',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='\nChIP-Seq Analysis Pipeline\nWritten by Tom Carroll,Suraj Menon and Rory Stark\nCRUK, CRI\n2012\n')
	 

group3 = parser.add_argument_group('Analysis Settings')
group = parser.add_argument_group('SLX and Project Management')
group1 = parser.add_argument_group('Coverage And Pileup')
group2 = parser.add_argument_group('Peak Calling')
ConfigArgs = parser.add_argument_group('Additional Config Arguments')
#argparse.RawDescriptionHelpFormatter


group3.add_argument("--genome",nargs=1,choices=["HG18","GRCh37","MM9","MM8"],dest="genome",default=None,help='Genome to use for analysis')
                    
group3.add_argument('--excludedRegions', nargs=1,type=file,dest="excludedregions",default=None,metavar="",
                    help='Bedfile of genomic regions to exclude from analysis')

group3.add_argument('--useExcludedRegions',action='store_true', default=True,dest="useexcludedregions",
                    help='Remove reads mapping to blacklist regions')

                   
group3.add_argument('--mapqFilter', nargs=1,dest="mapqfilter",metavar="",type=int,default=15,
                    help='MapQ quality filter setting')
                    
group3.add_argument('--removeDuplicates', action='store_true', default=False,
                    dest='removeduplicates',
                    help='Remove duplicates from analysis')

group1.add_argument('--bothStrands', action='store_true', default=False,
                    dest='BothStrands',
                    help='Generate seperate BedGraph files for either strand')

group.add_argument('--notProjectDirectory', action='store_false', default=True,
                    dest='usepresentdir',
                    help='Use directory basename for project assignment')

group.add_argument('--workingDirectory',nargs=1,metavar="",
		    #default=os.getcwd(),
                    dest='workingdirectory',
		    default=None,
                    help='Directory for project')


group.add_argument('--tempDirectory',nargs=1,metavar="",
		    #default=os.getcwd(),
                    dest='tempdirectory',
		    default=None,
                    help='Temporary Directory')
                    

group2.add_argument('--callMacsPeaks', nargs=1,choices=['Yes','No'],
                    dest='callmacspeaks',default="Yes",
                    help='Call MACS peaks')

group2.add_argument('--callSicerPeaks', nargs=1,choices=['Yes','No'],
                    dest='callsicerpeaks',default="No",
                    help='Call Sicer peaks')

group2.add_argument('--callTpicsPeaks', nargs=1,choices=['Yes','No'],
                    dest='calltpicspeaks',default="No",
                    help='Call T-PICS peaks')

group2.add_argument('--callMacsMotifs', nargs=1,choices=['Yes','No'],
                    dest='callmacsmotifs',default="No",
                    help='Call T-PICS peaks')

group.add_argument('--bamDirectory',nargs=1,metavar="",dest='bamdirectory',
#                    default=os.path.join(os.getcwd(),"bamFiles"),
		    default=None,
                    help='Directory for bam files (default: %(default)s)',
                    )
                    
group.add_argument('--fastqDirectory',nargs=1,metavar="",dest='fastqdirectory',
                    #default=os.path.join(os.getcwd(),"FQFiles"),
                    default=None,
                    help='Directory for fastq files',
                    )
                    
group.add_argument('--addSLXIDs',nargs="*",action='append', dest='SLXids',metavar="SLXID",
                    default=[],
                    help='SLXID/s to be added to the current project',
                    )

group.add_argument('--addProjects',nargs="*",action='append', dest='Projects',metavar="ProjectID",
                    default=[],
                    help='Project/s to be merged with the current project',
                    )

group.add_argument('--addMetadata',nargs="*",action='append', dest='metadata',metavar="SampleSheet.csv",
                    default=[],
                    help='SampleSheets containing metadata to be added to the current project',
                    )

parser.add_argument('--version', action='version', version='%(prog)s 0.1')
ConfigArgs.add_argument('Config Variables',nargs=argparse.REMAINDER,help="Overwrite or include additional variables into config.ini")

results = parser.parse_args()
CmdLineOptions = vars(results)
AllKeys = CmdLineOptions.keys()
if str(CmdLineOptions["tempdirectory"]) == "None":
	if str(CmdLineOptions["workingdirectory"]) == "None":
		Temptemp = os.path.join(os.getcwd(),"Temp")
	if str(CmdLineOptions["workingdirectory"]) != "None":
		Temptemp = os.path.join(ConfigOptions["workingdirectory"],"Temp")
if str(CmdLineOptions["tempdirectory"]) != "None":
	Temptemp = os.path.join(os.getcwd(),"Temp")
	Temptemp = os.path.join(ConfigOptions["workingdirectory"],"Temp")

config = ConfigParser.ConfigParser()
if os.path.exists(os.path.join(Temptemp,"config.ini")):
	config.read(os.path.join(Temptemp,"config.ini"))
	print "\nLocal config file found"
else:
	config.read("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/Config/Config.ini")
	print "\nUsing generic config\n"

ConfigOptions = {}
for section in config.sections():
    for option in config.options(section):
	ConfigOptions[option] = config.get(section, option)




for Key in AllKeys:
	if Key in ConfigOptions:
		#print Key+"\t"+ConfigOptions[Key]
		if str(ConfigOptions[Key]) != str(CmdLineOptions[Key]) and CmdLineOptions[Key] is not None:
			print "Overwriting config option for "+Key+" to "+str(CmdLineOptions[Key][0])+"\n"
			if Key != "genome":
				ConfigOptions[Key] = str(CmdLineOptions[Key])
			if Key == "genome":
				ConfigOptions[Key] = str(CmdLineOptions[Key][0])
			if Key == "callmacspeaks":
				ConfigOptions[Key] = str(CmdLineOptions[Key][0])
			if Key == "callsicerpeaks":
				ConfigOptions[Key] = str(CmdLineOptions[Key][0])
			if Key == "calltpicspeaks":
				ConfigOptions[Key] = str(CmdLineOptions[Key][0])
			if Key == "callmacsmotifs":
				ConfigOptions[Key] = str(CmdLineOptions[Key][0])

			#ConfigOptions[Key] = CmdLineOptions[Key]
			#print str(ConfigOptions[Key][0])+"\n"

if str(ConfigOptions["genome"]) == "None" and str(CmdLineOptions["genome"]) == "None":
	print "No Genome set in config or as commandline argument\nplease see usage with {ChipSeqPipeline --help}\n"
	sys.exit()

if str(ConfigOptions["workingdirectory"]) == "None" and str(CmdLineOptions["workingdirectory"]) == "None":
	print "No working directory set in config or as commandline argument\nworking directory as "+os.getcwd()+"\n"
	ConfigOptions["workingdirectory"] = os.getcwd()

if str(ConfigOptions["bamdirectory"]) == "None" and str(CmdLineOptions["bamdirectory"]) == "None":
	print "No Bam directory set in config or as commandline argument\nSetting Bam directory as "+os.path.join(ConfigOptions["workingdirectory"],"bamFiles\n")
	ConfigOptions["bamdirectory"] = os.path.join(ConfigOptions["workingdirectory"],"bamFiles")

if str(ConfigOptions["fastqdirectory"]) == "None" and str(CmdLineOptions["fastqdirectory"]) == "None":
	print "No FastQ directory set in config or as commandline argument\nSetting working directory as "+os.path.join(ConfigOptions["workingdirectory"],"FQFiles\n")
	ConfigOptions["fastqdirectory"] = os.path.join(ConfigOptions["workingdirectory"],"FQFiles")

if str(ConfigOptions["tempdirectory"]) == "None" and str(CmdLineOptions["tempdirectory"]) == "None":
	print "No Temp directory set in config or as commandline argument\nSetting temp directory as "+os.path.join(ConfigOptions["workingdirectory"],"Temp\n")
	ConfigOptions["tempdirectory"] = os.path.join(ConfigOptions["workingdirectory"],"Temp")

if not os.path.exists(ConfigOptions["tempdirectory"]):
    os.makedirs(ConfigOptions["tempdirectory"])

ExtraSLXids = []
if CmdLineOptions["metadata"]:
	metadata = CmdLineOptions["metadata"][0]
	metaFile = open(os.path.join(ConfigOptions["tempdirectory"],"metadata.txt"),"w")
	for meta in metadata:
		metaFile.write(str(meta)+"\n")
	metaFile.close()

ExtraMeta = []
if CmdLineOptions["SLXids"]:
	ExtraSLXids = CmdLineOptions["SLXids"][0]
	SLXFile = open(os.path.join(ConfigOptions["tempdirectory"],"Samples_SLXIDs.txt"),"w")
	for SLXid in ExtraSLXids:
		SLXFile.write(str(SLXid)+"\n")
	SLXFile.close()


ExtraProjects = []
if CmdLineOptions["Projects"] or CmdLineOptions["usepresentdir"]:
	ProjectFile = open(os.path.join(ConfigOptions["tempdirectory"],"Projects.txt"),"w")
	if CmdLineOptions["Projects"]:
		TempProjects=CmdLineOptions["Projects"][0]
		for aProj in TempProjects:
			ExtraProjects.append(aProj)
	if CmdLineOptions["usepresentdir"]:
		ExtraProjects.append(os.path.basename(os.getcwd()))
	for project in ExtraProjects:
		ProjectFile.write(str(project)+"\n")
	ProjectFile.close()

if not ExtraProjects and not ExtraSLXids:
	print "No Samples or Project specified!!! Can't do much"
	sys.exit()




subprocess.call(["bash", "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/BashScripts/GetLimsInfo.sh",ConfigOptions["tempdirectory"]])




if not os.path.exists(ConfigOptions["bamdirectory"]):
	os.makedirs(ConfigOptions["bamdirectory"])
	subprocess.call("lfs setstripe "+ConfigOptions["bamdirectory"],shell=True)
	


inifile = open(os.path.join(ConfigOptions["tempdirectory"],"config.ini"),'w')

OutConfig = ConfigParser.ConfigParser()
# add the settings to the structure of the file, and lets write it out...
OutConfig.add_section('Analysis Settings')
OutConfig.add_section('Peak Calling')
OutConfig.add_section('SLX and Project Management')
OutConfig.add_section('Executables')
OutConfig.add_section('Custom Scripts')
OutConfig.add_section('ExcludedRegions')
OutConfig.add_section('Genomes')

OutConfig.set('Analysis Settings','genome',str(ConfigOptions["genome"]))
OutConfig.set('Analysis Settings','excludedRegions',str(ConfigOptions["excludedregions"]))
OutConfig.set('Analysis Settings','mapQFilter',str(ConfigOptions["mapqfilter"]))
OutConfig.set('Analysis Settings','useExcludedRegionFilter',str(ConfigOptions["useexcludedregionfilter"]))
OutConfig.set('Analysis Settings','removeDuplicates',str(ConfigOptions["removeduplicates"]))


OutConfig.set('Peak Calling','callmacspeaks',str(ConfigOptions["callmacspeaks"]))
OutConfig.set('Peak Calling','callsicerpeaks',str(ConfigOptions["callsicerpeaks"]))
OutConfig.set('Peak Calling','calltpicspeaks',str(ConfigOptions["calltpicspeaks"]))
OutConfig.set('Peak Calling','callmacsmotifs',str(ConfigOptions["callmacsmotifs"]))


OutConfig.set('SLX and Project Management','workingdirectory',str(ConfigOptions["workingdirectory"]))
OutConfig.set('SLX and Project Management','bamdirectory',str(ConfigOptions["bamdirectory"]))
OutConfig.set('SLX and Project Management','fastqdirectory',str(ConfigOptions["fastqdirectory"]))
OutConfig.set('SLX and Project Management','tempdirectory',str(ConfigOptions["tempdirectory"]))

OutConfig.set('Executables','bwa',str(ConfigOptions["bwa"]))
OutConfig.set('Executables','python',str(ConfigOptions["python"]))
OutConfig.set('Executables','samtools',str(ConfigOptions["samtools"]))
OutConfig.set('Executables','picard',str(ConfigOptions["picard"]))
OutConfig.set('Executables','perl',str(ConfigOptions["perl"]))
OutConfig.set('Executables','rsync',str(ConfigOptions["rsync"]))
OutConfig.set('Executables','bedtools',str(ConfigOptions["bedtools"]))
OutConfig.set('Executables','java',str(ConfigOptions["java"]))

OutConfig.set('Custom Scripts','bam_processing_script',str(ConfigOptions["bam_processing_script"]))
OutConfig.set('Custom Scripts','metadata_script',str(ConfigOptions["metadata_script"]))
OutConfig.set('Custom Scripts','getgenome_script',str(ConfigOptions["getgenome_script"]))
OutConfig.set('Custom Scripts','bamlocations_script',str(ConfigOptions["bamlocations_script"]))
OutConfig.set('Custom Scripts','fastqlocations_script',str(ConfigOptions["fastqlocations_script"]))
OutConfig.set('Custom Scripts','sicer_cri_script',str(ConfigOptions["sicer_cri_script"]))
OutConfig.set('Custom Scripts','tpicszeta_cri_script',str(ConfigOptions["tpicszeta_cri_script"]))

OutConfig.set('ExcludedRegions','HG18',str(ConfigOptions["hg18"]))
OutConfig.set('ExcludedRegions','GRCh37',str(ConfigOptions["grch37"]))

OutConfig.set('Genomes','HG18',str(ConfigOptions["hg18"]))
OutConfig.set('Genomes','GRCh37',str(ConfigOptions["grch37"]))

OutConfig.write(inifile)
inifile.close()



subprocess.call(["/home/mib-cri/software/R-2.14.0/bin/Rscript","--vanilla","/lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/RMainPipeSetUp.r",str(ConfigOptions["tempdirectory"])])

#subprocess.call("bash /lustre/mib-cri/carrol09/Work/MyPipe/Process10/BashScripts/GetLimsInfoOld.sh",shell=True)


#subprocess.call(["bash", "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/BashScripts/GetLimsInfoOld.sh"])







#subprocess.call(["bash",ConfigOptions["lims_info_script"]])
#subprocess.call(["bash",str(ConfigOptions["lims_info_script"]),str(ConfigOptions["tempdirectory"])])
#subprocess.call(["bash", "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/BashScripts/GetLimsInfoOld.sh"])
#subprocess.call(["bash",ConfigOptions["lims_info_script"],ConfigOptions["tempdirectory"]])

#print(ConfigOptions["tempdirectory"])
#print(ConfigOptions["lims_info_script"])





