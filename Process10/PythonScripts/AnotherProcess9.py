import optparse
import logging as log
import os, sys
import ConfigParser
from collections import defaultdict
import simplejson, json
import urllib, urllib2
import time
import shutil
import csv
import re

SampleIDs = []
AllSamplenamesToGet = []
AllSampleSLXIDsToGet = []
AllSampleLocations = []
AllSampleBamNames = []
AllSampleFQNames = []
AllSampleBamLocations = []
AllSampleFQLocations = []
AllSampleRunsToGet = []
AllSampleLanesToGet = []
AnalysisState = []

ProjectNumber = 0
ProjectSampleNumber = 0
SLXIDSampleNumber = 0
SampleNameSampleNumber = 0
RetrievedBamLocations = 0
RetrievedFQLocations = 0

LocationsDir = sys.argv[1]
ProjectsFile = sys.argv[2]
SLXFile = sys.argv[3]
SampleNamesFile = sys.argv[4]

LocationsDir = "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/TestForTest/Locations"
ProjectsFile = "ProjectsFromLims.txt"
SLXFile = "SLXIDsFromLims.txt"
SampleNamesFile = "SamplenamesFromLims.txt"


for agr in sys.argv:
	print(agr)

#LocationsDir = '/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20120626_RobinsonJ_JC_ARchipFoxOE_ProperInput/Locations/'


from suds.client import Client
log.getLogger('suds').setLevel(log.INFO)
lims = Client("http://uk-cri-ldmz02.crnet.org/solexa-ws/SolexaExportBeanWS?wsdl")
log.debug(lims)


try:
   with open(os.path.join(LocationsDir,ProjectsFile)) as f:
   	for line in f:
   		ChompLine = line.rstrip("\n")
   		#print ChompLine
   		try:
   			ProjectSamples = lims.service.getSamplesByProjectName(ChompLine)
   			ProjectNumber = ProjectNumber + 1
   			for Sample in ProjectSamples.item:
   				SampleIDs.append(Sample.userId)
   				ProjectSampleNumber = ProjectSampleNumber + 1
   		except:
   			print "No Information for "+ChompLine
except IOError as e:
   print 'Error parsing or retrieving information for sample Project list'


try:
   with open(os.path.join(LocationsDir,SLXFile)) as f:
   	for line in f:
   		ChompLine = line.rstrip("\n")
   		#print ChompLine
   		SampleIDs.append(ChompLine)
   		SLXIDSampleNumber = SLXIDSampleNumber + 1
except IOError as e:
   print 'Error parsing or retrieving information for sample SLXids list'


try:
   with open(os.path.join(LocationsDir,SampleNamesFile)) as f:
   	for line in f:
   		ChompLine = line.rstrip("\n")
   		#print ChompLine
   		SampleIDs.append(ChompLine)
   		SampleNameSampleNumber = SampleNameSampleNumber + 1
except IOError as e:
   print 'Error parsing or retrieving information for sample samplename list'


for ASampleID in SampleIDs:
	try:
		sampleProcessIDs = lims.service.getAllSampleLanesForId(ASampleID)[0]
		for sampleProcessID in sampleProcessIDs:
			AllSamplenamesToGet.append(sampleProcessID.genomicsSampleId)
			AllSampleSLXIDsToGet.append(sampleProcessID.userSampleId)
			AllSampleRunsToGet.append(sampleProcessID.lane)
			AllSampleLanesToGet.append(sampleProcessID.solexaRun.runNumber)
			TempBam = lims.service.getFileLocations(sampleProcessID.sampleProcess_id,"HTTP","BAM")
			print ASampleID
			if TempBam == "":
				TempBam = lims.service.getFileLocations(sampleProcessID.sampleProcess_id,"FILE","BAM")
			if TempBam != "":
				AllSampleBamLocations.append(TempBam[0][0].host+":"+TempBam[0][0].path)
				AllSampleBamNames.append(TempBam[0][0].filename)
			if TempBam == "":
				AllSampleBamLocations.append("")
				AllSampleBamNames.append("")
			if TempBam != "":			
				RetrievedBamLocations = RetrievedBamLocations +1
			#AllSampleBamLocations.append("")
			#AllSampleBamNames.append("")
			TempFQ = lims.service.getFileLocations(sampleProcessID.sampleProcess_id,"HTTP","FASTQ")
			if TempFQ == "":
				TempFQ = lims.service.getFileLocations(sampleProcessID.sampleProcess_id,"FILE","FASTQ")			
			if TempFQ != "":
				AllSampleFQLocations.append(TempFQ[0][0].host+":"+TempFQ[0][0].path)
				AllSampleFQNames.append(TempFQ[0][0].filename)
			if TempFQ == "":
				AllSampleFQLocations.append("")
				AllSampleFQNames.append("")
			if TempFQ != "":			
				RetrievedFQLocations = RetrievedFQLocations +1	
			if sampleProcessID.solexaRun.multiplexed:
				AnalysisState.append("Parent_Of_Sample")
				limsSLXID = re.sub("SLX-","",sampleProcessID.genomicsSampleId)
				subSamples = lims.service.getSampleBySlxId(limsSLXID,"true")[0]
				for subSample in subSamples:
					AllSamplenamesToGet.append(sampleProcessID.genomicsSampleId+"."+subSample.fileSystemFriendlyUserId)
					AllSampleSLXIDsToGet.append(subSample.fileSystemFriendlyUserId)
					AllSampleRunsToGet.append(sampleProcessID.lane)
					AllSampleLanesToGet.append(sampleProcessID.solexaRun.runNumber)	
					subsSampleID = subSample.id
					TempSubFQ = lims.service.getFileLocationsForSubsample(sampleProcessID.sampleProcess_id,subSample.id)[0]
					AnalysisState.append("RunMe")
					if TempSubFQ != "":
						for SubLocation in TempSubFQ:
							if (SubLocation.scheme == "HTTP") &  (len(re.findall("fq",SubLocation.filename)) > 0):
								AllSampleFQLocations.append(SubLocation.host+":"+SubLocation.path)
								AllSampleFQNames.append(SubLocation.filename)
								RetrievedFQLocations = RetrievedFQLocations +1	
							if (SubLocation.scheme == "HTTP") &  (len(re.findall("bam",SubLocation.filename)) > 0):
								AllSampleBamLocations.append(SubLocation.host+":"+SubLocation.path)
								AllSampleBamNames.append(SubLocation.filename)
								RetrievedBamLocations = RetrievedBamLocations +1										
					if TempSubFQ == "":
						AllSampleBamLocations.append("")
						AllSampleBamNames.append("")							
						AllSampleFQLocations.append("")
						AllSampleFQNames.append("")
					print subSample.fileSystemFriendlyUserId
					#AllSampleLocations
			else:
				AnalysisState.append("RunMe") 	
	except:
		print "No Locations for "+ASampleID



with open(os.path.join(LocationsDir,'Lims_SampleLocations.txt'), 'wb') as csvfile:
	for i in range(len(AllSamplenamesToGet)):
		csvfile.write(AllSamplenamesToGet[i]+"\t"+AllSampleSLXIDsToGet[i]+"\t"+str(AllSampleRunsToGet[i])+"\t"+str(AllSampleLanesToGet[i])+"\t"+AllSampleBamLocations[i]+"\t"+AllSampleBamNames[i]+"\t"+AllSampleFQLocations[i]+"\t"+AllSampleFQNames[i]+"\t"+AnalysisState[i]+"\n")


print "Found "+str(ProjectNumber)+" Projects"
print "Found "+str(ProjectSampleNumber)+" samples associated with projects" 
print "Found "+str(SLXIDSampleNumber)+" samples by SLXID"
print "Found "+str(SampleNameSampleNumber)+" samples by samplename"
print "Found "+str(RetrievedBamLocations)+" BAM file locations" 
print "Found "+str(RetrievedFQLocations)+" FQ file locations" 
