###ChIP-seq HTML report

## Get R libraries

library(ggplot2)
library(hwriter)
library(googleVis)
library(reshape2)
library(Hmisc)
library(XML)
library(tractor.base)
library(raster)
library(scales)
library(gridSVG)
source("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/HTML_functionsV2.r")

##
ReadLengthActual <- 36


## Get arguments from the command line
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
Args <- commandArgs(trailingOnly = TRUE)
WkgDir <- Args[1]
WkgDir <- getwd()



source("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/Workflow_Functions3.r")

#ConfigFile <- readIniFile(file.path(WkgDir,"Config","Config.ini"))
#genome <- ConfigFile[ConfigFile[,2] %in% "genome",3]

PipeLineLocations <- GetImportantLocations(WkgDir,"Config")
genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config") 


system(paste("mkdir -p ",file.path(WkgDir,"HTML_Report","Plots")))
system(paste("mkdir -p ",file.path(WkgDir,"HTML_Report","IGV")))


## Read in SampleSheet
SampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),stringsAsFactors=F,sep=",",h=T)
NewGroup <- cbind(SampleSheet[,c("GenomicsID","SampleName")],paste(SampleSheet[,"Tissue"],SampleSheet[,"Factor"],sep="_"))
colnames(NewGroup) <- c("GenomicsID","SampleName","Tissue_And_Factor")



### Open HTML page for writing
p <- openPage(file.path(WkgDir,"HTML_Report","SampleSummary4.html"))
hwrite("ChIP-Seq HTML Report", heading=1,p,br=T)

AddDiffDivCode(p)


MakeIGVSampleMetadata(SampleSheet,file.path(WkgDir,"HTML_Report","IGV"))

####Make Images

##On/off Peak image
#Macs
tryCatch(
MakeOnOffPeakPics("Macs")
,error=function(e) ErrorForImageMaking(p),
finally=print("Macs On/Off peak Images done"))
#TPICs
tryCatch(
MakeOnOffPeakPics("TPICs")
,error=function(e) ErrorForImageMaking(p),
finally=print("TPICs On/Off peak Images done"))
#Sicer
tryCatch(
MakeOnOffPeakPics("Sicer")
,error=function(e) ErrorForImageMaking(p),
finally=print("Sicer On/Off peak Images done"))


tryCatch(
GetPeakProfilePlot(SampleSheet,"Macs")
,error=function(e) ErrorForImageMaking(p),
finally=print("Macs CoverageProfilePlots done"))

tryCatch(
GetPeakProfilePlot(SampleSheet,"Sicer")
,error=function(e) ErrorForImageMaking(p),
finally=print("Sicer CoverageProfilePlots done"))

tryCatch(
GetPeakProfilePlot(SampleSheet,"TPICs")
,error=function(e) ErrorForImageMaking(p),
finally=print("TPICs CoverageProfilePlots done"))


tryCatch(
GetGCInPeaks(SampleSheet,"Macs")
,error=function(e) ErrorForImageMaking(p),
finally=print("GC content in peaks done"))

tryCatch(
GetGCInPeaks(SampleSheet,"Sicer")
,error=function(e) ErrorForImageMaking(p),
finally=print("GC content in peaks done"))

tryCatch(
GetGCInPeaks(SampleSheet,"TPICs")
,error=function(e) ErrorForImageMaking(p),
finally=print("GC content in peaks done"))


tryCatch(
GetGenomicCov(SampleSheet)
,error=function(e) ErrorForImageMaking(p),
finally=print("On/Off peak Images done"))


tryCatch(
GetCoverageProfilePlot(SampleSheet)
,error=function(e) ErrorForImageMaking(p),
finally=print("CoverageProfilePlots done"))


tryCatch(
MakeCrossCor(SampleSheet)
,error=function(e) ErrorForImageMaking(p),
finally=print("Coverage shift plot done"))

tryCatch(
GetGiniSDPlots(SampleSheet)
,error=function(e) ErrorForImageMaking(p),
finally=print("Coverage shift plot done"))

tryCatch(
MakeFragmentLengthsPlots(SampleSheet)
,error=function(e) ErrorForImageMaking(p),
finally=print("Coverage shift plot done"))

#tryCatch(
#PlotPosAndNegInPeaks(SampleSheet,p)
#,error=function(e) ErrorForImageMaking(p),
#finally=print("Pos/Neg peak Images done")
#)

################################################################
#######   Generic Setup ######################################
##############################################################


hwrite("The Report",p,style='font-weight: bold',br=T)
hwrite("This document gives an overview of analysis performed as part of the CRI ChIP-seq analysis pipeline.",p,br=T)
hwrite("The organisation of this report follows the order of ChIP-seq analysis from generic and ChIP-seq specific QC through to the identification of genomic locations enriched for ChIP-seq signal, association with 
genomic features (i.e.genes,promoters) and any underlying enriched sequence motifs.",br=T,p)
hwrite("",p,br=T)

q <- openPage(file.path(WkgDir,"HTML_Report","Directory_Structure.html"))
OrganisationTable <- matrix(nrow=30,ncol=2)
OrganisationTable[1,1] <- basename(getwd())
OrganisationTable[1,2] <- "The project name"

OrganisationTable[2,1] <- "--HTML_Report"
OrganisationTable[2,2] <- "HTML report directory"

OrganisationTable[3,1] <- "&nbsp;&nbsp;|-SampleSummary.html" 
OrganisationTable[3,2] <- "Main HTML report (you are here!)" 

OrganisationTable[4,1] <- "&nbsp;&nbsp;|-IGV.html" 
OrganisationTable[4,2] <- "IGV session file directory"

OrganisationTable[5,1] <- "--bamFiles" 
OrganisationTable[5,2] <- "bamFile directory"

OrganisationTable[6,1] <- "&nbsp;&nbsp;|-*_Processed.bam" 
OrganisationTable[6,2] <- "Processed BamFiles for each sample"

OrganisationTable[7,1] <- "--Coverage" 
OrganisationTable[7,2] <- "Coverage and Pile-up directory"

OrganisationTable[8,1] <- "&nbsp;&nbsp;|-samplename_Processed.bw" 
OrganisationTable[8,2] <- "Processed BigWig (Coverage file) for each sample"

OrganisationTable[9,1] <- "&nbsp;&nbsp;|-samplename_Processed.bedGraph" 
OrganisationTable[9,2] <- "Processed BedGraph (Coverage file) for each sample" 

OrganisationTable[10,1] <- "--Peaks" 
OrganisationTable[10,2] <- "Main peak directory"

OrganisationTable[11,1] <- "&nbsp;&nbsp;|-Macs_Peaks" 
OrganisationTable[11,2] <- "MACS' peaks directory"

OrganisationTable[12,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;|-samplename_Processed_peaks.xls" 
OrganisationTable[12,2] <- "MACS output peak file per sample" 

OrganisationTable[13,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;|-samplename_Processed_peaks_Annotated.xls" 
OrganisationTable[13,2] <- "MACS output peak file  per sample annotated with overlapping and nearest genes"

OrganisationTable[14,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;|-samplename_Processed_GO_Results" 
OrganisationTable[14,2] <- "Gene Ontology terms enriched for peaks within each sample"

OrganisationTable[15,1] <- "&nbsp;&nbsp;|-Sicer_Peaks" 
OrganisationTable[15,2] <- "Sicer's peaks directory"

OrganisationTable[16,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;|-samplename" 
OrganisationTable[16,2] <- "Directory per sample containing Sicer results" 

OrganisationTable[17,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-samplename_Processed-islands-summary-FDR.01" 
OrganisationTable[17,2] <- "Sicer output peak file per sample" 

OrganisationTable[18,1] <- "&nbsp;&nbsp;|-TPICS_Peaks" 
OrganisationTable[18,2] <- "TPICS peaks directory"

OrganisationTable[19,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;|-samplename" 
OrganisationTable[19,2] <- "Directory per sample containing TPICS results" 

OrganisationTable[20,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-samplename_TPICS_Peaks.bed" 
OrganisationTable[20,2] <- "TPICS output peak file per sample" 

OrganisationTable[21,1] <- "--Motifs" 
OrganisationTable[21,2] <- "Main Motif directory"

OrganisationTable[22,1] <- "&nbsp;&nbsp;|--samplename" 
OrganisationTable[22,2] <- "Per sample motif directory"

OrganisationTable[23,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;|-Denovo" 
OrganisationTable[23,2] <- "Directory containin Denovo motif results per sample" 

OrganisationTable[24,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-index.html" 
OrganisationTable[24,2] <- "HTML file containing links to denvo motif finding results per sample" 

OrganisationTable[25,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;|-Known" 
OrganisationTable[25,2] <- "Directory containin Known motif results per sample" 

OrganisationTable[26,1] <- "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|-ame.html" 
OrganisationTable[26,2] <- "HTML file containing links to known motif finding results per sample" 

OrganisationTable[27,1] <- "&nbsp;&nbsp;|-Summary_Significance_Table_meme_Denovo_Motifs.txt" 
OrganisationTable[27,2] <- "Summary table of Meme Denovo motifs discovered in all sample" 

OrganisationTable[28,1] <- "&nbsp;&nbsp;|-Summary_Significance_Table_dreme_Denovo_Motifs.txt" 
OrganisationTable[28,2] <- "Summary table of Dreme Denovo motifs discovered in all sample" 

OrganisationTable[29,1] <- "&nbsp;&nbsp;|Summary_Significance_Table_KnownMotifs.txt" 
OrganisationTable[29,2] <- "Summary table of Ame Denovo motifs discovered in all sample" 




hwrite(OrganisationTable[1:29,1:2], q, border=0,style='font-family:monospace')
closePage(q)
hwrite("Directories and Files",p,style='font-weight: bold',br=T)
hwrite("This report contains only an overview of your ChIP-seq results and so much of the information can be found in accompanying files. References to the locations of accompanying files are provided both within the text and tables
describing file contents and a breakdown of relevant files and directory structure seen ",p) 
hwrite('here', p, link="Directory_Structure.html",br=T)
hwrite("",p,br=T)

hwrite("Visualisation of Sequence Files",p,style='font-weight: bold',br=T)
hwrite("Visual inspection of your dataset as a coverage Bedgraph or BigWig using tools such IGV, IGB, Deliance or UCSC allows for a rapid evaluation of known targets of ChIP enrichment and remains a standard for evaluation of ChIP enriched regions.",p,br=T) 

hwrite("This report contains links allowing for your data to be autoloaded in IGV but these require that an IGV session of the relevant genome be running. To use a webstart link to an IGV session for your genome
please use the below link",p,br=T)

#AddSessionLink(file.path(WkgDir,"HTML_Report","IGV"),file.path(WkgDir,"HTML_Report","TrialSampleSummary.html"),genome,p)
if(tolower(genome) == "hg18"){
hwrite(hmakeTag('a','Start New IGV Session for HG18',href='http://www.broadinstitute.org/igv/projects/current/igv.php?genome=hg18'),p,br=T)
}
if(tolower(genome) == "grch37"){
hwrite(hmakeTag('a','Start New IGV Session for GRCh37',href='http://www.broadinstitute.org/igv/projects/current/igv.php?genome=hg19'),p,br=T)
}
hwrite('',br=T)

hwrite("Once a blank session has been started, files may be loaded by simply selected a link for a new session (which will contain only the file selected) or by selected a link to the current session (which will load alongside other files in session)",p,br=T)
hwrite("For a detailed description of IGV genome browser and usage please visit their website ",p) 
hwrite('here', p, link='http://www.broadinstitute.org/igv/')





hwrite("Index", heading=2,p,br=T)
hwrite("Section 1 - Post Alignment QC", heading=3,p,style='font-weight: bold',br=T,link='#s1')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Read Filtering",p,br=T,style='font-weight: italics',link='#Read_Filtering')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Table 1 - Quality Control and Summary of Read Distributions of Samples Within Project",p,br=T,style='font-style: italic',link='#T1')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Read In Genes",p,br=T,style='font-weight: italics',link='#Read_In_Features')
hwrite("Section 2 - Coverage statistics", heading=3,style='font-weight: bold',p,br=T,link='#s2')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Coverage QC Metrics",p,br=T,link='#Coverage_QC')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Table 2 - Assessment of Enrichment Efficiency and Links to Genome Coverage Graphs",p,br=T,style='font-style: italic',link='#T2')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Genomic Depth of Coverage",p,br=T,link='#Depth_QC')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Cross-Coverage and Predicted Fragment Length",p,br=T,link='#FragmentLength_QC')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Inequality of Coverage",p,br=T,link='#Gini_QC')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Average Gene/TSS Coverage",p,br=T,link='#TSS_QC')
hwrite("Section 3 - Peak Calling", heading=3,style='font-weight: bold',p,br=T,link='#s3')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Correspondence Between Peak-Callers",p,br=T,link='#AcossPeakCallers_QC')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Table 3 - Percentage of overlapping peaks and the Jaccard Index score for pairwise comparisons between peak callers",p,br=T,style='font-style:italic',link='#T3')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Macs Peaks",p,br=T,link='#Macs_Peaks')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Calling Peaks With MACS",p,br=T,link='#Macs_Calling')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Annotation Of Peaks And Functional Enrichment",p,br=T,link='#AnnoAndGO_Calling')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Table 4 - Number of Peaks Calleds and Genes With Peaks using MACS peak calling algorithm",p,br=T,style='font-style:italic',link='#T4')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Reads and Coverage In Peaks",p,br=T,link='#OnAndOff_Calling')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("GC Content of Peaks",p,br=T,link='#GC_Calling')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Peaks in Genes",p,br=T,link='#PeaksGenes_Calling')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Correlation Between Sample Peak Sets",p,br=T,link='#GenomMetric_Calling')
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
hwrite("Table 5 - Pair-wise sample Jaccard Measures and links to GenoMetriCorr results",p,br=T,style='font-style:italic',link='#T4')
hwrite("Differential Binding Analysis (DiffBind)",p,br=T,link='#DiffBind_Calling')
#hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",p)
#hwrite("Quantitative Analysis of ChIP signal",p,br=T,link='#Diffbind_Calling')
####
######
###############################################################

################################################################
#######   Read QC Section ######################################
##############################################################


##Get Reads stats table
hwrite(hwrite("Section 1 - Post Alignment QC",name='s1'),p,heading=2)

hwrite("This first section discusses sequence read QC metrics and the distributions of filtered sequence reads within gene features. These metrics allow for a rapid identification of library of sequencing problems as well as identification of any sample group specific genomic feature enrichment",br=T,p)
hwrite(" ",br=T,p)
hwrite("Read Filtering",style='font-weight: bold',br=T,p,name='Read_Filtering')
hwrite("Prior to analysis of coverage and peak calling, ChIP-seq data is filtered to remove low quality reads and artifacts. Such artifacts and noise related signal may confound many downstream tools and the use of the following filters can 
help clean data to allow for a more successfull ChIP-seq analysis.",br=T,p)
hwrite("Inorder to remove noise related reads and sequencing artifacts, the aligned files are filtered successively by:",br=T,p)
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Removing/Counting reads mapped to main contigs/chromosomes (delRand) since downstream analysis does not include random or unassembled contigs.",br=T,p)
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Removing/Counting reads mapped to mapping outside Duke Excluded Regions (Excluded),regions previously identified as enriched across multiple unrelated ChIP-seq experiments.",br=T,p)
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Removing/Counting reads with low mapping quality ,MAPQ < 15, (Filtered) and hence less reliable in downstream analysis.",br=T,p)
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Marking/Counting reads which are mapped to the exact genomic locations of another mapped read (Duplicated),since such reads are often due to PCR amplification/bias.",br=T,p)
hwrite(" ",br=T,p)
hwrite("A Summary of the numbers of reads remaining after each successive filtering step can be seen in Table-1 and Figure 1 alongside the calculated duplication rate (Percentage of reads at final filtering step marked as Duplicated).",br=T,p)
hwrite(" ",br=T,p)


hwrite(hwrite("Table 1 - Quality Control and Summary of Read Distributions of Samples Within Project",name="T1"),heading=3,p,center=TRUE)

OptionsChart1 <- "{
height:500,
width:1000,
hAxis: {title: 'Number of Reads'},
chartArea:{left:100,top:0,width:\"50%\",height:\"80%\"},
isStacked:true
}"

OptionsChart2 <- "{
height:500,
width:1000,
hAxis: {title: 'Number of Reads'},
chartArea:{left:100,top:0,width:\"50%\",height:\"80%\"},
isStacked:false
}"

OptionsChart3 <- "{
height:500,
width:1000,
vAxis: {title: 'Number of Reads'},
isStacked:true
}"

chartOptions <- list(OptionsChart3,OptionsChart1,OptionsChart2)

Tommy <- ""

tryCatch(

SSRSMain <- MakeSSRS(SampleSheet,genome,chartOptions)

,error=function(e) ErrorForHTML(p),

finally=print("Read Stats HTML done"))

GoogleVisMain(SSRSMain,p)
addVisElement("table",p)
hwrite("Table shows sample information and the breakdown of read numbers and duplication rate after successive read quality filtering steps.",br=T,p)
hwrite("",br=T,p)
addVisElement("QCChart",p)
hwrite("Figure 1 - Barplot of reads remaining after successive filtering steps",br=T,p)
hwrite("",br=T,p)
#hwrite("Reads in Genes",style='font-weight: bold',br=T,p,name='Read_In_Features')
Divs1 <- GetDiffDivsAndSpan("Reads in Genes",'Read_In_Features')
cat(file=p,paste(Divs1[[1]][1],Divs1[[2]][1],Divs1[[3]][1],sep=""))
hwrite("",br=T,p)
hwrite("As well as assessment of the quality of reads within a dataset, it is often useful to examine the numbers of reads mapping to genomic locations of interest.
For many epigenetic marks the percentage of reads and hence signal at the promoter, TSS and along the gene body will reveal the ChIP success and
comparisons of genomic distribution of samples relating to the same ChIP antibody should share a highly similar distribution. Figure 2 and 3 show the total or percentage of reads which
map to genomic locations of gene upstream regions, TSS (+/- 500), gene body and intergenic regions
",br=T,p)
addVisElement("BarplotOfReads",p)
hwrite("Figure 2- Barplot of total reads and percentage of total reads within genomic features",br=T,p)

addVisElement("BarplotOfAbsolute",p)
hwrite("Figure 3 - Barplot of percentage of total reads within genomic features",br=T,p)

hwrite("",br=T,p)
cat(file=p,paste(Divs1[[3]][2],"\n",sep=""))
hwrite("",br=T,p)


################################################################
#######   Coverage QC Section ######################################
##############################################################

###############################################################
##Get Coverage table
hwrite(hwrite("Coverage statistics",name='s2'),heading=2,br=T,p)
hwrite("",br=T,p)
hwrite("This section investigates the distribution of reads along the genome and so provides information on the quality and effiency of ChIP for the samples under investigation
as well as links to files to visualise the read distributions along the genome.",br=T,p)
hwrite(" ",br=T,p)
hwrite("Coverage QC Metrics",style='font-weight: bold',br=T,p,name='Coverage_QC')
hwrite("For transcription factor ChIP experiments a more specific signal is expected along the genome whereas no IP input controls would have a diffuse non-specific signal through-out the genome. By the use of appropriate coverage related metrics
and visualisation in genome browsers the relative success of ChIP experiments can often be established ",p)
hwrite("and so two sets of metrics usefull to coverage QC are introduced within this section-",br=T,p)
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Measures of the uniformity of Coverage (SD, Gini and adjusted Gini scores)",br=T,p)
hwrite("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Measures of the correlation between coverage on positive or negative strands (Normalised Strand Cross-correlation (NSC) and Relative Strand Cross-correlation(RSC))",br=T,p)
hwrite(" ",br=T,p)
hwrite("Also introduced in this section are genomewide coverage graphs allowing for the visualation of read depth (pileup) at points along the genome. 
",br=T,p)

hwrite("Summary of NSC, RSC,SD and adjusted Gini scores as well as links to the relevant Bedgraph and BigWig coverage files are shown in Table 2.",br=T,p)

hwrite(hwrite("Table 2 - Assessment of Enrichment Efficiency and Links to Genome Coverage Graphs",name="T2"),heading=3,p,center=TRUE)

tryCatch(

MakeSSCT(SampleSheet,genome,p)

,error=function(e) ErrorForHTML(p),

finally=print("Coverage Stats HTML done"))

hwrite("Table shows sample information relating to coverage inequality within samples (SD and Gini scores) and links to BigWig coverage files",br=T,p)
hwrite("",br=T,p)
hwrite("Genomic Depth of Coverage",style='font-weight: bold',br=T,p,name='Depth_QC')
hwrite("A simple assessment of read coverage or pile-up for each sample is to plot the number of base pairs from the genome at a given read depth. 
A sample with a higher quality ChIP enrichment will have a greater proportion of the genome at higher read depths due to specific enrichment where as an input will show a far
smaller overall proportion of the genome at such depths. In figure 4 the log 2 base pairs of genome at differing read depths (including duplicated reads) are plotted and so 
a low number of base pairs covered at a high depth will often be observed for all samples due to some PCR bias.
",br=T,p)

hwrite("",br=T,p)



tryCatch(

addImage("CoveragePlot.png",file.path(WkgDir,"HTML_Report","Plots","CoveragePlot.png"),"","",p)

,error=function(e) ErrorForHTML(p),

finally=print("Coverage plot Image image inserted"))
hwrite("",br=T,p)
hwrite("Figure 4 - Plot of base pairs at given read depths.",br=T,p)
hwrite("",br=T,p)
 Divs1 <- GetDiffDivsAndSpan("Interactive Depth Plot",'interactive_Depth')
 cat(file=p,paste(Divs1[[1]][1],Divs1[[2]][1],Divs1[[3]][1],sep=""))


tryCatch(

#hwrite(paste(hmakeTag("object", data="CoveragePlot2.html",type="image/svg+xml",width="505px",height="505px"),"\n",sep=""),p)
InsertSVG(file.path(WkgDir,"HTML_Report","Plots","CoveragePlot2.html"),p)
,error=function(e) ErrorForHTML(p),

finally=print("SVG Coverage plot Image image inserted"))
cat(file=p,paste(Divs1[[3]][2],"\n",sep=""))
hwrite("",br=T,p)



hwrite("Cross-Correlation/Coverage and Predicted Fragment Length",style='font-weight: bold',br=T,p,name='FragmentLength_QC')
hwrite("The use of coverage profiles on the positive and negative strands to predict fragment length is an integral part of many peaks callers 
and provides a good QC measure for the quality of a ChIP-seq fragment.",br=T,p)
hwrite("Since ChIP-sequencing often involves short reads repesenting the end of a ChIPed fragment, reads from the positive and negative strand will be shifted around a true peak
or enriched region by half the fragment length. For a good ChIP-seq sample therefore this relationship between the reads on the negative and positive strand should be evident
and fragment length be close to that expected for this experiment. For ChIP samples, the distance at which reads from the positive or negative strands are most correlated provides
the measure of cross correlation. The ratio between maximum correlation and minimum observed correlation within 400 bp windows yields the NSC statisitc whereas the RSC statistic is
calculated as the maximum correlation within a 400bp window over the correlation observed at the read length distance. Cross-coverage similarly looks at the relationship between reads on either
strand by identifying the shift of the positive reads to the negative reads which results in the minimum fraction of the genome being covered by at least 1 read.  
",br=T,p)
hwrite("Figure 5 shows the resulting predicted fragment lengths by both cross-correlation and cross-coverage and figure 6 effects on shifting the genomic positions of positive and negative reads towards each other while measuring the proportion of genome covered at each
shift compared to proportion of genome covered wih no shifts. The shift at which the minimum proportion of the genome is covered represents the predicted fragment length.
.",br=T,p)

tryCatch(

addImage("FragmentLengthPlot.png",file.path(WkgDir,"HTML_Report","Plots","FragmentLengthPlot.png"),"","",p)

,error=function(e) ErrorForHTML(p),

finally=print("Coverage shift plot Image image inserted"))
 Divs1 <- GetDiffDivsAndSpan("Cross-Correlation plot",'CrossCorPlot')
 cat(file=p,paste(Divs1[[1]][1],Divs1[[2]][1],Divs1[[3]][1],sep=""))
 
tryCatch(

addImage("CoverageShiftPlot.png",file.path(WkgDir,"HTML_Report","Plots","CoverageShiftPlot.png"),"","",p)

,error=function(e) ErrorForHTML(p),

finally=print("Coverage shift plot Image image inserted"))


hwrite("Figure 6 - Plot of proportion of cross-correlation following shifting of reads on positive strand towards negative",br=T,p)
 cat(file=p,paste(Divs1[[3]][2],"\n",sep=""))

hwrite("",br=T,p)

hwrite("Inequality of Coverage",style='font-weight: bold',br=T,p,name='Gini_QC')
hwrite("A high quality TF or Histone ChIP-seq experiment will have genomic islands of high read depths and long spans of genomic regions devoid of reads. In contrast a good
input control will have few islands of high read depth and consist of small fluctuations in read depth along wider areas. The use of SD and Gini scores allow for the quantitative 
assessment of uneveneness in the read depth along a genome and so allow for a QC measure of how specifically enriched is a sample or how biased is an input control.
",br=T,p)
hwrite("Figures 7-9 show the SD, Gini and adjusted Gini scores respectively for samples and any input controls. These plots can be used to assess the specificity/bias in your sample/input
and hence give a measure of the specific relative enrichment of Sample-Input combinations",br=T,p)

Divs1 <- GetDiffDivsAndSpan("SD of Coverage between samples and Inputs",'SD')
cat(file=p,paste(Divs1[[1]][1],Divs1[[2]][1],Divs1[[3]][1],sep=""))

tryCatch(

addImage("SDsPlot.png",file.path(WkgDir,"HTML_Report","Plots","SDsPlot.png"),"","",p)

,error=function(e) ErrorForHTML(p),

finally=print("On/Off image inserted"))
hwrite("Figure 7 - Plots of coverage inequality as assessed by SD for samples",p)
 cat(file=p,paste(Divs1[[3]][2],"\n",sep=""))
 
Divs1 <- GetDiffDivsAndSpan("Gini of Coverage between samples and Inputs",'Gini')
cat(file=p,paste(Divs1[[1]][1],Divs1[[2]][1],Divs1[[3]][1],sep=""))


tryCatch(

addImage("GinisPlot.png",file.path(WkgDir,"HTML_Report","Plots","GinisPlot.png"),"","",p)

,error=function(e) ErrorForHTML(p),

finally=print("On/Off image inserted"))

hwrite("Figure 8 - Plots of coverage inequality as assessed by Gini scores for samples",p) 
cat(file=p,paste(Divs1[[3]][2],"\n",sep=""))

Divs1 <- GetDiffDivsAndSpan("Adj Gini of Coverage between samples and Inputs",'AdjGini')
cat(file=p,paste(Divs1[[1]][1],Divs1[[2]][1],Divs1[[3]][1],sep=""))
tryCatch(

addImage("AdjustedGinisPlot.png",file.path(WkgDir,"HTML_Report","Plots","AdjustedGinisPlot.png"),"","",p)

,error=function(e) ErrorForHTML(p),

finally=print("On/Off image inserted"))
 cat(file=p,paste(Divs1[[3]][2],"\n",sep=""))
hwrite("Figure 9 - Plots of coverage inequality as assessed by Adjusted Gini scores for samples",p) 
hwrite("(red)",p,style='font-family: monospace;color: #ff2233')
hwrite("and any corresponding inputs",p)
hwrite("(blue)",style='font-family: monospace;color:  #003F87',br=T,p)
hwrite("",br=T,p)
#hwrite("Average Gene/TSS Coverage",style='font-weight: bold',br=T,p,name='TSS_QC')
Divs1 <- GetDiffDivsAndSpan("Average Gene/TSS Coverage",'TSS_QC')
cat(file=p,paste(Divs1[[1]][1],Divs1[[2]][1],Divs1[[3]][1],sep=""))
hwrite("Further to the assessment of the number reads within certain genomic features seen in Section 1, many epigenetic marks are often characterised by their average distibution across the promoter and TSS.
Input samples conversely should show little enrichment across promoter,TSS or gene body and so any increase in average read depth across these regions may represent a feature specific bias. 
Figure 10 shows the average read coverage across all gene promoters and TSSs scaled to the total read length and from here any sample specific enrichment or input bias can be seen
",br=T,p)


tryCatch(

addImage("AverageTSSPlot.png",file.path(WkgDir,"HTML_Report","Plots","AverageTSSPlot.png"),"","",p)

,error=function(e) ErrorForHTML(p),

finally=print("Ave TSS image inserted"))
hwrite("Figure 10 - Average coverage across Genes scaled to total reads in sample",p) 

Divs2 <- GetDiffDivsAndSpan("Interactive TSS plot",'TSS_QCInter')
cat(file=p,paste(Divs2[[1]][1],Divs2[[2]][1],Divs2[[3]][1],sep=""))
tryCatch(

#hwrite(paste(hmakeTag("object", data="CoveragePlot2.html",type="image/svg+xml",width="505px",height="505px"),"\n",sep=""),p)
InsertSVG(file.path(WkgDir,"HTML_Report","Plots","AverageTSSPlot.html"),p)
,error=function(e) ErrorForHTML(p),

finally=print("SVG Coverage plot Image image inserted"))
cat(file=p,paste(Divs2[[3]][2],"\n",sep=""))
hwrite("",br=T,p)

cat(file=p,paste(Divs1[[3]][2],"\n",sep=""))


###############################################################
hwrite(hwrite("Peak Calling", name='s3'),heading=2,br=T,p)
hwrite("",p) 
hwrite("In this section, genomic regions of ChIP enrichment are identifed and profiled for association with other genomic features and their underlying genomic properties.",br=T,p)
hwrite("Identification of theses regions (Peaks) significantly enriched for the epigenetic mark under evaluations can be performed using many differing approaches and algorithms.",br=T,p)
hwrite("As part of the analysis pipeline, several peak callers are used to identify enriched ChIP regions. These differing peak callers allow for different classes of epigenetic marks to be identified as well as act as a further QC measure
. Samples showing little agreement between peak callers often contain higher degrees of non-specific signal or noise and so often lead to low quality peak calls",br=T,p)  

hwrite("",br=T,p)  
hwrite("Correspondence Between Peak-Callers",style='font-weight: bold',br=T,p,name='AcossPeakCallers_QC')
hwrite("In Table 3, the pair-wise agreement between peak calls is assessed. For each sample and pair-wise comparison the percentage of overlapping peaks and the Jaccard Index score for peak sets is shown. Samples showing a 
lower percentage or Jaccard score when compared to other samples within the group is often a sign of an outlying or poor sample.",br=T,p)


hwrite(hwrite("Table 3 - Percentage of overlapping peaks and the Jaccard Index score for pairwise comparisons between peak callers",name="T3"),heading=3,p,center=TRUE)

tryCatch(


MakeSSAcrossPeaks(SampleSheet,p)

,error=function(e) ErrorForHTML(p),

finally=print("Across Peaks HTML done"))

 
###############################################################
##Get Peaks tables (Macs)
hwrite("",br=T,p)  
hwrite("Macs Peaks", heading=2,style='font-weight: bold',p,name="Macs_Peaks",br=T)

hwrite("A well-known and popular peak caller is ",p)
hwrite("MACs",link="http://liulab.dfci.harvard.edu/MACS/index.html",p)
hwrite(" (Model-based Analysis of ChIP-Seq)",p,br=T)
hwrite("Calling Peaks With MACS",style='font-weight: bold',br=T,p,name='Macs_Calling')
hwrite("MACS's approach is to adjust read positions by the estimated fragment length (see section?) and to scan windows along the genome for enrichment over input as well as local area using a poisson distribution.
The output from MACs is the genomic locations of enriched regions, the position of maximum enrichment over input within enriched regions and associated score and p-value of this enrichment. For a more detailed
explanation of MACS algorithm and output please see their website ",p)
hwrite("here",link="http://liulab.dfci.harvard.edu/MACS/index.html",p)
hwrite("and their paper ",p)
hwrite("here",link="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pubmed&cmd=search&term=18798982[pmid]",p,br=T)
hwrite("Annotation Of Peaks And Functional Enrichment",style='font-weight: bold',br=T,p,name='AnnoAndGO_Calling')
hwrite("Following MACS peak calling, these peaks are associate with genes from the relevant genome and a Gene Ontology functional test is performed on every sample's peak calling results.",br=T,p)
hwrite("To annotate peaks to genes, all genes are first extended by 2000bp upstream from their TSS depending on their strand/orientation. Peaks are then first annotated to the extended genes they fall within. This can lead 
to a peak being associated with more than 1 gene if the extended gene boundaries cross eachother. Following this peaks which were not found to be within extended genes are then annotated to the most proximal gene irrespective of 
strand or distance. This results in every peak being associated with at least 1 gene. The number of peaks landing within 2000bp extended genes and the number of extended genes found to contain a peak are seen in Table4",br=T,p)
hwrite("Once peaks have been associated with gene, tests for an enrichment of peaks within genes in a specific functional term can be performed. Due to the differing lengths of genes and associated gene length bias of soe fucntional
terms and appropriate test accounting for differences in gene length is required. Here we apply a precision test, measuring the frequency of peaks in a functional group's genes while accounting for the proportion of genome 
related to this group. Such an approach is used in the functional enrichment tool ",p)
hwrite("GREAT",link="http://great.stanford.edu/public/html/splash.php",p,br=T)

MakePeaksHTMLPart(SampleSheet,p,"Macs")
MakePeaksHTMLPart(SampleSheet,p,"TPICs")
MakePeaksHTMLPart(SampleSheet,p,"Sicer")



#hwrite("Diffbind, written by Rory Stark and Gordon Brown, uses both co-occurence of peaks and pearson correlation of signal across samples to assess similarity between samples and provide a 
#metric of biological reproducibity within your experiment.",br=T,p)


#hwrite("",br=T,p)
#hwrite("Figure 20 and 21 show the results from a Diffbind occupancy and affinity analysis. Peaks were first merged across all samples and those peaks occuring in at least two samples are used for further analysis as shown in figure 19.
# A binary score for the occurence of peaks within samples is used to produce custering plot seen in figure 20 whereas normalised counts within peaks across samples used to produce heatmap in figure 21.",br=T,p)
#hwrite("",br=T,p)
#hwrite("Diffential binding analysis is performed by diffbind to quantitatively identify peaks whose read counts or enrichment is significantly different between conditions. All available pair-wise contrasts are 
# of tissure,factor and/or conditions are used for analysis and the results cna be seen in table 6",br=T,p)
#hwrite("",br=T,p)

#FDRCutOff <- 0.1
#PvalueCutOff <- 0.01

#if(file.exists(file.path(WkgDir,"DiffBind","Contrasts_macs.csv"))){

#   OccupancyHeatMap <-  file.path(WkgDir,"DiffBind","Occupancy_CorHeatmap_macs.png")
#   OverLapRate <- file.path(WkgDir,"DiffBind","Overlap_Rate_macs.png")
#   AffinityHeatMap <- file.path(WkgDir,"DiffBind","Affinity_CorHeatmap_macs.png")
#   MacsContrast <- read.delim(file.path(WkgDir,"DiffBind","Contrasts_macs.csv"),header=T,sep=",")
   
#   DiffMacsIGVLinks <-   vector("character",length=nrow(MacsContrast))
#   DiffMacsIGVLinksFresh <-  vector("character",length=nrow(MacsContrast))

   
#   AllContrasts <- dir(file.path(WkgDir,"DiffBind"),pattern="*.Full_Report_edgeR_macs.csv",full.names=T)
#   ContrastResults <- vector("list",length=nrow(MacsContrast))
#   Scores <- vector("list",length=nrow(MacsContrast))
#   MacsDiffMat <- matrix(nrow=nrow(MacsContrast),ncol=ncol(MacsContrast)+6)
#   for(i in 1:length(AllContrasts)){
#      ScoresTemp <- vector("numeric",length=5)
#     MacsContrastTemp <- read.delim(AllContrasts[[i]],header=T,sep=",")
#     ScoresTemp[1] <- sum(which(MacsContrastTemp$FDR < FDRCutOff))
 #    ScoresTemp[2] <- sum(which(MacsContrastTemp$FDR < FDRCutOff & MacsContrastTemp$Fold < 0))
#     ScoresTemp[3] <- sum(which(MacsContrastTemp$FDR < FDRCutOff & MacsContrastTemp$Fold > 0))
#     ScoresTemp[4] <- sum(which(MacsContrastTemp$p.value < PvalueCutOff))
#     ScoresTemp[5] <- sum(which(MacsContrastTemp$p.value < PvalueCutOff  & MacsContrastTemp$Fold < 0))
#     ScoresTemp[6] <- sum(which(MacsContrastTemp$p.value < PvalueCutOff  & MacsContrastTemp$Fold > 0))
#     ForIGV <- cbind(MacsContrastTemp[,c(1:3)],paste("Peak_",seq(1,nrow(MacsContrastTemp)),sep=""),MacsContrastTemp[,c(7,8)])
#     Links <- paste("=HYPERLINK(\"http://localhost:60151/goto?locus=",paste(MacsContrastTemp[,1],paste(MacsContrastTemp[,2],MacsContrastTemp[,3],sep="-"),sep=":"),"\"),",sep="")
#     WithLinks <- cbind(MacsContrastTemp,Links)
#     colnames(WithLinks)[ncol(WithLinks)] <- "Link"
     
#     DiffMacsIGVLinks[i] <- makeIGVLink(file.path(WkgDir,"DiffBind",paste("Contrast",i,".bed",sep="")),"MacsDiff",paste("Contrast",i,"_Macs",sep=""),file.path(WkgDir,"HTML_Report","IGV"),file.path(WkgDir,"HTML_Report","HTMLReport.html"),paste("MacsDiff",paste("Contrast",i,sep=""),sep=""),genomeName=genome,locusName="All")
#     DiffMacsIGVLinksFresh[i] <- makeIGVLinkFresh(file.path(WkgDir,"DiffBind",paste("Contrast",i,".bed",sep="")),"MacsDiff",paste("Contrast",i,"_Macs",sep=""),file.path(WkgDir,"HTML_Report","IGV"),file.path(WkgDir,"HTML_Report","HTMLReport.html"),paste("MacsDiff",paste("Contrast",i,sep=""),sep=""),genomeName=genome,locusName="All")

#     write.table(ForIGV,file.path(WkgDir,"DiffBind",paste("Contrast",i,".bed",sep="")),sep="\t",quote=F,row.names=F,col.names=F)
#     write.table(WithLinks,file.path(WkgDir,"DiffBind",paste("Contrast",i,"_WithLinks.csv",sep="")),sep=",",quote=F,row.names=F,col.names=F)
#     if(i == 1){
#      MacsDiffMat  <- cbind(MacsContrast[i,],matrix(ScoresTemp,nrow=1,ncol=6),paste("<a href=\"file:",relativePath(file.path(WkgDir,"DiffBind",paste("Contrast",i,"_WithLinks.csv",sep="")),file.path(WkgDir,"HTML_Report","HTMLReport.html")),"\">",paste("Contrast",i,"_WithLinks.csv",sep=""),"</a>",sep=""),DiffMacsIGVLinks,DiffMacsIGVLinksFresh)
#     }else{
#      MacsDiffMat  <- rbind(MacsDiffMat,cbind(MacsContrast[i,],matrix(ScoresTemp,nrow=1,ncol=6),relativePath(file.path(WkgDir,"DiffBind",paste("Contrast",i,"_WithLinks.csv",sep="")),file.path(WkgDir,"HTML_Report","HTMLReport.html")),DiffMacsIGVLinks,DiffMacsIGVLinksFresh))
#     }
     
#   }
#    colnames(MacsDiffMat)[6:14] <- c(paste("FDR <",FDRCutOff,sep=""),paste("FDR <",FDRCutOff,", Up in Group1",sep=""),paste("FDR <",FDRCutOff,", Up in Group2",sep=""),paste("P.Value <",PvalueCutOff,sep=""),paste("P.Value <",PvalueCutOff,", Up in Group1",sep=""),paste("P.Value <",PvalueCutOff,", Up in Group2",sep=""),"Diffential_Results","Current IGV session links","New IGV session links")
#    MacsDiffMat <- as.data.frame(rbind(as.matrix(MacsDiffMat),rep("",ncol(MacsDiffMat))))
#Gvis_DiffMacs <- gvisTable(MacsDiffMat, options=list(width=1750, height=35*nrow(MacsDiffMat)))
#Gvis_DiffMacs <- ReformatVisJS(Gvis_DiffMacs)


#   OccupancyHeatMap <-  file.path(WkgDir,"DiffBind","Occupancy_CorHeatmap_macs.png")
#   OverLapRate <- file.path(WkgDir,"DiffBind","Overlap_Rate_macs.png")
#   AffinityHeatMap <- file.path(WkgDir,"DiffBind","Affinity_CorHeatmap_macs.png")
#hwrite("",br=T,p)

#tryCatch(

#addImage("Overlap_Rate_macs.png",file.path(WkgDir,"DiffBind","Overlap_Rate_macs.png"),"","",p)

#,error=function(e) ErrorForHTML(p),

#finally=print("Ave Peak TSS image inserted"))
#hwrite("Figure 19 shows the number of peaks which overlap by at least 1bp across differing sample numbers and so provides an immediate reference to the 
#reproducibility of peaks across all samples. Marked in red is the number of overlapping peaks used for further diffbind analysis.",br=T,p)
#hwrite("",br=T,p)


#tryCatch(

#addImage("Occupancy_CorHeatmap_macs.png",file.path(WkgDir,"DiffBind","Occupancy_CorHeatmap_macs.png"),"","",p)

#,error=function(e) ErrorForHTML(p),

#finally=print("Ave Peak TSS image inserted"))

#hwrite("",br=T,p)
#hwrite("Figure 20 shows clustering of samples and heatmap of correlation for the co-occurence of merged peaks between samples.",br=T,p)
#hwrite("",br=T,p)



#tryCatch(

#addImage("Affinity_CorHeatmap_macs.png",file.path(WkgDir,"DiffBind","Affinity_CorHeatmap_macs.png"),"","",p)

#,error=function(e) ErrorForHTML(p),

#finally=print("Ave Peak TSS image inserted"))
#hwrite("",br=T,p)
#hwrite("Figure 21 shows clustering of samples and heatmap of correlation for the normalised counts within merged peaks between samples.",br=T,p)
#hwrite("",br=T,p)

#hwrite(hwrite("Table 6 - All-pairwise contrasts for differential binding between tissues, factors and/or conditions",name="T6"),heading=3,p,center=TRUE)
#cat(createGoogleGadget(Gvis_DiffMacs),file=p)
#hwrite("",br=T,p)


#}
#cat(file=p,paste(Divs1[[3]][2],"\n",sep=""))




