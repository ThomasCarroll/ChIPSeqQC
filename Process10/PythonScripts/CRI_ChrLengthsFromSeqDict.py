import pysam
import sys

BamFile = sys.argv[1]
OutFile = sys.argv[2]

#BamFile = "/lustre/mib-cri/carrol09/Work/AisleenNew/homo_sapiens.fa.sam"
#OutFile = "/lustre/mib-cri/carrol09/Work/AisleenNew/homo_sapiens.fa.txt"
SamFile = pysam.Samfile(BamFile,"r")
Header = SamFile.header
chrlengthsfile = open(OutFile,"w")
for ref in Header["SQ"]:
	chrlengthsfile.write(str(ref["SN"])+"\t"+str(ref["LN"])+"\n")

SamFile.close()
chrlengthsfile.close()

	