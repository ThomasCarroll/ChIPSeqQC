import sys
TrimTo = int(sys.argv[1])
from Bio.SeqIO.QualityIO import FastqGeneralIterator
trim = TrimTo
#handle = open(FQTrimed, "w")
write = sys.stdout.write
for title, seq, qual in FastqGeneralIterator(sys.stdin) :
    write("@%s\n%s\n+\n%s\n" % (title, seq[:trim], qual[:trim]))
