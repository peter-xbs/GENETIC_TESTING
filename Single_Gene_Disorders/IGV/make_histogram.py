import pysam
import sys

'''Input file is a BAM file with an index file'''

if len(sys.argv) > 3:
    chr = sys.argv[2]
    start = int(sys.argv[3])
    stop = int(sys.argv[4])
elif len(sys.argv) > 2:
    chr = sys.argv[2]
    start = None
    stop = None

samfile = pysam.Samfile(sys.argv[1], 'rb')

for pileupcolumn in samfile.pileup(chr, start, stop):
    print(' '.join([chr, str(pileupcolumn.pos), str(pileupcolumn.pos+1), str(pileupcolumn.n)]))
