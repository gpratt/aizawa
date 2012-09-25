import subprocess
import sys
#fa file to blat
fa = sys.argv[1]

print fa

for genome in open("/nas3/gpratt/ens/hg19.genome"):
    line = genome.strip().split()
    #wraps blat and blats against hg19 without getting memory errors

    call = ["blat","-out=pslx", "-t=dnax", "-q=prot", "/nas3/yeolab/Genome/ucsc/hg19/chromosomes/hg19.2bit:%s:1-%s" %(line[0], line[1]), fa, "hg.%s" % (line[0])]
    print call
    print subprocess.call(call)
