import sys
import collections


psl = collections.namedtuple('psl', 'matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, blockCount,blockSizes, qStarts, tStarts, qSeqs, tSeqs')
#parses out header


pad = 1000
blat_files = sys.argv[1:]
#print blat_files

names = open("names.txt", 'w')

for blat_file in blat_files:
    print >> sys.stderr,  blat_file
    chr = blat_file.strip().split(".")[1]

    blat_file = open(blat_file)

    blat_file.next()
    blat_file.next().strip().split("\t")
    blat_file.next().strip().split("\t")
    blat_file.next().strip().split("\t")
    blat_file.next().strip().split("\t")

    for line in blat_file:
        psl_result = psl._make(line.strip().split("\t"))
        
        #print only if there are no gaps and more than 90% of the target and query match
        start, stop =  [int(x) for x in psl_result.qName.split(",")[1].split("_")[0].split("-")]
        len = (stop - start) / 3 #convert to NT lengths
        aligned = int(psl_result.matches) + int(psl_result.misMatches) * 1.0
        
        percent_aligned = aligned / len
        
        #print len
        #print percent_aligned
        
        if percent_aligned > .90 and len > 30 and int(psl_result.tBaseInsert) < 100:
        #if int(psl_result.blockCount) == 1 and percent_aligned > .90: 
            #print psl_result
            start = int(psl_result.tStart) + 1
            print "\t".join([chr, str(start - pad), str(int(psl_result.tEnd) + pad), psl_result.qName, "0", psl_result.strand[1]])
            names.write(psl_result.qName + "\n") 
        #print "len", stop - start,
        #print "aligned",  int(psl_result.matches) + int(psl_result.misMatches) + int(psl_result.repMatches)
        
        
