import sys
import collections
from optparse import OptionParser

def cb(option, opt_str, value, parser):
    args=[]
    for arg in parser.rargs:
        if arg[0] != "-":
            args.append(arg)
        else:
            del parser.rargs[:len(args)]
            break
    if getattr(parser.values, option.dest):
        args.extend(getattr(parser.values, option.dest))
    setattr(parser.values, option.dest, args)
            

parser = OptionParser()

parser.add_option("--blat_files", action="callback", callback=cb, dest="blat_files", help="blat files to parse")
parser.add_option("--gap_size", dest="gap_size", type = 'int', help="only alignments with gap size or less gaps will be selected")
parser.add_option("--precent_alignment", '-p', dest="percent_aligned", type = 'float', help = "pairs of alignments higher than this will be selected",)
parser.add_option("--length", "-l", dest="length", type = 'int', help = "Only blatted sequences large than length will be selected")
parser.add_option("--padding", dest="padding", type = 'int', help = "Padding for the outputed bed file (so we can find an ORF)")
                    
(options, args) = parser.parse_args()

psl = collections.namedtuple('psl', 'matches, misMatches, repMatches, nCount, qNumInsert, qBaseInsert, tNumInsert, tBaseInsert, strand, qName, qSize, qStart, qEnd, tName, tSize, tStart, tEnd, blockCount,blockSizes, qStarts, tStarts, qSeqs, tSeqs')
#parses out header


pad = options.padding

#print blat_files

names = open("names.txt", 'w')

for blat_file in options.blat_files:
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
        
        if percent_aligned > options.percent_aligned and len > options.length and int(psl_result.tBaseInsert) < options.gap_size:

            start = int(psl_result.tStart) + 1
            print "\t".join([chr, str(start - pad), str(int(psl_result.tEnd) + pad), psl_result.qName, "0", psl_result.strand[1]])
            names.write(psl_result.qName + "\n") 

        
        
