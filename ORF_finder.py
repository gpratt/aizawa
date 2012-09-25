from Bio import SeqIO
import sys
from optparse import OptionParser

#need to take into account reverse complement, that will be in strand information...
def orf_finder(seq, padding, strand):
    reading_frame = padding % 3
    search_start = padding / 3
    #strand logic goes here...

    trans = seq.seq[reading_frame:].translate()
    #assume that negative direction is upstream and positive direction is downstream...

    #find start site
    cur = search_start
    while trans[cur] != 'M':
        cur -= 1
        if cur < 0:
            return None
        
    start = cur
    
    #find stop site
    while trans[cur] != "*":
        cur += 1
        if cur >= len(trans):
            return None
        
    stop = cur

    #    print start, stop
    if stop < search_start:
        return None #ORF is non-existant...

    
    else: return trans[start:stop]

if __name__ == "__main__":

    parser = OptionParser()
    
    parser.add_option("--fasta", dest="fasta", help="fasta (faa) file of sequences to search for ORFs, only searches one reading frame, the frame of the previously aligned sequence")
    parser.add_option("--padding", dest='padding', type= 'int', help="padding correction (so we know where to look for ORFs at the start, should be the same as the parse arg)")
    parser.add_option("--alignment_names", dest='alignment_names', help="names of alignments from parse_blat, lost in the fastafrombed conversion")
    (options, args) = parser.parse_args()

    names = open(options.alignment_names)
    result_list = []
    for seq in SeqIO.parse(open(options.fasta), 'fasta'):
        name = names.next().strip()
        
        result = orf_finder(seq, options.padding, "+")
        if result is not None:
            print ">%s\t%s" % (seq.id, name) 
            print result
