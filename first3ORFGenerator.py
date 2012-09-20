#!/usr/bin/python
#takes a bed file a number of orfs and the min length of orf outputs the first n orfs if that orf is more than the min length of the orf specificed

#finds all open reading frames in a specific RNA sequence (no reverse tranlsation) returns a list of ORFs
import optparse
from copy import deepcopy
from Bio import SeqIO
def openReadingFrameFinder(seq, numorfs):
    result_list = []
    matchCount = 0 
    for frame in range(3):

        trans = seq.seq[frame:].translate()
        aa_start = 0
        while aa_start < len(trans):

            if trans[aa_start] == 'M':
                aa_end = trans.find('*', aa_start) 
                
                #Calculate the Nucleotide start and stop locations
                start = frame + aa_start * 3
                end = min(len(seq), frame + aa_end * 3 + 3)
                
                #create a new sequence for the info
                storedSeq = deepcopy(seq)
                storedSeq.id = storedSeq.id + "," + str(start) + "-" + str(end)
                
                #add the result to matches
                matchCount += 1
                result_list.append((start,storedSeq[start:end]))
                
                #advance the pointer to the next possible place for an orf
                
#print "end: " + str(aa_end)
 #               print "start: " + str(aa_start)
                aa_start = aa_end
  
#              print "start after fixing: " + str(aa_start)
                if (len(result_list) > (numorfs * (frame + 1)) or
                    aa_end == -1): #if there is no last stop codon
                    break
            
            aa_start += 1
    return result_list

if __name__ == "__main__":
    
    parser = optparse.OptionParser(description='Find all ORFs in sequence given')

    parser.add_option('--file', '-f', type='string', help='The file to be parsed')

    parser.add_option('--numorf', '-n', default=3, type=int, help='The number or orfs after the first to analyze')
    parser.add_option('--minorf', '-m', default=90, type=int, help='The shortest uorf that will be printed on detection')

    parser.add_option('--out', '-o', type='string', default="out.faa", help="The output file")

    args, remainder = parser.parse_args()
    
    result_list = []
    count = 0 
    for seq in SeqIO.parse(open(args.file), 'genbank'):
        print "working on seq"
        count += 1

        ORFs = openReadingFrameFinder(seq, args.numorf)

    #gets the first n orfs and then filters based on k min nucleotide sequence

        firstThree = sorted(ORFs, key = lambda x : x[0])[:args.numorf]

        x= [ x[1] for x in firstThree if len(x[1]) > args.minorf ] 
        
        result_list += x

    SeqIO.write(result_list, open(args.out, 'w'), 'fasta')
