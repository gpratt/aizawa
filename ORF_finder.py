from Bio import SeqIO
import sys

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
    result_list = []
    names = open("names.txt")
    for seq in SeqIO.parse(open(sys.argv[1]), 'fasta'):
        result = orf_finder(seq, int(sys.argv[2]), "+")
        name = names.next().strip()
        if result is not None:
            print ">%s\t%s" %(seq.id, name) 
            print result
            #result_list.append(result)
            
    #SeqIO.write(result_list, open(sys.argv[3], 'w'), 'fasta')
