#Pairwise search and alignment of two different potentally protein coding genes in human and mouse genome, joined by location in genome
from Bio import SeqIO
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
matrix = matlist.blosum62
#argv[1] is fasta file (the mouse_aa file)
#argv[2] is fasta file

#Generate dict of seqs for fasta file 1
    
         
        
seqs = {}
trans = [] 
for seq in SeqIO.parse(open(sys.argv[1]), 'fasta'):
    key = seq.id

    #print key
    if key in seqs:
        pass
        #print "duplicate found, probably shouldn't happen"
        
    seqs[key] = seq

#dict created, on to the next one
for seq in SeqIO.parse(open(sys.argv[2]), 'fasta'):
    
    key = seq.description.split()[1]
    if key in seqs:
        #print seq.seq
        #print seqs[key].seq
        
        alignments =  pairwise2.align.globaldx(seq.seq[:-1], seqs[key].seq[:-1], matrix, one_alignment_only = True)

        if len(alignments) > 0:
            alignment = alignments[0]
            
            score = alignment[2]
            length = alignment[4]
            normalized_score = score / length
            #if normalized_score > 3.5:
            
            print key, score / length
            print format_alignment(*alignments[0])
                
            trans.append(seqs[key])
            
SeqIO.write(trans, open('transmembrane_domain_input.fna', 'w'), 'fasta')
        
