#Pairwise search and alignment of two different potentally protein coding genes in human and mouse genome, joined by location in genome
from Bio import SeqIO
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
matrix = matlist.blosum62
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--fasta1", dest="fasta1", help="first fasta file, must be the MouseRefGene_09132012_NR_not_overlapping_filtered_ORFs.fna or Human version of that because I parse the header")
parser.add_option("--fasta2", dest="fasta2", help="second fasta file, must be output from ORF_finder.py, again special parsing")
parser.add_option("--transmembrane_out", dest="transmembrane_out", help="file to write transmembrane stuff to")
parser.add_option("--cutoff_score", dest="cutoff_score", type = 'float', help= "the cutoff for alignments to print out")
(options, args) = parser.parse_args()


#argv[1] is fasta file (the mouse_aa file)
#argv[2] is fasta file

#Generate dict of seqs for fasta file 1
    
         
        
seqs = {}
trans1 = []
trans2 = [] 
for seq in SeqIO.parse(open(options.fasta1), 'fasta'):
    key = seq.id

    #print key
    if key in seqs:
        pass
        #print "duplicate found, probably shouldn't happen"
        
    seqs[key] = seq

#dict created, on to the next one
for seq in SeqIO.parse(open(options.fasta2), 'fasta'):
    
    key = seq.description.split()[1]
    if key in seqs:
        
        
        alignments =  pairwise2.align.globaldx(seq.seq[:-1], seqs[key].seq[:-1], matrix, one_alignment_only = True)

        if len(alignments) > 0:
            alignment = alignments[0]
            
            score = alignment[2]
            length = alignment[4]
            normalized_score = score / length
            if normalized_score > options.cutoff_score:
                
                print seqs[key].description, score / length
                print format_alignment(*alignments[0])
                
                trans1.append(seqs[key])
                trans2.append(seq)
            
SeqIO.write(trans1, open("fasta1_" + options.transmembrane_out, 'w'), 'fasta')
SeqIO.write(trans2, open("fasta2_" + options.transmembrane_out, 'w'), 'fasta')
        
