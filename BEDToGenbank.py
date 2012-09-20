#!/usr/bin/python
#Gets genbank fies from BED Formatted files
import sys
from Bio import SeqIO
from  Bio import Entrez
Entrez.email = "gpratt@ucsd.edu"

geneIds = []
for line in open(sys.argv[1]):
    line = line.split()
    geneIds.append((line[3], line[0], line[1], line[2], line[5]))
    
geneList = ""
count = 0
records = []
for gene, chromosome, start, stop, strand in geneIds:
    geneList += gene + ", "
    count += 1
    if count >= 50:

        handle = Entrez.efetch(db = "nucleotide", id=geneList[:-2], rettype='gb')
        
        records += SeqIO.parse(handle, 'genbank')
    
        count = 0
        geneList = ""

for t, record in zip(geneIds, records):

    gene, chromosome, start, stop, strand = t
    record.description = chromosome + ":" + start + "-" + stop + ":" + record.description + "," + strand

SeqIO.write(records, open(sys.argv[2], 'w'), 'genbank')
