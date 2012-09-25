#!/usr/bin/python

from Bio import SeqIO
import sys

resultList = []

for seq in SeqIO.parse(open(sys.argv[1]), 'gb'):
    description = seq.description.lower()
    if seq.name.startswith("NR_"):
        if not "pseudogene" in description:
            if not "transcript variant" in description:
                if not "transcribed rna" in description:
                    if not "small nucleolar rna" in description:
                        #hackity hack hack
                        if seq.description.lower()[:-2].endswith("non-coding rna."):
                            resultList.append(seq)
                            print seq.description

    
    

SeqIO.write(resultList, open(sys.argv[2], 'w'), 'genbank')
