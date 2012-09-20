#!/bin/bash
#/home/gabriel/bioinformatics/rawdata/RefGeneData
#Data Staging:

#Need Other NR
#Need human Not NR other not NR
#Need Mouse NR
#Mouse not NR other not NR
#Download RefGeneOther track from mouse and human UCSC genome browser
#as HumanRefGeneOther.bed and MouseRefGeneOther.bed
#Not sure what this means but downloading Other RefSeq (xenoRefGene) track
#got Whole Gene
#Download Mouse and Human RefGene tracks as
#MouseRefGene.bed and ${HUMAN}.bed
HUMAN=HumanRefGene_09132012
HUMAN_OTHER=HumanRefSeqOther
MOUSE=MouseRefGene_09132012
MOUSE_OTHER=MouseRefSeqOther

#Seperate NR From not NR
grep NR_ ${HUMAN_OTHER}.bed > ${HUMAN_OTHER}_NR.bed
grep NR_ ${MOUSE_OTHER}.bed > ${MOUSE_OTHER}_NR.bed
grep NR_ ${HUMAN}.bed > ${HUMAN}_NR.bed
grep NR_ ${MOUSE}.bed > ${MOUSE}_NR.bed

grep -v NR_ ${HUMAN_OTHER}.bed > ${HUMAN_OTHER}_not_NR.bed
grep -v NR_ ${MOUSE_OTHER}.bed > ${MOUSE_OTHER}_not_NR.bed
grep -v NR_ ${HUMAN}.bed > ${HUMAN}_not_NR.bed
grep -v NR_ ${MOUSE}.bed > ${MOUSE}_not_NR.bed

#Sort for later
sort -k4 ${HUMAN_OTHER}.bed > ${HUMAN_OTHER}.sorted.bed
sort -k4 ${MOUSE_OTHER}.bed > ${MOUSE_OTHER}.sorted.bed

#Combine all not NR transcripts
cat ${MOUSE_OTHER}_not_NR.bed ${MOUSE}_not_NR.bed > MouseAlignedRefGene_not_NR.bed
cat ${HUMAN_OTHER}_not_NR.bed ${HUMAN}_not_NR.bed > HumanAlignedRefGene_not_NR.bed

#find only NR transcripts that don't overlap with anything else
intersectBed -v -a ${MOUSE}_NR.bed -b MouseAlignedRefGene_not_NR.bed > ${MOUSE}_NR_not_overlapping.bed
intersectBed -v -a ${HUMAN}_NR.bed -b HumanAlignedRefGene_not_NR.bed > ${HUMAN}_NR_not_overlapping.bed

#Convert to genbank format to get sequence
./BEDToGenbank.py ${HUMAN}_NR_not_overlapping.bed ${HUMAN}_NR_not_overlapping.gbk

./BEDToGenbank.py ${MOUSE}_NR_not_overlapping.bed ${MOUSE}_NR_not_overlapping.gbk

#Filter by non-coding and pseudogenes (very hurestric)
./noncodingAndPseudogeneFilter.py ${HUMAN}_NR_not_overlapping.gbk ${HUMAN}_NR_not_overlapping_filtered.gbk
./noncodingAndPseudogeneFilter.py ${MOUSE}_NR_not_overlapping.gbk ${MOUSE}_NR_not_overlapping_filtered.gbk

#Find the first 3 ORFs in the remaining transcripts
../../scripts/../scripts/noncodingFiltering/first3ORFGenerator.py -n 100000 -m 30 -f ${HUMAN}_NR_not_overlapping_filtered.gbk -o ${HUMAN}_NR_not_overlapping_filtered_ORFs.faa
../../scripts/../scripts/noncodingFiltering/first3ORFGenerator.py -n 100000 -m 30 -f ${MOUSE}_NR_not_overlapping_filtered.gbk -o ${MOUSE}_NR_not_overlapping_filtered_ORFs.faa

#convert back to bed format
grep NR_ ${HUMAN}_NR_not_overlapping_filtered_ORFs.faa | awk '{print $1}' | sed 's/>\(.*\)\..*,\(.*\)-\(.*\)/\1\t\2\t\3/' > ${HUMAN}_NR_not_overlapping_filtered_ORFs.bed
grep NR_ ${MOUSE}_NR_not_overlapping_filtered_ORFs.faa | awk '{print $1}' | sed 's/>\(.*\)\..*,\(.*\)-\(.*\)/\1\t\2\t\3/' > ${MOUSE}_NR_not_overlapping_filtered_ORFs.bed

sort -k 1 ${HUMAN}_NR_not_overlapping_filtered_ORFs.bed > ${HUMAN}_NR_not_overlapping_filtered_ORFs.sorted.bed
sort -k 1 ${MOUSE}_NR_not_overlapping_filtered_ORFs.bed > ${MOUSE}_NR_not_overlapping_filtered_ORFs.sorted.bed

join -1 4 -2 1 ${HUMAN}.sorted.bed ${HUMAN}_NR_not_overlapping_filtered_ORFs.sorted.bed  > ${HUMAN}_NR_not_overlapping_filtered_ORFs.sorted.joined.bed

join -1 4 -2 1 ${MOUSE}.sorted.bed ${MOUSE}_NR_not_overlapping_filtered_ORFs.sorted.bed > ${MOUSE}_NR_not_overlapping_filtered_ORFs.sorted.joined.bed

transeq --sequence ${HUMAN}_NR_not_overlapping_filtered_ORFs.faa --outseq ${HUMAN}_NR_not_overlapping_filtered_ORFs.fna

transeq --sequence ${MOUSE}_NR_not_overlapping_filtered_ORFs.faa --outseq ${MOUSE}_NR_not_overlapping_filtered_ORFs.fna