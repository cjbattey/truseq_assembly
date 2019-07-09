#!/bin/bash
# Scaffold alignment pipeline for avian whole genomes (for Mac OSX)
# run time is ~12hrs. 

#install homebrew if needed:
#/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
cd ~/Desktop/
brew install genomeTools
brew install mummer

#download repeat-masked zebra finch genome and split fasta entries. 
curl "http://hgdownload.cse.ucsc.edu/goldenPath/taeGut2/bigZips/taeGut2.fa.masked.gz" -o "tgut2.fa.gz"
gunzip tgut2.fa.gz
mkdir tgut2_split
gt splitfasta -splitdesc tgut2_split/ tgut2.fa #don't open splitfasta/ in finder
rm tgut2.*  
cd ./tgut2_split/ #remove short contigs, leave the chromosome-level scaffolds. (?) 
find . -name "*_*" -delete
cd ~/Desktop/

#download anna's hummingbird genome. unzip failed on last dl for unknown reasons...
#curl "ftp://climb.genomics.cn/pub/10.5524/101001_102000/101004/Calypte_anna.fa.gz" -o "Canna.fa.gz"
#gunzip Canna.fa.gz

#output directories
mkdir mummer_out
mkdir ./mummer_out/coords
mkdir ./mummer_out/alignments

#Run mummer alignment against each chromosome. 
refSeqs=tgut2_split/*

for chr in $refSeqs

	do
	echo "aligning to $chr"
	
	nucmer $chr ./Canna.fa
	
	show-coords ./out.delta > ~/Desktop/mummer_out/coords/${chr##*/}.coords #Check regex if changing paths.
	
	mv out.delta ~/Desktop/mummer_out/alignments/${chr##*/}.delta
	
done

