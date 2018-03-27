#truseq resequencing SNP calling pipeline
setwd("/media/burke/bigMac/cj/")
library(magrittr);library(ggplot2);library(foreach);library(doMC)
registerDoMC(cores=8)

#######################################################################################
# setup: you'll need R, java (v1.8), samtools, picard tools, gatk, and AdapterRemoval #
# Budget around 5x the initial download size for hard drive space

######################################################################################################
#run AdapterRemoval to trim adapters, merge overlapping paired reads, and drop low-quality base calls#
R1 <- list.files("demultiplexed_reads",full.names=T) %>% grep("R1_001",.,value=T)
R2 <- list.files("demultiplexed_reads",full.names=T) %>% grep("R2_001",.,value=T)
commands <- c()
for(i in 1:40){
  commands[i] <- paste0("AdapterRemoval --file1 ",R1[i],
                        " --file2 ",R2[i],
                        " --basename trimmed/", R1[i] %>% basename() %>% strsplit("_") %>% unlist() %>% .[1],
                        " --trimns --trimqualities --collapse --threads 30")
}

for(i in commands){
  system(i)
}

#foreach(i=commands) %dopar% system(i) #untested parallel version

################################################################
# align trimmed reads to the C. anna reference using bowtie2 ###
# ANHU genome: https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=MUGM01#
R1 <- list.files("selasphorus_assembly/trimmed",full.names=T) %>% grep("pair1",.,value=T)
R2 <- list.files("selasphorus_assembly/trimmed",full.names=T) %>% grep("pair2",.,value=T)
commands <- c()
for(i in 1:40){
  sampleID <- basename(R1[i]) %>% strsplit("\\.") %>% unlist() %>% .[1]
  commands[i] <- paste0("bowtie2",
                        " -p 30",
                        " -x /media/burke/bigMac/cj/anhu_alignment/Canna",
                        " -1 ",R1[i],
                        " -2 ",R2[i],
                        " -S ","/media/burke/bigMac/cj/selasphorus_assembly/bam/",sampleID,".sam;",
                        
                        "samtools view",
                        " -S -b",
                        " /media/burke/bigMac/cj/selasphorus_assembly/bam/",sampleID,".sam",
                        " >",
                        " /media/burke/bigMac/cj/selasphorus_assembly/bam/",sampleID,".bam;",

                        " rm ","/media/burke/bigMac/cj/selasphorus_assembly/bam/",sampleID,".sam;")
}
for(i in 1:40){
  system(commands[i])
}
#foreach(i=commands) %dopar% system(i) #untested parallel version

################################################################
# mark duplicates, sort, and index BAM files ###################
setwd("/media/burke/bigMac/cj/selasphorus_assembly/")
bam <- list.files("bam",full.names=T)
commands <- c()
for(i in 1:40){
  sampleID <- basename(bam[i]) %>% strsplit("\\.") %>% unlist() %>% .[1]
  commands[i] <- paste0("samtools sort -@ 30 ",bam[i]," tmp;",
                        " mv tmp.bam ",bam[i],"; ",
                        " samtools index ",bam[i],";",
                        " java -Xmx16g -jar /media/burke/bigMac/cj/picard.jar MarkDuplicates",
                        " I=",bam[i],
                        " O=dedup/",sampleID,".bam",
                        " M=dedup_metrics/",sampleID,"_dup_metrics.txt;"
                        )
}
for(i in 1:40){
  system(commands[i])
}
#foreach(i=commands) %dopar% system(i) #untested parallel version

##############################################################################################
# add readgroup info to BAM files (this is gatk bullshit - switch to bwa to avoid? may not be necessary for angsd)  #########
setwd("/media/burke/bigMac/cj/selasphorus_assembly/")
bam <- list.files("dedup",full.names=T)
commands <- c()
for(i in 1:40){
  sampleID <- basename(bam[i]) %>% strsplit("\\.") %>% unlist() %>% .[1]
  commands[i] <- paste0("java -Xmx10g -jar /media/burke/bigMac/Dropbox/tools/picard.jar AddOrReplaceReadGroups",
                        " I=",bam[i],
                        " O=tmp.bam",
                        " R=/media/burke/bigMac/Dropbox/selasphorus/anhu_alignment/Canna_MUGM01_haploid.fa",
                        " RGID=1",
                        " RGLB=selasphorus",
                        " RGPL=illumina",
                        " RGPU=pool1",
                        " RGSM=",sampleID,";",
                        
                        "mv tmp.bam ",bam[i],";",
                        
                        "samtools index ",bam[i]
  )
}
for(i in 2:40){
  system(commands[i])
}
#foreach(i=commands) %dopar% system(i) #untested parallel version

################################################################
# call SNP's from BAM files ####################################
#attempt command for unifiedgenotyper (8-10gb RAM & ~8hrs for 7 samples). Note GATK will throw an error if there are any spaces after the "\".
system("java -Xmx40g -jar /media/burke/bigMac/Dropbox/tools/GenomeAnalysisTK.jar -T UnifiedGenotyper \
-R /media/burke/bigMac/Dropbox/selasphorus/anhu_alignment/Canna_MUGM01_haploid.fa \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/AK1.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/AK3.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/CA10.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/CA15.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/CA19.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/CA22.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/CA33.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/CA5.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/CA6.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/NM1.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/NM2.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/NM3.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/NM4.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/NM5.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/NM6.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/NM7.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/OR10.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/OR2.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Scal1.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Scal2.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Scal3.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Scal4.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Scal5.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Scal6.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Scal7.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Ssasin10.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Ssasin11.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Ssasin2.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Ssasin4.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Ssasin5.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Ssasin6.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Ssasin7.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/Ssasin9.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/UT1.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/UT2.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/UT3.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/WA22.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/WA31.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/WA35.bam \
-I /media/burke/bigMac/cj/selasphorus_assembly/dedup/WA36.bam \
-o selasphorus.vcf
")

#note I abandoned gatk and switched to angsd here, but the above should still work to generate static calls with no filters
################################################################
# apply additional filters in vcftools #########################



################################################################
# calculate & plot summary stats for final VCF #################



