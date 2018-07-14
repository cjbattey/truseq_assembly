#selasphorus Fst scan visualization
#get angsd windowed fst output
setwd("~/Dropbox/selasphorus/")
library(data.table);library(plyr);library(ggplot2);library(magrittr);library(foreach);library(ggrepel)
fst <- fread("ruf_sas_fstwindow.txt",sep="\t")[,-1]
colnames(fst) <- c("contig","window_start","snps","Fst")
#tajD <- fread("rufous.theta_windows.gz.pestPG",sep="\t")[,-1]
gff <- read.table("~/Dropbox/augustus/Canna_MUGM01_fwd.gff",comment.char = "#",skip = 19,sep="\t",stringsAsFactors = F)
colnames(gff) <- c("contig","whytho","feature","start","end","score","strand","frame","attribute")

#summarize mummer output
#read in and clean up formatting for coords files
files <- list.files("~/Dropbox/mummer_alignments/mummer_out/Canna_MUGM01/coords/",full.names = T)
files <- files[grepl("\\.coords",files)]
i <- 1
for(f in files){
  if(i==1){
    dat <- fread(f,sep=" ",data.table=F)
    i=i+1
  } else {
    tmp <- fread(f,sep=" ",data.table=F)
    dat <- rbind(tmp,dat)
    i=i+1
  }
}
colnames(dat) <- c("refStart","refStop","sep1","qStart","qStop","sep2","refLength","qLength","sep3",
                   "p.identity","sep4","names")
dat$refName <- strsplit(dat$names,"\t") %>% sapply(function(e) unlist(e)[1])
dat$qName <- strsplit(dat$names,"\t") %>% sapply(function(e) unlist(e)[2])
dat <- arrange(dat,refName,refStart)
sum <- ddply(dat,
             .(qName,refName),
             summarize,
             totalMatch=sum(qLength),
             refStart=min(refStart))
sum <- arrange(sum,refName,refStart,totalMatch)
sum <- subset(sum,totalMatch>1000)                                                 
sum <- ddply(sum,.(qName),function(e){                                             #get contigs with 1 hit > 1000bp
  a <- subset(e,refName==e$refName[e$totalMatch==max(totalMatch)])
  if(nrow(a)==1){
    a
  }
})

#chromosome order for pretty plots
chr_order <- c("1","1A","1B","2","3","4","4A",as.character(5:28),"Z","M","NA")

#merge mummer info with angsd windowed Fst's
fst <- merge(fst,sum,by.x="contig",by.y="qName",all.x=T,all.y=F)
fst$chr <- gsub("chr","",fst$refName)
fst$chr[!fst$chr %in% chr_order] <- "NA"
fst$chr <- factor(fst$chr,levels=chr_order)
fst <- arrange(fst,chr,refStart)
fst$row <- 1:nrow(fst)

# tajD <- merge(tajD,sum3,by.x="Chr",by.y="qName",all.x=T,all.y=F)
# tajD$chr <- gsub("chr","",tajD$refName)
# tajD$chr <- factor(tajD$chr,levels=chr_order)
# tajD <- arrange(tajD,chr,refStart)
# tajD$row <- 1:nrow(tajD)

#build second dataset for zebra finch chromosomes
chr_labels <- ddply(fst,.(chr),summarize,mid=median(row),start=min(row),stop=max(row))
chr_labels$chr <- as.character(chr_labels$chr)
chr_labels$chr[chr_labels$chr %in% as.character(21:27)] <- "21-27"
chr_labels$mid[chr_labels$chr=="21-27"] <- median(chr_labels$mid[chr_labels$chr=="21-27"],na.rm=T)
chr_labels$start[chr_labels$chr=="21-27"] <- min(chr_labels$start[chr_labels$chr=="21-27"],na.rm=T)
chr_labels$stop[chr_labels$chr=="21-27"] <- max(chr_labels$stop[chr_labels$chr=="21-27"],na.rm=T)
chr_labels$start[chr_labels$chr=="1B"] <- 0
chr_labels$stop[chr_labels$chr=="1B"] <- 0
chr_labels$Fst <- min(fst$Fst)-.25*min(fst$Fst)
chr_labels <- subset(chr_labels,!is.na(chr) & !duplicated(chr))

# chr_labels_tajD <- ddply(tajD,.(chr),summarize,mid=median(row),start=min(row),stop=max(row))
# chr_labels_tajD$chr <- as.character(chr_labels_tajD$chr)
# chr_labels_tajD$chr[chr_labels_tajD$chr %in% as.character(21:27)] <- "21-27"
# chr_labels_tajD$mid[chr_labels_tajD$chr=="21-27"] <- median(chr_labels_tajD$mid[chr_labels_tajD$chr=="21-27"],na.rm=T)
# chr_labels_tajD$start[chr_labels_tajD$chr=="21-27"] <- min(chr_labels_tajD$start[chr_labels_tajD$chr=="21-27"],na.rm=T)
# chr_labels_tajD$stop[chr_labels_tajD$chr=="21-27"] <- max(chr_labels_tajD$stop[chr_labels_tajD$chr=="21-27"],na.rm=T)
# chr_labels_tajD$start[chr_labels_tajD$chr=="1B"] <- 0
# chr_labels_tajD$stop[chr_labels_tajD$chr=="1B"] <- 0
# chr_labels_tajD$Tajima <- min(tajD$Tajima)
# chr_labels_tajD <- subset(chr_labels_tajD,!is.na(chr) & !duplicated(chr))

#rolling means
fst$rollmean <- zoo::rollmean(fst$Fst,50)

#Fst ~ distance to a gene
cds <- subset(gff,feature=="CDS")
cdfst <- subset(fst,contig %in% cds$contig)
cdfst$dist_to_cds <- pbapply(cdfst,1,function(e){ #~5 min
  tmp <- subset(cds,contig==e[1])
  min(abs(as.numeric(e[2])-tmp$start))
})
ggplot(data=cdfst,aes(x=dist_to_cds,y=Fst))+theme_minimal()+
  geom_point(alpha=0.1,size=0.5)

#closest coding sequences to top 100 Fst windows
top <- subset(cdfst,Fst>quantile(fst$Fst,0.999) & dist_to_cds<10000) %>% arrange(desc(Fst))
top_cds_hits <- pbapply(top,1,function(e){
  cds$attribute[which.min(abs(cds$start-as.numeric(e[2])))]
})
top_cds_hits <- top_cds_hits %>% strsplit(";") %>% sapply(function(e) e[2]) %>% gsub(" gene_id ","",.)

genes <- readLines("~/Dropbox/augustus/Canna_MUGM01_fwd.gff")
startrows <- grep("start gene",genes) %>% genes[.] %>% grep(top_cds_hits)
endrows <- genes[grep("end gene",genes)]


# #plot
# png(width=6,height=1.5,units="in",res=600,file="~/Desktop/rufous_tajD.png")
# ggplot(data=tajD,aes(x=row,y=Tajima,col=Chr))+theme_bw()+
#   theme(legend.position="none",
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         panel.grid = element_blank(),
#         panel.grid.major.y=element_line(color="grey20",size=0.3),
#         panel.grid.minor.y=element_line(color="grey20",size=0.1),
#         axis.title = element_text(size=7),
#         axis.text = element_text(size=7),
#         axis.ticks.x=element_blank())+
#   scale_color_manual(values=rep(c("grey60","grey80"),length(levels(factor(tajD$Chr)))/2+1))+
#   scale_y_continuous(limits=c(-3,max(tajD$Tajima)),breaks=c(-1,0,1,2),minor_breaks = c(-1.5,-0.5,0.5,1.5))+
#   geom_point(size=0.01)+
#   geom_point(data=subset(tajD,Tajima>=quantile(tajD$Tajima,0.99)),size=0.01,col="red")+
#   geom_segment(data=chr_labels_tajD,aes(x=start+25,xend=stop-25,y=Tajima,yend=Tajima,col=NA),col="black")+
#   geom_text_repel(data=chr_labels_tajD,aes(label=chr,x=mid,y=Tajima,col=NA),force=2,ylim=c(-3,-1.75),
#                 col="black",size=2,angle=0,direction="y",box.padding = 0.15,
#                 segment.size=0.2)
# dev.off()

png(width=6,height=1.5,units="in",res=600,file="~/Desktop/west_east_fst.png")
ggplot(data=fst,aes(x=row,y=Fst,col=chr))+theme_bw()+
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y=element_line(color="grey60",size=0.2),
        panel.grid.minor.y=element_line(color="grey60",size=0.1),
        axis.title.y = element_text(size=7,hjust=.75),
        axis.text = element_text(size=7),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2))+
  scale_color_manual(values=rep(c("grey60","grey80"),length(levels(factor(fst$chr)))/2+1))+
  scale_y_continuous(limits=c(-max(fst$Fst),max(fst$Fst)),breaks=seq(0,max(fst$Fst),0.1),minor_breaks = NULL)+
  geom_point(size=0.4,shape=21,stroke=0.2)+
  geom_point(data=subset(fst,Fst>=quantile(fst$Fst,0.999)),shape=21,stroke=0.4,size=0.5,col="red")+
  geom_line(aes(y=rollmean),lwd=0.15,col="black")+
  geom_segment(data=chr_labels,aes(x=start+25,xend=stop-25,y=Fst,yend=Fst,col=NA),col="black")+
  geom_text_repel(data=chr_labels,aes(label=chr,x=mid,y=Fst,col=NA),force=2,ylim=c(-max(fst$Fst),0),
                  col="black",size=2,angle=0,direction="y",box.padding = 0.15,
                  segment.size=0.2)
dev.off()
