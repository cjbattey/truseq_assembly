#angsd sfs plots

#sfs=output of realSFS
#nind1/2=number of diploid individuals in population 1/2
#pop1/2=character string labels for each population

plot_angsd_sfs_2d <- function(sfs,nind1,nind2,pop1,pop2,log=T){
  library(data.table);library(ggplot2)
  sfs <- unlist(fread(sfs,data.table = F)[1,])
  sfs <- matrix(sfs,nrow=2*nind1+1,ncol=2*nind2+1)
  meltsfs <- melt(sfs)
  meltsfs$Var1 <- meltsfs$Var1-1 #first category is 0 
  meltsfs$Var2 <- meltsfs$Var2-1
  meltsfs$logvalue <- log(meltsfs$value)
  meltsfs$logvalue[meltsfs$logvalue == -Inf] <- NA
  meltsfs <- meltsfs[-1,]
  if(log){
    ggplot(data=meltsfs,aes(x=Var1,y=Var2,fill=logvalue))+
      theme_classic()+
      scale_x_continuous(breaks=0:(2*nind1))+scale_y_continuous(breaks=0:(nind2*2))+
      scale_fill_distiller(palette = "RdYlBu",na.value = "white",name="log(sites)")+
      xlab(paste0("n derived: ",pop1))+ylab(paste0("n derived: ",pop2))+
      geom_tile(alpha=0.9)
  } else {
    ggplot(data=meltsfs,aes(x=Var1,y=Var2,fill=value))+
      theme_classic()+
      scale_x_continuous(breaks=0:(2*nind1))+scale_y_continuous(breaks=0:(nind2*2))+
      scale_fill_distiller(palette = "RdYlBu",na.value = "white",name="log(sites)")+
      xlab(paste0("n derived: ",pop1))+ylab(paste0("n derived: ",pop2))+
      geom_tile(alpha=0.9)
  }
}
