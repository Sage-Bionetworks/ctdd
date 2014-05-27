hit_call_molecular_copy<-function(inputfile,inputfile1,inputfile2,outfolder){
  library(gplots)
  M0<-read.delim(inputfile2)
  
  q0<-subset(M0,M0[,"Group"] == "Library")
  ee<-which(q0$Gene_symbol =="EMPTY" | q0$Gene_symbol == "Control" | q0$Gene_symbol == "Contr")
  if(length(ee)>0){
    q0<-q0[-ee,]  
  }
  
  INPUT<-read.delim(inputfile1)
  load(paste(inputfile,"COPY.Rdata",sep = ""))
  COPY.Sanger<-copy.Sanger
  COPY.TCGA<-copy.TCGA
  
  #### GSEA KS test should be applied here
  source("~/Sage-Analysis-Pipeline/PathwayAnalysis/myPathwayAnalysis.R")
  geneset<-as.character(as.matrix(INPUT$Gene_symbol)) # hits only 
  Geneset<-as.character(as.matrix(q0$Gene_symbol)) # all libraries
  
  # Cell Line
  REF<- COPY.Sanger$Mean
    
  # hits only
  pathwayAnalysis<-myPathwayAnalysis$new()
  pathwayAnalysis$gsea(REF,geneset,np=1000,w =1)
  pathwayAnalysis$gseaResult$p.value
  png(paste(outfolder,"/COPY/GSEA_Sanger_COPY_hits.png",sep = ""),width = 800,height = 800)
  pathwayAnalysis$gseaPlot(REF,pathwayAnalysis$gseaResult$geneset)
  dev.off()
  # all libraries
  pathwayAnalysis<-myPathwayAnalysis$new()
  pathwayAnalysis$gsea(REF,Geneset,np=1000,w =1)
  pathwayAnalysis$gseaResult$p.value
  png(paste(outfolder,"/COPY/GSEA_Sanger_COPY_libraries.png",sep = ""),width = 800,height = 800)
  pathwayAnalysis$gseaPlot(REF,pathwayAnalysis$gseaResult$geneset)
  dev.off()
  
  # TCGA
  REF1<- COPY.TCGA$Mean
  
  pathwayAnalysis1<-myPathwayAnalysis$new()
  pathwayAnalysis1$gsea(REF1,geneset,np=1000,w =1)
  pathwayAnalysis1$gseaResult$p.value
  png(paste(outfolder,"/COPY/GSEA_TCGA_COPY_hits.png",sep = ""),width = 800,height = 800)
  pathwayAnalysis1$gseaPlot(REF1,pathwayAnalysis1$gseaResult$geneset)
  dev.off()
  
  pathwayAnalysis1<-myPathwayAnalysis$new()
  pathwayAnalysis1$gsea(REF1,Geneset,np=1000,w =1)
  pathwayAnalysis1$gseaResult$p.value
  png(paste(outfolder,"/COPY/GSEA_TCGA_COPY_libraries.png",sep = ""),width = 800,height = 800)
  pathwayAnalysis1$gseaPlot(REF1,pathwayAnalysis1$gseaResult$geneset)
  dev.off()
  
  
  # Heatmap Sanger
  pos.CG<-which(COPY.Sanger$Mean >= mean(COPY.Sanger$Mean)+sd(COPY.Sanger$Mean))
  MAT<-COPY.Sanger$copy
  hit.copy.sanger<-match(INPUT$Gene_symbol,rownames(COPY.Sanger$copy))
  library.copy.sanger<-match(q0$Gene_symbol,rownames(COPY.Sanger$copy))
  # Hit with DOE(differentially overexpressed)
  pos.hit.copy.sanger<-intersect(hit.copy.sanger[which(!is.na(hit.copy.sanger))],pos.CG)
  pos.library.copy.sanger<-library.copy.sanger[which(!is.na(library.copy.sanger))]
  
  Hit.copy.sanger<-rownames(COPY.Sanger$copy)[pos.hit.copy.sanger]
  write.table(Hit.copy.sanger,file = paste(outfolder,"/COPY/hit_Copynumbergained_Sanger.txt",sep = ""),row.names=F,col.names=F,quote=F)
  
  # because of huge matrix dimensionality, let's only concentrate libraries and hits when drawing heatmap
  MAT.library<-MAT[pos.library.copy.sanger,]
  D2<-apply(MAT.library,1,rank)
  
  pos<-match(Hit.copy.sanger,rownames(MAT.library))
  A1<-rep("green",nrow(MAT.library))
  A1[pos]<-"red"
  
  mycol <- colorpanel(n=99,low="green",mid="white",high="red")
  ColColorBar2<-rep("grey",ncol(MAT.library))
  ColColorBar2[1:ncol(COPY.Sanger$copy)]<-"darkmagenta"
  
  png(paste(outfolder,"/COPY/Heatmap_Sanger_overexpressed.png",sep = ""),width = 1600,height = 3200)
  aaa<-heatmap.2(t(D2), col=mycol,Rowv = T,Colv=T,
                 xlab = NULL,
                 ylab = NULL,
                 scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
                 trace="none", cexRow=1, cexCol=1, mar = c(6,10),symkey=F,breaks=seq(min(D2),max(D2),length = length(mycol)+1),
                 ColSideColors = ColColorBar2,
                 RowSideColors = A1)
  dev.off()
  
  
  
  # TCGA Heatmap
  pos.CG1<-which(COPY.TCGA$Mean >= mean(COPY.TCGA$Mean)+sd(COPY.TCGA$Mean))
  MAT<-COPY.TCGA$copy
  hit.copy.tcga<-match(INPUT$Gene_symbol,rownames(COPY.TCGA$copy))
  library.copy.tcga<-match(q0$Gene_symbol,rownames(COPY.TCGA$copy))
  # Hit with DOE(differentially overexpressed)
  pos.hit.copy.tcga<-intersect(hit.copy.tcga[which(!is.na(hit.copy.tcga))],pos.CG1)
  pos.library.copy.tcga<-library.copy.tcga[which(!is.na(library.copy.tcga))]
  
  Hit.copy.tcga<-rownames(COPY.TCGA$copy)[pos.hit.copy.tcga]
  write.table(Hit.copy.tcga,file = paste(outfolder,"/COPY/hit_Copynumbergained_TCGA.txt",sep = ""),row.names=F,col.names=F,quote=F)
  
  # because of huge matrix dimensionality, let's only concentrate libraries and hits when drawing heatmap
  MAT.library<-MAT[pos.library.copy.tcga,]
  D2<-apply(MAT.library,1,rank)
  
  pos<-match(Hit.copy.tcga,rownames(MAT.library))
  A1<-rep("green",nrow(MAT.library))
  A1[pos]<-"red"
  
  mycol <- colorpanel(n=99,low="green",mid="white",high="red")
  ColColorBar2<-rep("grey",ncol(MAT.library))
  ColColorBar2[1:ncol(COPY.TCGA$copy)]<-"darkmagenta"
  
  png(paste(outfolder,"/COPY/Heatmap_TCGA_Copynumbergained.png",sep = ""),width = 1600,height = 3200)
  aaa<-heatmap.2(t(D2), col=mycol,Rowv = T,Colv=T,
                 xlab = NULL,
                 ylab = NULL,
                 scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
                 trace="none", cexRow=1, cexCol=1, mar = c(6,10),symkey=F,breaks=seq(min(D2),max(D2),length = length(mycol)+1),
                 ColSideColors = ColColorBar2,
                 RowSideColors = A1)
  dev.off()
  
}