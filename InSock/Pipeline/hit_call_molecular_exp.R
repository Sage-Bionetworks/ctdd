hit_call_molecule_exp<-function(inputfile,inputfile1,inputfile2,outfolder){
  library(gplots)
  M0<-read.delim(inputfile2)
  
  q0<-subset(M0,M0[,"Group"] == "Library")
  ee<-which(q0$Gene_symbol =="EMPTY" | q0$Gene_symbol == "Control" | q0$Gene_symbol == "Contr")
  if(length(ee)>0){
    q0<-q0[-ee,]  
  }
  
  INPUT<-read.delim(inputfile1)
  load(paste(inputfile,"EXP.Rdata",sep = ""))
  
  #### GSEA KS test should be applied here
  source("~/Sage-Analysis-Pipeline/PathwayAnalysis/myPathwayAnalysis.R")
  geneset<-as.character(as.matrix(INPUT$Gene_symbol)) # hits only 
  Geneset<-as.character(as.matrix(q0$Gene_symbol)) # all libraries
  
  # Cell Line
  REF<- -log10(EXP.Sanger$pval)
  names(REF)<-rownames(EXP.Sanger$cancer)
  REF[which(EXP.Sanger$pval==0)]<-300
  
  # hits only
  pathwayAnalysis<-myPathwayAnalysis$new()
  pathwayAnalysis$gsea(REF,geneset,np=1000,w =1)
  pathwayAnalysis$gseaResult$p.value
  png(paste(outfolder,"/EXP/GSEA_Sanger_EXP_hits.png",sep = ""),width = 800,height = 800)
  pathwayAnalysis$gseaPlot(REF,pathwayAnalysis$gseaResult$geneset)
  dev.off()
  # all libraries
  pathwayAnalysis<-myPathwayAnalysis$new()
  pathwayAnalysis$gsea(REF,Geneset,np=1000,w =1)
  pathwayAnalysis$gseaResult$p.value
  png(paste(outfolder,"/EXP/GSEA_Sanger_EXP_libraries.png",sep = ""),width = 800,height = 800)
  pathwayAnalysis$gseaPlot(REF,pathwayAnalysis$gseaResult$geneset)
  dev.off()
  
  # TCGA
  REF1<- -log10(EXP.TCGA$pval)
  names(REF1)<-rownames(EXP.TCGA$cancer)
  REF1[which(EXP.TCGA$pval==0)]<-300
  
  pathwayAnalysis1<-myPathwayAnalysis$new()
  pathwayAnalysis1$gsea(REF1,geneset,np=1000,w =1)
  pathwayAnalysis1$gseaResult$p.value
  png(paste(outfolder,"/EXP/GSEA_TCGA_EXP_hits.png",sep = ""),width = 800,height = 800)
  pathwayAnalysis1$gseaPlot(REF1,pathwayAnalysis1$gseaResult$geneset)
  dev.off()
  
  pathwayAnalysis1<-myPathwayAnalysis$new()
  pathwayAnalysis1$gsea(REF1,Geneset,np=1000,w =1)
  pathwayAnalysis1$gseaResult$p.value
  png(paste(outfolder,"/EXP/GSEA_TCGA_EXP_libraries.png",sep = ""),width = 800,height = 800)
  pathwayAnalysis1$gseaPlot(REF1,pathwayAnalysis1$gseaResult$geneset)
  dev.off()
  
  
  # Heatmap Sanger
  pos.DE<-which(EXP.Sanger$pval <=0.05)
  MAT<-cbind(EXP.Sanger$cancer,EXP.Sanger$other)
  hit.exp.sanger<-match(INPUT$Gene_symbol,rownames(EXP.Sanger$cancer))
  library.exp.sanger<-match(q0$Gene_symbol,rownames(EXP.Sanger$cancer))
  # Hit with DOE(differentially overexpressed)
  pos.hit.exp.sanger<-intersect(hit.exp.sanger[which(!is.na(hit.exp.sanger))],pos.DE)
  pos.library.exp.sanger<-library.exp.sanger[which(!is.na(library.exp.sanger))]
  
  Hit.exp.sanger<-rownames(EXP.Sanger$cancer)[pos.hit.exp.sanger]
  write.table(Hit.exp.sanger,file = paste(outfolder,"/EXP/hit_Overexpressed_Sanger.txt",sep = ""),row.names=F,col.names=F,quote=F)
  
  # because of huge matrix dimensionality, let's only concentrate libraries and hits when drawing heatmap
  MAT.library<-MAT[pos.library.exp.sanger,]
  D2<-apply(MAT.library,1,rank)
  
  pos<-match(Hit.exp.sanger,rownames(MAT.library))
  A1<-rep("green",nrow(MAT.library))
  A1[pos]<-"red"
  
  mycol <- colorpanel(n=99,low="green",mid="white",high="red")
  ColColorBar2<-rep("grey",ncol(MAT.library))
  ColColorBar2[1:ncol(EXP.Sanger$cancer)]<-"darkmagenta"
  
  png(paste(outfolder,"/EXP/Heatmap_Sanger_overexpressed.png",sep = ""),width = 1600,height = 3200)
  aaa<-heatmap.2(t(D2), col=mycol,Rowv = T,Colv=T,
                 xlab = NULL,
                 ylab = NULL,
                 scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
                 trace="none", cexRow=1, cexCol=1, mar = c(6,10),symkey=F,breaks=seq(min(D2),max(D2),length = length(mycol)+1),
                 ColSideColors = ColColorBar2,
                 RowSideColors = A1)
  dev.off()
  
  
  
  # TCGA Heatmap
  pos.DE1<-which(EXP.TCGA$pval <=0.05)
  MAT<-cbind(EXP.TCGA$cancer,EXP.TCGA$other)
  hit.exp.tcga<-match(INPUT$Gene_symbol,rownames(EXP.TCGA$cancer))
  library.exp.tcga<-match(q0$Gene_symbol,rownames(EXP.TCGA$cancer))
  # Hit with DOE(differentially overexpressed)
  pos.hit.exp.tcga<-intersect(hit.exp.tcga[which(!is.na(hit.exp.tcga))],pos.DE1)
  pos.library.exp.tcga<-library.exp.tcga[which(!is.na(library.exp.tcga))]
  
  Hit.exp.tcga<-rownames(EXP.TCGA$cancer)[pos.hit.exp.tcga]
  write.table(Hit.exp.tcga,file = paste(outfolder,"/EXP/hit_Overexpressed_TCGA.txt",sep = ""),row.names=F,col.names=F,quote=F)
  
  # because of huge matrix dimensionality, let's only concentrate libraries and hits when drawing heatmap
  MAT.library<-MAT[pos.library.exp.tcga,]
  D2<-apply(MAT.library,1,rank)
  
  pos<-match(Hit.exp.tcga,rownames(MAT.library))
  A1<-rep("green",nrow(MAT.library))
  A1[pos]<-"red"
  
  mycol <- colorpanel(n=99,low="green",mid="white",high="red")
  ColColorBar2<-rep("grey",ncol(MAT.library))
  ColColorBar2[1:ncol(EXP.TCGA$cancer)]<-"darkmagenta"
  
  png(paste(outfolder,"/EXP/Heatmap_TCGA_overexpressed.png",sep = ""),width = 1600,height = 3200)
  aaa<-heatmap.2(t(D2), col=mycol,Rowv = T,Colv=T,
                 xlab = NULL,
                 ylab = NULL,
                 scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
                 trace="none", cexRow=1, cexCol=1, mar = c(6,10),symkey=F,breaks=seq(min(D2),max(D2),length = length(mycol)+1),
                 ColSideColors = ColColorBar2,
                 RowSideColors = A1)
  dev.off()
  
}