inputfile = "~/DATA/siRNA_Screens/insock_processed/HNSCC_druggable/CellLines_E_uM_compound_0.txt"
inputfile1 = "~/DATA/siRNA_Screens/insock_processed/HNSCC_druggable/CellLines_E_uM_compound_0.249.txt"

outfolder = "~/RNAi_Analysis/HNSCC/z_transformation/Untreated/Patient/Druggable/Hit_Call/Pathway_Relevant/"
outfolder1 = "~/RNAi_Analysis/HNSCC/z_transformation/Treated/Patient/Druggable/Hit_Call/Pathway_Relevant/"

outfolder = "~/RNAi_Analysis/HNSCC/pct_normalization_control/Untreated/Patient/Druggable/Hit_Call/Pathway_Relevant/"
outfolder1 = "~/RNAi_Analysis/HNSCC/pct_normalization_control//Treated/Patient/Druggable/Hit_Call/Pathway_Relevant/"

hit_call_modifier_pathway<-function(inputfile,inputfile1,outfolder,outfolder1,norm_method){
  
  # kinome or druggable
  A1<-read.delim(inputfile)
  A2<-read.delim(inputfile1)
  
  AA1<-subset(A1,A1$Group == "Library")
  AA2<-subset(A2,A2$Group == "Library")
  
  AAA1.1<-subset(A1,A1$Group == "Blank")
  AAA1.2<-subset(A1,A1$Group == "Mock")
  
  AAA2.1<-subset(A2,A2$Group == "Blank")
  AAA2.2<-subset(A2,A2$Group == "Mock")
  
  AAA1<-rbind(AAA1.1,AAA1.2)
  AAA2<-rbind(AAA2.1,AAA2.2)
  
  b<-match(AA1$Well_ID,AA2$Well_ID)
    
  if(norm_method == "pct"){
    bb1<-apply(AA1[,c(12:14)],1,median) -100
    bb2<-apply(AA2[,c(12:14)],1,median) -100
        
  }else{
    bb1<-apply(AA1[,15:17],1,median)  
    bb2<-apply(AA2[,15:17],1,median)  
  }
  
  names(bb1)<-AA1$Gene_symbol
  names(bb2)<-AA2$Gene_symbol
  
  
  # Untreated
  a1<-read.table(file = paste(outfolder,"/Canonical/KEGG.txt",sep = ""))
  a2<-read.table(file = paste(outfolder,"/Canonical/BIOCARTA.txt",sep = ""))
  a3<-read.table(file = paste(outfolder,"/Canonical/REACTOME.txt",sep = ""))
  a4<-read.table(file = paste(outfolder,"/GO/GO_BP.txt",sep = ""))
  a5<-read.table(file = paste(outfolder,"/GO/GO_MF.txt",sep = ""))    
  # Treated
  b1<-read.table(file = paste(outfolder1,"/Canonical/KEGG.txt",sep = ""))
  b2<-read.table(file = paste(outfolder1,"/Canonical/BIOCARTA.txt",sep = ""))
  b3<-read.table(file = paste(outfolder1,"/Canonical/REACTOME.txt",sep = ""))
  b4<-read.table(file = paste(outfolder1,"/GO/GO_BP.txt",sep = ""))
  b5<-read.table(file = paste(outfolder1,"/GO/GO_MF.txt",sep = ""))
  
  commonLIST.1<-which(a1[,2]<=0.05 & b1[,2]<=0.05)
  commonLIST.2<-which(a2[,2]<=0.05 & b2[,2]<=0.05)
  commonLIST.3<-which(a3[,2]<=0.05 & b3[,2]<=0.05)
  commonLIST.4<-which(a4[,2]<=0.05 & b4[,2]<=0.05)
  commonLIST.5<-which(a5[,2]<=0.05 & b5[,2]<=0.05)
  
  length(which(a1[,2]<=0.05))
  length(which(a2[,2]<=0.05))
  length(which(a3[,2]<=0.05))
  length(which(a4[,2]<=0.05))
  length(which(a5[,2]<=0.05))
  
  length(which(b1[,2]<=0.05))
  length(which(b2[,2]<=0.05))
  length(which(b3[,2]<=0.05))
  length(which(b4[,2]<=0.05))
  length(which(b5[,2]<=0.05))
  source("~/Sage-Analysis-Pipeline/PathwayAnalysis/leadingEdgeFind.R")
  findLeadingEdgeGenes<-function(LIST,results,bb){
    genes<-c()
    for(k in 1:length(LIST)){
      genes<-union(genes,leadingEdgeFind(bb,results[[LIST[k]]]$gseaResult$geneset))
    }
    return(genes)
  }
  
  load(paste(outfolder,"/Canonical/GSEA.Rdata",sep =""))
  load(paste(outfolder,"/GO/GSEA.Rdata",sep =""))
  
  LEGenes1.kegg      <- findLeadingEdgeGenes(commonLIST.1,results1.gsea.kegg,bb1)
  LEGenes1.biocarta  <- findLeadingEdgeGenes(commonLIST.2,results1.gsea.biocarta,bb1)
  LEGenes1.reactome  <- findLeadingEdgeGenes(commonLIST.3,results1.gsea.reactome,bb1)
  LEGenes1.go_bp     <- findLeadingEdgeGenes(commonLIST.4,results1.gsea.go_bp,bb1)
  LEGenes1.go_mf     <- findLeadingEdgeGenes(commonLIST.5,results1.gsea.go_mf,bb1)
  
  load(paste(outfolder1,"/Canonical/GSEA.Rdata",sep =""))
  load(paste(outfolder1,"/GO/GSEA.Rdata",sep =""))
  
  LEGenes2.kegg      <- findLeadingEdgeGenes(commonLIST.1,results1.gsea.kegg,bb2)
  LEGenes2.biocarta  <- findLeadingEdgeGenes(commonLIST.2,results1.gsea.biocarta,bb2)
  LEGenes2.reactome  <- findLeadingEdgeGenes(commonLIST.3,results1.gsea.reactome,bb2)
  LEGenes2.go_bp     <- findLeadingEdgeGenes(commonLIST.4,results1.gsea.go_bp,bb2)
  LEGenes2.go_mf     <- findLeadingEdgeGenes(commonLIST.5,results1.gsea.go_mf,bb2)
  
  length(intersect(LEGenes1.kegg,LEGenes2.kegg))
  length(intersect(LEGenes1.biocarta,LEGenes2.biocarta))
  length(intersect(LEGenes1.reactome,LEGenes2.reactome))
  length(intersect(LEGenes1.go_bp,LEGenes2.go_bp))
  length(intersect(LEGenes1.go_mf,LEGenes2.go_mf))
  
  write.table(LEGenes1.kegg,    file =paste(outfolder,"Canonical/commonLeadingEdgeGenes_Modifier_KEGG.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes1.biocarta,file =paste(outfolder,"Canonical/commonLeadingEdgeGenes_Modifier_BIOCARTA.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes1.reactome,file =paste(outfolder,"Canonical/commonLeadingEdgeGenes_Modifier_REACTOME.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes1.go_bp,   file =paste(outfolder,"GO/commonLeadingEdgeGenes_Modifier_GO_BP.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes1.go_mf,   file =paste(outfolder,"GO/commonLeadingEdgeGenes_Modifier_GO_MF.txt",sep = ""),row.names=F,col.names=F,quote=F)
  
  write.table(LEGenes2.kegg,    file =paste(outfolder1,"Canonical/commonLeadingEdgeGenes_Modifier_KEGG.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes2.biocarta,file =paste(outfolder1,"Canonical/commonLeadingEdgeGenes_Modifier_BIOCARTA.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes2.reactome,file =paste(outfolder1,"Canonical/commonLeadingEdgeGenes_Modifier_REACTOME.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes2.go_bp,   file =paste(outfolder1,"GO/commonLeadingEdgeGenes_Modifier_GO_BP.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes2.go_mf,   file =paste(outfolder1,"GO/commonLeadingEdgeGenes_Modifier_GO_MF.txt",sep = ""),row.names=F,col.names=F,quote=F)
  
  if(norm_method != "pct"){
    write.table(intersect(LEGenes1.kegg,LEGenes2.kegg),        file ="~/RNAi_Analysis/HNSCC/z_transformation/commonLeadingEdgeGenes_Both_KEGG.txt",row.names=F,col.names=F,quote=F)
    write.table(intersect(LEGenes1.biocarta,LEGenes2.biocarta),file ="~/RNAi_Analysis/HNSCC/z_transformation/commonLeadingEdgeGenes_Both_BIOCARTA.txt",row.names=F,col.names=F,quote=F)
    write.table(intersect(LEGenes1.reactome,LEGenes2.reactome),file ="~/RNAi_Analysis/HNSCC/z_transformation/commonLeadingEdgeGenes_Both_REACTOME.txt",row.names=F,col.names=F,quote=F)
    write.table(intersect(LEGenes1.go_bp,LEGenes2.go_bp),      file ="~/RNAi_Analysis/HNSCC/z_transformation/commonLeadingEdgeGenes_Both_GO_BP.txt",row.names=F,col.names=F,quote=F)
    write.table(intersect(LEGenes1.go_mf,LEGenes2.go_mf),      file ="~/RNAi_Analysis/HNSCC/z_transformation/commonLeadingEdgeGenes_Both_GO_MF.txt",row.names=F,col.names=F,quote=F)
    
    genes.CANONICAL <-union(intersect(LEGenes1.kegg,LEGenes2.kegg),union(intersect(LEGenes1.biocarta,LEGenes2.biocarta),intersect(LEGenes1.reactome,LEGenes2.reactome)))
    genes.GO <-union(intersect(LEGenes1.go_bp,LEGenes2.go_bp),intersect(LEGenes1.go_mf,LEGenes2.go_mf))
    
    write.table(genes.CANONICAL,        file ="~/RNAi_Analysis/HNSCC/z_transformation/commonLeadingEdgeGenes_Canonical.txt",row.names=F,col.names=F,quote=F)
    write.table(genes.GO,               file ="~/RNAi_Analysis/HNSCC/z_transformation/commonLeadingEdgeGenes_GO.txt",row.names=F,col.names=F,quote=F)
    
    bb<-match(AAA1$Well_ID,AAA2$Well_ID)
    qq1<-apply(AAA1[,c(15:17)],1,median)
    qq2<-apply(AAA2[bb,c(15:17)],1,median)
    names(qq1)<-AAA1$Gene_symbol
    names(qq2)<-AAA2$Gene_symbol[bb]    
    f <- function (x, y, a) sqrt(sum((y - a*x)^2))
    xmin <- optimize(f, c(-6, 6), tol = 0.00000001, x=qq1,y=qq2)
    slope<-xmin$minimum    
    
    plot(qq1,qq2,pch = 19,cex =0.5,col = "darkgrey",xlab = "Untreated",ylab = "Treated")
    abline(v=0,h=0,col = "blue",lty = 2)
    abline(0,1,col = "blue",lty = 2)
    abline(0,slope,col = "red",lty = 2)
    
    
    p.CANONICAL<-match(genes.CANONICAL,names(bb1))
    p.GO<-match(genes.GO,names(bb1))
    p.SYNERGY<-which(slope*bb1>bb2[b])
    
    pp.CANONICAL<-intersect(p.SYNERGY,p.CANONICAL)
    pp.GO<-intersect(p.SYNERGY,p.GO)
    
    
    write.table(names(bb1)[pp.CANONICAL],        file ="~/RNAi_Analysis/HNSCC/z_transformation/sensitizer_commonLeadingEdgeGenes_Canonical.txt",row.names=F,col.names=F,quote=F)
    write.table(names(bb1)[pp.GO],               file ="~/RNAi_Analysis/HNSCC/z_transformation/sensitizer_commonLeadingEdgeGenes_GO.txt",row.names=F,col.names=F,quote=F)
    
    png("~/RNAi_Analysis/HNSCC/z_transformation/sensitizer_commonLeadingEdgeGenes.png",width = 1280,height = 640)
    par(mfrow = c(1,2))
    par(pty="s")
    plot(bb1,bb2[b],pch = 19,cex =0.5,col = "darkgrey",main = "Pathway Relevant Sensitizer Hits", xlab = "Untreated (z-score)", ylab = "Treated (z-score)")
    points(bb1[pp.CANONICAL],bb2[b[pp.CANONICAL]],pch = 19,cex = 0.5,col = "orange")
    abline(0,slope,lty = 2)
    legend("topleft",legend = c("libraries","Hits with canonical pathway"),pch = 19,col = c("darkgrey","orange"),bty = "n")
    
    par(pty="s")
    plot(bb1,bb2[b],pch = 19,cex =0.5,col = "darkgrey",main = "Pathway Relevant Sensitizer Hits", xlab = "Untreated (z-score)", ylab = "Treated (z-score)")
    points(bb1[pp.GO],bb2[b[pp.GO]],pch = 19,cex = 0.5,col = "orange")
    abline(0,slope,lty = 2)
    legend("topleft",legend = c("libraries","Hits with GO-term"),pch = 19,col = c("darkgrey","orange"),bty = "n")
    dev.off()
    
  }else{
    write.table(intersect(LEGenes1.kegg,LEGenes2.kegg),        file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/commonLeadingEdgeGenes_Both_KEGG.txt",row.names=F,col.names=F,quote=F)
    write.table(intersect(LEGenes1.biocarta,LEGenes2.biocarta),file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/commonLeadingEdgeGenes_Both_BIOCARTA.txt",row.names=F,col.names=F,quote=F)
    write.table(intersect(LEGenes1.reactome,LEGenes2.reactome),file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/commonLeadingEdgeGenes_Both_REACTOME.txt",row.names=F,col.names=F,quote=F)
    write.table(intersect(LEGenes1.go_bp,LEGenes2.go_bp),      file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/commonLeadingEdgeGenes_Both_GO_BP.txt",row.names=F,col.names=F,quote=F)
    write.table(intersect(LEGenes1.go_mf,LEGenes2.go_mf),      file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/commonLeadingEdgeGenes_Both_GO_MF.txt",row.names=F,col.names=F,quote=F)
    
    genes.CANONICAL <-union(intersect(LEGenes1.kegg,LEGenes2.kegg),union(intersect(LEGenes1.biocarta,LEGenes2.biocarta),intersect(LEGenes1.reactome,LEGenes2.reactome)))
    genes.GO <-union(intersect(LEGenes1.go_bp,LEGenes2.go_bp),intersect(LEGenes1.go_mf,LEGenes2.go_mf))
    
    write.table(genes.CANONICAL,        file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/commonLeadingEdgeGenes_Canonical.txt",row.names=F,col.names=F,quote=F)
    write.table(genes.GO,               file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/commonLeadingEdgeGenes_GO.txt",row.names=F,col.names=F,quote=F)
    
    slope<-1
    
    
    p.CANONICAL<-match(genes.CANONICAL,names(bb1))
    p.GO<-match(genes.GO,names(bb1))
    p.SYNERGY<-which(slope*bb1>bb2[b])
    
    pp.CANONICAL<-intersect(p.SYNERGY,p.CANONICAL)
    pp.GO<-intersect(p.SYNERGY,p.GO)
    
    
    write.table(names(bb1)[pp.CANONICAL],        file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/sensitizer_commonLeadingEdgeGenes_Canonical.txt",row.names=F,col.names=F,quote=F)
    write.table(names(bb1)[pp.GO],               file ="~/RNAi_Analysis/HNSCC/pct_normalization_control/sensitizer_commonLeadingEdgeGenes_GO.txt",row.names=F,col.names=F,quote=F)
    
    png("~/RNAi_Analysis/HNSCC/pct_normalization_control/sensitizer_commonLeadingEdgeGenes.png",width = 1280,height = 640)
    par(mfrow = c(1,2))
    par(pty="s")
    plot(bb1,bb2[b],pch = 19,cex =0.5,col = "darkgrey",main = "Pathway Relevant Sensitizer Hits", xlab = "Untreated (z-score)", ylab = "Treated (z-score)")
    points(bb1[pp.CANONICAL],bb2[b[pp.CANONICAL]],pch = 19,cex = 0.5,col = "orange")
    abline(0,slope,lty = 2)
    legend("topleft",legend = c("libraries","Hits with canonical pathway"),pch = 19,col = c("darkgrey","orange"),bty = "n")
    
    par(pty="s")
    plot(bb1,bb2[b],pch = 19,cex =0.5,col = "darkgrey",main = "Pathway Relevant Sensitizer Hits", xlab = "Untreated (z-score)", ylab = "Treated (z-score)")
    points(bb1[pp.GO],bb2[b[pp.GO]],pch = 19,cex = 0.5,col = "orange")
    abline(0,slope,lty = 2)
    legend("topleft",legend = c("libraries","Hits with GO-term"),pch = 19,col = c("darkgrey","orange"),bty = "n")
    dev.off()
    
  }
}