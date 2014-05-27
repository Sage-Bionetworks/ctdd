hit_call_pathway<-function(inputfile,outfolder,norm_method){
  
  # kinome or druggable
  A1<-read.delim(inputfile)
  
  AA1<-subset(A1,A1$Group == "Library")
  if(norm_method == "pct"){
    bb1<-apply(AA1[,12:14],1,median) - 100
  }else{
    bb1<-apply(AA1[,15:17],1,median)  
  }
  
  names(bb1)<-AA1$Gene_symbol
  
  
  AAA1.1<-subset(A1,A1$Group == "Blank")
  AAA1.2<-subset(A1,A1$Group == "Mock")
  
  
  AAA1<-rbind(AAA1.1,AAA1.2)
  
  source("~/Sage-Analysis-Pipeline/PathwayAnalysis/pathwayAnalysis.R")
  
  results1.gsea.kegg<-pathwayAnalysis(synID="syn1681370",pathwayName = "KEGG",Reference = bb1,Test.method = "GSEA",cores = 14)
  results1.gsea.biocarta<-pathwayAnalysis(synID="syn1681370",pathwayName = "BIOCARTA",Reference = bb1,Test.method = "GSEA",cores = 14)
  results1.gsea.reactome<-pathwayAnalysis(synID="syn1681370",pathwayName = "REACTOME",Reference = bb1,Test.method = "GSEA",cores = 14)
  results1.gsea.go_bp<-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_BP",Reference = bb1,Test.method = "GSEA",cores = 14)
  results1.gsea.go_mf<-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_MF",Reference = bb1,Test.method = "GSEA",cores = 14)
  
  
  save(results1.gsea.kegg,results1.gsea.biocarta,results1.gsea.reactome,file = paste(outfolder,"/Canonical/GSEA.Rdata",sep =""))
  save(results1.gsea.go_bp,results1.gsea.go_mf,file = paste(outfolder,"/GO/GSEA.Rdata",sep =""))
  
  
  pvalFunction1<-function(results.fet){
    Pval.fet<-c()
    for(k in 1:length(results.fet)){
      if(results.fet[[k]]$gseaResult$nes > 0 | is.nan(results.fet[[k]]$gseaResult$nes)){
        Pval.fet<-c(Pval.fet,NA)
      }else{
        
        Pval.fet<-c(Pval.fet,results.fet[[k]]$gseaResult$p.value)
      }
    }
    names(Pval.fet)<-names(results.fet)
    #     Pval.sort<-sort(Pval.fet)
    #     return(Pval.sort)
    return(Pval.fet)
  }
  
  Pval1.gsea.kegg     <-pvalFunction1(results1.gsea.kegg)
  Pval1.gsea.biocarta <-pvalFunction1(results1.gsea.biocarta)
  Pval1.gsea.reactome <-pvalFunction1(results1.gsea.reactome)
  Pval1.gsea.go_bp    <-pvalFunction1(results1.gsea.go_bp)
  Pval1.gsea.go_mf    <-pvalFunction1(results1.gsea.go_mf)
  
  write.table(Pval1.gsea.kegg,file = paste(outfolder,"/Canonical/KEGG.txt",sep = ""),row.names=T,col.names=F,quote=F)
  write.table(Pval1.gsea.biocarta,file = paste(outfolder,"/Canonical/BIOCARTA.txt",sep = ""),row.names=T,col.names=F,quote=F)
  write.table(Pval1.gsea.reactome,file = paste(outfolder,"/Canonical/REACTOME.txt",sep = ""),row.names=T,col.names=F,quote=F)
  write.table(Pval1.gsea.go_bp,file = paste(outfolder,"/GO/GO_BP.txt",sep = ""),row.names=T,col.names=F,quote=F)
  write.table(Pval1.gsea.go_mf,file = paste(outfolder,"/GO/GO_MF.txt",sep = ""),row.names=T,col.names=F,quote=F)
  
  
  findGenes<-function(LIST,results){
    genes<-c()
    for(k in 1:length(LIST)){
      genes<-union(genes,results[[LIST[k]]]$gseaResult$geneset)
    }
    return(genes)
  }
  
  LIST1.1<-which(Pval1.gsea.kegg<=0.05)
  LIST1.2<-which(Pval1.gsea.biocarta<=0.05)
  LIST1.3<-which(Pval1.gsea.reactome<=0.05)
  LIST1.4<-which(Pval1.gsea.go_bp<=0.05)
  LIST1.5<-which(Pval1.gsea.go_mf<=0.05)
  
  
  Genes1.kegg      <-findGenes(LIST1.1,results1.gsea.kegg)
  Genes1.biocarta  <-findGenes(LIST1.2,results1.gsea.biocarta)
  Genes1.reactome  <-findGenes(LIST1.3,results1.gsea.reactome)
  Genes1.go_bp     <-findGenes(LIST1.4,results1.gsea.go_bp)
  Genes1.go_mf     <-findGenes(LIST1.5,results1.gsea.go_mf)
  
  
  source("~/Sage-Analysis-Pipeline/PathwayAnalysis/leadingEdgeFind.R")
  findLeadingEdgeGenes<-function(LIST,results,bb){
    genes<-c()
    for(k in 1:length(LIST)){
      genes<-union(genes,leadingEdgeFind(bb,results[[LIST[k]]]$gseaResult$geneset))
    }
    return(genes)
  }
  
  LEGenes1.kegg      <- findLeadingEdgeGenes(LIST1.1,results1.gsea.kegg,bb1)
  LEGenes1.biocarta  <- findLeadingEdgeGenes(LIST1.2,results1.gsea.biocarta,bb1)
  LEGenes1.reactome  <- findLeadingEdgeGenes(LIST1.3,results1.gsea.reactome,bb1)
  LEGenes1.go_bp     <- findLeadingEdgeGenes(LIST1.4,results1.gsea.go_bp,bb1)
  LEGenes1.go_mf     <- findLeadingEdgeGenes(LIST1.5,results1.gsea.go_mf,bb1)
  
  
  write.table(LEGenes1.kegg,    file = paste(outfolder,"/Canonical/LeadingEdgeGenes_KEGG.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes1.biocarta,file = paste(outfolder,"/Canonical/LeadingEdgeGenes_BIOCARTA.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes1.reactome,file = paste(outfolder,"/Canonical/LeadingEdgeGenes_REACTOME.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes1.go_bp,   file = paste(outfolder,"/GO/LeadingEdgeGenes_GO_BP.txt",sep = ""),row.names=F,col.names=F,quote=F)
  write.table(LEGenes1.go_mf,   file = paste(outfolder,"/GO/LeadingEdgeGenes_GO_MF.txt",sep = ""),row.names=F,col.names=F,quote=F)
  
  
  # genes.CANONICAL <-union(intersect(LEGenes1.kegg,LEGenes2.kegg),union(intersect(LEGenes1.biocarta,LEGenes2.biocarta),intersect(LEGenes1.reactome,LEGenes2.reactome)))
  # genes.GO <-union(intersect(LEGenes1.go_bp,LEGenes2.go_bp),intersect(LEGenes1.go_mf,LEGenes2.go_mf))
  # 
  # genes.1<-genes.CANONICAL
  # genes.2<-genes.GO
  # genes<-union(genes.CANONICAL,genes.GO)
  # 
  # # b<-match(genes,AA1$Gene_symbol)
  # 
  # 
  # # batch effect correction between vehicle and modifier screen 
  # b<-match(AAA1$Well_ID,AAA2$Well_ID)
  # qq2<-apply(AAA2[b,c(15:17)],1,median)
  # qq1<-apply(AAA1[,c(15:17)],1,median)
  # 
  # names(qq2)<-AAA2$Gene_symbol[b]
  # names(qq1)<-AAA1$Gene_symbol
  # 
  # f <- function (x, y, a) sqrt(sum((y - a*x)^2))
  # xmin <- optimize(f, c(-6, 6), tol = 0.00000001, x=qq2,y=qq1)
  # slope<-xmin$minimum
  # 
  # plot(qq2,qq1,pch = 19,cex =0.5,col = "darkgrey",xlab = "Untreated",ylab = "Treated")
  # abline(v=0,h=0,col = "blue",lty = 2)
  # abline(0,1,col = "blue",lty = 2)
  # abline(0,slope,col = "red",lty = 2)
  # 
  # 
  # a<-match(AA1$Well_ID,AA2$Well_ID)
  # q1<-apply(AA1[,c(15:17)],1,median)
  # q2<-apply(AA2[a,c(15:17)],1,median)
  # names(q1)<-AA1$Gene_symbol
  # names(q2)<-AA2$Gene_symbol[a]
  # 
  # 
  # p.CANONICAL<-match(genes.1,names(q1))
  # p.GO<-match(genes.2,names(q1))
  # p.ALL<-match(genes,names(q1))
  # p.SYNERGY<-which(slope*q2>q1)
  # 
  # pp.CANONICAL<-intersect(p.SYNERGY,p.CANONICAL)
  # pp.GO<-intersect(p.SYNERGY,p.GO)
  # pp.ALL<-intersect(p.SYNERGY,p.ALL)
  # 
  # plot(q2,q1,pch = 19,cex =0.5,col = "gold",main = "Pathway Relevant Genes")
  # points(q2[pp.CANONICAL],q1[pp.CANONICAL],pch = 3,cex = 0.5,col = "blue")
  # points(q2[pp.GO],q1[pp.GO],pch = 4,cex = 0.5,col = "blue")
  # abline(0,slope,lty = 2)
  # abline(v=0,h=0,col = "blue",lty = 2)
  # 
  # 
  # QQ<-cbind(AA1[,-c(10:12)],q1,AA2[a,c(15:17)],q2)
  # 
  # names(QQ)[13]<-c("Treated.norm.Median")
  # names(QQ)[17]<-c("Untreated.norm.Median")
  # names(QQ)[c(10:12)]<-c("Treated.norm.Replicate.1","Treated.norm.Replicate.2","Treated.norm.Replicate.3")
  # names(QQ)[c(14:16)]<-c("Untreated.norm.Replicate.1","Untreated.norm.Replicate.2","Untreated.norm.Replicate.3")
  # 
  # outfolder = "~/ScreenRNAi/HNSCC/Mock_normalization/hit_call_pathway/"
  # write.table(QQ[pp.ALL,],file = paste(outfolder, "/table_hit_All.txt",sep = ""),row.names = F,quote = F,sep = "\t")
  # write.table(QQ[pp.CANONICAL,],file = paste(outfolder, "/table_hit_Canonical.txt",sep = ""),row.names = F,quote = F,sep = "\t")
  # write.table(QQ[pp.GO,],file = paste(outfolder, "/table_hit_Go.txt",sep = ""),row.names = F,quote = F,sep = "\t")
  # 
  # 
  # AA<-cbind(AA1[,c(15:17)],AA2[a,c(15:17)])
  # colnames(AA)<-c("modifier_1","modifier_2","modifier_3","vehicle_1","vehicle_2","vehicle_3")
  # rownames(AA)<-AA1$Gene_symbol
  # 
  # D<-apply(AA,2,rank)
  # D2<-apply(AA[b,],2,rank)
  # rownames(D2)<-genes
  # mycol <- colorpanel(n=99,low="red",mid="white",high="green")
  # 
  # png("~/ScreenRNAi/HNSCC/Mock_normalization/hit_call_pathway/hits_heatmap_pathways.png",width = 640,height = 1960)
  # aaa<-heatmap.2(D2, col=mycol,Rowv = T,Colv=F,
  #                scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
  #                trace="none", cexRow=0.5, cexCol=1, mar = c(6,10),symkey=F,breaks=seq(min(D2),max(D2),length = length(mycol)+1))
  # dev.off()
  # 
  # plot(aaa$rowDendrogram)
}