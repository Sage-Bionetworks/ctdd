A<-read.delim("~/HNSCC//screen.txt",sep = "",header = F)
folds<-as.numeric(A$V2)
names(folds)<-as.character(as.matrix(A$V1))
pvals<-as.numeric(A$V3)
names(pvals)<-as.character(as.matrix(A$V1))

KK<-(intersect(union(which(folds>=4),which(folds<=-4)),which(pvals<=0.05)))
source("~/Sage-Analysis-Pipeline/PathwayAnalysis/pathwayAnalysis.R")

# preparation : myReference list for GSEA
myReference<-abs(folds)

# preparation : my test set for FET(fisher exact test)
B<-names(folds)[KK]

# MSigDB (synID="syn2135029" for GraphiteDB)
# run test with multicore with FET (fisher exact test)
results.fet.kegg     <-pathwayAnalysis(synID="syn1681370",pathwayName = "KEGG",Reference = B,Test.method = "FET",cores = 8)
results.fet.biocarta <-pathwayAnalysis(synID="syn1681370",pathwayName = "Biocarta",Reference = B,Test.method = "FET",cores = 8)
results.fet.reactome <-pathwayAnalysis(synID="syn1681370",pathwayName = "Reactome",Reference = B,Test.method = "FET",cores = 8)
results.fet.go_bp    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_BP",Reference = B,Test.method = "FET",cores = 8)
results.fet.go_mf    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_MF",Reference = B,Test.method = "FET",cores = 8)


# combine all p-value from FET
summaryFET<-function(results.fet){
  Pval.fet<-c()
  for(k in 1:length(results.fet)){
    Pval.fet<-c(Pval.fet,results.fet[[k]]$fetResult$p.value)
  }
  names(Pval.fet)<-names(results.fet)
  return(Pval.fet)
}

fet.kegg<-summaryFET(results.fet.kegg)
fet.biocarta<-summaryFET(results.fet.biocarta)
fet.reactome<-summaryFET(results.fet.reactome)
fet.go_bp<-summaryFET(results.fet.go_bp)
fet.go_mf<-summaryFET(results.fet.go_mf)

# run test with multicore with GSEA
results.gsea.kegg     <-pathwayAnalysis(synID="syn1681370",pathwayName = "KEGG",Reference = myReference,Test.method = "GSEA",cores = 8)
results.gsea.biocarta <-pathwayAnalysis(synID="syn1681370",pathwayName = "Biocarta",Reference = myReference,Test.method = "GSEA",cores = 8)
results.gsea.reactome <-pathwayAnalysis(synID="syn1681370",pathwayName = "Reactome",Reference = myReference,Test.method = "GSEA",cores = 8)
results.gsea.go_bp    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_BP",Reference = myReference,Test.method = "GSEA",cores = 8)
results.gsea.go_mf    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_MF",Reference = myReference,Test.method = "GSEA",cores = 8)



summaryGSEA<-function(results.fet){
  Pval.fet<-c()
  NES.fet<-c()
  
  for(k in 1:length(results.fet)){
    Pval.fet<-c(Pval.fet,results.fet[[k]]$gseaResult$p.value)
    NES.fet<-c(NES.fet,results.fet[[k]]$gseaResult$nes)
  }
  names(Pval.fet)<-names(results.fet)
  names(NES.fet)<-names(results.fet)
  
  return(cbind(Pval.fet,NES.fet))
}



gsea.kegg<-summaryGSEA(results.gsea.kegg)
gsea.biocarta<-summaryGSEA(results.gsea.biocarta)
gsea.reactome<-summaryGSEA(results.gsea.reactome)
gsea.go_bp<-summaryGSEA(results.gsea.go_bp)
gsea.go_mf<-summaryGSEA(results.gsea.go_mf)


intersect(intersect(which(gsea.kegg[,1]<=0.05),which(gsea.kegg[,2]>0)),which(fet.kegg<=0.05))

mat.make<-function(fet,gsea,threshold = 0.05){
  mat<-matrix(0,nrow = length(fet),ncol = 2)
  gsea[which(gsea[,2]<0),1]<-1
  mat[,1]<-gsea[,1]
  mat[,2]<-fet
  mat[which(mat[,1]>=threshold),1]<-1
  mat[which(mat[,2]>=threshold),2]<-1
  
  rownames(mat)<-rownames(gsea)
  Mat<-mat[which(apply(mat,1,sum)!=2),]
  X<-sort(Mat[,1],index.return=T)
  
  return(Mat[X$ix,])  
}

mat.kegg<-mat.make(fet.kegg,gsea.kegg,0.05)
mat.biocarta<-mat.make(fet.biocarta,gsea.biocarta,0.05)
mat.reactome<-mat.make(fet.reactome,gsea.reactome,0.05)
mat.go_bp<-mat.make(fet.go_bp,gsea.go_bp,0.05)
mat.go_mf<-mat.make(fet.go_mf,gsea.go_mf,0.05)

write.table(mat.kegg,file = "~/HNSCC/pathway_analysis/cisplatin_modifier//kegg.txt",row.names = T,col.names = F, sep = "\t",quote=F)
write.table(mat.biocarta,file = "~/HNSCC/pathway_analysis/cisplatin_modifier/biocarta.txt",row.names = T,col.names = F, sep = "\t",quote=F)
write.table(mat.reactome,file = "~/HNSCC/pathway_analysis/cisplatin_modifier/reactome.txt",row.names = T,col.names = F, sep = "\t",quote=F)
write.table(mat.go_bp,file = "~/HNSCC/pathway_analysis/cisplatin_modifier/go_bp.txt",row.names = T,col.names = F, sep = "\t",quote=F)
write.table(mat.go_mf,file = "~/HNSCC/pathway_analysis/cisplatin_modifier/go_mf.txt",row.names = T,col.names = F, sep = "\t",quote=F)








# combine all p-value from FET
summaryFET<-function(results.fet){
  Pval.fet<-c()
  for(k in 1:length(results.fet)){
    Pval.fet<-c(Pval.fet,results.fet[[k]]$fetResult$p.value)
  }
  names(Pval.fet)<-names(results.fet)
  return(Pval.fet)
}
summaryGSEA<-function(results.fet){
  Pval.fet<-c()
  NES.fet<-c()
  
  for(k in 1:length(results.fet)){
    Pval.fet<-c(Pval.fet,results.fet[[k]]$gseaResult$p.value)
    NES.fet<-c(NES.fet,results.fet[[k]]$gseaResult$nes)
  }
  names(Pval.fet)<-names(results.fet)
  names(NES.fet)<-names(results.fet)
  
  return(cbind(Pval.fet,NES.fet))
}



gsea.kegg<-summaryGSEA(results.gsea.kegg)
gsea.biocarta<-summaryGSEA(results.gsea.biocarta)
gsea.reactome<-summaryGSEA(results.gsea.reactome)
gsea.go_bp<-summaryGSEA(results.gsea.go_bp)
gsea.go_mf<-summaryGSEA(results.gsea.go_mf)
gsea.go_cc<-summaryGSEA(results.gsea.go_cc)

fet.kegg<-summaryFET(results.fet.kegg)
fet.biocarta<-summaryFET(results.fet.biocarta)
fet.reactome<-summaryFET(results.fet.reactome)
fet.go_bp<-summaryFET(results.fet.go_bp)
fet.go_mf<-summaryFET(results.fet.go_mf)
fet.go_cc<-summaryFET(results.fet.go_cc)


dim(gsea.kegg)
dim(gsea.biocarta)
dim(gsea.reactome)
dim(gsea.go_bp)
dim(gsea.go_mf)
dim(gsea.go_cc)

length(intersect(which(gsea.kegg[,1]<=0.05),which(gsea.kegg[,2]>=0)))
length(intersect(which(gsea.biocarta[,1]<=0.05),which(gsea.biocarta[,2]>=0)))
length(intersect(which(gsea.reactome[,1]<=0.05),which(gsea.reactome[,2]>=0)))
length(intersect(which(gsea.go_bp[,1]<=0.05),which(gsea.go_bp[,2]>=0)))
length(intersect(which(gsea.go_mf[,1]<=0.05),which(gsea.go_mf[,2]>=0)))
length(intersect(which(gsea.go_cc[,1]<=0.05),which(gsea.go_cc[,2]>=0)))



length(fet.kegg)
length(fet.biocarta)
length(fet.reactome)
length(fet.go_bp)
length(fet.go_mf)
length(fet.go_cc)

length(which(fet.kegg<=0.05))
length(which(fet.biocarta<=0.05))
length(which(fet.reactome<=0.05))
length(which(fet.go_bp<=0.05))
length(which(fet.go_mf<=0.05))
length(which(fet.go_cc<=0.05))

# difference between gsea vs. fet
length(intersect(intersect(which(gsea.kegg[,1]<=0.05),which(gsea.kegg[,2]>=0)),which(fet.kegg<=0.05)))
length(intersect(intersect(which(gsea.biocarta[,1]<=0.05),which(gsea.biocarta[,2]>=0)),which(fet.biocarta<=0.05)))
length(intersect(intersect(which(gsea.reactome[,1]<=0.05),which(gsea.reactome[,2]>=0)),which(fet.reactome<=0.05)))
length(intersect(intersect(which(gsea.go_bp[,1]<=0.05),which(gsea.go_bp[,2]>=0)),which(fet.go_bp<=0.05)))
length(intersect(intersect(which(gsea.go_mf[,1]<=0.05),which(gsea.go_mf[,2]>=0)),which(fet.go_mf<=0.05)))
length(intersect(intersect(which(gsea.go_cc[,1]<=0.05),which(gsea.go_cc[,2]>=0)),which(fet.go_cc<=0.05)))



###################################################
# preparation : myReference list for GSEA
myReference<-(folds)

# preparation : my test set for FET(fisher exact test)
B1<-names(sort(myReference,decreasing = T))[1:100]
B2<-names(sort(myReference,decreasing = F))[1:100]

# MSigDB (synID="syn2135029" for GraphiteDB)
# run test with multicore with FET (fisher exact test)
results1.fet.kegg     <-pathwayAnalysis(synID="syn1681370",pathwayName = "KEGG",Reference = B1,Test.method = "FET",cores = 8)
results1.fet.biocarta <-pathwayAnalysis(synID="syn1681370",pathwayName = "Biocarta",Reference = B1,Test.method = "FET",cores = 8)
results1.fet.reactome <-pathwayAnalysis(synID="syn1681370",pathwayName = "Reactome",Reference = B1,Test.method = "FET",cores = 8)
results1.fet.go_bp    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_BP",Reference = B1,Test.method = "FET",cores = 8)
results1.fet.go_mf    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_MF",Reference = B1,Test.method = "FET",cores = 8)
results1.fet.go_cc    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_CC",Reference = B1,Test.method = "FET",cores = 8)

results2.fet.kegg     <-pathwayAnalysis(synID="syn1681370",pathwayName = "KEGG",Reference = B2,Test.method = "FET",cores = 8)
results2.fet.biocarta <-pathwayAnalysis(synID="syn1681370",pathwayName = "Biocarta",Reference = B2,Test.method = "FET",cores = 8)
results2.fet.reactome <-pathwayAnalysis(synID="syn1681370",pathwayName = "Reactome",Reference = B2,Test.method = "FET",cores = 8)
results2.fet.go_bp    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_BP",Reference = B2,Test.method = "FET",cores = 8)
results2.fet.go_mf    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_MF",Reference = B2,Test.method = "FET",cores = 8)
results2.fet.go_cc    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_CC",Reference = B2,Test.method = "FET",cores = 8)

# run test with multicore with GSEA
result.gsea.kegg     <-pathwayAnalysis(synID="syn1681370",pathwayName = "KEGG",Reference = myReference,Test.method = "GSEA",cores = 8)
result.gsea.biocarta <-pathwayAnalysis(synID="syn1681370",pathwayName = "Biocarta",Reference = myReference,Test.method = "GSEA",cores = 8)
result.gsea.reactome <-pathwayAnalysis(synID="syn1681370",pathwayName = "Reactome",Reference = myReference,Test.method = "GSEA",cores = 8)
result.gsea.go_bp    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_BP",Reference = myReference,Test.method = "GSEA",cores = 8)
result.gsea.go_mf    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_MF",Reference = myReference,Test.method = "GSEA",cores = 8)
result.gsea.go_cc    <-pathwayAnalysis(synID="syn1681370",pathwayName = "GO_CC",Reference = myReference,Test.method = "GSEA",cores = 8)

gsea1.kegg     <-summaryGSEA(result.gsea.kegg)
gsea1.biocarta <-summaryGSEA(result.gsea.biocarta)
gsea1.reactome <-summaryGSEA(result.gsea.reactome)
gsea1.go_bp    <-summaryGSEA(result.gsea.go_bp)
gsea1.go_mf    <-summaryGSEA(result.gsea.go_mf)
gsea1.go_cc    <-summaryGSEA(result.gsea.go_cc)

fet1.kegg      <-summaryFET(results1.fet.kegg)
fet1.biocarta  <-summaryFET(results1.fet.biocarta)
fet1.reactome  <-summaryFET(results1.fet.reactome)
fet1.go_bp     <-summaryFET(results1.fet.go_bp)
fet1.go_mf     <-summaryFET(results1.fet.go_mf)
fet1.go_cc     <-summaryFET(results1.fet.go_cc)

fet2.kegg      <-summaryFET(results2.fet.kegg)
fet2.biocarta  <-summaryFET(results2.fet.biocarta)
fet2.reactome  <-summaryFET(results2.fet.reactome)
fet2.go_bp     <-summaryFET(results2.fet.go_bp)
fet2.go_mf     <-summaryFET(results2.fet.go_mf)
fet2.go_cc     <-summaryFET(results2.fet.go_cc)

# gsea pos
length(intersect(which(gsea1.kegg[,1]<=0.05),which(gsea1.kegg[,2]>=0)))
length(intersect(which(gsea1.biocarta[,1]<=0.05),which(gsea1.biocarta[,2]>=0)))
length(intersect(which(gsea1.reactome[,1]<=0.05),which(gsea1.reactome[,2]>=0)))
length(intersect(which(gsea1.go_bp[,1]<=0.05),which(gsea1.go_bp[,2]>=0)))
length(intersect(which(gsea1.go_mf[,1]<=0.05),which(gsea1.go_mf[,2]>=0)))
length(intersect(which(gsea1.go_cc[,1]<=0.05),which(gsea1.go_cc[,2]>=0)))
#gsea neg
length(intersect(which(gsea1.kegg[,1]<=0.05),which(gsea1.kegg[,2]<0)))
length(intersect(which(gsea1.biocarta[,1]<=0.05),which(gsea1.biocarta[,2]<0)))
length(intersect(which(gsea1.reactome[,1]<=0.05),which(gsea1.reactome[,2]<0)))
length(intersect(which(gsea1.go_bp[,1]<=0.05),which(gsea1.go_bp[,2]<0)))
length(intersect(which(gsea1.go_mf[,1]<=0.05),which(gsea1.go_mf[,2]<0)))
length(intersect(which(gsea1.go_cc[,1]<=0.05),which(gsea1.go_cc[,2]<0)))



length(fet.kegg)
length(fet.biocarta)
length(fet.reactome)
length(fet.go_bp)
length(fet.go_mf)
length(fet.go_cc)

length(which(fet1.kegg<=0.05))
length(which(fet1.biocarta<=0.05))
length(which(fet1.reactome<=0.05))
length(which(fet1.go_bp<=0.05))
length(which(fet1.go_mf<=0.05))
length(which(fet1.go_cc<=0.05))

length(which(fet2.kegg<=0.05))
length(which(fet2.biocarta<=0.05))
length(which(fet2.reactome<=0.05))
length(which(fet2.go_bp<=0.05))
length(which(fet2.go_mf<=0.05))
length(which(fet2.go_cc<=0.05))

# difference between gsea vs. fet
length(intersect(intersect(which(gsea1.kegg[,1]<=0.05),which(gsea1.kegg[,2]>=0)),which(fet1.kegg<=0.05)))
length(intersect(intersect(which(gsea1.biocarta[,1]<=0.05),which(gsea1.biocarta[,2]>=0)),which(fet1.biocarta<=0.05)))
length(intersect(intersect(which(gsea1.reactome[,1]<=0.05),which(gsea1.reactome[,2]>=0)),which(fet1.reactome<=0.05)))
length(intersect(intersect(which(gsea1.go_bp[,1]<=0.05),which(gsea1.go_bp[,2]>=0)),which(fet1.go_bp<=0.05)))
length(intersect(intersect(which(gsea1.go_mf[,1]<=0.05),which(gsea1.go_mf[,2]>=0)),which(fet1.go_mf<=0.05)))
length(intersect(intersect(which(gsea1.go_cc[,1]<=0.05),which(gsea1.go_cc[,2]>=0)),which(fet1.go_cc<=0.05)))

length(intersect(intersect(which(gsea1.kegg[,1]<=0.05),which(gsea1.kegg[,2]<0)),which(fet2.kegg<=0.05)))
length(intersect(intersect(which(gsea1.biocarta[,1]<=0.05),which(gsea1.biocarta[,2]<0)),which(fet2.biocarta<=0.05)))
length(intersect(intersect(which(gsea1.reactome[,1]<=0.05),which(gsea1.reactome[,2]<0)),which(fet2.reactome<=0.05)))
length(intersect(intersect(which(gsea1.go_bp[,1]<=0.05),which(gsea1.go_bp[,2]<0)),which(fet2.go_bp<=0.05)))
length(intersect(intersect(which(gsea1.go_mf[,1]<=0.05),which(gsea1.go_mf[,2]<0)),which(fet2.go_mf<=0.05)))
length(intersect(intersect(which(gsea1.go_cc[,1]<=0.05),which(gsea1.go_cc[,2]<0)),which(fet2.go_cc<=0.05)))










a1<-which(Pval<=0.05)
b1<-which(ES<0)
b2<-which(ES>0)
d1<-intersect(a1,b1) # cisplatin suppressor
d2<-intersect(a1,b2) # cisplatin enhancer
# Do you want to plot all Enrichment score? then do following:
Pval.suppress<-sort(Pval[b1])
Pval.enhance<-sort(Pval[b2])
write.table(Pval.suppress,file = paste("~/Documents/Data/HNSCC/",pathwayName,"_Cisplatin_Suppressor.txt",sep = ""),row.names = TRUE, col.names = FALSE, quote = F,)
write.table(Pval.enhance,file = paste("~/Documents/Data/HNSCC/",pathwayName,"_Cisplatin_Enhancer.txt",sep = ""),row.names = TRUE, col.names = FALSE, quote = F,)



D1<-d1[which.min(Pval[d1])]
D2<-d2[which.min(Pval[d2])]
results[[D1]]$gseaPlot(myReference,results[[D1]]$gseaResult$geneset)
title(names(results)[D1])
results[[D2]]$gseaPlot(myReference,results[[D2]]$gseaResult$geneset)
title(names(results)[D2])
# word cloud form

pathwayDB<-function(pathwayName){
  pathwayName = toupper(pathwayName)
  require(wordcloud)
  MSIGDB<-synGet("syn1681370",load = TRUE)
  
  if(is.element(pathwayName,"BIOCARTA")){
    allPathways <- MSIGDB$objects$C2.CP.BIOCARTA
  }
  if(is.element(pathwayName,"KEGG")){
    allPathways <- MSIGDB$objects$C2.CP.KEGG
  }
  if(is.element(pathwayName,"REACTOME")){
    allPathways <- MSIGDB$objects$C2.CP.REACTOME
  }
  if(is.element(pathwayName,"GO_BP")){
    allPathways <- MSIGDB$objects$C5.GO_BP
  }
  if(is.element(pathwayName,"GO_CC")){
    allPathways <- MSIGDB$objects$C5.GO_CC
  }
  if(is.element(pathwayName,"GO_MF")){
    allPathways <- MSIGDB$objects$C5.GO_MF
  }
  return(allPathways)
}

allPathways<-pathwayDB("Kegg")

gene.suppress<-c()
for(k in 1:length(d1)){
  gene.suppress<-c(gene.suppress,allPathways[[d1[k]]])
}
K.suppress<-table(gene.suppress)


gene.enhance<-c()
for(k in 1:length(d2)){
  gene.enhance<-c(gene.enhance,allPathways[[d2[k]]])
}
K.enhance<-table(gene.enhance)

par(mfrow = c(1,2))
wordcloud(names(K.suppress),K.suppress,min.freq=2)
wordcloud(names(K.enhance),K.enhance,min.freq=2)