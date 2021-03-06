rm(list = ls())
A1<-read.delim("~/DATA/siRNA_Screens/insock_processed_with_mock//HNSCC_kinome/CellLines_JHU019_uM_compound_0.txt")
A2<-read.delim("~/DATA/siRNA_Screens/insock_processed_with_mock/HNSCC_kinome/CellLines_UM-SCC-14A_uM_compound_0.txt")
A3<-read.delim("~/DATA/siRNA_Screens/insock_processed_with_mock/HNSCC_kinome/CellLines_UM-SCC-14C_uM_compound_0.txt")
A4<-read.delim("~/DATA/siRNA_Screens/insock_processed_with_mock/HNSCC_kinome/CellLines_UM-SCC-15A_uM_compound_0.txt")
A5<-read.delim("~/DATA/siRNA_Screens/insock_processed_with_mock/HNSCC_kinome/CellLines_UM-SCC-15B_uM_compound_0.txt")
A0<-read.delim("~/DATA/siRNA_Screens/insock_processed_with_mock/HNSCC_druggable/CellLines_E_uM_compound_0.txt")

a1<-subset(A1,A1$Group == "Library" & A1$Gene_symbol != "Control" & A1$Gene_symbol != "EMPTY" & A1$Gene_symbol != "Contr")
a2<-subset(A2,A2$Group == "Library" & A2$Gene_symbol != "Control" & A2$Gene_symbol != "EMPTY" & A2$Gene_symbol != "Contr")
a3<-subset(A3,A3$Group == "Library" & A3$Gene_symbol != "Control" & A3$Gene_symbol != "EMPTY" & A3$Gene_symbol != "Contr")
a4<-subset(A4,A4$Group == "Library" & A4$Gene_symbol != "Control" & A4$Gene_symbol != "EMPTY" & A4$Gene_symbol != "Contr")
a5<-subset(A5,A5$Group == "Library" & A5$Gene_symbol != "Control" & A5$Gene_symbol != "EMPTY" & A5$Gene_symbol != "Contr")
a0<-subset(A0,A0$Group == "Library" & A0$Gene_symbol != "Control" & A0$Gene_symbol != "EMPTY" & A0$Gene_symbol != "Contr")

# find consensus iinome hits
B1<-read.delim("~/RNAi_Analysis/HNSCC_Mock//z_transformation/Untreated/JHU019/Kinome/Hit_Call/General/table_general_hits.txt")
B2<-read.delim("~/RNAi_Analysis/HNSCC_Mock/z_transformation/Untreated/SCC-14A/Kinome/Hit_Call/General/table_general_hits.txt")
B3<-read.delim("~/RNAi_Analysis/HNSCC_Mock/z_transformation/Untreated/SCC-14C/Kinome/Hit_Call/General/table_general_hits.txt")
B4<-read.delim("~/RNAi_Analysis/HNSCC_Mock/z_transformation/Untreated/SCC-15A/Kinome/Hit_Call/General/table_general_hits.txt")
B5<-read.delim("~/RNAi_Analysis/HNSCC_Mock/z_transformation/Untreated/SCC-15B/Kinome/Hit_Call/General/table_general_hits.txt")
B0<-read.delim("~/RNAi_Analysis/HNSCC_Mock/z_transformation/Untreated/Patient/Druggable/Hit_Call/General/table_general_hits.txt")

hit.consensus<-intersect(B2$Gene_symbol,intersect(B3$Gene_symbol,intersect(B4$Gene_symbol,intersect(B5$Gene_symbol,B0$Gene_symbol))))
# intersect(B2$Gene_symbol,intersect(B3$Gene_symbol,intersect(B4$Gene_symbol,B5$Gene_symbol)))

write.table(hit.consensus,file = "~/RNAi_Analysis/HNSCC/z_transformation/Untreated/Hit_consensus_kinomes_MedianSummarization.txt",row.names=F,col.names=F,quote=F)
# sanity check
sanityCheck<-function(K1,K2){
  p<-match(K1$Accession_number,K2$Accession_number)
  if(length(which(is.na(p))>0)){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

sanityCheck(a1,a2)
sanityCheck(a1,a3)
sanityCheck(a1,a4)
sanityCheck(a1,a5)
sanityCheck(a1,a0)

p0<-match(a1$Accession_number,a0$Accession_number)

ALL<-cbind(a0[p0[-which(is.na(p0))],c(15:17)],
           a1[-which(is.na(p0)),c(15:17)],
           a2[-which(is.na(p0)),c(15:17)],
           a3[-which(is.na(p0)),c(15:17)],
           a4[-which(is.na(p0)),c(15:17)],
           a5[-which(is.na(p0)),c(15:17)]
)
ALL<-as.matrix(ALL)
rownames(ALL)<-a1$Gene_symbol[-which(is.na(p0))]
colnames(ALL)<-rep(c("HNSCC_druggable","HNSCC_JHU019","HNSCC_14A","HNSCC_14C","HNSCC_15A","HNSCC_15B"
),each=3)

COL<-as.character(as.matrix(read.delim("~/TEMP/color.txt",header = F)))

cc<-names(table(colnames(ALL)))
ColColorBar1<-rep("grey",ncol(ALL))
for(k in 1:length(cc)){
  a1<-which(colnames(ALL)==cc[k])
  ColColorBar1[a1]<-COL[k]
}


mycol <- colorpanel(n=99,low="red",mid="white",high="green")

D3<-apply(ALL,2,rank)
png("~/RNAi_Analysis/HNSCC/z_transformation/Untreated/heatmap_hit_consensus_kinomes_Clustering.png",width = 360,height = 1800)
aaa<-heatmap.2((D3), col=mycol,Rowv = T,Colv=T,
               scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
               trace="none", cexRow=0.5, cexCol=1.5, mar = c(15,10),symkey=F,breaks=seq(min(D3),max(D3),length = length(mycol)+1),
               ColSideColors = ColColorBar1)
dev.off()

# lrll
aaa1<-cut(aaa$rowDendrogram,2624)
aaa2<-cut(aaa1$lower[[1]],2042)
aaa3<-cut(aaa2$lower[[2]],1785)
aaa4<-cut(aaa3$lower[[1]],1563)
gene.cluster<-labels(aaa4$lower[[1]])

ColColorBar2<-rep("grey",nrow(ALL))
ppp<-match(gene.cluster, rownames(ALL))
ColColorBar2[ppp]<-"gold"

png("~/RNAi_Analysis/HNSCC//z_transformation/Untreated/heatmap_hit_consensus_kinomes_Clustering.png",width = 360,height = 1800)
aaa<-heatmap.2((D3), col=mycol,Rowv = T,Colv=T,
               scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
               trace="none", cexRow=0.5, cexCol=1.5, mar = c(15,10),symkey=F,breaks=seq(min(D3),max(D3),length = length(mycol)+1),
               ColSideColors = ColColorBar1,
               RowSideColors = ColColorBar2)
dev.off()




# exclude JHU 
ALL.exclude<-cbind(a0[p0[-which(is.na(p0))],c(15:17)],           
           a2[-which(is.na(p0)),c(15:17)],
           a3[-which(is.na(p0)),c(15:17)],
           a4[-which(is.na(p0)),c(15:17)],
           a5[-which(is.na(p0)),c(15:17)]
)
ALL.exclude<-as.matrix(ALL.exclude)
rownames(ALL.exclude)<-a2$Gene_symbol[-which(is.na(p0))]
colnames(ALL.exclude)<-rep(c("HNSCC_druggable","HNSCC_14A","HNSCC_14C","HNSCC_15A","HNSCC_15B"
),each=3)

COL<-as.character(as.matrix(read.delim("~/TEMP/color.txt",header = F)))

cc<-names(table(colnames(ALL.exclude)))
ColColorBar1<-rep("grey",ncol(ALL.exclude))
for(k in 1:length(cc)){
  a1<-which(colnames(ALL.exclude)==cc[k])
  ColColorBar1[a1]<-COL[k]
}


mycol <- colorpanel(n=99,low="red",mid="white",high="green")

D3.exclude<-apply(ALL.exclude,2,rank)
png("~/RNAi_Analysis/HNSCC/z_transformation/Untreated/heatmap_hit_consensus_kinomes_Clustering_exclude.png",width = 360,height = 1800)
aaa<-heatmap.2((D3.exclude), col=mycol,Rowv = T,Colv=T,
               scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
               trace="none", cexRow=0.5, cexCol=1.5, mar = c(15,10),symkey=F,breaks=seq(min(D3),max(D3),length = length(mycol)+1),
               ColSideColors = ColColorBar1)
dev.off()

# lrrrrl
aaa.1<-cut(aaa$rowDendrogram,2398)
aaa.2<-cut(aaa.1$lower[[1]],1872)
aaa.3<-cut(aaa.2$lower[[2]],1633)
aaa.4<-cut(aaa.3$lower[[2]],1377)
aaa.5<-cut(aaa.4$lower[[2]],1345)
aaa.6<-cut(aaa.5$lower[[2]],1207)

gene.cluster<-labels(aaa.6$lower[[1]])
write.table(gene.cluster,file = "~/RNAi_Analysis/HNSCC/z_transformation/Untreated/Hit_consensus_kinomes_Cluster_exclude.txt",row.names=F,col.names=F,quote=F)



ColColorBar2<-rep("grey",nrow(ALL.exclude))
ppp<-match(gene.cluster, rownames(ALL.exclude))
ColColorBar2[ppp]<-"gold"

png("~/RNAi_Analysis/HNSCC/z_transformation/Untreated/heatmap_hit_consensus_kinomes_Clustering_exclude.png",width = 360,height = 1800)
aaa<-heatmap.2(D3.exclude, col=mycol,Rowv = T,Colv=T,
               scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
               trace="none", cexRow=0.5, cexCol=1.5, mar = c(15,10),symkey=F,breaks=seq(min(D3),max(D3),length = length(mycol)+1),
               ColSideColors = ColColorBar1,
               RowSideColors = ColColorBar2)
dev.off()



####### not rank clustering
mycol <- colorpanel(n=99,low="red",mid="white",high="green")

D3<-ALL.exclude

png("~/RNAi_Analysis/HNSCC/z_transformation/Untreated/heatmap_hit_consensus_kinomes_Clustering_z.png",width = 480,height = 1800)
aaa<-heatmap.2((D3), col=mycol,Rowv = T,Colv=T,
               scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
               trace="none", cexRow=0.5, cexCol=1.5, mar = c(15,10),symkey=F,breaks=seq(-4,4,length = length(mycol)+1),
               ColSideColors = ColColorBar1)
dev.off()

#lllr
aaa1<-cut(aaa$rowDendrogram,30.2)
aaa2<-cut(aaa1$lower[[1]],19.8)
aaa3<-cut(aaa2$lower[[1]],17.2)
aaa4<-cut(aaa3$lower[[1]],14.4)
gene.cluster<-labels(aaa4$lower[[2]])

write.table(gene.cluster,file = "~/RNAi_Analysis/HNSCC/z_transformation/Untreated/Hit_consensus_kinomes_Cluster_z.txt",row.names=F,col.names=F,quote=F)

ColColorBar2<-rep("grey",nrow(ALL.exclude))
ppp<-match(gene.cluster, rownames(ALL.exclude))
ColColorBar2[ppp]<-"gold"

png("~/RNAi_Analysis/HNSCC//z_transformation/Untreated/heatmap_hit_consensus_kinomes_Clustering_z.png",width = 480,height = 1800)
aaa<-heatmap.2((D3), col=mycol,Rowv = T,Colv=T,
               scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
               trace="none", cexRow=0.5, cexCol=1.5, mar = c(15,10),symkey=F,breaks=seq(-4,4,length = length(mycol)+1),
               ColSideColors = ColColorBar1,
               RowSideColors = ColColorBar2)
dev.off()


source("~/RNAi_Analysis/work_flow/rankProd.R")
PP<-rankProd(ALL.exclude)
p1<-which(PP$pval<=0.05)
q1<-which(PP$zscore<=0)
P<-intersect(p1,q1)
write.table(rownames(ALL.exclude)[P],file = "~/RNAi_Analysis/HNSCC_Mock/z_transformation/Untreated/Hit_consensus_rankProd.txt",row.names=F,col.names=F,quote=F)


gene.cluster<-rownames(ALL.exclude)[P]

ColColorBar2<-rep("grey",nrow(ALL.exclude))
ppp<-match(gene.cluster, rownames(ALL.exclude))
ColColorBar2[ppp]<-"gold"

mycol <- colorpanel(n=99,low="red",mid="white",high="green")

png("~/RNAi_Analysis/HNSCC//z_transformation/Untreated/heatmap_hit_consensus_kinomes_rankProduct.png",width = 480,height = 1800)
aaa<-heatmap.2((D3.exclude), col=mycol,Rowv = T,Colv=T,
               scale="none", key=TRUE, density.info="none",keysize = 0.75, #breaks = pairs.breaks, 
               trace="none", cexRow=0.5, cexCol=1.5, mar = c(15,10),symkey=F,breaks=seq(min(D3.exclude),max(D3.exclude),length = length(mycol)+1),
               ColSideColors = ColColorBar1,
               RowSideColors = ColColorBar2)
dev.off()
