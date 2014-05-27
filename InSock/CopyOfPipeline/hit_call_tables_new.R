rm(list = ls())

P<-read.delim("~/DATA/siRNA_Screens/insock_processed/HNSCC_druggable/CellLines_E_uM_compound_0.txt")
P.lib<-subset(P,P$Group == "Library")
median.Untreated.z.uni<-apply(P.lib[,c(15:17)],1,median)

S<-read.delim("~/DATA/siRNA_Screens/insock_processed_with_mock/HNSCC_druggable/CellLines_E_uM_compound_0.txt")
S.lib<-subset(S,S$Group == "Library")

median.Untreated.pct<-apply(S.lib[,c(12:14)],1,median)
median.Untreated.z.mock<-apply(S.lib[,c(15:17)],1,median)

Q<-cbind(P.lib[,c(1:11)],median.Untreated.pct,median.Untreated.z.uni,median.Untreated.z.mock)

P1<-read.delim("~/DATA/siRNA_Screens/insock_processed/HNSCC_druggable/CellLines_E_uM_compound_0.249.txt")
P1.lib<-subset(P1,P1$Group == "Library")
median.Treated.z.uni<-apply(P1.lib[,c(15:17)],1,median)

S1<-read.delim("~/DATA/siRNA_Screens/insock_processed_with_mock/HNSCC_druggable/CellLines_E_uM_compound_0.249.txt")
S1.lib<-subset(S1,S1$Group == "Library")

median.Treated.pct<-apply(S1.lib[,c(12:14)],1,median)
median.Treated.z.mock<-apply(S1.lib[,c(15:17)],1,median)

Q1<-cbind(P1.lib[,c(1:11)],median.Treated.pct,median.Treated.z.uni,median.Treated.z.mock)

q1<-paste(Q1$Row,Q1$Column,Q1$Gene_symbol,Q1$Well_ID,sep="-")
q<-paste(Q$Row,Q$Column,Q$Gene_symbol,Q$Well_ID,sep="-")

aa<-match(q,q1)

QQ<-cbind(Q,Q1[aa,c(12:14)])[,-c(1,4,5,7,11)]


outfold = "~/RNAi_Analysis/HNSCC_Mock/pct_normalization_control///"
#QQQ<-QQ[,-c(4,5,6,8,9)]
QQQ<-QQ[,-c(8,9,11,12)]
X<-sort(QQQ$median.Untreated.pct,index.return=T)

# kinome hits from rank cluster
A1<-read.delim(paste(outfold,"/Untreated/Hit_consensus_kinomes_Cluster_exclude.txt",sep = ""),header = F)
hit.kinome.cluster<-as.character(as.matrix(A1$V1))

# kinome hits using median summarization
A2<-read.delim(paste(outfold,"/Untreated/Hit_consensus_kinomes_MedianSummarization.txt",sep = ""),header = F)
# hit.kinome.median<-as.character(as.matrix(A2$V1))
hit.kinome.median<-c()

# kinome hits from z-cluster
A3<-read.delim(paste(outfold,"/Untreated/Hit_consensus_kinomes_Cluster_z.txt",sep = ""),header = F)
# hit.kinome.zcluster<-as.character(as.matrix(A3$V1))
hit.kinome.zcluster<-c()

# kinome hits from rank product
A4<-read.delim(paste(outfold,"/Untreated/Hit_consensus_rankProd.txt",sep = ""),header = F)
hit.kinome.rankproduct<-as.character(as.matrix(A4$V1))

load("~/RNAi_Analysis/Data_Molecular_Chracterization/HNSCC/EXP.Rdata")
A3<-EXP.Sanger$pval
names(A3)<-rownames(EXP.Sanger$cancer)
hit.exp.sanger<-names(which(A3<=0.05))
A4<-EXP.TCGA$pval
names(A4)<-rownames(EXP.TCGA$cancer)
hit.exp.tcga<-names(which(A4<=0.05))

load("~/RNAi_Analysis/Data_Molecular_Chracterization/HNSCC/COPY.Rdata")
A5<-copy.Sanger$pval
names(A5)<-rownames(copy.Sanger$copy)
hit.copy.sanger<-names(which(A5<=0.05 & copy.Sanger$Mean>0))
A6<-copy.TCGA$pval
names(A6)<-rownames(copy.TCGA$copy)
hit.copy.tcga<-names(which(A6<=0.05 & copy.TCGA$Mean >0))

# pathway relevant hits from untreated hits
A7.1<-read.delim(paste(outfold,"/Untreated/Patient/Druggable/Hit_Call//Pathway_Relevant/Canonical/LeadingEdgeGenes_KEGG.txt",sep = ""),header = F)
A7.2<-read.delim(paste(outfold,"/Untreated/Patient/Druggable/Hit_Call//Pathway_Relevant/Canonical/LeadingEdgeGenes_BIOCARTA.txt",sep = ""),header = F)
A7.3<-read.delim(paste(outfold,"/Untreated/Patient/Druggable/Hit_Call//Pathway_Relevant/Canonical/LeadingEdgeGenes_REACTOME.txt",sep = ""),header = F)
hit.untreated.pathway_canonical <-as.character(as.matrix(union(A7.1$V1,union(A7.2$V1,A7.3$V1))))

A8.1<-read.delim(paste(outfold,"/Untreated/Patient/Druggable/Hit_Call//Pathway_Relevant/GO/LeadingEdgeGenes_GO_BP.txt",sep = ""),header = F)
A8.2<-read.delim(paste(outfold,"/Untreated/Patient/Druggable/Hit_Call//Pathway_Relevant/GO/LeadingEdgeGenes_GO_MF.txt",sep = ""),header = F)
hit.untreated.pathway_go <-as.character(as.matrix(union(A8.1$V1,A8.2$V1)))

# pathway relevant hits from treated hits
A9.1<-read.delim(paste(outfold,"/Treated/Patient/Druggable/Hit_Call//Pathway_Relevant/Canonical/LeadingEdgeGenes_KEGG.txt",sep = ""),header = F)
A9.2<-read.delim(paste(outfold,"/Treated/Patient/Druggable/Hit_Call//Pathway_Relevant/Canonical/LeadingEdgeGenes_BIOCARTA.txt",sep = ""),header = F)
A9.3<-read.delim(paste(outfold,"/Treated/Patient/Druggable/Hit_Call//Pathway_Relevant/Canonical/LeadingEdgeGenes_REACTOME.txt",sep = ""),header = F)
hit.treated.pathway_canonical <-as.character(as.matrix(union(A9.1$V1,union(A9.2$V1,A9.3$V1))))

A10.1<-read.delim(paste(outfold,"/Treated/Patient/Druggable/Hit_Call//Pathway_Relevant/GO/LeadingEdgeGenes_GO_BP.txt",sep = ""),header = F)
A10.2<-read.delim(paste(outfold,"/Treated/Patient/Druggable/Hit_Call//Pathway_Relevant/GO/LeadingEdgeGenes_GO_MF.txt",sep = ""),header = F)
hit.treated.pathway_go <-as.character(as.matrix(union(A10.1$V1,A10.2$V1)))


# leading edge common pathways + 45 degree line application
A.11<-read.delim(paste(outfold,"/sensitizer_commonLeadingEdgeGenes_Canonical.txt",sep = ""),header = F)
hit.sensitizer.pathway_canonical<-as.character(as.matrix(A.11$V1))
A.12<-read.delim(paste(outfold,"/sensitizer_commonLeadingEdgeGenes_GO.txt",sep = ""),header = F)
hit.sensitizer.pathway_go<-as.character(as.matrix(A.12$V1))


Libraries<-as.character(as.matrix(QQQ$Gene_symbol[X$ix]))
kinome.cluster<-match(Libraries,hit.kinome.cluster)
kinome.median<-match(Libraries,hit.kinome.median)
kinome.zcluster<-match(Libraries,hit.kinome.zcluster)
kinome.rankproduct<-match(Libraries,hit.kinome.rankproduct)
exp.sanger<-match(Libraries,hit.exp.sanger)
exp.tcga<-match(Libraries,hit.exp.tcga)
copy.sanger<-match(Libraries,hit.copy.sanger)
copy.tcga<-match(Libraries,hit.copy.tcga)
untreated.pathway_canonical<-match(Libraries,hit.untreated.pathway_canonical)
untreated.pathway_go<-match(Libraries,hit.untreated.pathway_go)
treated.pathway_canonical<-match(Libraries,hit.treated.pathway_canonical)
treated.pathway_go<-match(Libraries,hit.treated.pathway_go)
sensitizer.pathway_canonical<-match(Libraries,hit.sensitizer.pathway_canonical)
sensitizer.pathway_go<-match(Libraries,hit.sensitizer.pathway_go)


vec1<-rep(1,length(Libraries))
vec2<-rep(1,length(Libraries))
vec3<-rep(1,length(Libraries))
vec4<-rep(1,length(Libraries))
vec5<-rep(0,length(Libraries))
vec6<-rep(0,length(Libraries))
vec7<-rep(0,length(Libraries))
vec8<-rep(0,length(Libraries))
vec9<-rep(1,length(Libraries))
vec0<-rep(1,length(Libraries))
vec11<-rep(1,length(Libraries))
vec12<-rep(1,length(Libraries))
vec13<-rep(1,length(Libraries))
vec14<-rep(1,length(Libraries))


vec1[which(is.na(kinome.cluster))]<-0
vec2[which(is.na(kinome.median))]<-0
vec3[which(is.na(kinome.zcluster))]<-0
vec4[which(is.na(kinome.rankproduct))]<-0
vec5[exp.sanger[which(!is.na(exp.sanger))]]<-1
vec6[exp.tcga[which(!is.na(exp.tcga))]]<-1
vec7[copy.sanger[which(!is.na(copy.sanger))]]<-1
vec8[copy.tcga[which(!is.na(copy.tcga))]]<-1
vec9[which(is.na(untreated.pathway_canonical))]<-0
vec0[which(is.na(untreated.pathway_go))]<-0
vec11[which(is.na(treated.pathway_canonical))]<-0
vec12[which(is.na(treated.pathway_go))]<-0
vec13[which(is.na(sensitizer.pathway_canonical))]<-0
vec14[which(is.na(sensitizer.pathway_go))]<-0

library(gdata)
IN<-read.xls("~/RNAi_Analysis/Rosetta Sigma Druggable Genome Library Drug Gene Compilation v2.xlsx")
b<-match(Libraries,IN$Gene.Symbol)
drug<-IN$Drug.s.[b]
Drug<-gsub("--","",drug)
Drug[which(is.na(b))]<-""



MAT<-cbind(as.matrix(QQQ[X$ix,]),vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8,vec9,vec0,vec11,vec12,vec13,vec14,Drug)
colnames(MAT)<-c(names(QQQ),"hit.kinome.cluster","hit.kinome.median","hit.kinome.zcluster","hit.kinome.rankproduct",
                 "hit.exp.sanger","hit.exp.tcga","hit.copy.sanger","hit.copy.tcga",
                 "hit.untreated.pathway.canonical","hit.untreated.pathway.go",
                 "hit.treated.pathway.canonical","hit.treated.pathway.go",
                 "hit.sensitizer.pathway.canonical","hit.sensitizer.pathway.go",
                 "Inhibitor")

write.table(MAT,   file =paste(outfold,"/new_Summary_Hits_ver01.txt",sep = ""),row.names=F,col.names=T,quote=F,sep = "\t")



load("~/Result_priorIncorporateStepwiseRegression_filterVar02/Analysis/drugInfo_new/Sanger/drugInfo_E_BIOCARTA.Rdata")
genes.biocarta<-Genes$Cisplatin
load("~/Result_priorIncorporateStepwiseRegression_filterVar02/Analysis/drugInfo_new/Sanger/drugInfo_E_KEGG.Rdata")
genes.kegg<-Genes$Cisplatin
load("~/Result_priorIncorporateStepwiseRegression_filterVar02/Analysis/drugInfo_new/Sanger/drugInfo_E_NCI.Rdata")
genes.nci<-Genes$Cisplatin
load("~/Result_priorIncorporateStepwiseRegression_filterVar02/Analysis/drugInfo_new/Sanger/drugInfo_E_GO_BP.Rdata")
genes.gobp<-Genes$Cisplatin
load("~/Result_priorIncorporateStepwiseRegression_filterVar02/Analysis/drugInfo_new/Sanger/drugInfo_E_GO_MF.Rdata")
genes.gomf<-Genes$Cisplatin



X1<-sort(genes.biocarta[,2],decreasing = T, index.return=T)
X2<-sort(genes.kegg[,2],decreasing = T, index.return=T)
X3<-sort(genes.nci[,2],decreasing = T, index.return=T)
X4<-sort(genes.gobp[,2],decreasing = T, index.return=T)
X5<-sort(genes.gomf[,2],decreasing = T, index.return=T)

x1<-genes.biocarta[X1$ix,]
x2<-genes.kegg[X2$ix,]
x3<-genes.nci[X3$ix,]
x4<-genes.gobp[X4$ix,]
x5<-genes.gomf[X5$ix,]


Y<-(union(x1[,1],union(x2[,1],union(x3[,1],union(x4[,1],x5[,1])))))
y1<-match(x1[,1],Y)
y2<-match(x2[,1],Y)
y3<-match(x3[,1],Y)
y4<-match(x4[,1],Y)
y5<-match(x5[,1],Y)

MAT<-matrix(0,nrow = length(Y),ncol = 5)
rownames(MAT)<-Y
colnames(MAT)<-c("BIOCARTA","KEGG","NCI","GO_BP","GO_MF")
MAT[y1,1]<-as.numeric(x1[,2])
MAT[y2,2]<-as.numeric(x2[,2])
MAT[y3,3]<-as.numeric(x3[,2])
MAT[y4,4]<-as.numeric(x4[,2])
MAT[y5,5]<-as.numeric(x5[,2])


M<-apply(MAT,1,mean)
plot(M)
M1<-apply(MAT[,c(1:3)],1,mean)
plot(M1)


load("~/SGLR/graphite_pathways_structure.Rdata")
Z<-(union(rownames(structure.BIOCARTA),union(rownames(structure.GO_BP),union(rownames(structure.KEGG),union(rownames(structure.NCI),rownames(structure.GO_MF))))))
MAT1<-matrix(0, nrow = length(Z),ncol = 5)
rownames(MAT1)<-Z
colnames(MAT1)<-c("BIOCARTA","KEGG","NCI","GO_BP","GO_MF")
z1<-match(rownames(structure.BIOCARTA),Z)
z2<-match(rownames(structure.KEGG),Z)
z3<-match(rownames(structure.NCI),Z)
z4<-match(rownames(structure.GO_BP),Z)
z5<-match(rownames(structure.GO_MF),Z)
MAT1[z1,1]<-1
MAT1[z2,2]<-1
MAT1[z3,3]<-1
MAT1[z4,4]<-1
MAT1[z5,5]<-1

N<-apply(MAT1,1,sum)


NN<-N[match(names(M),names(N))]

P<-cbind(M,NN)
PP<-P[,1]/P[,2]
plot(PP)

Q<-sort(PP,decreasing = T)



BB<-read.xls("~/ctdd2013/HNSCC/HNSCC_interestingGenes.xls")
as.character(as.matrix(BB$GeneSymbol))

length(intersect(toupper(as.character(as.matrix(BB$GeneSymbol))),names(M)))
intersect(toupper(as.character(as.matrix(BB$GeneSymbol))),names(M)[which(M>=1.5)])

C1<- read.delim("~/ctdd2013/merged_ch.on_PCT_Pkill.ins_PNC.tsv")
C2<- read.delim("~/ctdd2013/merged_ch.on_PCT_Psen.ins_PNC.tsv")


R.30<-as.character(as.matrix(C1$Symbol[(which(C1$Resist_30 == 1))]))
S.30<-as.character(as.matrix(C1$Symbol[(which(C1$Sen_30 == 1))]))
length(intersect(toupper(as.character(as.matrix(BB$GeneSymbol))),R.30))
length(intersect(toupper(as.character(as.matrix(BB$GeneSymbol))),S.30))


RR.30<-as.character(as.matrix(C1$Symbol[(which(C2$Resist_30 == 1))]))
SS.30<-as.character(as.matrix(C1$Symbol[(which(C2$Sen_30 == 1))]))
length(intersect(toupper(as.character(as.matrix(BB$GeneSymbol))),RR.30))
length(intersect(toupper(as.character(as.matrix(BB$GeneSymbol))),SS.30))

A<-read.delim("~/RNAi_Analysis/HNSCC_Mock/pct_normalization_control/new_Summary_Hits_ver01.txt")
names(M)[which(M>=1.5)]
cc1<-match(A$Gene_symbol,names(M)[which(M>=1.5)])

B<-names(M)[which(M>=1.5)]
cc1<-match(A$Gene_symbol,B)


Pred.Features<-matrix(0,nrow = nrow(A),ncol = 1)
Pred.Features[which(!is.na(cc1))]<-1


load("~/ctdd2013/OV/mut_freq_TCGAlive_OV.Rdata")
mutfreq<-function(x){
  length(which(!is.na(x)))
}
mm<-apply(MUT,1,mutfreq)

mut.freq<-mm/ncol(MUT)


match(A$Gene_symbol,names(mut.freq))
cc<-match(A$Gene_symbol,names(mut.freq))
TCGA.Mut.freq<-mut.freq[cc]

AA<-cbind(A,TCGA.Mut.freq,Pred.Features)

sum(AA$Pred.Features)

write.table(AA,   file ="~/RNAi_Analysis/OV_Mock/pct_normalization_control/new_Summary_Hits_ver02.txt",row.names=F,col.names=T,quote=F,sep = "\t")



