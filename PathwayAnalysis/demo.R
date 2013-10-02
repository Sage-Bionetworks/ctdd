source("~/Sage-Analysis-Pipeline/PathwayAnalysis/pathwayAnalysis.R")

pathway.type = c("KEGG","Biocarta","Reactome","GO_BP","GO_CC","GO_MF")

# preparation : myReference list for GSEA
A<-read.table("~/Sage-Analysis-Pipeline/PathwayAnalysis/reference_demo.txt")
myReference<-A[,2]
names(myReference)<-A[,1]

# preparation : my test set for FET(fisher exact test)
B<-sample(names(myReference),300)

# MSigDB (synID="syn2135029" for GraphiteDB)
# run test with multicore with FET (fisher exact test)
results.fet<-pathwayAnalysis(synID="syn1681370",pathwayName = "KEGG",Reference = B,Test.method = "FET",cores = 8)

# combine all p-value from FET
Pval.fet<-c()
for(k in 1:length(results.fet)){
  Pval.fet<-c(Pval.fet,results.fet[[k]]$fetResult$p.value)
}
names(Pval.fet)<-names(results.fet)


# run test with multicore with GSEA
results<-pathwayAnalysis(synID="syn1681370",pathwayName = "KEGG",Reference = myReference,Test.method = "GSEA",cores = 8)

# combine all p-value from GSEA
Pval<-c()
for(k in 1:length(results)){
  Pval<-c(Pval,results[[k]]$gseaResult$p.value)
}
names(Pval)<-names(results)


# Do you want to plot all Enrichment score? then do following:
results[[135]]$gseaPlot(myReference,results[[135]]$gseaResult$geneset)