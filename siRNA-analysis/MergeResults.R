# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# Quick results merge

require(gdata) 

bimodalFiles <- list("PAAD_druggable.tsv")
results <- list("0002_PAAD_Druggable_Vehicle_TGENQ3_2012_hits.tsv")

lapply(seq(1:length(results)), function(idx){
  mergeResults(bimodalFiles[[idx]], results[[idx]])
})

mergeResults <- function(bimodalFile,
                         resultsFile){

  pathBiM <- "/home/onikolov/projects/RNAi/ctdd/results/bimodal_zscore_InSock/"
  pathRes <- "/home/onikolov/projects/RNAi/ctdd/"
  
  d <- read.delim(paste(pathBiM, bimodalFile, sep=""), header=TRUE, sep="\t")
  
  cols <- list("Well","Gene.Symbol", "newZscore")
  idxCols <- unlist(lapply(cols, function(colName){
    return(pmatch(colName, names(d)))
  }))
  
  dd <- d[,idxCols]
  names(dd) <- c("Well","Gene_symbol","Zscore_Bimodal")
  
  az <- read.table(paste(pathRes, resultsFile, sep=""), header=TRUE, sep="\t")
  
  x <- merge(az, dd, by=c("Well", "Gene_symbol"))
  
  write.table(x, file=paste(pathRes, "merged_", resultsFile, sep=""), sep="\t", row.names=FALSE)
}
