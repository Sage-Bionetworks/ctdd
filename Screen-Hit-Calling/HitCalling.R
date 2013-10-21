# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# Comparison of hit calling procedures
# Produces hit calls according to different protocols
# Input:   
# synID:          synapse id of an siRNA screen
# outputFileName: output file name
# Output: 
# .tsv file with all computed statistics and final calls
# Usage: getHits(<synid>[, <outputfilename>])
# =======================================
require(synapseClient)
source("RankProduct.R")

path <- "/home/onikolov/projects/RNAi/ctdd/"

synIDs <- list("syn2275557"),
                "syn2275535",
                "syn2275534",
                "syn2275556",
                "syn2275538")

files <- list("0002_PAAD_Druggable_Vehicle_TGENQ3_2012_hits.tsv"),
                "0001_HNSCC_Druggable_Cisplatin_Mendez_hits.tsv",
                "0001_BRCA_ext_Kinome_Doxorubicin_Kemp_2012_hits.tsv",
                "0001_PAAD_Kinome_Vehicle_Moser_Feb2013_hits.tsv",
                "0002_HNSCC_Kinome_Vehicle_Mendez_2010_hits.tsv"
)

lapply(seq(1:length(files)), function(idx){
  getHits(synIDs[[idx]], paste(path, files[[idx]], sep=""))
})


getHits <- function(synID,      # Syanpse ID of the entity to analyze
                    outputFileName=NULL # output file name
){
  obj <- synGet(synID)  
  all <- read.delim(obj@filePath)
  
  # extract library only
  all <- all[all$Group=='Library',]
  
  # if treatment - extract vehicle
  if(length(grep('None', all$Compound)) < nrow(all)){  
    if(length(grep('Vehicle', all$Compound)) > 0){
      all <- all[all$Compound=='Vehicle',]
    }else if(length(grep('Vehicle', all$Compound)) > 0){
      stop("Error: Cannot analyze input file; data contains treatment but no vehicle...\n")
    }
  }
  
  # split by replicate
  splitData <-split(all, all$Replicate)
  
  # match replicates 
  idxValue <- grep("Pct_normalization_control", colnames(all))
  idxWell <- pmatch("Well", names(all))
  
  removeCols <- list("Barcode","Signal", "Replicate")
  idxOtherCols <- unlist(lapply(removeCols, function(colName){
    return(pmatch(colName, names(all)))
  }))
  
  if(sum(is.na(idxOtherCols)) > 0){
    stop("Error: siRNA input file is not in the expected format; 
         missing one of the following columns: 'Barcode','Signal', 'Replicate'...\n")
  }
  
  reps0 <- cbind(splitData[[1]][,-c(idxValue, idxOtherCols)])
  reps <- reps0[order(reps0$Gene_symbol,reps0$Well),]
  
  for(i in 1:length(splitData)){
    tmp <- splitData[[i]][order(splitData[[i]]$Gene_symbol, splitData[[i]]$Well),]    
    reps <- cbind.data.frame(reps, tmp[,idxValue])
    colnames(reps)[ncol(reps)] <- sprintf("replicate.%s", i)
  }
  
  if(sum(is.na(reps)) > 0){
    stop("Error: Some replicates are incomplete...\n")
  }
  
  # gene names
  gnames <- reps$Gene_symbol
  
  # generate unique gene names as duplicates present
  gnames.new <- lapply(as.list(seq(1:length(gnames))), function(idx){
    return(paste(gnames[idx],"_",idx, sep=""))
  })
  
  #################################
  # Compute Z-scores
  #################################
  
  # mean of the respective replicates for each gene
  my.mean <- apply(reps[,grep("replicate", names(reps))],1, mean)  
  #my.med <- apply(reps[,grep("replicate", names(reps))],1, median)
  
  # mean and median of the mean of three replicates
  my.mean.mean <- mean(my.mean)
  my.mean.med <- median(my.mean)
  
  #################################
  # Z-score N1: standard
  #################################
  my.zscore.mean <- (my.mean - my.mean.mean)/sd(my.mean)
  names(my.zscore.mean) <- gnames.new
  
  #################################
  # Z-score N2: robust 
  #################################
  MADCONST <- 1.4826
  my.mad <- MADCONST*median(abs(my.mean - my.mean.med )) # 21.0568
  my.zscore.robust.mean <- (my.mean - my.mean.med)/my.mad
  names(my.zscore.robust.mean) <- gnames.new
  
  #################################
  # Score N3: rank product
  #################################
  my.replicates <- lapply(grep("replicate", names(reps), value=TRUE), function(colName){
    return(reps[,colName])
  })
  
  rp.new <- rankProduct(my.replicates, gnames.new)
  
  # run permutation test
  numgenes <- length(reps$Gene_symbol)
  # generate p-value for rank product test
  pv_10000 <- permTest(ng = numgenes, # number of elements (genes)
                       nr = 3, # number of replicates
                       np = 10000,  # number of permutations
                       erp = rp.new # experimental rank product result; named!
  )
  
  #################################
  # Compute population-based hom. 
  # t-test p-value as recommended by James
  #################################
  pop <- do.call('c', reps[,grep("replicate", names(reps))])
  ttestPvalue <- unlist(lapply(seq(1:length(gnames.new)), function(rowIdx){
    return(t.test(reps[,grep("replicate", names(reps))][rowIdx,], pop, var.equal=TRUE)$p.value)
  }))
  
  #####################################
  # Summarize computed stats
  #####################################
  # add unique identifiers to distinguish duplicates
  reps$Gene_symbol_un <- unlist(gnames.new)
  # add mean
  reps$Mean <- my.mean
  # add z-score
  reps$Zscore_Mean <- my.zscore.mean
  # add robust z-score
  reps$Zscore_Robust <- my.zscore.robust.mean
  # add t-test p-value
  reps$TTest_PValue <- ttestPvalue
  # add rank product's p-value
  reps$RankProduct_PValue <- pv_10000
  
  ######################################
  # Call hits using different criteria
  #####################################
  # Significant by Z-score (Carla)
  reps$sig.Zscore_Mean <- rep(0, length(gnames.new))
  reps$sig.Zscore_Mean[reps$Zscore_Mean <= -0.6 & reps$TTest_PValue <= 0.1] <- 1
  
  reps$sig.Zscore_Robust <- rep(0, length(gnames.new))
  reps$sig.Zscore_Robust[reps$Zscore_Robust <= -0.6 & reps$TTest_PValue <= 0.1] <- 1
  
  reps$sig.RankProduct <- rep(0, length(gnames.new))
  reps$sig.RankProduct[reps$RankProduct_PValue <= 0.02738 ] <- 1 
  #  0.02738 gives exactly length = 805 in PAAD as the z-score with p-value
  
  # write to file here!
  if(is.null(outputFileName)){
    outputFileName=paste(synID,".tsv", sep="")
  }
  
  write.table(reps, file=outputFileName, sep="\t", row.names=FALSE)
  }

