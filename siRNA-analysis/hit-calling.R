# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# Comparison of hit calling procedures
# Produces hit calls according to different protocols
# Input:   
# synID:          synapse id of an siRNA screen
# outputFileName: output file name
# modifier:           modifier of nalysis: 0 or 1
#                 FALSE: single group (None, Vechicle-only, Treatment-only)
#                 TRUE: 2 groups: Vechicle and some treatment (non-"None")
# Output: 
# .tsv file with all computed statistics and final calls
# Usage: getHits(<synid>[, <outputfilename>])
# =======================================
require(synapseClient)
source("RankProduct.R")

getHits <- function(synID,  # Syanpse ID of the entity to analyze
                    inputFileName=NULL, # output file name
                    modifier=0 # modifier of screen analysis: 
                    #       0: 1 group (None, Vechicle only, Treatment only)
                    #       1: 2 groups: Vechicle and the other treatment (non-"None")
                    ){
  # declare constants
  MADCONST <- 1.4826
  
  # read in data
  obj <- synGet(synID)  
  d <- read.delim(obj@filePath)
  
  # extract library only
  d <- d[d$Group=='Library',]
  
  # extract all cell lines
  cls <- unique(as.vector(d$Cell_line))
  
  # extract all treatments
  comps <- unique(as.vector(d$Compound))
  
for(cl in cls){
  if(modifier){
    all.treatment <- NULL
  }
  for(comp in comps){
    
    all <- d[d$Cell_line == cl & d$Compound == comp,]
    
    #Exclude misc non-library rows: Empty, Treated Mock, Blank, and Controls (PLK1, UNI
    exclude <- grep( paste(c("EMPTY","MOCK","BLANK","UNI","PLK1"), collapse="|"), 
                     toupper(all$Gene_symbol))
    
    if(length(exclude) > 0){
      all <- all[-exclude,]
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
    my.mad <- MADCONST*median(abs(my.mean - my.mean.med )) # 21.0568
    my.zscore.robust.mean <- (my.mean - my.mean.med)/my.mad
    names(my.zscore.robust.mean) <- gnames.new
    
    #################################
    # Score N3: rank product
    #################################
    my.replicates <- lapply(grep("replicate", names(reps), value=TRUE), function(colName){
      return(reps[,colName])
    })
    
    # compute rank product
    rp <- rankProduct(my.replicates, gnames.new)
    
    # generate p-value for rank product test
    numgenes <- length(reps$Gene_symbol)
    pv_10000 <- permTest(ng = numgenes, # number of elements (genes)
                         nr = 3,        # number of replicates
                         np = 10000,    # number of permutations
                         erp = rp   # experimental rank product result; named!
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
    
    #######################################
    # Record for across group analysis
    #######################################
    if(modifier){
      # get the following columns for each treatment
      cols <- unlist(lapply(list("Gene_symbol","Gene_symbol_un","Gene_ID","Mean"), function(colName){
        return(pmatch(colName, names(reps)))
      }))
      
      cols <- c(cols, grep("REPLICATE", toupper(names(reps))))
      
      if(is.null(all.treatment)){
        
        all.treatment <- cbind(reps[,cols])
        
        names(all.treatment) <- c("Gene_symbol","Gene_symbol_un","Gene_ID",paste("Mean",comp, sep="_"), 
                                  "replicate.1","replicate.2","replicate.3")
        
        
      } else {
        newnames <- c(names(all.treatment),"Gene_symbol","Gene_symbol_un","Gene_ID",
                                  paste("Mean",comp, sep="_"), "replicate.1","replicate.2","replicate.3")

        all.treatment <- cbind(all.treatment, reps[,cols])
        names(all.treatment) <- newnames
        
      }
      
    }#_end recording for across groups
    
    # write to file within group hits
    outputFileName <- ""
    
    if(is.null(inputFileName)){
      outputFileName=paste(synID, cl, comp, "hits.tsv", sep="_")
    } else {
      base <- strsplit(inputFileName, ".tsv")
      outputFileName=paste(base, cl, comp, "hits.tsv", sep="_")
    }
    
    write.table(reps, file=outputFileName, sep="\t", row.names=FALSE)
  } #_end for comps
  
  if(modifier){
    
    # get the treatment and vehicle rows
    vecmeanidx <- grep("MEAN_VEHICLE", toupper(names(all.treatment)))
    if(is.na(vecmeanidx)){
      stop("Error: No mean for vehicle was found for across groups analysis...\n")
      }
      
    trtname <- ""
    for(comp in comps){
      if(toupper(comp) != "VEHICLE"){
        trtname <- comp
      }
    }
    trtmeanidx <- grep(toupper(paste("Mean", trtname, sep="_")), toupper(names(all.treatment)))
    if(is.na(trtmeanidx)){
      stop("Error: No mean for treatment was found for across groups analysis...\n")
    }
    
    # ratio treatment to vehicle
    fc <- all.treatment[,trtmeanidx] / all.treatment[,vecmeanidx]
    
    fc.mean <- mean(fc)
    fc.med <- median(fc)
    
    # Z-score
    fc.zscore <- (fc - fc.mean)/sd(fc)
    
    # Robust Z-score
    fc.mad <- MADCONST*median(abs(fc - fc.med )) # 21.0568
    fc.zscore.robust <- (fc - fc.med)/fc.mad
    
    # T-tes p-value
    ttestpv <- unlist(lapply(1:nrow(all.treatment), function(row){
      return(t.test(all.treatment[row, c(5,6,7)], all.treatment[row, c(12,13,14)])$p.value) # change to grepping
    }))
    
    # record statistics
    all.treatment$Z_score <- fc.zscore
    all.treatment$Z_score_robust <- fc.zscore.robust
    all.treatment$TTestPV <- ttestpv
      
    # write file modifier to vechicle hits
    outputFileName <- ""
    
    if(is.null(inputFileName)){
      outputFileName=paste(synID, cl, trtname, "Vehicle_hits.tsv", sep="_")
    } else {
      base <- strsplit(inputFileName, ".tsv")
      outputFileName=paste(base, cl, trtname, "Vehicle_hits.tsv", sep="_")
    }
    
    write.table(all.treatment, file=outputFileName, sep="\t", row.names=FALSE)
    
  }
  
} #_end of for cls
} #_end get hits
