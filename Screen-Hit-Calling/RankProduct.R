# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# Utility functions for Rank Product
# =======================================

#####################################
# Computes rank product of a list of
# replicates of equal lentgths
# Computes by Geometric d-n
#####################################
rankProduct <- function(replicates = inreplicates, # list of replicates of equal length
                        gnames = innames,
                        ranks = FALSE
){
  
  if(!is.list(replicates)){
    stop("Error: Input is not a list of replicates ...\n")
  }
  
  l <- unlist(lapply(replicates, length))
  names(l) <- rep("", length(l))
  
  if(!isTRUE(all.equal(l, rep(l[1], length(l))))){
    stop("Error: Need to use alternative rank product test with different lengths of replicates...\n")
  }  
  
  r.means <- lapply(replicates, mean)
  
  # fold changes
  r.fcs <- lapply(as.list(seq(1:length(replicates))), function(idx){
    x <- replicates[[idx]]/r.means[[idx]]
    names(x) <- gnames 
    return(x)
  })
  
  # ranks of the fold changes
  r.ranks <- lapply(replicates, function(vect){
    x <- as.numeric(factor(vect))
    names(x) <- gnames 
    return(x)
  })
  
  # matrix form
  r.m.ranks <- t(do.call("rbind", r.ranks))
  
  # rank product
  rp <- apply(r.m.ranks, 1, function(row){
    (prod(row))^(1/length(row))
  })
  
  if(ranks == FALSE){
    # return the original rank product
    return(rp)
  }
  else{
    # return the rank of the rank product
    r <- as.numeric(factor(rp))
    names(r) <- gnames
    return(r)
  }
} #_end function rank product

permTest <- function(ng = in_ng, # number of elements (genes)
                     nr = in_nr, # number of replicates
                     np = 1000,  # number of permutations
                     erp = in_erp # experimental rank product result; named!
                     #gnames = in_names # names of elements
){
  # names of the elements
  gnames <- names(erp)
  
  set.seed(124)
  
  # create a list of np permuted rank products
  # of an experiments with ng elements and nr replicates
  perms <- lapply(as.list(seq(1:np)), function(x){
    # create a list of nr replicates
    y <- lapply(as.list(seq(1:nr)), function(idx){
      return(sample.int(ng, replace=FALSE))
    })
    
    # name them
    lapply(y, function(rep){names(rep) <- gnames})
    
    # matrix
    y.m <- t(do.call("rbind", y))
    
    # compute permuted rank product
    rp <- apply(y.m, 1, function(row){
      (prod(row))^(1/length(row))
    })
    
    return(rp)
  })
  
  # compute empirical p-value
  perms.m <- do.call("rbind", perms)
  my.pv <- sapply(seq(1:ng), function(idx){
    # how many times the permuted rank product was
    # lower than the experimentally determined one
    # for each gene
    cnt <- length(perms.m[,idx][perms.m[,idx] < erp[idx]])
    pv <- cnt/(np+1)
    return(pv)
  })
  
  names(my.pv) <- gnames
  return(my.pv)
}# end compute p-value
