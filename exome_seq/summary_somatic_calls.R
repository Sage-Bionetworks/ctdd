library(synapseClient)
options(stringsAsFactors=FALSE)
setwd("projects/CTDD_Exome_Seq_201405/results/")

read<-function(SYNID){
  return(read.delim(file=synGet(SYNID)@filePath, sep="\t", header=T, stringsAsFactors=F, skip=1))  
}

QUERY <- synapseQuery("select name, id from entity where parentId == 'syn2497971'")

DATA <- lapply(QUERY$entity.id, read)
names(DATA) <- QUERY$entity.name

DATA.Hugo <- lapply(DATA, function(df){split(df,df$Hugo_Symbol)})

RES <- lapply(1:length(DATA.Hugo), function(i){
  
  x.Hugo<-DATA.Hugo[[i]]
  res.list <-lapply(x.Hugo, function(a){
    HUGO <- a$Hugo_Symbol[1]
    CLASSIFICATION <- paste(a$Variant_Classification, collapse=";")
    VARIANT <- paste(a$Variant_Type, collapse=";")
    TDEPTH <- paste(a$t_depth, collapse=";")
    return(c(HUGO, CLASSIFICATION,VARIANT,TDEPTH))
  })
  
  res.mat <- do.call("rbind", res.list)
  colnames(res.mat) <- paste(names(DATA)[i], c("Hugo_Symbol","Variant_Classification","Variant_Type","t_depth"), sep="_")
  
  OUTPUTFILE <- gsub(".maf",".tsv",names(DATA)[i])
  write.table(res.mat, file=OUTPUTFILE, col.names=T, row.names=F, sep="\t", quote=F)
  return(as.data.frame(res.mat))
})
