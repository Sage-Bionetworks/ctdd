hit_call_general<-function(inputfile,outfolder,norm_method){  
  # kinome or druggable
  A1<-read.delim(inputfile)
  
  AA1<-subset(A1,A1$Group == "Library" & A1$Gene_symbol != "Control" & A1$Gene_symbol != "EMPTY" & A1$Gene_symbol != "Contr")
  if(norm_method == "pct"){
    bb1<-apply(AA1[,12:14],1,median)
    a<-which(bb1<= 65)
    AA2<-AA1[a,]
    AA2<-cbind(AA1[a,c(1:14)],bb1[a])
    names(AA2)[15]<-"pct_Median"
    write.table(AA2,file = paste(outfolder,"/table_general_hits.txt",sep = ""),row.names=F,col.names=T,quote=F,sep = "\t")
  }else{
    bb1<-apply(AA1[,15:17],1,median)  
    a<-which(bb1<= -1)
    AA2<-cbind(AA1[a,c(1:11,15:17)],bb1[a])
    names(AA2)[15]<-"norm_Median"
    write.table(AA2,file = paste(outfolder,"/table_general_hits.txt",sep = ""),row.names=F,col.names=T,quote=F,sep = "\t")
  }  
}