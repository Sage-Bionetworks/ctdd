# inputTable can be downloaded from Synapse directly
# outfile must be temporary folder to save the normalized table
inputTable<-"~/DATA/siRNA_Screens/real/MendezheadAndNeckKinome2010.tsv"
outfile<- "~/DATA/siRNA_Screens/insock_processed/HNSCC_kinome//"
system(paste("mkdir ",outfile,"/Distribution",sep = ""))

INSOCKsummaryTable(inputTable,outfile)

INSOCKsummaryTable<-function(inputTable,outfile){
  # standard input table from James
  INPUT<-read.delim(inputTable)
  
  # how many cell lines?
  cellLines<-names(table(INPUT$Cell_line))
  for(k1 in 1:length(cellLines)){
    a1<-which(INPUT$Cell_line == cellLines[k1])
    INPUT.1 <- INPUT[a1,]
    # is it modifier screen?
    uM_compound<-names(table(INPUT.1$uM_compound))
    for(k2 in 1:length(uM_compound)){
      a2<-which(INPUT.1[,"uM_compound"] == uM_compound[k2])
      INPUT.2 <- INPUT.1[a2,] 
      # how many replicates?
      
      Replicate = sort(unique(INPUT.2$Replicate))
      
      b1<-which(INPUT.2[,"Replicate"] == Replicate[1])
      b2<-which(INPUT.2[,"Replicate"] == Replicate[2])
      b3<-which(INPUT.2[,"Replicate"] == Replicate[3])
      INPUT.3.1 <- INPUT.2[b1,] 
      INPUT.3.2 <- INPUT.2[b2,] 
      INPUT.3.3 <- INPUT.2[b3,] 
      
      p1<-match(INPUT.3.1$Well_ID,INPUT.3.2$Well_ID)
      p2<-match(INPUT.3.1$Well_ID,INPUT.3.3$Well_ID)
      
      c1<-sum(INPUT.3.1$Accession_number == INPUT.3.2$Accession_number[p1])
      c2<-sum(INPUT.3.1$Accession_number == INPUT.3.3$Accession_number[p2])
      c3<-sum(INPUT.3.3$Accession_number[p2] == INPUT.3.2$Accession_number[p1])
      
      if(c1 != nrow(INPUT.3.1) | c2 != nrow(INPUT.3.1) | c3 != nrow(INPUT.3.1)){
        print("error")
        print(k1)
        print(k2)
        print(k3)
      }
      ccc1<-sum(INPUT.3.1$Well_ID == INPUT.3.2$Well_ID[p1])
      ccc2<-sum(INPUT.3.1$Well_ID == INPUT.3.3$Well_ID[p2])
      ccc3<-sum(INPUT.3.3$Well_ID[p2] == INPUT.3.2$Well_ID[p1])
      if(ccc1 != nrow(INPUT.3.1) | ccc2 != nrow(INPUT.3.1) | ccc3 != nrow(INPUT.3.1)){
        print("error")
        print(k1)
        print(k2)
        print(k3)
      }
      
      OUTPUT<-INPUT.3.1[,c(1,4,5,8,11,20,21,22,23,24,25)]
      OUTPUT$Pct_Replicate.1<-INPUT.3.1$Pct_normalization_control
      OUTPUT$Pct_Replicate.2<-INPUT.3.2$Pct_normalization_control[p1]
      OUTPUT$Pct_Replicate.3<-INPUT.3.3$Pct_normalization_control[p2]
      
      
      sortTable<-function(INPUT.3){
        Barcode<-names(table(as.character(as.matrix(INPUT.3$Barcode))))
        origSignal<-c()
        normSignals<-c()      
        for(k4 in 1:length(Barcode)){
          a4<-which(INPUT.3[,"Barcode"] == Barcode[k4])
          b4<-as.character(as.matrix(INPUT.3$Well_ID[a4]))
          
          INPUT.4<-INPUT.3[a4,] 
          
          # You have to change this part to select which normalization you would like to do
          ### Pct_normalization
          ### Z_transformation with Mock
          ### Z_transformation with UNI
          
          # Also, due to mislabeling by James, column and Group name should be manually modified
#           e1.0<-which(INPUT.4$Group=="Kif11" | INPUT.4$Group=="KIF11")

          e1.0<-which(INPUT.4$Group=="Kif11-Control")          
          #           e1.1<-which(INPUT.4$Column == 23 & INPUT.4$Group == "Mock") # Mock
          #           e1.2<-which(INPUT.4$Column == 24 & INPUT.4$Group == "Blank") # Blank             
          e1.1<-which(INPUT.4$Group == "Universal Control") # Mock
          e1.4<-which(INPUT.4$Column == 23) # Mock
          e1.2<-which(INPUT.4$Column == 24) # Blank             
          e1.3<-which(INPUT.4$Group == "Spike In") # Blank  
          
          if(length(e1.0)==0 & length(e1.1)==0 & length(e1.2)==0){
            next
          }
          
          cellLines.1<-as.character(unique(INPUT.4$Cell_line))
          uM_compound.1<-as.character(unique(INPUT.4$uM_compound))
          Replicate.1<-as.character(unique(INPUT.4$Replicate))
          
          png(paste(outfile,"/Distribution/Cellline_",cellLines.1,"_compound_",uM_compound.1,"_replicate_",Replicate.1,"_plateBarcode_",Barcode[k4],".png",sep = ""),width = 480)
          stripchart(list(Kif11=INPUT.4$Signal[e1.0],Mock = INPUT.4$Signal[e1.4],Blank = INPUT.4$Signal[e1.2],UNI = INPUT.4$Signal[e1.1]), vertical=T, pch=19,method="jitter", cex=1,col=c(2:5))
          dev.off()
          
          
          NE1<-(INPUT.4$Signal - mean(INPUT.4$Signal[e1.2]))        
          normSignal<-((NE1-mean(NE1[e1.1])) / sd(NE1[e1.1]))  
          
          
          origSignal<-rbind(origSignal,
                            c(cellLines.1,uM_compound.1,Replicate.1,Barcode[k4],
                              mean(INPUT.4$Signal[e1.0]),
                              mean(INPUT.4$Signal[e1.1]),
                              mean(INPUT.4$Signal[e1.2]),
                              mean(INPUT.4$Signal[e1.3]),
                              mean(INPUT.4$Signal[e1.4]),
                              sd(INPUT.4$Signal[e1.0]),
                              sd(INPUT.4$Signal[e1.1]),
                              sd(INPUT.4$Signal[e1.2]),
                              sd(INPUT.4$Signal[e1.3]),
                              sd(INPUT.4$Signal[e1.4])))
          
          
          names(normSignal)<-b4
          d1<-match(names(normSignal),OUTPUT$Well_ID)
          normSignals[d1]<-normSignal
        }
        return(list(N=normSignals,O=origSignal))
      }
      
      
      out.1<-sortTable(INPUT.3.1)
      out.2<-sortTable(INPUT.3.2[p1,])
      out.3<-sortTable(INPUT.3.3[p2,])
      
      OUTPUT$norm_Replicate.1<-out.1$N
      OUTPUT$norm_Replicate.2<-out.2$N
      OUTPUT$norm_Replicate.3<-out.3$N
      
      V1<-rbind(out.1$O,out.2$O,out.3$O)
      colnames(V1)<-c("cellLines","uM_compound","Replicate","Barcode","Mean.Kif11","Mean.UNI","Mean.Blank","Mean.SpikeIn","Mean.Mock","SD.Kif11","SD.UNI","SD.Blank","SD.SpikeIn","SD.Mock")
      write.table(OUTPUT,file = paste(outfile,"/CellLines_",cellLines[k1],"_uM_compound_",uM_compound[k2],".txt",sep = ""),row.names=F,col.names=T,quote=F,sep = "\t")
      write.table(V1,file = paste(outfile,"/Distribution/Distribution_CellLines_",cellLines[k1],"_uM_compound_",uM_compound[k2],".txt",sep = ""),row.names=F,col.names=T,quote=F,sep = "\t")
      
    }
    
  }
  
}
