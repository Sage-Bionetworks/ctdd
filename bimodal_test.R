library(gdata)
A0<-read.xls("~/HNSCC/CisplatinModifierSummaryDataFinalForExport.xlsx")
inputScreenData1<-A0$Test.Mean[which(A0$Compound == "Cisplatin")]
inputScreenData2<-A0$Cognate.Mean[which(A0$Compound == "Cisplatin")]
filename<-"~/TEMP/test_folder/HNSCC_vehicle_bimodal.png"


distributionPlot_bimodal<-function(filename,inputScreenData){
  
  require(mixtools)
  png(filename)
  Y<-Y0.vehicle<-InputScreenData
  wait0.vehicle <- normalmixEM(as.numeric(Y), lambda = .5, mu = c(10,200), maxit=100000)
    
  plot(wait0.vehicle, which = 2,density = TRUE, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8, main2 = "", xlab2 = "Percent Normalized Cell Viability")
  curve(dnorm(x, mean=mean(Y0.vehicle), sd=sd(Y0.vehicle)), add=TRUE,col = "blue")
  dev.off()
  
}

# multi-variate analysis is needed implemented
# two distribution is considered simultaneously to compute P-value to identify hits 
