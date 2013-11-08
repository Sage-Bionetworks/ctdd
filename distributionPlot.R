P<-read.xls("/Volumes/ijang//HNSCC//CisplatinModifierSummaryDataFinalForExport.xlsx")

P1<-P[which(P$Compound == "Cisplatin"),]
P2<-P[which(P$Compound == "Vehicle"),]

A.1 <- subset(P1, P1[,"Group"] == "PLK1-Control")
A.2 <- subset(P1, P1[,"Group"] == "Library")
A.3 <- subset(P1, P1[,"Group"] == "Kif11")
A.4 <- subset(P1, P1[,"Group"] == "Spike In")
A.5 <- subset(P1, P1[,"Group"] == "UNI")
A.6 <- subset(P1, P1[,"Group"] == "Mock")
A.7 <- subset(P1, P1[,"Group"] == "Treated Mock")

NAME<-c("PLK1","Library","Kif11","Spike In","UNI","Mock","Treated Mock")
COLOR<-c("blue","brown","cyan","green","black","red","magenta")

B<-list(A.1$Cognate.Mean,
        A.2$Cognate.Mean,
        A.3$Cognate.Mean,
        A.4$Cognate.Mean,
        A.5$Cognate.Mean,
        A.6$Cognate.Mean,
        A.7$Cognate.Mean)
names(B)<-NAME

png("~/Documents/Data//HNSCC/Cisplatin/eachDistribution.png",width = 640, height = 480)
par(mar=c(5,5,5,2))
stripchart(B, vertical=T, pch=15,method="jitter", cex=0.35,col=COLOR,xlab = "Group", ylab = "Viability (%)",cex.axis = 1.25,cex.lab = 1.5,add=T)
abline(h = 100, lty = 2,col = "darkgrey")
rect(1.89, min(A.2$Cognate.Mean)-0.1, 2.11, max(A.2$Cognate.Mean)+0.1, density = NULL, angle = 45,
     col = NA, border = COLOR[2], lty = par("lty"), lwd = par("lwd"))
dev.off()

require(vioplot)
violinPlot<-function(A.1,maxY = 0.75,...){  
  plot(c(1,length(A.1)),c(0,maxY),type = "n",xlab = "", ylab = "Correlation",xaxt = "n",...)
  axis(1,at = c(1:length(A.1)),labels = F)
  for(k in 1:length(A.1)){
    AA.1<-A.1[[k]]
    
    vioplot(AA.1[!is.na(AA.1)],at = k,col = "gold",add = T,drawRect = F,border = "gold")
  }
}

violinPlot(B,270)
stripchart(B, vertical=T, pch=15,method="jitter", cex=0.35,col=COLOR,xlab = "Group", ylab = "Viability (%)",cex.axis = 1.25,cex.lab = 1.5,add=T)

# T-test for synergy in library only
C.1<-cbind(A.2$Test.1,A.2$Test.2,A.2$Test.3)
C.2<-cbind(A.2$Cognate.1,A.2$Cognate.2,A.2$Cognate.3)

pval<-c()
tstats<-c()
for(k in 1:nrow(C.1)){
  aa<-t.test(C.1[k,],C.2[k,])
  pval<-c(pval,aa$p.value)
  tstats<-c(tstats,aa$statistic)  
}

e.1<-which(pval<=0.005)
e.2<-which(tstats<0)
e.3<-which(A.3$Gene.Symbol!="EMPTY")
e<-intersect(e.3,intersect(e.1,e.2))
length(e)


X1<-sort(A.2$Test.Mean[e],decreasing = T, index.return = T)

a1<-as.character(as.matrix(A.2$Gene.Symbol[e[X1$ix]]))



png("~/Documents/Data//HNSCC/Cisplatin/temp_new05_pval_005_Cisplatin_fig2.png",width = 640, height = 640)
par(pty="s")
plot(P2$Cognate.Mean,P2$Test.Mean,pch = 15,cex =0.45,col = "black",xlim = c(0,340),ylim = c(0,340),ylab = "Viability (%) : Cisplatin treated", xlab = "Viability (%) : Untreated",cex.lab = 1.5)
points(A.2$Cognate.Mean,A.2$Test.Mean,pch = 15,cex =0.45,col = COLOR[2],xlim = c(0,340),ylim = c(0,340))
points(A.3$Cognate.Mean,A.3$Test.Mean,pch = 15,cex =0.45,col = COLOR[3],xlim = c(0,340),ylim = c(0,340))
points(A.4$Cognate.Mean,A.4$Test.Mean,pch = 15,cex =0.45,col = COLOR[4],xlim = c(0,340),ylim = c(0,340))
points(A.5$Cognate.Mean,A.5$Test.Mean,pch = 15,cex =0.45,col = COLOR[5],xlim = c(0,340),ylim = c(0,340))
points(A.1$Cognate.Mean,A.1$Test.Mean,pch = 15,cex =0.45,col = COLOR[1],xlim = c(0,340),ylim = c(0,340))
legend("topleft",NAME,pch = 15,col = COLOR,bty="n",cex = 1.25)



x1<-  A.2$Cognate.Mean[e[X1$ix]]
y1<-  A.2$Test.Mean[e[X1$ix]]
points(x1,y1,pch = 16,col = "darkgreen",cex = 1.25)
x<-300
y<-seq(345,1,len = length(e))
for(k in 1:length(a1)){
  text(x,y[k],a1[k],pos=4,srt = 0,font = 3,col = "darkgreen",cex = 1)
  arrows(x=x, y = y[k], x1=x1[k], y1=y1[k], length = 0.15, angle = 10,code = 2, col = "darkgrey", lty = par("lty"),lwd = par("lwd"))
}

dev.off()
