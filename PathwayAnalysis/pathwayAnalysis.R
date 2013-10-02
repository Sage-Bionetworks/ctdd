# work flow or Demo using MSigDB or Graphite DB
require(synapseClient)
require(multicore)
source("~/Sage-Analysis-Pipeline/PathwayAnalysis/myPathwayAnalysis.R")

pathwayAnalysis<-function(synID=NULL,pathwayName = NULL,Reference = NULL,Test.method = c("FET","GSEA"),cores = 1){
  
  if(is.null(synID)==1){
    error("Please select which database you would like to use in your analysis: MSigDB(syn1681370) or Graphite(syn2135025)")
    break;
  }
  if(is.null(pathwayName)==1){
    error("Please select which pathway you would like to use in your analysis: KEGG, BioCarta, ...")
    break;
  }
  if(is.null(Reference)==1){
    error("Please provide your reference list to be tested in pathway analysis: list of genes or two column reference list (gene name and its statistics")
    break;
  }
  
  pathwayName = toupper(pathwayName)
  
  if(synID == "syn1681370"){
    MSIGDB<-synGet("syn2227979")
    load(MSIGDB@filePath)
    
    if(is.element(pathwayName,"BIOCARTA")){
      allPathways <- MSigDB$C2.CP.BIOCARTA
    }
    if(is.element(pathwayName,"KEGG")){
      allPathways <- MSigDB$C2.CP.KEGG
    }
    if(is.element(pathwayName,"REACTOME")){
      allPathways <- MSigDB$C2.CP.REACTOME
    }
    if(is.element(pathwayName,"GO_BP")){
      allPathways <- MSigDB$C5.GO_BP
    }
    if(is.element(pathwayName,"GO_CC")){
      allPathways <- MSigDB$C5.GO_CC
    }
    if(is.element(pathwayName,"GO_MF")){
      allPathways <- MSigDB$C5.GO_MF
    }
  }
  if(synID == "syn2135029"){
    GRAPHITE<-synGet("syn2135029",load = TRUE)
    
    if(is.element(pathwayName,"BIOCARTA")){
      allPathways <- GRAPHITE@objects$BIOCARTA
    }
    if(is.element(pathwayName,"KEGG")){
      allPathways <- GRAPHITE@objects$KEGG
    }
    if(is.element(pathwayName,"REACTOME")){
      allPathways <- GRAPHITE@objects$REACTOME
    }
    if(is.element(pathwayName,"NCI")){
      allPathways <- GRAPHITE@objects$NCI
    }    
  }
    
  # FET test
  if(is.element(Test.method,"FET")){
    allGenes <-c()
    for (i in 1:length(allPathways)){
      allGenes<-union(allGenes,allPathways[[i]])
    }      
    testSet = Reference
    AllGenes<-union(testSet,allGenes)
    
    pathwayTest <-function(x){    
      curPathwayGenes <- allPathways[[x]]        
      pathwayAnalysis<-myPathwayAnalysis$new()
      pathwayAnalysis$fet(AllGenes,curPathwayGenes,testSet)
      return(pathwayAnalysis)      
    }    
    results<-mclapply(1:length(allPathways), function(x)pathwayTest(x),mc.cores= cores)     
    
  }else{ # GSEA test
    referenceSet<-sort(Reference, decreasing =TRUE, index.return =TRUE)
    pathwayTest <-function(x){    
      curPathwayGenes <- allPathways[[x]]        
      pathwayAnalysis<-myPathwayAnalysis$new()
      pathwayAnalysis$gsea(referenceSet$x,curPathwayGenes,np=1000,w =1)
      return(pathwayAnalysis)      
    }     
    results<-mclapply(1:length(allPathways), function(x)pathwayTest(x),mc.cores= cores)               
  }
  
  names(results)<-names(allPathways)  
  return(results)  
}
