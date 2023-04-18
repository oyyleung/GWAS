
"
Created for specific use at MS & PD Tissue Bank Database, Imperial College London
Copyright (c) [2023] [Yeung Yeung Leung]
"
runCOLOC <- function(eN,#eQTL sample N number
                     sdY,#standard deviation of phenotype
                     inFileName,#filename of csv file containing eQTL and GWAS statistics (from prepData4COLOC.py)
                     outFileName){
  
  library(coloc)
  data <- read.csv(file=inFileName,header=T)
  data=data[!duplicated(data$SNP),] 
  genes<-unique(data$gene_id)
  print(length(genes))
  summary<- data.frame()
  for (gene in genes){
    try({
      input<-data[data$gene_id==gene,]
      input<-input[input$maf>0,]    
      result <- coloc.abf(dataset1=list(beta=input$BETA,
                                        varbeta='^'(input$SE,2),
                                        pvalues=input$P,
                                        snp=input$SNP,
                                        position=input$POS,
                                        type="quant",sdY=sdY),
                          dataset2=list(beta=input$slope,
                                        varbeta='^'(input$slope_se,2),
                                        snp=input$SNP,
                                        pvalues=input$pval_nominal, 
                                        type="quant", N=eN),
                          MAF=input$maf)
      summary_ <- data.frame(t(result$summary))
      summary_$gene_id<-gene
      summary<- rbind(summary, summary_)
      rm(result, summary_) 
    })                
  }
  write.csv(summary, outFileName, na = "", row.names = F, sep=",")

}