
"
Created for specific use at MS & PD Tissue Bank Database, Imperial College London
Copyright (c) [2023] [Yeung Yeung Leung]
"
runCOLOCevalSNP <- function(
    ci=0.99,#confidence interval for the selected SNP being the shared variant for both GWAS and eWTL significance
    PPH4thres=0.5,#PPH4 threshold for gene selection
    eN,#eQTL sample N number
    sdY,#standard deviation of phenotype
    inFileName,#filename of csv file containing eQTL and GWAS statistics (from prepData4COLOC.py)
    colocFileName,#saved COLOC result (csv file) from runCOLOC.R
    outFileName
){
  
  library(coloc)
  coloc.res=read.csv(colocFileName)
  coloc.res<-coloc.res[coloc.res$PP.H4.abf>PPH4thres,] 
  summary<- data.frame()
  
  
  data <- read.csv(file=inFileName,header=T)
  #print(dim(data))
  data=data[!duplicated(data$SNP),]
  data<-data[data$maf>0,] 
  #print(dim(data))
  for (gene in coloc.res$gene_id){
    try({
      
      
      input<-data[data$gene_id==gene,]
      input <- input[order(input$POS),]
      #print(dim(input))
      
      
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
      
      
      o <- order(result$results$SNP.PP.H4,decreasing=TRUE)
      cs <- cumsum(result$results$SNP.PP.H4[o])
      w <- which(cs > ci)[1]
      summary_<-result$results[o,][1:w,]
      summary_<-dplyr::rename(summary_, SNP=snp)
      #    input<-input[input$P<0.05,]
      summary_<- merge(summary_,input, by="SNP", all=FALSE)
      summary_$gene_id<-gene
      summary_$gene_name<-coloc.res$gene_name[coloc.res$gene_id==gene]
      summary_$chr<-coloc.res$chr[coloc.res$gene_id==gene]
      
      #outputFile = paste(GTEXdir,'susie_',filename,pheno,norm,'_',tissue,'_',gene_name,'.csv', sep='')
      #write.csv(summary_, outputFile, na = "", row.names = F, sep=",")
      summary<- rbind(summary, summary_)
      rm(result, summary_) 
    })
  }
  
  write.csv(summary, outputFile, na = "", row.names = F, sep=",")
  
}
