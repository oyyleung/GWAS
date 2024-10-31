
"
Created for specific use at MS & PD Tissue Bank Database, Imperial College London
Copyright (c) [2023] [Yeung Yeung Leung]
"
#######################################################################
# function for running and saving coloc results
# assumed input file processed by preproCOLOC.py 
# [fastGWA summary statistics output & GTEX eQTL results format]
#######################################################################

library(coloc)
runCOLOC <- function(eN,#eQTL sample N number
                     sdY,#standard deviation of phenotype
                     inFileName,#filename of csv file containing eQTL and GWAS statistics (from prepData4COLOC.py)
                     outFileName){ 
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


#######################################################################
# run colo
#######################################################################

GTEXdir = "E://GTEXv7/"
filename='MSf3MAF_Het_imputed_rsid_'
tissues=c(
'Brain_Substantia_nigra',    
'Brain_Spinal_cord_cervical_c-1',      
 'Brain_Anterior_cingulate_cortex_BA24',
 'Brain_Caudate_basal_ganglia',
 'Brain_Cerebellar_Hemisphere',
 'Brain_Cerebellum',
 'Brain_Cortex',
 'Brain_Frontal_Cortex_BA9',
 'Brain_Hippocampus',
 'Brain_Hypothalamus',
 'Brain_Nucleus_accumbens_basal_ganglia',
 'Brain_Putamen_basal_ganglia',
 'Cells_EBV-transformed_lymphocytes',
 'Muscle_Skeletal',
 'Nerve_Tibial',
 'Pituitary',
 'Thyroid',
 'Whole_Blood'
)
features= c('nSFG', 'nCG', 'nThal', 'nOcc', 'nPons', 'nCbDN')
# load file containing sample size of eQTL study
tissueN=read.csv(paste0(GTEXdir,'tissueN.csv'), header=TRUE)
# load file containing standard deviation of phenotype feature
sdYs=read.csv(paste0(GTEXdir,'sdYs.csv'), header=TRUE)

# loop through phentypes and tissues
for (pheno in features){
  for (tissue in tissues){
    eN<-tissueN$N[tissues==tissue]
    sdY<-sdYs$sdYs[features==pheno]
    for (i in c(1:22)){
      inFileName = paste(GTEXdir,filename,pheno,'_',tissue,'_chr',i,'.txt.gz', sep='')
      outFileName = paste(GTEXdir,'coloc_',filename,pheno,'_',tissue,'_chr',i,'.csv', sep='')
      if (file.exists(outFileName)) next
      runCOLOC(eN,#eQTL sample N number
              sdY,#standard deviation of phenotype
              inFileName,#filename of csv file containing eQTL and GWAS statistics (from prepData4COLOC.py)
              outFileName)
    }
  }
}