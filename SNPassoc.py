"""
Created for specific use at MS & PD Tissue Bank Database, Imperial College London
Copyright (c) [2023] [Yeung Yeung Leung]
"""

import os

class SNPassoc:
    def __init__(self, genotype, sampleID, covariate,
                plinkDir='', gtcaDir=''):
        #provide program file location to ensure program being recognized
        self.plinkDir = plinkDir
        self.gtcaDir = gtcaDir
        self.genotype = genotype
        self.sampleID = sampleID
        self.covariate = covariate

        
    def plink_association(self, phenotype, plink_assoc_outName):    
         
        """Perform standard GWAS using Plink (for clumping)"""
        os.system(self.plinkDir+'plink --bfile '+self.genotype+' \
        --keep '+self.sampleID+' \
        --pheno '+phenotype+' \
        --covar '+self.covariate+' \
        --assoc --out '+plink_assoc_outName)
        # Rename log file
        os.rename(plink_assoc_outName+'.log', plink_assoc_outName+'.assoc.log')
        
    def clumping(self, plink_assoc_outName, clumping_outName):    
        # Clumping (remove SNPs with linkage disequilibrium r2 > 0.9)
        os.system(self.plinkDir+'plink --bfile '+self.genotype+' \
        --keep '+self.sampleID+' \
        --clump-p1 1 \
        --clump-r2 0.9 \
        --clump-kb 1000 \
        --clump '+plink_assoc_outName+'.qassoc \
        --clump-snp-field SNP \
        --clump-field P \
        --out '+clumping_outName)  
        # Rename log file
        os.rename(clumping_outName+'.log', clumping_outName+'.clumping.log')
        # Extract clumped SNPs ID (at 3rd output column)
        os.system("awk 'NR!=1{print $3}' "+clumping_outName+".clumped >  \
                  "+clumping_outName+".snp")     
        
    def fastGWAS(self, phenotype, GRM,
                 fastGWAS_outName,clumping_outName=None):       
        # Perform fastGWAS
        if clumping_outName==None: extractSNP=''
        else:  extractSNP='--extract '+clumping_outName+'.snp'              
        os.system(self.gtcaDir+'gcta64 --bfile '+self.genotype+' \
                '+extractSNP+' \
                --keep '+self.sampleID+' \
                --grm-sparse '+GRM+' \
                --qcovar '+self.covariate+' \
                --pheno '+phenotype+' \
                --fastGWA-mlm \
                --save-fastGWA-mlm-residual \
                --out '+fastGWAS_outName)       
        # Rename log file
        os.rename(fastGWAS_outName+'.log', fastGWAS_outName+'.fastGWAS.log')
        
    def SNP_level_association(self, phenotype, GRM, SNPassoc_outName):
        self.plink_association(phenotype, plink_assoc_outName=SNPassoc_outName)
        self.clumping(plink_assoc_outName=SNPassoc_outName,
                 clumping_outName=SNPassoc_outName)
        self.fastGWAS(phenotype, GRM,                  
                        fastGWAS_outName=SNPassoc_outName,
                        clumping_outName=SNPassoc_outName)