"""
Created for specific use at MS & PD Tissue Bank Database, Imperial College London
Copyright (c) [2023] [Yeung Yeung Leung]
"""

import os
import pandas as pd

def prepare_SNPpval(SNP_outName, SNPpval_outName):
    # extract fastGWAS pvalue for MAGMA Gene analysis
        pd.read_csv(SNP_outName+'.fastGWA',
                delimiter = '\s+',
                usecols=["SNP","P","N"],
                low_memory=False).to_csv(SNPpval_outName+'.fastGWA.P.N',index=None, sep=' ')

class GENEassoc:
    def __init__(self, genotype=''):
        self.genotype = genotype
                                 
    def GENE_level_association(self, SNPpval_outName, geneAnno_outName, 
                               MAGMAgene_outName,
                               genemodel='SNPmean.top', 
                               phenotype='',sampleID='', covariate=''):
           
        if genemodel=='SNPmean': 
            # Mean SNP association model, snp-wise=mean(default)  
            MAGMApara = ''
                 
        if genemodel=='SNPmean.top': 
            # Aggregate of snp-wise=mean(default) and snp-wise=top p-value
            MAGMApara = ' --gene-model multi '
        
        if genemodel=='SNPmean.top.PCA': 
            # Aggregate of linreg (gene SNPs PCA), snp-wise=mean(default) and snp-wise=top p-value

            if not os.path.isfile(covariate+'.keepsample'):
                # since --keep option for sample selection is not available in MAGMA
                # create covariate file that only contains sampleID information
                # to ensure that only selected sample go into analysis
                dfco=pd.read_csv(covariate, 
                                 header=None,
                                delimiter = '\s+')
                dfco=dfco.rename(columns={0:"FID",1: "IID"})
                dfsam=pd.read_csv(sampleID,
                                delimiter = '\s+')
                dfco=pd.merge( dfco, dfsam['IID'],on=['IID'])
                dfco.to_csv(covariate+'.keepsample',index=None, header=None, sep=' ')

            MAGMApara = '--pheno file='+phenotype+' \
                        --covar file='+covariate+'.keepsample \
                        --gene-model multi '   
            
        # Run MAGMA Gene analysis     
        os.system('magma --bfile '+self.genotype+MAGMApara+' \
                    --pval '+SNPpval_outName+'.fastGWA.P.N ncol=N \
                    --gene-annot '+geneAnno_outName +'.genes.annot \
                    --out '+MAGMAgene_outName+'.'+genemodel)