"""
Created for specific use at MS & PD Tissue Bank Database, Imperial College London
Copyright (c) [2023] [Yeung Yeung Leung]
"""

import os
import pandas as pd

GTEXdir = "E://GTEXv7/"  # eQTL file directory
GWASdir = "E://GWAS/"  #  gwas file directory
filename = 'MSf3MAF_Het_imputed_rsid_' # gwas filename prefix
# gwas phenotypes
features = ['nSFG', 'nCG', 'nThal', 'nOcc', 'nPons', 'nCbDN']
# eQTL tissues
tissues = ['Brain_Substantia_nigra',    
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
'Whole_Blood']

# extract required eQTL summary data columns, separate files by chromosome
for tissue in tissues:
    dfe = pd.read_csv(GTEXdir+tissue+'.allpairs.txt.gz',
                    usecols=['gene_id', 'variant_id', 'tss_distance',
                    'maf', 'pval_nominal', 'slope','slope_se'],
                    delimiter = '\t',
                    low_memory=False)
    for i in range(1,23):
        dfe.loc[dfe['variant_id'].str.startswith(str(i)+'_'),:].to_csv(
            GTEXdir+tissue+'.allpairs_chr'+str(i)+'.txt.gz',index=None)
        
# curate and merge eQTL & GWAS file     
for pheno in features:
    for tissue in tissues:
        for i in range(1,23):
            dfe = pd.read_csv(GTEXdir+tissue+'.allpairs_chr'+str(i)+'.txt.gz')
            dfref = pd.read_csv(GTEXdir+'GTEx_v7_lookup_table_chr'+str(i)+'.txt.gz')
            dfref= dfref.rename(columns =  {'rs_id_dbSNP147_GRCh37p13':'SNP'})
            dfm = pd.merge(dfe, dfref, on=['variant_id'])
            dfm = dfm.sort_values(by=['pval_nominal'])  
            dfm = dfm.drop_duplicates(subset=['SNP'])
            dfgwas = pd.read_csv(GWASdir+filename+pheno+'.fastGWA', 
                                 delimiter = '\t',
                                 low_memory = False)
            dfgwas = dfgwas.sort_values(by=['P'])  
            dfgwas = dfgwas.drop_duplicates(subset=['SNP'])
            dfm = pd.merge(dfgwas, dfm, on=['SNP'])
            dfm = dfm.drop_duplicates(subset=['SNP'])
            dfm.to_csv(GTEXdir+filename+pheno+'_'+tissue+'_chr'+str(i)+'.txt.gz',index=None)

        df=pd.DataFrame()
        for i in range(1,23):
            dfm = pd.read_csv(GTEXdir+filename+pheno+'_'+tissue+'_chr'+str(i)+'.txt.gz') 
            print(dfm.shape)
            df=df.append(dfm)
        df.to_csv(GTEXdir+filename+pheno+'_'+tissue+'.txt.gz',index=None)
        