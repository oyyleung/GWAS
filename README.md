# GWAS
Functions and code examples made for the analysis performed in the manuscript "Genome wide screening of causal genes underlying Multiple Sclerosis clinical and pathological severity"


## Usage Example

```python
import SNPassoc
import GENEassoc

"""
SNP level association
"""

# Input files in plink format
plinkDir='./plink_mac_20201019/'#plink package directory
gtcaDir='./gcta_1.93.2beta_mac/'#fastGWAS package directory
fileDir='./' # datafile directory
genotype=fileDir+'MSf3MAF' # genotype [QC filtered] file prefix [plink format ie .bed, .bim .fam]
sampleID=fileDir+'MSHet.sample' # sample ID to be used in analysis
covariate=fileDir+'MS.PCASexAgePmi' # covariates to be accounted for in analysis
GRM=fileDir+'MSspGRM'# file suffix for sparse genetic relationship matrix [fastGWA format ie .grm.id .grm.sp]
outDir='./out/'# directory for output results

GWAS_MSneuropath=SNPassoc.SNPassoc(genotype, sampleID, covariate, plinkDir, gtcaDir)
for pheno in ['nSFG', 'nCG', 'nThal', 'nOcc', 'nPons', 'nCbDN']:
    print(pheno)    
    phenotype=fileDir+'MS.'+pheno    
    SNPassoc_outName=outDir+'MS_'+pheno
    GWAS_MSneuropath.SNP_level_association(phenotype, GRM, SNPassoc_outName)
```

## Contact

mail@leungyeungyeung.com
