Imputed Genotypes
===================

This folder contains plink2 format files and variant summary statistics (MAF + 
INFO scores) for the three versions of imputation from genotypes that have been
performed by UK Biobank. Data has been converted to plink format as it is much 
faster to work with than the BGEN format the imputed genotype data is dispensed 
in, which can be found under 'Bulk/Imputed Genotypes/'. 

The HRCUK10K (GRCh37)/ folder contains the version of the imputed genotype data
that has historically and most widely used - it is the same as the imputed 
genotype data we used to have stored on CSD3.

More recently, UK Biobank have also imputed their genotype data to the TopMed
and Genomics England reference panels, which are aligned to the genome build 
GRCh38 rather than the much older GRCh37 genome build used by the HRC/UK10K 
imputed genotype data.

In addition to being aligned to a much more modern genome build, the TopMed and
Genomics England imputations provide better coverage of rarer variants - 
Genomics England imputation is the most accurate for people of white british 
and south asian ancestries, while the TopMed imputation is more accurate for 
other ancestry groups.

For further details, see the technical report comparing the three different 
imputation panels: https://biobank.ndph.ox.ac.uk/ukb/ukb/docs/GEL_imputation.pdf

Importantly, please also be aware that whole genome sequencing data is now also
available for all UK Biobank participants; plink2 format files for the WGS data
have been provided by UK Biobank under 
'Bulk/DRAGEN WGS/DRAGEN population level WGS variants, PLINK format [500k release]/'

Note that sex chromosome imputation has not been performed for the Genomics 
England imputation panel, and only for chromosome X for the TopMed imputation
panel, whereas HRC/UK10K have imputed chrosomes X and XY (with one PAR1 SNP in 
the XY chromosome file).

Variant summary statistics the Genomics England imputation panel have been 
computed de-novo from the pfiles using plink2 --freq as unlike the HRC/UK10K 
and TopMed imputation this information has not been precomputed by UK Biobank. 
The output of plink2 --freq have been post-processed to match the information 
provided in the varstats files for the HRC/UK10K imputation.
