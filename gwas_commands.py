import pandas as pd
import numpy as np
# GWAS quality control functions

## uppercase alleles for easy matching
def upper_alleles(gwas_df, ref, alt):
    gwas_df[f'{ref}'] = gwas_df[f'{ref}'].str.upper()
    gwas_df[f'{alt}'] = gwas_df[f'{alt}'].str.upper()
    return gwas_df

# create var_id that matches our database
def make_var_id(gwas_df, chr, bp, ref, alt):
    gwas_df.dropna(inplace = True)
    gwas_df = upper_alleles(gwas_df, ref, alt)
    var_id = gwas_df[f'{chr}'].astype(str) + ':' + gwas_df[f'{bp}'].astype(str) + ':' + gwas_df[f'{ref}'].astype(str) + ':' + gwas_df[f'{alt}'].astype(str)
    return var_id

# filter maf based on chosen value
def filter_maf(gwas_df,maf_col,maf_val):
    gwas_df[f'{maf_col}'] = pd.to_numeric(gwas_df[f'{maf_col}'])
    maf_filtered = gwas_df[gwas_df[f'{maf_col}'] > maf_val]
    return maf_filtered

# filter out duplicated snps
def filter_duplicates(gwas_df,snp_col):
    gwas_df = gwas_df[gwas_df[f'{snp_col}'].duplicated() == False]
    return gwas_df

# filter out ambiguous snps
def filter_ambiguous(gwas_df, ref, alt):
    cond1 = ((gwas_df[f'{ref}'] == 'A') & (gwas_df[f'{alt}'] == 'T'))
    cond2 = ((gwas_df[f'{ref}'] == 'T') & (gwas_df[f'{alt}'] == 'A'))
    cond3 = ((gwas_df[f'{ref}'] == 'G') & (gwas_df[f'{alt}'] == 'C'))
    cond4 = ((gwas_df[f'{ref}'] == 'C') & (gwas_df[f'{alt}'] == 'G'))
    all_cond = cond1 | cond2 | cond3 | cond4
    gwas_df = gwas_df.loc[~all_cond,:]
    return gwas_df

# Adjust columns for downstream comparison with plink files
def adjust_cols(gwas_df, snp_col, chr, bp, ref, alt, maf_col):
    gwas_df = gwas_df.astype(
        {f'{snp_col}':'str',f'{chr}':'str', f'{bp}':'int32',f'{ref}':'str',f'{alt}':'str',f'{maf_col}':'float64'}
        )
    gwas_df.rename(
        columns={f'{snp_col}':'SNP',f'{chr}':'CHR', f'{bp}':'BP',f'{ref}':'A2',f'{alt}':'A1',f'{maf_col}':'MAF'},
        inplace=True
    )
    return gwas_df

# do all steps at once
def gwas_qc(gwas_file, snp_col, chr, bp, ref, alt, maf_col, maf_val):
    print('Reading in GWAS file...')
    gwas_df = pd.read_table(f'{gwas_file}',sep = '\s+')
    print('Done. \n')

    print('Creating Var_ID to match database markers...')
    gwas_df[f'{snp_col}'] = make_var_id(gwas_df, chr, bp, ref, alt)
    print('Done. \n')

    print(f'Filtering out snps with MAF < {maf_val}...')
    gwas_df = filter_maf(gwas_df,maf_col,maf_val)
    print('Done. \n')

    print('Filtering out duplicate snps...')
    gwas_df = filter_duplicates(gwas_df,snp_col)
    print('Done. \n')

    print('Filtering out ambiguous snps...')
    gwas_df = filter_ambiguous(gwas_df,ref,alt)
    print('Done. \n')

    print('Adjusting columns for downstream quality control...')
    gwas_df = adjust_cols(gwas_df, snp_col, chr, bp, ref, alt, maf_col)
    print('Done. \n')

    return gwas_df
