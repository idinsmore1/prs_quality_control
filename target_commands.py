import pandas as pd
import numpy as np
import gwas_commands as gc
import subprocess
import os
import gzip

# Plink commands as subprocesses
def backup_plink_files(prefix, threads, indir = './', backup_dir = './', in_path = False):
    if in_path:
        cmd = 'plink '
    else:
        cmd = './plink '
    cmd += f'--bfile {indir}{prefix} --make-bed --out {backup_dir}{prefix}.bk --threads {threads}'
    print('Writing backup of plink files...')
    subprocess.check_call(cmd,shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    print('Done. \n')

def make_var_id(bim, chr, bp, ref, alt):
    bim.dropna(inplace = True)
    bim = gc.upper_alleles(bim, ref, alt)
    var_id = bim[f'{chr}'].astype(str) + ':' + bim[f'{bp}'].astype(str) + ':' + bim[f'{ref}'].astype(str) + ':' + bim[f'{alt}'].astype(str)
    return var_id

def update_bim_var_id(prefix,indir = './', outdir = './', write_map = True):
    bim = pd.read_table(f'{indir}{prefix}.bim', header = None,dtype={0:'str'})
    original_ids = list(bim.iloc[:,1])
    bim.columns = ['CHR','SNP','CM','BP','A1','A2']

    print('Creating var_id in bim file...')
    bim['SNP'] = make_var_id(bim,'CHR','BP','A2','A1')
    print('Done. \n')
    var_ids = list(bim['SNP'])
    print('Writing a bim file with var_id to match GWAS...')
    bim.to_csv(f'{indir}{prefix}.bim', header = False, index = False, sep = '\t')
    print('Done. \n')

    if write_map:
        id_map = pd.DataFrame({'original_id':original_ids, 'var_id':var_ids})
        id_map.to_csv(f'{outdir}{prefix}.original_to_var_id_map.tsv',index=False,sep='\t')
        print(f'Original ID to Var_ID map written to {outdir}{prefix}.original_id_to_var_id.tsv. \n')

def gwas_qc_on_target(prefix, threads, maf = 0.01, hwe = '1e-6', geno = 0.01, mind = 0.01, indir = './',outdir = './', in_path = False):
    if in_path:
        cmd = 'plink '
    else:
        cmd = './plink '
    
    cmd += f'--bfile {indir}{prefix} \
    --maf {str(maf)} \
    --hwe {hwe} \
    --geno {geno} \
    --mind {str(mind)} \
    --write-snplist \
    --make-just-fam \
    --threads {threads} \
    --out {outdir}{prefix}.QC'

    print(f'Performing filtering on plink files using \n \
        MAF : {maf} \n \
        HWE : {hwe} \n \
        SNP Missingness Rate : {geno} \n \
        Genotype Missingness Rate : {mind} \n')
    subprocess.check_call(cmd,shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    print('Done. \n')

def get_snps_in_LD(prefix, threads, indir = './', outdir = './', in_path = False):
    if in_path:
        cmd = 'plink '
    else:
        cmd = './plink '
    cmd += f'--bfile {indir}{prefix} \
    --keep {outdir}{prefix}.QC.fam \
    --extract {outdir}{prefix}.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --threads {threads} \
    --out {outdir}{prefix}.QC'
    print('Computing SNPs in LD...')
    subprocess.check_call(cmd,shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    print(f'SNPs not in LD written to {outdir}{prefix}.QC.prune.in. \n')

def backup_and_symlink_bim(prefix, indir = './', backup_dir = './'):
    backup_cmd = f'mv {indir}{prefix}.bim {backup_dir}{prefix}.bim.backup'
    symlink_cmd = f'ln -s {indir}{prefix}.QC.adj.bim {prefix}.bim'
    # make backup of bim file 
    print(f'Writing backup of bim file to {backup_dir}{prefix}.bim.backup \n')
    subprocess.check_call(backup_cmd, shell = True, executable='/bin/bash')

    print(f'Symlinking {prefix}.QC.adj.bim to {prefix}.bim \n')
    subprocess.check_call(symlink_cmd, shell = True, executable='/bin/bash')

def get_het_rates(prefix, threads, snp_set,indir = './', outdir = './', in_path = False):
    if in_path:
        cmd = 'plink '
    else:
        cmd = './plink '
    cmd += f'--bfile {indir}{prefix} \
    --extract {snp_set} \
    --keep {outdir}{prefix}.QC.fam \
    --het \
    --threads {threads} \
    --out {outdir}{prefix}.QC'
    print('Computing heterozygosity rates...')
    subprocess.check_call(cmd,shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    print(f'Heterozygosity rates written to {outdir}{prefix}.QC')

# Function for removing outliers based on heterozygosity rates
def remove_het_outliers(prefix, outdir = './'):
    dat = pd.read_csv(f'{outdir}{prefix}.QC.het', sep = '\s+')
    m = dat['F'].mean()
    s = dat['F'].std()

    print('Removing outliers outside 3 standard deviations from mean F stat...')
    # filter outliers that are outside 3 SD of the mean
    valid_dat = dat.loc[(dat['F'] <= m + 3 * s) & (dat['F'] >= m-3*s),:]
    print('Done. \n')

    print(f'Writing non-outliers to {outdir}{prefix}.valid.sample...')
    valid_dat[['FID','IID']].to_csv(f'{outdir}{prefix}.valid.sample',sep='\t', index = False)
    print('Done. \n')

# Functions to compare the plink and gwas files
def create_qc_df(bim, snp_set):
    qc = pd.read_table(f'{snp_set}', header = None)
    qc = qc.merge(bim,right_on='SNP',left_on=0)
    return qc

def create_info_df(bim, qc, gwas_df):
    print(f'Merging bim and gwas dataframes...')
    info = bim.merge(gwas_df, on = ['CHR','BP','SNP'])
    print('Done. \n')

    print("Filtering out SNPs that didn't pass quality control..." )
    info = info[info['SNP'].isin(qc['SNP'])]
    print('Done. \n')
    return info

# Functions to compare snps across gwas and bim
def get_matching_df(info):
    info_match = info.loc[(info['A1'] == info['B.A1']) & (info['A2'] == info['B.A2']),:]
    return info_match

def complement(x):
    if x == 'A':
        return('T')
    elif x == 'T':
        return('A')
    elif x == 'C':
        return('G')
    elif x == 'G':
        return('C')
    else: 
        return(np.nan)

def create_complementary_allele_columns(info):
    info['C.A1'] = info['B.A1'].apply(lambda x: complement(x))
    info['C.A2'] = info['B.A2'].apply(lambda x: complement(x))
    return info

def get_complementary_alleles(info,bim):
    info_complement = info.loc[(info['A1'] == info['C.A1']) & (info['A2'] == info['C.A2']),:]
    complement_snps = bim['SNP'].isin(info_complement['SNP'])
    return complement_snps, info_complement

def update_complementary_alleles(bim, complement_snps):
    bim.loc[complement_snps,'B.A1'] = bim.loc[complement_snps,'B.A1'].apply(lambda x: complement(x))
    bim.loc[complement_snps,'B.A2'] = bim.loc[complement_snps,'B.A2'].apply(lambda x: complement(x))
    return bim

def get_recode_alleles(bim,info):
    info_recode = info_recode = info.loc[(info['A1'] == info['B.A2']) & (info['A2'] == info['B.A1']),:]
    recode_snps = bim['SNP'].isin(info_recode['SNP'])
    return recode_snps, info_recode

def update_recode_alleles(bim, recode_snps):
    temp = bim.loc[recode_snps,'B.A1']
    bim.loc[recode_snps,'B.A1'] = bim.loc[recode_snps,'B.A2']
    bim.loc[recode_snps,'B.A2'] = temp
    return bim

def get_recode_and_compl_alleles(bim, info):
    info_crecode = info.loc[(info['A1'] == info['C.A2']) & (info['A2'] == info['C.A1']),:]
    cre_snps = bim['SNP'].isin(info_crecode['SNP'])
    return cre_snps, info_crecode

def update_crecode_alleles(bim, cre_snps):
    #Update the snps
    temp = bim.loc[cre_snps,'B.A1']
    bim.loc[cre_snps,'B.A1'] = bim.loc[cre_snps,'B.A2'].apply(lambda x: complement(x))
    bim.loc[cre_snps,'B.A2'] = temp.apply(lambda x: complement(x))
    return bim
    # Write bim after this in main function 

#### MISMATCH
def get_mismatching_alleles(prefix, bim, info_match, info_complement, info_recode, info_crecode, outdir = './'):
    cond1 = bim['SNP'].isin(info_match['SNP'])
    cond2 = bim['SNP'].isin(info_complement['SNP'])
    cond3 = bim['SNP'].isin(info_recode['SNP'])
    cond4 = bim['SNP'].isin(info_crecode['SNP'])
    all_cond = cond1 | cond2 | cond3 | cond4
    mismatch = bim.loc[~all_cond,'SNP']
    print(f'Writing mismatching snps to {outdir}{prefix}.mismatch...')
    mismatch.to_csv(f'{outdir}{prefix}.mismatch', index=False, header = False)
    print('Done. \n')

def adjust_target_alleles(prefix, gwas_df, snp_set, indir = './', outdir = './'):
    # read in files
    print(f'Reading in {indir}{prefix}.bim... \n')
    bim = pd.read_table(f'{indir}{prefix}.bim', header = None,dtype={0:'str'})

    # create backup qc dataframe
    bim.columns = ['CHR','SNP','CM','BP','B.A1','B.A2']
    bim = gc.upper_alleles(bim,'B.A1','B.A2')
    qc = create_qc_df(bim, snp_set)

    gwas_df = gwas_df.astype({'CHR': str,'BP':str,'SNP':str})
    bim = bim.astype({'CHR': str,'BP':str,'SNP':str})
    qc = qc.astype({'CHR': str,'BP':str,'SNP':str})

    # merge bim and gwas to get info dataframe
    info = create_info_df(bim, qc, gwas_df)
    info_match = get_matching_df(info)

    # complementary alleles
    print('Updating complementary alleles...')
    info = create_complementary_allele_columns(info)
    complementary_snps, info_complement  = get_complementary_alleles(info,bim)
    bim = update_complementary_alleles(bim,complementary_snps)
    print('Done. \n')
    # recode alleles
    print('Updating alleles that need recoding...')
    recode_snps, info_recode = get_recode_alleles(bim, info)
    bim = update_recode_alleles(bim,recode_snps)
    print('Done. \n')

    # crecode alleles
    print('Updating alleles that need complement and recoding...')
    cre_snps, info_crecode = get_recode_and_compl_alleles(bim, info)
    bim = update_crecode_alleles(bim,cre_snps)
    print('Done. \n')

    # Write out mismatching alleles
    get_mismatching_alleles(prefix, bim, info_match, info_complement, info_recode, info_crecode, outdir)
    # write updated bim file out
    print(f'Writing updated bim file to {indir}{prefix}.QC.adj.bim')
    bim.to_csv(f'{indir}{prefix}.QC.adj.bim', index = False, header = False, sep = '\t')
    print('Done. \n')



def sex_check(prefix, threads, indir = './', outdir='./',in_path = False):
    if in_path:
        cmd = 'plink '
    else:
        cmd = './plink '
    
    cmd += f'--bfile {indir}{prefix} \
    --extract {outdir}{prefix}.QC.prune.in \
    --keep {outdir}{prefix}.valid.sample \
    --check-sex --threads {threads} \
    --out {outdir}{prefix}.QC'

    print('Performing sex check on plink data...')
    subprocess.check_call(cmd,shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    print('Done. \n')

    print('Reading in samples from {outdir}{prefix}.valid.sample...')
    valid = pd.read_table(f'{prefix}.valid.sample')

    print('Reading in results from plink sexcheck')
    dat = pd.read_table(f'{prefix}.QC.sexcheck',sep='\s+')

    print('Comparing sexcheck results to sex column in fam file...')
    for i in range(len(dat.index)):
        if (dat.loc[i,'SNPSEX'] == 1) & (dat.loc[i,'F'] > 0.8):
            dat.loc[i,'STATUS'] = 'OK'
        elif (dat.loc[i,'SNPSEX'] == 2) & (dat.loc[i,'F'] < 0.2):
            dat.loc[i,'STATUS'] = 'OK'

    valid = dat.loc[(dat['STATUS'] == 'OK') & (dat['FID'].isin(valid['FID'])),:]
    # write the valid file
    print(f'Writing updated valid sample to {outdir}{prefix}.QC.valid...')
    valid[['FID','IID']].to_csv(f'{outdir}{prefix}.QC.valid', index = False, header = False, sep = '\t')
    print('Done. \n')

def get_related_patients(prefix, threads, indir = './', outdir = './', in_path = False, sex_check = False):
    if in_path:
        cmd = 'plink '
    else:
        cmd = './plink '
    
    if sex_check:
        cmd += f'--bfile {indir}{prefix} \
                --extract {outdir}{prefix}.QC.prune.in \
                --keep {outdir}{prefix}.QC.valid \
                --rel-cutoff 0.125 \
                --threads {threads} \
                --out {outdir}{prefix}.QC'
    else:
        cmd += f'--bfile {indir}{prefix} \
                --extract {outdir}{prefix}.QC.prune.in \
                --keep {outdir}{prefix}.valid.sample \
                --rel-cutoff 0.125 \
                --threads {threads} \
                --out {outdir}{prefix}.QC'
    print('Calculating relationship matrix...')
    subprocess.check_call(cmd,shell=True, executable='/bin/bash')
    print(f'List of unrelated patients written to {outdir}{prefix}.QC.rel.id \n')

def generate_final_files(prefix, threads, unrelated_ids, snp_set, indir = './', outdir = './',in_path = False):
    if in_path:
        cmd = 'plink '
    else:
        cmd = './plink '
    cmd += f'--bfile {indir}{prefix} \
    --make-bed \
    --threads {threads} \
    --keep {unrelated_ids} \
    --out {outdir}{prefix}.QC \
    --extract {snp_set} \
    --exclude {outdir}{prefix}.mismatch'
    print('Creating quality controlled plink files...')
    subprocess.check_call(cmd,shell=True, executable='/bin/bash')
    print(f'Plink files for PRS written to {outdir}{prefix}.QC.bed/bim/fam')