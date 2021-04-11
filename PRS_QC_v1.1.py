import gwas_commands as gc
import target_commands as tc
import pandas as pd
import numpy as np
import argparse
import os

# Establish command line arguments and their meanings
parser = argparse.ArgumentParser()
parser.add_argument('--gwas', help = 'Base GWAS Filepath, in gzip format.')
parser.add_argument('--prefix', help = 'Prefix of Plink files (excluding .bed/.bim/.fam).')
parser.add_argument('--pheno', help = 'Phenotype of GWAS.')
parser.add_argument('--input_dir', default='./', help = 'Input directory. Default current directory.')
parser.add_argument('--output_dir', default='PRSiceQC_output/', help='Output directory. Default PRSiceQC_output/.')
parser.add_argument('--backup_dir', default='plink_backup/', help='Backup directory. Default plink_backup/.')
parser.add_argument('--threads', default = 16, help = 'Number of threads to use with plink. Default 16')
parser.add_argument('--unrelated_ids', default = 'NA', help = 'List of unrelated patient ids in plink format. Will calculate if not given.')
parser.add_argument('--snp_set', default = 'NA', help = 'List of quality controlled snps in plink format. Will calculate if not given.')
parser.add_argument('--plink_in_path', default=False, action='store_true', help='Use this switch if plink is in your PATH variable. Otherwise plink should be in current directory.')
parser.add_argument('--no_plink_backup', default=True, action='store_false', help = 'Do not make backups of original plink files.')
parser.add_argument('--no_gwas_qc', default=True, action='store_false', help = 'Do not perform GWAS QC.')
parser.add_argument('--maf_val', default = 0.01, help = 'Minor Allele Frequency threshold for filtering datasets. Default 0.01')
parser.add_argument('--hwe', default = '1e-6', help = 'Hardy Weinberg Equilibrium threshold. Default 1e-6')
parser.add_argument('--mind', default = 0.01, help = 'Genotype Missingness Rate Threshold for Target Data. Default 0.01')
parser.add_argument('--geno', default = 0.01, help = 'SNP Missingness threshold. Default 0.01')
parser.add_argument('--snp_id', default = 'SNP', help = 'SNP Column name in GWAS. Default SNP')
parser.add_argument('--ref', default = 'A1', help = 'Reference Allele Column name in GWAS. Default A1')
parser.add_argument('--alt', default = 'A2', help = 'Alt Allele Column name in GWAS. Default A2')
parser.add_argument('--maf_col', default = 'MAF', help = 'MAF Column name in GWAS. Default MAF')
parser.add_argument('--chr', default = 'CHR', help = 'Chromosome Column name in GWAS. Default CHR')
parser.add_argument('--bp', default = 'BP', help = 'Base Position Column name in GWAS. Default BP')
parser.add_argument('--no_gwas_out', default=False, action= 'store_true', help = 'Do not write QC GWAS file.')
parser.add_argument('--no_id_map', default=True, action= 'store_false', help = 'Do not write var_id map')
parser.add_argument('--perform_sex_check', default=False, action= 'store_true', help = 'Perform sex check on plink data (not recommended with MyCode data). Sex column in fam file must not be empty.')
args = parser.parse_args()

# localize parser arguments as variables
gwas = args.gwas
prefix = args.prefix
pheno = args.pheno
input_dir = args.input_dir
output_dir = args.output_dir
backup_dir = args.backup_dir
plink_in_path = args.plink_in_path
plink_backup = args.no_plink_backup
maf_val = args.maf_val
hwe = args.hwe
mind = args.mind
geno = args.geno
snp_id = args.snp_id
ref = args.ref
alt = args.alt
maf_col = args.maf_col
chr = args.chr
bp = args.bp
no_gwas_out = args.no_gwas_out
gwas_qc = args.no_gwas_qc
id_map = args.no_id_map
sex_check = args.perform_sex_check
threads = args.threads
unrelated_ids = args.unrelated_ids
snp_set = args.snp_set

# Make output directories if they do not exist.
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
if not os.path.exists(args.backup_dir):
    os.makedirs(args.backup_dir)

# GWAS QC
if gwas:
    if gwas_qc:
        gwas_df = gc.gwas_qc(gwas, snp_id, chr, bp, ref, alt, maf_col, maf_val)
        if not no_gwas_out:
            if pheno:
                print(f'Writing Quality Controlled GWAS to {output_dir}{pheno}.gwas.QC.gz...')
                gwas_df.to_csv(f'{output_dir}{pheno}_gwas.QC.gz', index = False, header = True, sep = '\t',compression = 'gzip')
                print('Done. \n')
            else:
                print(f'No phenotype given. Please use --pheno flag to specify GWAS phenotype. \n')
                exit()
    else:
        print("Reading in QC'd GWAS file...")
        gwas_df = pd.read_table(f'{gwas}',sep = '\s+')
        print('Done. \n')
else:
    print("NO GWAS file given. Please give filepath to GWAS.")
    exit()

# Target QC
if prefix:
    if plink_backup:
        # Backup plink files
        tc.backup_plink_files(prefix, threads, input_dir, backup_dir,in_path=plink_in_path)

    # Get var_id to match our database
    tc.update_bim_var_id(prefix, input_dir, output_dir,id_map)

    # Perform GWAS level QC on plink files
    tc.gwas_qc_on_target(prefix, threads, maf=maf_val, hwe=hwe, geno=geno, mind=mind, indir =input_dir, outdir=output_dir, in_path = plink_in_path)

    if snp_set == 'NA':
        # Calculate SNPs in LD
        tc.get_snps_in_LD(prefix,threads,input_dir,output_dir,plink_in_path)
        snp_set_file = f'{output_dir}{prefix}.QC.prune.in'
    else:
        snp_set_file = args.snp_set
    # Get Heterozygosity rates and remove outliers
    tc.get_het_rates(prefix, threads, snp_set_file, input_dir,output_dir, plink_in_path)
    tc.remove_het_outliers(prefix, output_dir)

    # Adjust alleles in Plink bim file to match the gwas
    tc.adjust_target_alleles(prefix, gwas_df, snp_set_file, input_dir, output_dir)
    tc.backup_and_symlink_bim(prefix, input_dir,backup_dir)

    if sex_check:
        tc.sex_check(prefix, threads, input_dir,output_dir, plink_in_path)

    if unrelated_ids == 'NA':
        tc.get_related_patients(prefix, threads, input_dir,output_dir, plink_in_path, sex_check)
        unrelated_id_file = f'{output_dir}{prefix}.QC.rel.id'
    else:
        unrelated_id_file = unrelated_ids

    tc.generate_final_files(prefix, threads, unrelated_id_file, snp_set_file, input_dir,output_dir, plink_in_path)
else:
    print('No plink prefix given, use --prefix to specify the name of your files.')
