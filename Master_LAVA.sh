##########################################
#      Master file for LAVA pipeline
#      Date: 24/01/2025
##########################################


##########################################
# Creating a new environment
##########################################

# conda create -n LAVA_2024
## Package Plan ##
#  environment location: /home/dkramarenk/.conda/envs/LAVA_2024
#
# To activate this environment, use
conda activate LAVA_2024
# To deactivate an active environment, use
#     $ conda deactivate

# Installing module
conda search -c bioconda plink
conda install -c bioconda plink=1.90b6.21
conda search -c bioconda bcftools
conda install -c bioconda bcftools=1.17
conda search -c bioconda r
conda install -c bioconda r-base=3.6.0
conda install -c bioconda r-base=4.2.0
conda install -c bioconda r-base=4.4

conda search -c bioconda plink2
conda install -c bioconda plink2=2.00a2.3
# newer plink2 available but some conflicts
conda search -c bioconda zstd
conda install -c bioconda zstd=1.5.5 

install.packages("BiocManager")
BiocManager::install("snpStats")
install.packages("remotes")
remotes::install_github("josefin-werme/LAVA")
# /gpfs/home2/dkramarenk/.conda/envs/LAVA_2024/lib/R/library/LAVA

export WD_PROJECT="/home/dkramarenk/projects/LAVA/DCM_HCM"

# /home/dkramarenk/projects/LAVA/DCM_HCM -> ${WD_PROJECT}
##########################################
#      0.  Downloading the data          #
##########################################
cd ${WD_PROJECT}

# HCM - build 37
cp /archive/expcard/Lisa_Projects/HCM-GWAS/HCM-2022/hcm.gwama.250621.gel.nofin.txt .

# DCM
cp /archive/expcard/DCM_GWAS_2023/data/summary_stats/MTAG/Garnier_Meder_Amsterdam_FinnGen_UKB_MGB__DCM__META1_chr1_22_MAF0.005.tsv.gz .

#  1000G
cp /home/dkramarenk/projects/Project_DCM_GWAS/POP/1000_genomes_ALEX/1KGautos_fil_biallelic* .

##########################################
#      How file should look like        #
##########################################

# # SumStats
# CHR BP SNPID_UKB A1 A2 MAF BETA SE STAT P NMISS INFO_UKB
# 1 113418061 rs2360008 G A 0.2261 -0.01392 0.01282 -1.086 0.2775 379831 1
# 1 113418144 rs9429498 T C 0.001218 0.2492 0.1538 1.62 0.1051 379682 0.951417
# 1 113418415 rs1237670 G A 0.2242 -0.01563 0.01286 -1.215 0.2244 379470 0.998238

# # Input
# phenotype	cases	controls	filename
# asthma	44301	341521	vignettes/data/asthma.sumstats.txt
# bmi	NA	NA	vignettes/data/bmi.sumstats.txt
# depression	170756	329443	vignettes/data/depression.sumstats.txt
# diabetes	18483	366937	vignettes/data/diabetes.sumstats.txt
# hypothyroidism	13043	231847	vignettes/data/hypothyroidism.sumstats.txt
# neuro	NA	NA	vignettes/data/neuro.sumstats.txt
# rheuma	14361	43923	vignettes/data/rheuma.sumstats.txt

##########################################
#     1.1    Files processing            #
##########################################

# Same as FLAMES
sbatch ${WD_PROJECT}/scripts/DCM_MTAG_process.sh
sbatch ${WD_PROJECT}/scripts/HCM_MTAG_process.sh
sbatch ${WD_PROJECT}/scripts/CC_GWAS_process.sh
sbatch ${WD_PROJECT}/scripts/CC_MTAG_process.sh
sbatch ${WD_PROJECT}/scripts/DCM_GWAS_process.sh
sbatch ${WD_PROJECT}/scripts/HCM_GWAS_process.sh

####################################################################################

##########################################
# 1.2  Check so all SNPs on the same ref #
##########################################
sbatch ${WD_PROJECT}/scripts/compare_SNPs.sh

##########################################
#   2          Make infofile             #
##########################################

cp data/sumstats/DCM_stats_LAVA.txt  data/sumstats/dcm.sumstats.txt
cp data/sumstats/HCM_stats_LAVA.txt  data/sumstats/hcm.sumstats.txt

nano data/input.info.txt
##########################################
phenotype	cases	controls	filename
DCM     37083   37083   /home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/DCM_GWAS_37_exclMYBPC3reg.txt
HCM     21708   21708   /home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/HCM_GWAS_37_exclMYBPC3reg.txt

# DCM	9365	946368	/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/DCM_GWAS_37_exclMYBPC3reg.txt
# HCM	5900	68359	/home/dkramarenk/projects/LAVA/DCM_HCM/data/sumst_processed/HCM_GWAS_37_exclMYBPC3reg.txt

# 9,365 DCM cases and 946,368
# 5,900 HCM cases, 68,359 controls
##########################################
scp /Users/drkramarenko/Library/CloudStorage/OneDrive-Personal/Computation/LAVA/WCPG_LAVA_tutorial/cluster_setup/data/refdat.tar.gz  dkramarenk@snellius.surf.nl:/home/dkramarenk/projects/LAVA/DCM_HCM/data

diff DCM_stats_LAVA.txt  data/sumst_processed/DCM_GWAS_37_exclMYBPC3reg.txt > dcm_changes.txt
##########################################
#      Make ref file (with chr:bp)       #
##########################################

plink2 --bfile g1000_eur.maf01 \
--set-all-var-ids @:# \
--make-bed \
--out g1000_eur.maf01.ids


cd /home/dkramarenk/projects/LAVA/DCM_HCM



for i in {0..24}; do
  start=$((i * 100 + 1))
  end=$((start + 99))
  if [ $end -gt 2495 ]; then
    end=2495
  fi
  sbatch --export=ALL,START=$start,END=$end lava_par.job
done
