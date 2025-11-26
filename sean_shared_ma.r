######## Running shared-effects meta-analysis ############

################################
# DCM vs HCM; constain rg at 1 
# (fixed-effects meta-analysis taking into account sample overlap)
################################

mtag/mtag.py \
--sumstats sum_stats/DCM_for_ccMTAG.txt,sum_stats/HCM_for_ccMTAG.txt \
--out ./DCM__HCM__ccMTAG_rg1constrained \
--beta_name beta \
--snp_name snpid \
--se_name se \
--z_name z \
--n_name n \
--eaf_name freq \
--a1_name a1 \
--a2_name a2 \
--p_name pval \
--stream_stdout \
--equal_h2 --perfect_gencov

R
library(data.table)
dat <- fread('DCM__HCM__ccMTAG_rg1constrained_mtag_meta.txt', stringsAsFactors=F, data.table=F)
dat <- dat[order(dat$mtag_pval), ]
### check for SEs

dcm <- fread(paste0('sum_stats/Jurgens_DCM_sumstats.tsv.gz'), stringsAsFactors=F, data.table=F, select=c("rsID", "EA", "NEA", "BETA", "SE", "P"))
colnames(dcm) <- c("SNP", "EA", "NEA", "BETA_DCM", "SE_DCM", "P_DCM")
dcm1 <- dcm2 <- dcm
dcm1$SNP_EA <- paste0(dcm1$SNP, "_", dcm1$EA)
dcm2$SNP_EA <- paste0(dcm2$SNP, "_", dcm2$NEA)
dcm2$BETA_DCM <- -1*dcm2$BETA_DCM 
dcm <- rbind(dcm1, dcm2)
dcm <- dcm[,c("SNP_EA", "BETA_DCM", "SE_DCM", "P_DCM")]

hcm <- fread('sum_stats/Tadros_HCM_sumstats.txt', stringsAsFactors=F, data.table=F, select=c("rsid", "effect_allele", "noneffect_allele", "beta", "se", "pvalue"))
colnames(hcm) <- c("SNP", "EA", "NEA", "BETA_HCM", "SE_HCM", "P_HCM")
hcm1 <- hcm2 <- hcm
hcm1$SNP_EA <- paste0(hcm1$SNP, "_", hcm1$EA)
hcm2$SNP_EA <- paste0(hcm2$SNP, "_", hcm2$NEA)
hcm2$BETA_HCM <- -1*hcm2$BETA_HCM 
hcm <- rbind(hcm1, hcm2)
hcm <- hcm[,c("SNP_EA", "BETA_HCM", "SE_HCM", "P_HCM")]

dat$SNP_EA <- paste0(dat$SNP, "_", dat$A1)
dat <- merge(dat, dcm, by="SNP_EA", all.x=T, all.y=F)
dat <- merge(dat, hcm, by="SNP_EA", all.x=T, all.y=F)

length(which(is.na(dat$BETA_DCM)))
length(which(is.na(dat$BETA_HCM)))
nrow(dat)

dat$Directions <- paste0(ifelse(dat$BETA_DCM>=0, "+", "-"), ifelse(dat$BETA_HCM>=0, "+", "-"))
summary(lm(dat$SE_DCM ~ dat$SE_HCM))
###Estimate Std. Error t value Pr(>|t|)
###(Intercept) 7.312e-05  4.236e-06   17.26   <2e-16 ***
###  dat$SE_HCM  7.833e-01  1.011e-04 7747.47   <2e-16 ***
  
summary(dat$SE_DCM / dat$SE_HCM)
###Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
###0.3257  0.7494  0.7808  0.7842  0.8134 11.9739

#### Fix some issues with Standard error ratios (outliers have disproportiante power in one GWAS over the other)
rm <- which((dat$SE_DCM / dat$SE_HCM)<0.65 | (dat$SE_DCM / dat$SE_HCM)>1.2)
if(length(rm)>0){dat <- dat[-rm, ]}
nrow(dat)
###[1] 4843297

dat <- dat[order(dat$mtag_pval), ]
nrow(dat[dat$mtag_pval<5e-8 &
           dat$P_DCM<0.001 & dat$P_HCM<0.001 & 
           sign(dat$BETA_DCM)==sign(dat$BETA_HCM), ])
#8
dat[dat$mtag_pval<5e-8 &
           dat$P_DCM<0.001 & dat$P_HCM<0.001 & 
           sign(dat$BETA_DCM)==sign(dat$BETA_HCM), ]

dat <- dat[order(dat$CHR, dat$BP), ]

###write.table(dat, file='DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCfull.txt', col.names=T, row.names=F, quote=F, sep='\t')
###write.table(dat[,c("SNP", "CHR", "BP", "A1", "A2", "meta_freq", "mtag_beta", "mtag_se", "mtag_pval")], file='DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCshort.txt', col.names=T, row.names=F, quote=F, sep='\t')
###system("gzip DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCfull.txt")
###system("gzip DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCshort.txt")


###########################################
# Run secondary random-effects meta-analysis
# (using the meta R package)
############################################

library(meta)
dat <- fread('DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCfull.txt.gz', stringsAsFactors = F, data.table=F)
nrow(dat)
#[1] 4843297
dat[dat$mtag_pval<0.0001, c("metagen_res")] <- apply(X=dat[dat$mtag_pval<0.0001, c("BETA_DCM", "SE_DCM", "BETA_HCM", "SE_HCM")], MARGIN=1, FUN=function(line){TEs <- c(line[1], line[3]); seTEs <- c(line[2], line[4]); met <- metagen(TE=TEs, seTE=seTEs); Q <- met$Q; p.Q <- met$pval.Q; beta <- met$TE.random; se <- met$seTE.random; p <- met$pval.random; return(paste0(c(Q, p.Q, beta, se, p), collapse=";"))})
dat2 <- dat[dat$mtag_pval<0.0001, ]
dat2 <- tidyr::separate(dat2, col="metagen_res", into=c("Q", "Q.pvalue", "random_beta", "random_se", "random_pval"), sep=";")
dat <- dat[dat$mtag_pval>=0.0001, -(which(colnames(dat)=="metagen_res"))]
dat[,c("Q", "Q.pvalue", "random_beta", "random_se", "random_pval")] <- NA
dat <- rbind(dat2, dat)
class(dat$Q) <- class(dat$Q.pvalue) <- class(dat$random_beta) <- class(dat$random_se) <- class(dat$random_pval) <- "numeric"
head(dat)
nrow(dat)
##[1] 4843297
dat[dat$mtag_pval<5e-8 & dat$random_pval<1e-4, ]

# Take most conservative P-value from fixed-effects and random-effects approach
dat[, c('harmonized_beta', 'harmonized_se', 'harmonized_pval')] <- dat[,c('mtag_beta', 'mtag_se', 'mtag_pval')]
dat[which(!is.na(dat$mtag_pval) & dat$random_pval>dat$mtag_pval), c('harmonized_beta', 'harmonized_se', 'harmonized_pval')] <- dat[which(!is.na(dat$mtag_pval) & dat$random_pval>dat$mtag_pval), c('random_beta', 'random_se', 'random_pval')]

dcm <- fread(paste0('sum_stats/Jurgens_DCM_sumstats.tsv.gz'), stringsAsFactors=F, data.table=F, select=c("rsID", "CHRBP_B37"))
dcm <- dcm[-which(duplicated(dcm$rsID)), ]
dat_new <- merge(dcm, dat, by.y="SNP", by.x="rsID", all=F)
nrow(dat_new)
head(dat_new)
dat_new <- tidyr::separate(dat_new, col="CHRBP_B37", into=c("CHR_B37", "BP_B37"), sep=":")

dat_new <- dat_new[order(dat_new$CHR_B37, dat_new$BP_B37), ]
write.table(dat_new, file='DCM__HCM__ccMTAG_rg1constrained_mtagandrandom_meta_QCfull.txt', col.names=T, row.names=F, quote=F, sep='\t')
write.table(dat_new[,c("rsID", "CHR_B37", "BP_B37", "A1", "A2", "meta_freq", "mtag_beta", "mtag_se", "mtag_pval")], file='DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCshort.txt', col.names=T, row.names=F, quote=F, sep='\t')
write.table(dat_new[,c("rsID", "CHR_B37", "BP_B37", "A1", "A2", "meta_freq", "harmonized_beta", "harmonized_se", "harmonized_pval")], file='DCM__HCM__ccMTAG_rg1constrained_harmonized_meta_QCshort.txt', col.names=T, row.names=F, quote=F, sep='\t')
system("gzip DCM__HCM__ccMTAG_rg1constrained_mtagandrandom_meta_QCfull.txt")
system("gzip DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCshort.txt")
system("gzip DCM__HCM__ccMTAG_rg1constrained_harmonized_meta_QCshort.txt")

dat_new <- fread("DCM__HCM__ccMTAG_rg1constrained_mtagandrandom_meta_QCfull.txt.gz", stringsAsFactors=F, data.table=F)
####rsID CHR_B37    BP_B37       SNP_EA CHR       BP A1 A2 meta_freq
####1 rs10875231       1 100000012 rs10875231_T   1 99534456  T  G 0.2536205
####2  rs6678176       1 100000827  rs6678176_T   1 99535271  T  C 0.3162560
####3 rs78286437       1 100000843 rs78286437_C   1 99535287  C  T 0.0648655
####4 rs76909621       1 100001201 rs76909621_T   1 99535645  T  G 0.0929655
####5 rs78642210       1 100002490 rs78642210_T   1 99536934  T  C 0.0638320
####6 rs77140576       1 100002713 rs77140576_T   1 99537157  T  C 0.0972180
####mtag_beta    mtag_se      mtag_z mtag_pval BETA_DCM SE_DCM   P_DCM
####1 -0.0192938773 0.01603962 -1.20288864 0.2290194  -0.0327 0.0194 0.09126
####2 -0.0151021517 0.01466993 -1.02946304 0.3032621  -0.0244 0.0175 0.16300
####3 -0.0009875843 0.02858159 -0.03455316 0.9724361  -0.0020 0.0334 0.95200
####4 -0.0078659641 0.02442266 -0.32207653 0.7473947  -0.0154 0.0287 0.59050
####5  0.0031976286 0.02876149  0.11117741 0.9114757   0.0030 0.0336 0.92830
####6 -0.0109460637 0.02387229 -0.45852593 0.6465746  -0.0139 0.0283 0.62270
####BETA_HCM   SE_HCM    P_HCM Directions  Q Q.pvalue random_beta random_se
####1  0.002520 0.025759 0.922057         -+ NA       NA          NA        NA
####2  0.001363 0.024159 0.954988         -+ NA       NA          NA        NA
####3  0.000773 0.047747 0.987072         -+ NA       NA          NA        NA
####4  0.005916 0.040788 0.884665         -+ NA       NA          NA        NA
####5  0.003445 0.048049 0.942825         ++ NA       NA          NA        NA
####6 -0.005558 0.039401 0.887806         -- NA       NA          NA        NA
####random_pval harmonized_beta harmonized_se harmonized_pval
####1          NA   -0.0192938773    0.01603962       0.2290194
####2          NA   -0.0151021517    0.01466993       0.3032621
####3          NA   -0.0009875843    0.02858159       0.9724361
####4          NA   -0.0078659641    0.02442266       0.7473947
####5          NA    0.0031976286    0.02876149       0.9114757
####6          NA   -0.0109460637    0.02387229       0.6465746
## excl MYBPC3 region
dat_new <- dat_new[-(which(dat_new$CHR_B37==11 & dat_new$BP_B37>=30000000 & dat_new$BP_B37<=80000000)), ]
nrow(dat_new)
##[1] 4769682

write.table(dat_new, file='DCM__HCM__ccMTAG_rg1constrained_mtagandrandom_meta_QCfull_exclMYBPC3reg.txt', col.names=T, row.names=F, quote=F, sep='\t')
write.table(dat_new[,c("rsID", "CHR_B37", "BP_B37", "A1", "A2", "meta_freq", "mtag_beta", "mtag_se", "mtag_pval")], file='DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCshort_exclMYBPC3reg.txt', col.names=T, row.names=F, quote=F, sep='\t')
write.table(dat_new[,c("rsID", "CHR_B37", "BP_B37", "A1", "A2", "meta_freq", "harmonized_beta", "harmonized_se", "harmonized_pval")], file='DCM__HCM__ccMTAG_rg1constrained_harmonized_meta_QCshort_exclMYBPC3reg.txt', col.names=T, row.names=F, quote=F, sep='\t')
system("gzip DCM__HCM__ccMTAG_rg1constrained_mtagandrandom_meta_QCfull_exclMYBPC3reg.txt")
system("gzip DCM__HCM__ccMTAG_rg1constrained_mtag_meta_QCshort_exclMYBPC3reg.txt")
system("gzip DCM__HCM__ccMTAG_rg1constrained_harmonized_meta_QCshort_exclMYBPC3reg.txt")