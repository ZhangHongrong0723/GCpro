
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")



library(TwoSampleMR)
setwd(path)  
exposure_dat<-read_exposure_data(
  filename = "input.txt",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)
exposure_dat<- exposure_dat[exposure_dat$pval.exposure < 5e-6,]
exposure_dat <- clump_data(exposure_dat,clump_kb = 10000,clump_r2 = 0.001)
write.table(exposure_dat, "exposure.txt",row.names = F,sep = "\t",quote = F)
outcomeID="out-id"     


# exposure_dat <- extract_instruments(exposureID, p1=5e-05, clump=TRUE)

# out_rt<-read_outcome_data(snn="hm_rsid")
# 
# out_rt<-read_outcome_data(
#   filename = "34594039-GCST90018629-EFO_0000178.txt",
#   sep = "\t",
#   # snn="hm_rsid"
#   snp_col = "variant_id",
#   
#   )

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcomeID)

dat <- harmonise_data(exposure_dat, outcome_dat)


outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)


mrResult=mr(dat)

#mr_method_list()$obj
#mrResult=mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode"))

mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)


heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)


pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)


pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

res_single=mr_singlesnp(dat)      
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()


pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()


pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()



