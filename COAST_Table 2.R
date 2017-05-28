#COAST results ----
#Run this for true_d= 0.2, 0.4, and 0.6, setting lower_bound to 20,  40, or 60 to reproduce Table 2

#Simulation settings
n_sim<-20000 #numbber of simulated studies
true_d<-0.6 #True effect size when sd = 1 (for between) or 1/sqrt2) (for within)
true_sd<-1/sqrt(2) #Set True standard deviation (set to 1 to make D cohen's d.
lower_bound<-40 #After how any participants do you start to look? Frick calls this the lower bound 

#Sequential design settings NHST (alpha and # looks)
n<-1000 #total number of datapoints (per condition) you are willing to collect
sided <- 2 #Set whether tests are 1-sided or 2-sided
paired_test<- TRUE

pmat<-readRDS(paste("pmat_d_",true_d,"_n_",n,"_n_sim_",n_sim,"_sided_",sided,"_paired_",paired_test,".Rdata", sep=""))

#COAST requires a minimum number of observations - Frick recommends not testing first 20. So these are dropped.  
pmat[,1:lower_bound-1]<-NA

#Find at which N first p<0.01 and first p > 0.36 is observed 
#Results are not exactly the same. Not sure why, maybe because Frick writes:
#"The value of t needed for p = .36 was computed from a linear average" - so he did not use p-values directly?
#I'm assuming the difference in simulations is because of the t-value cut-off he used. 
n_first_sig<-apply(pmat<0.01, 1, function(x) match(TRUE, x))
n_first_nonsig<-apply(pmat>0.36, 1, function(x) match(TRUE, x))
#if NA no p-value below 0.01 or above 0.36) code as n
n_first_sig[is.na(n_first_sig)] <- n
n_first_nonsig[is.na(n_first_nonsig)] <- n

coast_n_stop<-pmin(n_first_sig,n_first_nonsig) #use pmin function to take lowest value (earliest stop)
coast_mean_n<-mean(coast_n_stop) #calculate mean sample size for COAST procedure

coast_stop_high_before_low<-sum(n_first_sig>n_first_nonsig)/n_sim
coast_stop_low_before_high<-sum(n_first_sig<n_first_nonsig)/n_sim
coast_inconclusive<-1-coast_stop_high_before_low-coast_stop_low_before_high
coast_mean_n #How quickly would we stop?
coast_stop_high_before_low #How often is p < 0.01 passed first
coast_stop_low_before_high #How often is p > 0.36 passed first
coast_inconclusive #How often does the p-value 0.01 < p < 0.036?

cat("COAST would stop after",coast_mean_n,"participants, have",coast_stop_low_before_high,"power, and a Type 2 error rate of",coast_stop_high_before_low,"with",coast_inconclusive,"inconclusive results.")


result_all <- as.data.frame(matrix(0, ncol = 1, nrow = 1))

result_all[1,1]<-n_sim
result_all[1,2]<-true_d
result_all[1,3]<-true_sd
result_all[1,4]<-lower_bound
result_all[1,5]<-n
result_all[1,6]<-sided
result_all[1,7]<-paired_test
result_all[1,8]<-coast_mean_n #How quickly would we stop?
result_all[1,9]<-coast_stop_high_before_low #How often is p < 0.01 passed first
result_all[1,10]<-coast_stop_low_before_high #How often is p > 0.36 passed first
result_all[1,11]<-coast_inconclusive #How often does the p-value 0.01 < p < 0.036?

colnames(result_all)<-c("n_sim","cohens_d", "sd", "lower_bound", "n_max", "sided", "paired_test", "coast_mean_n", "coast_stop_high_before_low", "coast_stop_low_before_high", "coast_inconclusive")

write.table(result_all, "COAST_results.csv", sep = ",", append = T, row.names = FALSE, col.names = FALSE)