#TO DO: TOST lower bound - 99999 for one-sidedness
# Create option to switch between one-sided and two-sided tests

set.seed(2)

if(!require(BayesFactor)){install.packages('BayesFactor')}
library(BayesFactor) 
if(!require(gsDesign)){install.packages('gsDesign')}
library(gsDesign)

options(scipen=20) #disable scientific notation for numbers

#TOST sim function, adapted from TOSTER package, removed plot and output
TOSTsim<-function (m1, m2, sd1, sd2, n1, n2, low_eqbound_d, high_eqbound_d, 
                   alpha){
  sdpooled <- sqrt((((n1 - 1) * (sd1^2)) + (n2 - 1) * (sd2^2))/((n1 + 
                                                                   n2) - 2))
  low_eqbound <- low_eqbound_d * sdpooled
  high_eqbound <- high_eqbound_d * sdpooled
  degree_f <- n1 + n2 - 2
  t1 <- ((m1 - m2) - low_eqbound)/(sdpooled * sqrt(1/n1 + 
                                                     1/n2))
  p1 <- pt(t1, degree_f, lower.tail = FALSE)
  t2 <- ((m1 - m2) - high_eqbound)/(sdpooled * sqrt(1/n1 + 
                                                      1/n2))
  p2 <- pt(t2, degree_f, lower.tail = TRUE)
  ptost <- max(p1, p2)
  ttost <- ifelse(abs(t1) < abs(t2), t1, t2)
  dif <- (m1 - m2)
  TOSTresults <- data.frame(t1, p1, t2, p2, degree_f)
  colnames(TOSTresults) <- c("t-value 1", "p-value 1", "t-value 2", 
                             "p-value 2", "df")
  invisible(list(TOST_t1 = t1, TOST_p1 = p1, TOST_t2 = t2, 
                 TOST_p2 = p2, TOST_df = degree_f, alpha = alpha, low_eqbound = low_eqbound, 
                 high_eqbound = high_eqbound, low_eqbound_d = low_eqbound_d, 
                 high_eqbound_d = high_eqbound_d, ptost = ptost, ttost = ttost))
}

#Create function to run script below, so it can easily be looped
#anything you want to look, block out setting line below, add variable to function.
analyze_showdown<-function(true_d,n, looks){

#Simulation settings (n_sim, true effect size D)
n_sim<-20000 #numbber of simulated studies
#true_d<-0 #True effect size (Keep SD below to 1, otherwise, this is just mean dif, not d)
true_sd<-1 #Set True standard deviation (set to 1 to make D cohen's d.

#Sequential design settings NHST (alpha and # looks)
#n<-150 #total number of datapoints (per condition) you are willing to collect
#looks<-3 #set number of looks at the data
sided <- 2 #Set whether tests are 1-sided or 2-sided

#Bayes Factor settings
bf_cutoff_high<-6
bf_cutoff_low<-1/6
rscale_bf<-sqrt(2)/2 #Rscale matters little for sample size. Would not vary, but still explore for 0.5 as well

#TOST settings
high_eqbound_d <- 0.5 #Set lower equivalence bound for TOST (postive value!)
low_eqbound_d <- -0.5 #Set upper equivalence bound for TOST (negative value!)

seqdesign <- gsDesign(k=looks, test.type=1, alpha=0.025, sfu="Pocock")
#get alpha boundaries for each look
alphalook<-(pnorm(-abs(seqdesign$upper$bound)))*2
#get N at each look
look_n<-(ceiling(n*seqdesign$timing))

#define some variables to store simulations after each n
x<-numeric(n) #store all data group x
y<-numeric(n) #store all data group y

#define some variables to store simulations after each Look
pmat_looks<-matrix(NA, nrow=n_sim, ncol=looks) #Matrix for p-values at sequential tests
dmat_looks<-matrix(NA, nrow=n_sim, ncol=looks) #Matrix for d at sequential tests
bfmat_looks<-matrix(NA, nrow=n_sim, ncol=looks) #Matrix for BayesFactors
ptostmat_looks<-matrix(NA, nrow=n_sim, ncol=looks) #Matrix for tost

xmat<-readRDS(paste("xmat_d_",0,"_true_sd_",true_sd,"_n_",1000,"_n_sim_",n_sim,".Rdata", sep=""))
ymat<-readRDS(paste("ymat_d_",true_d,"_true_sd_",true_sd,"_n_",1000,"_n_sim_",n_sim,".Rdata", sep=""))


#### RUN SIMULATION ----

# create progress bar
pb <- winProgressBar(title = "progress bar", min = 0, max = n_sim, width = 300)
#Run simulation - create 1 participant at a time, store all data in a list, so I can do tests after any number of participants I want
for (j in 1:n_sim){
  setWinProgressBar(pb, j, title=paste( round(j/n_sim*100, 2), "% done"))
  x<-xmat[j,]
  y<-ymat[j,]
  #Calculate initial 10 datapoints before performing first t-test, cannot do test on n = 1 and small samples vary widely.
  #Perform tests at each look, store p-value
  #Note if main interest is in some nuber of looks, not each n, tweaking code below to calc p and bf etc is much more efficient
  for (i in 1:looks){
    mean_x<-mean(x[1:look_n[i]])  
    mean_y<-mean(y[1:look_n[i]])
    sd_x<-sd(x[1:look_n[i]])
    sd_y<-sd(y[1:look_n[i]])
    if (sided == 2){
      z<-t.test(x[1:look_n[i]],y[1:look_n[i]], alternative = "two.sided", var.equal=TRUE, alpha = alphalook[i]) #perform the t-test
      bfmat_looks[j,i] <- exp(ttest.tstat(z$statistic, n1=length(x[1:look_n[i]]), n2=length(y[1:look_n[i]]), rscale = rscale_bf)[['bf']])
      ptostmat_looks[j,i] <- TOSTsim(m1=mean_x, m2=mean_y, sd1=sd_x,sd2=sd_y,n1=length(x[1:look_n[i]]),n2=length(y[1:look_n[i]]),low_eqbound_d=low_eqbound_d,high_eqbound_d=high_eqbound_d, alpha=alphalook[i])$ptost
      pmat_looks[j,i]<-z$p.value #store all p-values for all simulations
    }
    if (sided == 1){
      z<-t.test(x[1:look_n[i]],y[1:look_n[i]], alternative = "less", var.equal=TRUE, alpha = alphalook[i]) #perform the t-test
      bfmat_looks[j,i] <- exp(ttest.tstat(z$statistic, n1=length(x[1:look_n[i]]), n2=length(y[1:look_n[i]]), nullInterval = c(0, -Inf), rscale = rscale_bf)[['bf']])
      ptostmat_looks[j,i] <- TOSTsim(m1=mean_x, m2=mean_y, sd1=sd_x,sd2=sd_y,n1=length(x[1:look_n[i]]),n2=length(y[1:look_n[i]]),low_eqbound_d=low_eqbound_d,high_eqbound_d=99999, alpha=alphalook[i])$ptost
      pmat_looks[j,i]<-z$p.value #store all p-values for all simulations
    }
    obs_d<-(mean_y-mean_x)/(sqrt((((length(x[1:look_n[i]]) - 1)*((sd_x^2))) + (length(y[1:look_n[i]]) - 1)*((sd_y^2)))/((length(x[1:look_n[i]])+length(y[1:look_n[i]])-2))))
    dmat_looks[j,i]<-obs_d #store all effect sizes for all simulations
  }
}
#close progress bar
close(pb)


#### INTERPRET SIMULATION RESULTS ----

result_all <- as.data.frame(matrix(0, ncol = 22, nrow = 1))

result_all[1,1]<-true_d
result_all[1,2]<-looks
result_all[1,3]<-n_sim
result_all[1,4]<-n
result_all[1,5]<-alphalook[1]
result_all[1,6]<-bf_cutoff_high
result_all[1,7]<-bf_cutoff_low
result_all[1,8]<-rscale_bf
result_all[1,9]<-high_eqbound_d
result_all[1,10]<-low_eqbound_d
result_all[1,11]<-alphalook[1]

##FREQUENTIST NHST

#Find at which N first significant p-value is observed based on looks
n_first_sig_looks<-apply(pmat_looks<alphalook[1], 1, function(x) match(TRUE, x))
result_all[1,12]<-(n_sim-sum(is.na(n_first_sig_looks)))/n_sim #power for sequential analyses based on looks
n_first_sig_looks[is.na(n_first_sig_looks)] <- looks
result_all[1,13]<-mean(n_first_sig_looks*look_n[1]) #Average sample size for sequential analyses based on Max N

#BAYES FACTORS

#BF stops when BF > X or BF < Y. 
bf_n_stop<-apply(bfmat_looks<bf_cutoff_low|bfmat_looks>bf_cutoff_high, 1, function(x) match(TRUE, x))
bf_n_stop[is.na(bf_n_stop)] <- looks
bf_n_stop<-mean(bf_n_stop*look_n[1])
result_all[1,14]<-bf_n_stop

#Find at which N first BF>X is observed based on looks 
bf_first_high_looks<-apply(bfmat_looks>bf_cutoff_low, 1, function(x) match(TRUE, x))
bf_first_high_looks[is.na(bf_first_high_looks)] <- looks
result_all[1,16]<-mean(bf_first_high_looks*look_n[1]) #Average sample size Bayes factor for null after each Look

#Find at which N first BF<X is observed based on looks 
bf_first_low_looks<-apply(bfmat_looks<bf_cutoff_low, 1, function(x) match(TRUE, x))
bf_first_low_looks[is.na(bf_first_low_looks)] <- looks
result_all[1,17]<-mean(bf_first_low_looks*look_n[1]) #Average sample size Bayes factor for null after each Look

#Count how many BF stop for lower or higher bound for each look
bf_first_low_looks<-apply(bfmat_looks<bf_cutoff_low, 1, function(x) match(TRUE, x))
bf_first_low_looks[is.na(bf_first_low_looks)] <- looks+1 #replace NA with 1 value higher than max look to count
bf_first_high_looks<-apply(bfmat_looks>bf_cutoff_high, 1, function(x) match(TRUE, x))
bf_first_high_looks[is.na(bf_first_high_looks)] <- looks+1 #replace NA with 1 value higher than max look to count
#How many analyses are inconclusive?
# If we add two matrices, highest values are inconclusive trials (have value 2*looks+1)
inconclusive_bf<-bf_first_low_looks+bf_first_high_looks
inconclusive_bf<-sum(inconclusive_bf==(looks+1)*2)/n_sim
result_all[1,18]<-inconclusive_bf #inconclusive bf based on looks

#Calculate Decision error
#Only a decision error if we cross low BF before high BF (when D <> 0)
bf_stop_high_before_low_looks<-sum(bf_first_low_looks>bf_first_high_looks)/n_sim
bf_stop_low_before_high_looks<-sum(bf_first_low_looks<bf_first_high_looks)/n_sim
result_all[1,19]<-bf_stop_high_before_low_looks
result_all[1,20]<-bf_stop_low_before_high_looks

#### TOST Results ----
#Find at which N first significant p-value is observed based on 4 looks
n_first_sig_looks<-apply(ptostmat_looks<alphalook[1], 1, function(x) match(TRUE, x))
result_all[1,21]<-(n_sim-sum(is.na(n_first_sig_looks)))/n_sim #power for sequential analyses
n_first_sig<-apply(ptostmat_looks<alphalook[1], 1, function(x) match(TRUE, x))
n_first_sig[is.na(n_first_sig)] <- looks
result_all[1,22]<-mean(n_first_sig*look_n[1]) #Mean sample size for TOST looking at each Look

#If true effect, does seq or BF have more power?
#If no true effect, does seq or BF have higher error rate?
result_all[12] #NHST
result_all[19] #BF

#If true effect, does TOST or BF have more power?
#If no true effect, does TOST or BF have higher error rate?
result_all[21] #TOST
result_all[20] #BF

#Which has lower ASN?
result_all[13] #NHST
result_all[14] #BF
result_all[22] #TOST

colnames(result_all)<-c("cohens_d", "looks", "number_of_simulations", "n", "alpha_each_look", "bf_cutoff_high", "bf_cutoff_low", "rscale_bf", "high_eqbound_d", "low_eqbound_d", "alpha_each_n", "power_sequential", "asn_sequential", "asn_bf", "asn_coast", "asn_bf_low", "asn_bf_low_look", "inconclusive_bf", "bf_stop_high_before_low_looks", "bf_stop_low_before_high_looks", "power_tost", "asn_tost")

write.table(result_all, "showdown_results.csv", sep = ",", append = T, row.names = FALSE, col.names = FALSE)

}


#Note: d_list needs to consist of effect sizes simulated using the sim_x_mat_y_mat.R script
d_list<-seq(0,0.8,0.1)
#N_list max is 1000 (or whatever max looks were generated using the sim_x_mat_y_mat.R)
n_list<-c(50,100,150,200,250,300,400,500,750,1000)
looks_list<-seq(2,8,1)

for (loop_d in 1:length(d_list)){
  for (loop_n in 1:length(n_list)){
    analyze_showdown(true_d=d_list[loop_d],n=n_list[loop_n])
  }
}

for (loop_d in 1:length(d_list)){
  for (loop_n in 1:length(n_list)){
    for (loop_looks in 1:length(looks_list)){
      analyze_showdown(true_d=d_list[loop_d],n=n_list[loop_n],looks=looks_list[loop_n])
    }
  }
}
