#TO DO: TOST lower bound - 99999 for one-sidedness
# Create option to switch between one-sided and two-sided tests

#Script to simulate studies and apply the Frequentist and Bayesian sequential sampling recommendation.

#set.seed(2) #Set seed 2 to recreate simulation in manuscript
if(!require(pscl)){install.packages('pscl')}
library(pscl) 
if(!require(MBESS)){install.packages('MBESS')}
library(MBESS) 
if(!require(pwr)){install.packages('pwr')}
library(pwr) 
if(!require(BayesFactor)){install.packages('BayesFactor')}
library(BayesFactor) 
if(!require(gsDesign)){install.packages('gsDesign')}
library(gsDesign)

options(scipen=999) #disable scientific notation for numbers

#Simulation settings (n_sim, true effect size D)
n_sim<-10000 #numbber of simulated studies
true_d<-0.6 #True effect size (Keep SD below to 1, otherwise, this is just mean dif, not d)
true_sd<-1/sqrt(2) #Set True standard deviation (set to 1 to make D cohen's d.
lower_bound<-20 #After how any participants do you start to look? Frick calls this the lower bound 

#Sequential design settings NHST (alpha and # looks)
n<-1000 #total number of datapoints (per condition) you are willing to collect
sided <- 2 #Set whether tests are 1-sided or 2-sided

#define some variables to store simulations after each n
x<-numeric(n) #store all data group x
y<-numeric(n) #store all data group y
xmat<-matrix(NA, nrow=n_sim, ncol=n) #matrix for all x values 
ymat<-matrix(NA, nrow=n_sim, ncol=n) #matrix for all y values
pmat<-matrix(NA, nrow=n_sim, ncol=n) #matrix for all p-values (for plot)
dmat<-matrix(NA, nrow=n_sim, ncol=n) #matrix for all effect sizes d (for plot) 
cormat<-matrix(NA, nrow=n_sim, ncol=n) #matrix for all cors 
#### RUN SIMULATION ----

# create progress bar
pb <- winProgressBar(title = "progress bar", min = 0, max = n_sim, width = 300)
#Run simulation - create 1 participant at a time, store all data in a list, so I can do tests after any number of participants I want
for (j in 1:n_sim){
  setWinProgressBar(pb, j, title=paste( round(j/n_sim*100, 2), "% done"))
  #Calculate initial 10 datapoints before performing first t-test, cannot do test on n = 1 and small samples vary widely.
  for(i in 1:10){ #for each simulated experiment
    x[i]<-rnorm(n = 1, mean = 0, sd = true_sd)
    y[i]<-rnorm(n = 1, mean = true_d, sd = true_sd) #Use D as a difference score, because sd=1 equals Cohen's D
  }  
  
  #After the first 10, I perform one-sided statistical tests after every participant. 
  for(i in 11:n){ #for each simulated participants after the first 10
    x[i]<-rnorm(n = 1, mean = 0, sd = true_sd)
    y[i]<-rnorm(n = 1, mean = true_d, sd = true_sd)
    mean_x<-mean(x[1:i])
    mean_y<-mean(y[1:i])
    sd_x<-sd(x[1:i])
    sd_y<-sd(y[1:i])
    cormat[j,i]<-cor(x,y)
    if (sided == 2){
      z<-t.test(x[1:i],y[1:i], paired = TRUE, alternative = "two.sided", var.equal=TRUE, alpha = alpha_each_n) #perform the t-test
      pmat[j,i]<-z$p.value #store all p-values for all simulations
    }
    if (sided == 1){
      z<-t.test(x[1:i],y[1:i], paired = TRUE, alternative = "less", var.equal=TRUE, alpha = alpha_each_n) #perform the t-test
      pmat[j,i]<-z$p.value #store all p-values for all simulations
    }
    #d<-(mean_y-mean_x)/(sqrt((((length(x[1:i]) - 1)*((sd_x^2))) + (length(y[1:i]) - 1)*((sd_y^2)))/((length(x[1:i])+length(y[1:i])-2)))) #Cohen;s d for between designs
    d<-z$statistic/sqrt(length(x[1:i])) #Cohen's d for within design
    dmat[j,i]<-d #store all effect sizes for all simulations
  }
  xmat[j,]<-x
  ymat[j,]<-y
}
#close progress bar
close(pb)


#### INTERPRET SIMULATION RESULTS ----

#COAST
#COAST requires a minimum number of observations - Frick recommends not testing first 20. So these are dropped.  
pmat_coast<-pmat
pmat_coast[,1:lower_bound-1]<-NA

#Find at which N first p<0.01 and first p > 0.36 is observed 
#Results are not exactly the same. Not sure why, maybe because Frick writes:
#"The value of t needed for p = .36 was computed from a linear average" - so he did not use p-values directly?
n_first_sig<-apply(pmat_coast<0.01, 1, function(x) match(TRUE, x))
n_first_nonsig<-apply(pmat_coast>0.36, 1, function(x) match(TRUE, x))
#if NA no p-value below 0.01 or above 0.37) code as n
n_first_sig[is.na(n_first_sig)] <- n
n_first_nonsig[is.na(n_first_nonsig)] <- n

coast_n_stop<-pmin(n_first_sig,n_first_nonsig)
mean(coast_n_stop)

coast_stop_high_before_low<-sum(n_first_sig>n_first_nonsig)/n_sim
coast_stop_low_before_high<-sum(n_first_sig<n_first_nonsig)/n_sim
coast_inconclusive<-1-coast_stop_high_before_low-coast_stop_low_before_high
coast_stop_high_before_low #How often is p < 0.01 passed first
coast_stop_low_before_high #How often is p > 0.36 passed first
coast_inconclusive #How often does the p-value 0.01 < p < 0.036?


#### Check if simulation returned the expected results----

#check if SD's are 1
mean(apply(xmat, 1, sd))
mean(apply(ymat, 1, sd))

#check if means are 0 and 0.4
mean(apply(xmat, 1, mean))
mean(apply(ymat, 1, mean))
#check differences are 0-0.4
difmat<-xmat-ymat
mean(apply(difmat, 1, mean))
#check sd of dif scores (when sd = 1, sd_dif = sqrt(2) in within design)
mean(apply(difmat, 1, sd))
#Check if correlation is 0
cormat<-drop(cormat[,11:n])
mean(cormat)     
#check if effect size is correct
dmat<-drop(dmat[,11:n])
mean(dmat)
