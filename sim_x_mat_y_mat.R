#Create x_mat and y_mat

options(scipen=999) #disable scientific notation for numbers

set.seed(1)

#Simulation settings
n_sim<-20000 #numbber of simulated studies
true_d<-0 #True effect size when sd = 1 (for between) or 1/sqrt2) (for within)
true_sd<-1 #Set True standard deviation (set to 1 to make D cohen's d.
n<-1000 #total number of datapoints (per condition) you are willing to collect

#define some variables to store simulations after each n
x<-numeric(n) #store all data group x
y<-numeric(n) #store all data group y
xmat<-matrix(NA, nrow=n_sim, ncol=n) #matrix for all x values 
ymat<-matrix(NA, nrow=n_sim, ncol=n) #matrix for all y values
#### RUN SIMULATION ----

# create progress bar
pb <- winProgressBar(title = "progress bar", min = 0, max = n_sim, width = 300)
#Run simulation - create 1 participant at a time, store all data in a list, so I can do tests after any number of participants I want
for (j in 1:n_sim){
  setWinProgressBar(pb, j, title=paste( round(j/n_sim*100, 2), "% done"))
  for(i in 1:n){ #for each simulated experiment
    #x[i]<-rnorm(n = 1, mean = 0, sd = true_sd)
    y[i]<-rnorm(n = 1, mean = true_d, sd = true_sd) #Use D as a difference score, because sd=1 equals Cohen's D
  }  
  #xmat[j,]<-x
  ymat[j,]<-y
}
#close progress bar
close(pb)

#Store data
#saveRDS(xmat, paste("xmat_d_",true_d,"_true_sd_",true_sd,"_n_",n,"_n_sim_",n_sim,".Rdata", sep=""))
saveRDS(ymat, paste("ymat_d_",true_d,"_true_sd_",true_sd,"_n_",n,"_n_sim_",n_sim,".Rdata", sep=""))
