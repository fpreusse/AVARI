####
# Simulate data
####
# In this script, the anytime-valid simultaneous lower confidence bounds for the TDP
# are computed based on simulated data as described in Section 5 of the manuscript entitled:
# "Anytime-valid simultaneous lower confidence bounds for the true discovery proportion"

# The considered scenarios are run in parallel, using 12 cores.
# The number of cores can be adapted at the makeCluster()- function
# Note that running this script will take several hours.

# If only one scenario is to be run, the function "sim_tdp_bounds_fun" (see line 49) 
# can be run outside the "foreach" loop with the corresponding parameters of the
# considered scenario. 
# This function outputs a list; to get the results in a dataframe, the
# code from line 63 onward needs to be adapted manually.

# Set the working directory to the folder "simulation_study"
require(MASS)
require(safestats)
require(hommel)
require(foreach)
require(doParallel)

source("Basis_functions_sim.R")

# Simulation study  ------------------------------
## 1) Set the parameters as in Section 5 of the manuscript
rhos <- c(0.2, 0.6)
effects <- c(0.5,1,1.5)
size_Rs <- c(500, 50)

## 2) Define the design object needed to compute the e-processes
design_obj <- designSafeT(nPlan=100, deltaMin = 0.5, alternative = "greater", pb=F)

## 3) Initialize new environment used in foreach loop
env.tdp <- new.env()
env.tdp$design_obj <- design_obj
env.tdp$scenarios <- expand.grid(rhos, effects, size_Rs)

## 4) Run the function sim_e_vals_fun in parallel
cl <- makeCluster(12)
registerDoParallel(cl)
clusterExport(cl,
              c("scenarios", "design_obj"), envir = env.tdp)
clusterExport(cl, c("seq_TDP_e",  
                    "seq_TDP_hom", 
                    "m0_bound_average"))
TDP_bounds_full <- foreach(i= 1:12, .packages=c("MASS", "hommel", "safestats"))%dopar%{
  sim_tdp_bounds_fun(rho=scenarios[i,1],
                     size_R=scenarios[i,3], 
                     effect=scenarios[i,2], 
                     design_obj=design_obj, 
                     seed=14*i+i, 
                     nplan=100, 
                     B=1000, 
                     alpha=0.2)
}
stopCluster(cl)

# The foreach loop returns a list with 12 elements. 
# Save the entries in a dataframe instead to facilitate evaluation of the 
# simulated bounds.
df_full <- expand.grid(c(0.9,0.5,0.1),   # pi_1
                       c("mom", "ARI"),  # methods
                       c(0.5,1,1.5),     # effect size
                       c(0.2,0.6),       # rho
                       c(50, 500),       # size_R
                       c(11:100))        # stopping times
colnames(df_full) <- c("pi_1", "method", "effect", "rho", "size_R", "n")
df_full <- cbind(df_full, matrix(nrow=nrow(df_full), ncol=1000))


for(i in 1:12){
  if(grepl("0.5",names(TDP_bounds_full[[i]]),fixed=T)){effect <- 0.5}
  if(grepl("1",names(TDP_bounds_full[[i]]),fixed=T)){effect <- 1}
  if(grepl("1.5",names(TDP_bounds_full[[i]]),fixed=T)){effect <- 1.5}
  
  if(grepl("0.2",names(TDP_bounds_full[[i]]),fixed=T)){rho <- 0.2}else{rho <- 0.6}
  
  if(grepl("500",names(TDP_bounds_full[[i]]),fixed=T)){size_R <- 500}else{size_R <- 50}
  
  df_full[df_full[,1]==0.9&df_full[,2]=="mom" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom`$`0.9`[11:100,]
  df_full[df_full[,1]==0.9&df_full[,2]=="ARI" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`ARI`$`0.9`[11:100,]
  
  df_full[df_full[,1]==0.5&df_full[,2]=="mom" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom`$`0.5`[11:100,]
  df_full[df_full[,1]==0.5&df_full[,2]=="ARI" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`ARI`$`0.5`[11:100,]
  
  df_full[df_full[,1]==0.1&df_full[,2]=="mom" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom`$`0.1`[11:100,]
  df_full[df_full[,1]==0.1&df_full[,2]=="ARI" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`ARI`$`0.1`[11:100,]
}

file_out <- file.path("intermediate_results", "Result_simulations_df.RData")
save(df_full,file=file_out)

####
# Evaluate the simulated data
####
# In this script, the anytime-valid simultaneous lower confidence bounds for the TDP
# computed based on simulated data are evaluated as described in Section 5 of the manuscript entitled:
# "Anytime-valid simultaneous lower confidence bounds for the true discovery proportion"

# This includes the code to generate Figures 1 and 2 in the manuscript.

# Set the working directory to the folder "simulation_study".
require(ggplot2)
source("Basis_functions_sim.R")

load(file.path("intermediate_results","Results_simulations_df.RData"))

# Summary statistics and corresponding dataframe------------------------------------------------------
## Compute summary statistics to evaluate the simulation results
# These include the mean, standard deviation, range and non-coverage rate
df_sum <- expand.grid(c(0.9,0.5,0.1),   # pi_1
                      c("mom", "ARI"),  # methods
                      c(0.5,1,1.5),     # effect size
                      c(0.2,0.6),       # rho
                      c(50, 500),       # size_R
                      c(11:100))        # observation time points
colnames(df_sum) <- c("pi_1", "method", "effect", "rho", "size_R", "n")
df_sum$mean <- rowMeans(df_full[,-c(1:6)], na.rm=T)
df_sum$sd <- apply(df_full[,-c(1:6)],1,sd, na.rm=T)
df_sum$max <- apply(df_full[,-c(1:6)],1,max, na.rm=T)
df_sum$min <- apply(df_full[,-c(1:6)],1,min, na.rm=T)
df_sum$non_cov <- NA
df_sum$non_cov[df_sum$pi_1==0.9] <- apply(df_full[df_full$pi_1==0.9,-c(1:6)],1,non_cov_rate, true_pi=0.9)
df_sum$non_cov[df_sum$pi_1==0.5] <- apply(df_full[df_full$pi_1==0.5,-c(1:6)],1,non_cov_rate, true_pi=0.5)
df_sum$non_cov[df_sum$pi_1==0.1] <- apply(df_full[df_full$pi_1==0.1,-c(1:6)],1,non_cov_rate, true_pi=0.1)

# Add the true values for comparsion
df_helper <- expand.grid(c(0.9,0.5,0.1), 
                         c("pi_1"),      
                         c(0.5,1,1.5),
                         c(0.2,0.6),
                         c(50, 500),
                         c(11:100))
colnames(df_helper) <- c("pi_1", "method", "effect", "rho", "size_R", "n")
df_helper$mean <- df_helper$max <- df_helper$sd <-df_helper$min <- df_helper$pi_1
df_helper$non_cov <- 0.2

df_sum <- rbind(df_sum, df_helper)

# Based on the dataframe containing the summary statistics,
# a dataframe which is used to generate the figures using ggplot is defined.

df_work <- df_sum[df_sum$pi_1%in%c(0.9,0.1),]
df_work$Methods <- factor(df_work$method,
                          levels =c("pi_1","ARI", "mom"), ordered=T) 
df_work$Pi_1 <- factor(df_work$pi_1, levels=c(0.1,0.9),
                       ordered=T, labels=c(expression(paste(pi[1],"(R)=0.1")), 
                                           expression(paste(pi[1],"(R)=0.9"))))
df_work$Rho <- factor(df_work$rho, levels=c(0.2, 0.6),
                      ordered=T, labels=c(expression(paste(rho,"=0.2")), 
                                          expression(paste(rho,"=0.6"))))
df_work$Effect <- factor(df_work$effect, levels=c(0.5,1, 1.5),
                         ordered=T, labels=c(expression(paste(mu,"=0.5")), 
                                             expression(paste(mu,"=1")),
                                             expression(paste(mu,"=1.5"))))
df_work$Size_R <- factor(df_work$size_R, levels=c(50, 500),
                         ordered=T, labels=c(expression(paste(abs(R),"=50")), 
                                             expression(paste(abs(R),"=500"))))


