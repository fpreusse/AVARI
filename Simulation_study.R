# Simulation study (is run in parallel) ------------------------------
## Set- up
source("Basis_functions.R")
library(MASS)
library(safestats)
library(hommel)
library(foreach)
library(doParallel)

## 1) define the design object needed to compute the e-processes, is the same for all scenarios
design_obj <- designSafeT(nPlan=100, deltaMin = 0.5, alternative = "greater", pb=F)

## 2) initialize the different scenarios:
rhos <- c(0.2, 0.6)
effects <- c(0.5,1,1.5)
size_Rs <- c(500, 50)

## 3) initialize new environment used in foreach loop
env.tdp <- new.env()
env.tdp$design_obj <- design_obj
env.tdp$scenarios <- expand.grid(rhos, effects, size_Rs)

## 4) run the function sim_e_vals_fun in parallel
cl <- makeCluster(12)
registerDoParallel(cl)
clusterExport(cl,
              c("scenarios", "design_obj"), envir = env.tdp)
clusterExport(cl, c("seq_TDP_p_s", "seq_TDP_mom_s", "seq_TDP_ARC_mom_s", 
                    "seq_TDP_hom", "m0_bound_average"))
TDP_bounds_full <- foreach(i= 1:12, .packages=c("MASS", "hommel", "safestats"))%dopar%{
  sim_e_vals_fun(rho=scenarios[i,1],size_R=scenarios[i,3], effect=scenarios[i,2], design_obj=design_obj, seed=12*i+i, nplan=100, B=1000, alpha=0.2)
}
stopCluster(cl)

# Evaluation of simulation study ------------------------------------------
## 1) restructure output TDP_bounds_full to be in a data.frame
df_full <- expand.grid(c(0.9,0.5,0.1), #pi_1
                       c("e-to-p", "mom", "mom_ARD","ARI"), #methods
                       c(0.5,1,1.5), #effect size
                       c(0.2,0.6), #rho
                       c(50, 500),#size_R
                       c(11:100))#time points at which we have computed the TDP bounds
colnames(df_full) <- c("pi_1", "method", "effect", "rho", "size_R", "n")
df_full <- cbind(df_full, matrix(nrow=nrow(df_full), ncol=1000))

for(i in 1:12){
  # what effect size, rho and size of discovery set did we simulate?
  if(grepl("0.5",names(TDP_bounds_full[[i]]),fixed=T)){effect <- 0.5}
  if(grepl("1",names(TDP_bounds_full[[i]]),fixed=T)){effect <- 1}
  if(grepl("1.5",names(TDP_bounds_full[[i]]),fixed=T)){effect <- 1.5}
  
  if(grepl("0.2",names(TDP_bounds_full[[i]]),fixed=T)){rho <- 0.2}else{rho <- 0.6}
  
  if(grepl("500",names(TDP_bounds_full[[i]]),fixed=T)){size_R <- 500}else{size_R <- 50}
  
  df_full[df_full[,1]==0.9&df_full[,2]=="e-to-p" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`e-to-p`$`0.9`[11:100,]
  df_full[df_full[,1]==0.9&df_full[,2]=="mom" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom`$`0.9`[11:100,]
  df_full[df_full[,1]==0.9&df_full[,2]=="mom_ARD" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom_ARD`$`0.9`[11:100,]
  df_full[df_full[,1]==0.9&df_full[,2]=="ARI" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`ARI`$`0.9`[11:100,]
  
  df_full[df_full[,1]==0.5&df_full[,2]=="e-to-p" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`e-to-p`$`0.5`[11:100,]
  df_full[df_full[,1]==0.5&df_full[,2]=="mom" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom`$`0.5`[11:100,]
  df_full[df_full[,1]==0.5&df_full[,2]=="mom_ARD" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom_ARD`$`0.5`[11:100,]
  df_full[df_full[,1]==0.5&df_full[,2]=="ARI" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`ARI`$`0.5`[11:100,]
  
  df_full[df_full[,1]==0.1&df_full[,2]=="e-to-p" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`e-to-p`$`0.1`[11:100,]
  df_full[df_full[,1]==0.1&df_full[,2]=="mom" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom`$`0.1`[11:100,]
  df_full[df_full[,1]==0.1&df_full[,2]=="mom_ARD" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`mom_ARD`$`0.1`[11:100,]
  df_full[df_full[,1]==0.1&df_full[,2]=="ARI" & df_full[,3]==effect &df_full[,4]==rho& df_full[,5]==size_R,-c(1:6)] <- TDP_bounds_full[[i]][[1]]$`ARI`$`0.1`[11:100,]
}
## 2) Compute summary statistics based on observed bounds
df_sum <- expand.grid(c(0.9,0.5,0.1), #pi_1
                      c("e-to-p", "mom", "mom_ARD", "ARI"),#methods
                      c(0.5,1,1.5), #effect size
                      c(0.2,0.6), #rho
                      c(50, 500),#size_R
                      c(11:100)) #time points of observations
colnames(df_sum) <- c("pi_1", "method", "effect", "rho", "size_R", "n")
df_sum$mean <- rowMeans(df_full[,-c(1:6)], na.rm=T)
df_sum$sd <- apply(df_full[,-c(1:6)],1,sd, na.rm=T)
df_sum$max <- apply(df_full[,-c(1:6)],1,max, na.rm=T)
df_sum$min <- apply(df_full[,-c(1:6)],1,min, na.rm=T)
df_sum$non_cov <- NA
df_sum$non_cov[df_sum$pi_1==0.9] <- apply(df_full[df_full$pi_1==0.9,-c(1:6)],1,non_cov_rate, true_pi=0.9)
df_sum$non_cov[df_sum$pi_1==0.5] <- apply(df_full[df_full$pi_1==0.5,-c(1:6)],1,non_cov_rate, true_pi=0.5)
df_sum$non_cov[df_sum$pi_1==0.1] <- apply(df_full[df_full$pi_1==0.1,-c(1:6)],1,non_cov_rate, true_pi=0.1)


# Graphics ----------------------------------------------------------------
library(ggplot2)
## 1) Initialize data.frame to be used for graphics
df_helper <- expand.grid(c(0.9,0.5,0.1), #pi_1
                         c("pi_1"), #methods
                         c(0.5,1,1.5), #effect size
                         c(0.2,0.6), #rho
                         c(50, 500),#size_R
                         c(11:100)) #time points of observations
colnames(df_helper) <- c("pi_1", "method", "effect", "rho", "size_R", "n")
df_helper$mean <- df_helper$max <- df_helper$sd <-df_helper$min <- df_helper$pi_1
df_helper$non_cov <- 0.2

df_sum <- rbind(df_sum, df_helper)

# since we are only interested in comparing a true TDP of 0.9 and 0.1:
df_work <- df_sum[df_sum$pi_1%in%c(0.9,0.1),]
df_work$Methods <- factor(df_work$method,
                          levels =c("pi_1","ARI", "e-to-p", "mom", "mom_ARD"), ordered=T) 

df_work$Pi_1 <- factor(df_work$pi_1, levels=c(0.1,0.9),
                       ordered=T, labels=c(expression(paste(pi[1],"(R)=0.1")), expression(paste(pi[1],"(R)=0.9"))))
df_work$Rho <- factor(df_work$rho, levels=c(0.2, 0.6),
                      ordered=T, labels=c(expression(paste(rho,"=0.2")), expression(paste(rho,"=0.6"))))
df_work$Effect <- factor(df_work$effect, levels=c(0.5,1, 1.5),
                         ordered=T, labels=c(expression(paste(mu,"=0.5")), expression(paste(mu,"=1")),expression(paste(mu,"=1.5"))))
df_work$Size_R <- factor(df_work$size_R, levels=c(50, 500),
                         ordered=T, labels=c(expression(paste(abs(R),"=50")), expression(paste(abs(R),"=500"))))

## 2) Graphics displaying the non-coverage rate
# Compare pi_1 and mean for rho=0.2, size_R=500
ggplot(df_work[df_work$size_R==500 & df_work$rho==0.2,], aes(x=n, y=non_cov, group = Methods))+
  theme_bw()+
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5))+
  geom_line(aes(linetype=Methods,color=Methods), linewidth=1.25)+
  scale_color_manual(breaks= c("e-to-p", "mom", "mom_ARD", "ARI"), #,"ARI")
                     values=c("pi_1"="black", "e-to-p"="cyan2", "mom"="cyan4", "mom_ARD"="aquamarine3", "ARI"="grey70"))+ #legend in color, do not show m_1
  scale_linetype_manual(breaks= c("e-to-p", "mom", "mom_ARD", "ARI"),
                        values=c("pi_1"=1, "e-to-p"=2, "mom"=5, "mom_ARD"=4, "ARI"=1))+
  ylab("Empirical non-coverage rate")+
  xlab("Sample size n")+
  scale_x_continuous(breaks=c(20,40,60,80,100), limits = c(11, 100))+
  facet_grid(Pi_1~Effect, scales="free_y", labeller=label_parsed)

# Compare rho and size_R for pi_1=0.1, mean=0.5
ggplot(df_work[df_work$effect==0.5 & df_work$pi_1==0.1,], aes(x=n, y=non_cov, group = Methods))+
  theme_bw()+
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5))+
  geom_line(aes(linetype=Methods,color=Methods), linewidth=1.25)+
  scale_color_manual(breaks= c("e-to-p", "mom", "mom_ARD", "ARI"), #,"ARI")
                     values=c("pi_1"="black", "e-to-p"="cyan2", "mom"="cyan4", "mom_ARD"="aquamarine3", "ARI"="grey70"))+ #legend in color, do not show m_1
  scale_linetype_manual(breaks= c("e-to-p", "mom", "mom_ARD", "ARI"),
                        values=c("pi_1"=1, "e-to-p"=2, "mom"=5, "mom_ARD"=4, "ARI"=1))+
  ylab("Empirical non-coverage rate")+
  xlab("Sample size n")+
  scale_x_continuous(breaks=c(20,40,60,80,100), limits = c(11, 100))+
  facet_grid(Size_R~Rho, scales="free_y", labeller=label_parsed)

## 3) Graphics displaying the power
# Compare pi_1 and mean for rho=0.2, size_R=500
ggplot(df_work[df_work$size_R==500 & df_work$rho==0.2,], aes(x=n, y=mean, group = Methods))+
  theme_bw()+
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5))+
  geom_line(aes(linetype=Methods,color=Methods), linewidth=1.25)+
  scale_color_manual(breaks= c("e-to-p", "mom", "mom_ARD", "ARI"),
                     values=c("pi_1"="black", "e-to-p"="cyan2", "mom"="cyan4", "mom_ARD"="aquamarine3", "ARI"="grey70"))+ #legend in color, do not show m_1
  scale_linetype_manual(breaks= c("e-to-p", "mom", "mom_ARD", "ARI"),
                        values=c("pi_1"=1, "e-to-p"=2, "mom"=5, "mom_ARD"=4, "ARI"=1))+
  ylab(expression(bar(paste(hat(pi)[1]^n,"(R)", sep=""))))+
  xlab("Sample size n")+
  scale_x_continuous(breaks=c(20,40,60,80,100), limits = c(11, 100))+
  facet_grid(Pi_1~Effect, scales="free_y", labeller=label_parsed)

# Compare rho and size_R for pi_1=0.1, mean=0.5
ggplot(df_work[df_work$effect==0.5 & df_work$pi_1==0.1,], aes(x=n, y=mean, group = Methods))+
  theme_bw()+
  theme(text=element_text(size=20), plot.title = element_text(hjust = 0.5))+
  geom_line(aes(linetype=Methods,color=Methods), linewidth=1.25)+
  scale_color_manual(breaks= c("e-to-p", "mom", "mom_ARD", "ARI"), #,"ARI")
                     values=c("pi_1"="black", "e-to-p"="cyan2", "mom"="cyan4", "mom_ARD"="aquamarine3", "ARI"="grey70"))+ #legend in color, do not show m_1
  scale_linetype_manual(breaks= c("e-to-p", "mom", "mom_ARD", "ARI"),
                        values=c("pi_1"=1, "e-to-p"=2, "mom"=5, "mom_ARD"=4, "ARI"=1))+
  ylab(expression(bar(paste(hat(pi)[1]^n,"(R)", sep=""))))+
  xlab("Sample size n")+
  scale_x_continuous(breaks=c(20,40,60,80,100), limits = c(11, 100))+
  scale_y_continuous(breaks=c(0,0.025,0.05,0.075,0.1))+
  facet_grid(Size_R~Rho, scales="free_y", labeller=label_parsed)