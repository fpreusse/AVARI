require(safestats)
require(hommel)

###
# Compute anytime-valid simultaneous upper bound for the number of false discoveries based on e-processes at time n
###

## Input:
# vec_e: vector of values of the e-processes at time n corresponding to the elementary hypotheses
# setS: indices of elementary hypotheses that are part of the discovery set
# alpha: confidence level
# ell_rej: for non-increasing bounds (i.e., ARD guarantee): upper bound from previous sequences? 

## Output:
# upper (1-\alpha) confidence bound for the number of false discoveries at time n

## Notes:
# This function implements the shortcut introduced in Lemma 1.
# Therefore, the average is used as e-merging function for e-processes.
# If ell_rej is given, ARD is ensured through bookkeeping (i.e., not testing hypotheses that have already been rejected)

m0_bound_average <- function(vec_e, setS, alpha, ell_rej=NULL){
  # Sort e-process values corresponding to discoveries
  S_e <- sort(vec_e[setS])
  # Define maximum value of h (see Lemma 1)
  if(!is.null(ell_rej)){
    h <- ell_rej
  }else{
    h <- length(S_e)
  }
  c_alpha <- 0
  
  ## If all elementary hypotheses are discoveries:
  # check at each level if the intersection hypothesis with the smallest corresponding e-process value is rejected
  if(length(setS)==length(vec_e)){
    # Iterative approach: decrease h 
    while(c_alpha==0 & h >=1){
      e_h <- sum(S_e[1:h])
      if(e_h<(h/alpha)){
        c_alpha <- h
      }else{
        h <- h-1
      }
    }
  }else{
    ## Compute discovery set specific critical value as in Eq. 6 (right side of inequation) 
    
    # Since the average is used as e-merging function: 
    # only merge e-processes corresponding to hypotheses in H_K_h with e-process values that are SMALLER than the e-process value corresponding to H_K_h
    # otherwise, the e-process value of the intersection is larger than the e-process value of H_K_h.
    BarS_e <- sort(vec_e[-setS])
    if(sum(BarS_e>mean(S_e))!=0& sum(BarS_e>mean(S_e))!=length(BarS_e)){
      BarS_e <- BarS_e[-which(BarS_e>mean(S_e))]
    }
    max_Cumsum_e <- max((1:length(BarS_e))/alpha-cumsum(BarS_e))
    
    ## Apply shortcut: iteratively decrease size of the intersection hypothesis H_K_h until H_K_h is not rejected by closed testing.
    while(c_alpha== 0 & h>=1){
      e_h <- sum(S_e[1:h])
      if(e_h<(h/alpha)){
        # if H_K_H is rejected by local level-alpha test it is also rejected by closed testing.
        c_alpha <- h
      }else{
        if((e_h-h/alpha)<max_Cumsum_e){ 
          # If true, then hypothesis is not rejected by closed testing
          c_alpha <- h
        }else{
          # If H_K_h is rejected by closed testing: decrease h
          h <- h-1
        }
      }
    }
  }
  return(c_alpha)
}


# Anytime-valid simultaneous lower confidence bounds for the TDP ------------------------------------
# These function compute the anytime-valid simultaneous confidence bounds at multiple time points for multiple discovery sets
# Different anytime-valid local level-alpha tests are used to define the rejection set of closed testing

###
# Anytime-valid simultaneous lower confidence bounds  based on the e-to-p process
###

## Input:
# mat_e_seq: matrix containing the e-processes for each of the m elementary hypotheses. 
# alpha: confidence level
# setSs: list of discovery sets to be considered at each time point (fixed over time).
#         Each list element contains a vector of indices indicating the hypotheses in the discovery set
# n_min: minimal sample size based on which the bounds are computed, has to be the same as the 
# N: maximal sample size used for computation of the TDP bounds.

## Output:
# seq_TDP: matrix containing the TDP bounds at each considered time point (rows) for each discovery set (columns)
# seq_p: matrix containing the p-processes, constructed as mat_e_seq

## Notes:
# The matrix mat_e_seq has N-n_min+1 rows and m columns
# Computation of the bounds utilizes the computational shortcut by Goeman et al (2019), implemented in the "hommel" package
# Bounds are computed at time n=n_min,n_min+1,...,N

seq_TDP_p_s <- function(mat_e_seq,alpha, setSs, n_min,N){
  m <- ncol(mat_e_seq)
  n_sets <- length(setSs)
  length_sets <- sapply(setSs, length)
  
  # initialize matrix for sequential p-values
  mat_p_seq <- matrix(nrow=(N-n_min+1), ncol=m)
  # initialize matrix for sequential TDP bounds
  mat_TDP <- mat_ell_rej <- matrix(nrow=(N-n_min+1), ncol=n_sets)
  
  ## compute bounds for every observation time point, starting with n=n_min
  for(n in n_min:N){
    # at each time point:
    # 1) compute the e-to-p-process
    for(j in 1:m){
      mat_p_seq[(n-n_min+1),j] <- min(1/mat_e_seq[c(1:((n-n_min+1))),j])
    }
    # 2) compute the TDP bounds for each considered discovery set
    hom <- hommel(mat_p_seq[(n-n_min+1),], simes=T)
    for(s in 1:n_sets){
      mat_TDP[(n-n_min+1),s] <- tdp(hom, ix=setSs[[s]], alpha=alpha)
    }
  }
  return(list(seq_TDP=mat_TDP,
              seq_p = mat_p_seq))
}

###
# Anytime-valid simultaneous lower confidence bounds based on an e-process
###

## Input:
# mat_e_seq: matrix containing the e-processes corresponding to each of the m elementary hypotheses. 
# alpha: confidence level
# setSs: list of discovery sets to be considered at each time point (fixed over time).
#         Each list element contains a vector of indices indicating the hypotheses in the discovery set
# n_min: minimal sample size based on which the bounds are computed.
# N: maximal sample size used for computation of the TDP bounds.

## Output:
# seq_TDP: matrix containing the TDP bounds at each considered time point (rows) for each discovery set (columns)
# seq_e: mat_e_seq

## Notes:
# The matrix mat_e_seq has N-n_min+1 rows and m columns
# The bounds may be decreasing
# The average is used as e-merging function
# Bounds are computed at time n=n_min,n_min+1,...,N

seq_TDP_mom_s <- function(mat_e_seq,alpha, setSs,  n_min, N){
  m <- ncol(mat_e_seq)
  n_sets <- length(setSs)
  length_sets <- sapply(setSs, length)
  # initialize matrix for sequential TDP bounds
  mat_TDP <- matrix(nrow=(N-n_min+1), ncol=n_sets)

  ## compute bounds for every observation time point, starting with n=n_min
  for(n in n_min:N){
    # at each time point:
    # compute the TDP bounds for each considered discovery set 
    for(s in 1:n_sets){
      mat_TDP[(n-n_min+1),s] <- 1-(m0_bound_average(vec=mat_e_seq[(n-n_min+1),], setS=setSs[[s]], alpha=alpha))/length_sets[s]
    }
  }
  return(list(seq_TDP=mat_TDP,
              seq_e = mat_e_seq))
}

###
# Anytime-valid simultaneous lower confidence bounds based on an e-process with ARD guarantee
###

## Input:
# mat_e_seq: matrix containing the e-processes corresponding to each of the m elementary hypotheses. 
# alpha: confidence level
# setSs: list of discovery sets to be considered at each time point (fixed over time).
#         Each list element contains a vector of indices indicating the hypotheses in the discovery set
# n_min: minimal sample size based on which the bounds are computed.
# N: maximal sample size used for computation of the TDP bounds.

## Output:
# seq_TDP: matrix containing the TDP bounds at each considered time point (rows) for each discovery set (columns)
# seq_e: mat_e_seq

## Notes:
# The matrix mat_e_seq has N-n_min+1 rows and m columns
# The bounds are non-decreasing, this is ensured by bookkeeping (i.e., hypotheses rejected at an earlier time point are not tested again)
# The average is used as e-merging function
# Bounds are computed at time n=n_min,n_min+1,...,N

seq_TDP_ARC_mom_s <- function(mat_e_seq,alpha, setSs,  n_min,N ){
  m <- ncol(mat_e_seq)
  n_sets <- length(setSs)
  length_sets <- sapply(setSs, length)
  
  # initialize matrix for sequential TDP bounds and for upper bounds for the number of false rejections (ell_rej)
  mat_TDP <- mat_ell_rej <- matrix(nrow=(N-n_min+1), ncol=n_sets)
  
  ## compute bounds for every observation time point, starting with n=n_min
  for(n in n_min:N){
    # at each time point:
    # compute the (non-decreasing) TDP bounds for each considered discovery set 
     if(n==n_min){
       # at the first time point, no previous rejections are known, therefore compute bounds the same as without ARD guarantee
      for(s in 1:n_sets){
        mat_ell_rej[1,s] <- m0_bound_average(vec_e=mat_e_seq[1,], setS=setSs[[s]], alpha=alpha)
        mat_TDP[1,s] <- 1-mat_ell_rej[1,s]/length_sets[s]
      }
    }else{
      # after the first time point, utilize knowledge about bounds at previous observation time points to ensure ARD
      for(s in 1:n_sets){
        mat_ell_rej[(n-n_min+1),s] <- m0_bound_average(vec_e=mat_e_seq[(n-n_min+1),], setS=setSs[[s]], alpha=alpha, ell_rej =mat_ell_rej[(n-n_min),s])
        mat_TDP[(n-n_min+1),s] <- 1-mat_ell_rej[(n-n_min+1),s]/length_sets[s]
      }
    }
  }
  return(list(seq_TDP=mat_TDP,
              seq_e = mat_e_seq))
}

###
# Simultaneous lower confidence bounds based on work by Goeman and Solari (2011)
###

## Input
# mat_x_seq: matrix containing the observed values, with ncol=m, nrow=n_obs, which are the imput for a one-sample t-test
# alpha: confidence level
# setSs: list of discovery sets to be considered at each time point (fixed over time).
#         Each list element contains a vector of indices indicating the hypotheses in the discovery set
# n_min: minimal sample size based on which the bounds are computed
# N: maximal sample size used for computation of the TDP bounds.

## Output
# matrix containing the TDP bounds at each considered time point (rows) for each discovery set (columns)

## Notes
# These bounds are not anytime-valid
# Computation of the bounds utilizes the computational shortcut by Goeman et al (2019), implemented in the "hommel" package

seq_TDP_hom <- function(mat_x_seq, alpha, setSs, n_min=11){ 
  n_obs <- nrow(mat_x_seq)
  m <- ncol(mat_x_seq) 
  n_sets <- length(setSs)
  
  # Initialize matrix for p-values at each time point
  mat_p_seq <- matrix(nrow=(n_obs-(n_min-1)), ncol=m)
  # initialize matrix for sequential TDP bounds
  mat_TDP <- matrix(nrow=(n_obs-(n_min-1)), ncol=n_sets)
  
  ## compute bounds for every observation time point, starting with n=n_min
  for(n in n_min:n_obs){
    # at each time point:
    # 1) compute the p-values corresponding to the m elementary null hypotheses
    for(j in 1:m){
      mat_p_seq[(n-(n_min-1)),j] <- t.test(x= mat_x_seq[c(1:n),j],
                                           alternative="greater")$p.value
      # t-test based on the observations in mat_x_seq
    }
    # 2) compute the TDP bounds for each considered discovery set
    hom <- hommel(mat_p_seq[(n-(n_min-1)),], simes=T)
    for(s in 1:n_sets){
      mat_TDP[(n-(n_min-1)),s] <- tdp(hom, ix=setSs[[s]], alpha=alpha)
    }
  }
  return(mat_TDP)
}

# Functions used in the simulation study -----------------------------------------
# This first function can be used to simulate data and compute (anytime-valid) simultaneous lower confidence bounds of the TDP
# as described in Section 4.
# The second function is used to compute the empirical non-coverage rate

###
# Simulating one scenario
###

## Input
# rho: level of dependency between the elementary hypotheses
# size_R: size of the rejection sets
# effect: what is the mean of the data if the alternative hypothesis is true
# design_obj: Design object needed to compute e-processes. Does not depent on the observed data
# seed: set a seed to make study reproducible
# nplan: maximal sample size used for computation of the TDP bounds.
# B: number of Monte-Carlo iterations
# alpha: confidence level
# n_min: minimal sample size based on which the bounds are computed.

## Output
# List containing the (anytime-valid) simultaneous lower confidence bounds based on the e-to-p process, the mom process, mom_ARD process and ARI
# The list is structured as follows: [[Procedure]][[true TDP]][[TDP_bounds]]
# TDP_bounds are saved in a matrix with N rows and B columns

## Notes:
# This function automatically saves the output as well as the e-processes for each elementary hypothesis
# The file is named "mu=[effect]_rho=[rho]_size_R=[size_R].Rdata"
# Only rejection sets with a true TDP of 0.1, 0.5 and 0.9 are computed, 
# should this be adapted: change m0, list_sets and for-loops over the discovery sets (indexed by i) in the function

sim_e_vals_fun <- function(rho,size_R, effect, design_obj, seed=12, nplan=100, B=1000, alpha=0.2, n_min=10){
  set.seed(seed)
  # Initialize parameters for simulation of the data
  m <- 1000
  pi0 <- 0.5 # half of the m elementary hypotheses are true.
  V <- matrix(rho,ncol=m, nrow=m)
  diag(V) <- 1
  n_obs <- nplan
  
  ## Initialize discovery sets of size size_R with true TDP of 0.1, 0.5 and 0.9
  m0 <- c(size_R*0.1, size_R*0.5, size_R*0.9)
  list_sets <- list(pi0_1=c(sample(c(1:500),m0[1]),sample(c(501:m),m0[3])),
                    pi0_5=c(sample(c(1:500),m0[2]),sample(c(501:m),m0[2])),
                    pi0_9=c(sample(c(1:500),m0[3]),sample(c(501:m),m0[1])))
  
  # Function to initialize lists where output is saved
  make_list_matrix <- function(i,nplan, B){
    pi_1_bounds<- lapply(c(1:3), matrix, data=0, nrow=nplan, ncol=B)
    names(pi_1_bounds) <- c("0.9", "0.5", "0.1")
    return(pi_1_bounds)
  }
  
  ## Initialize lists where output is saved
  TDP_bounds <- lapply(c(1:4), make_list_matrix, nplan=nplan, B=B)
  names(TDP_bounds)<- c("e-to-p", "mom", "mom_ARD", "ARI")
  E_vals<- lapply(rep(B,3), vector, mode="list")
  names(E_vals) <- c("e-to-p", "mom", "mom_ARD")
  
  ## Simulate the data and compute the e-processes and TDP bounds
  # do this for B independent Monte Carlo repetitions
  for(b in 1:B){
    set.seed(seed*b)
    ## Randomly draw data from multivariate normal distribution
    mat_x_seq <- mvrnorm(n=nplan, mu=c(rep(0,m*pi0), rep(effect,m*(1-pi0))), Sigma=V) 
    
    ## Compute the mom e-processes corresponding to the elementary hypotheses mu_j<=0.
    # use functions from package "safestats" to do so
    mat_e_seq <- matrix(nrow=(nplan-n_min+1), ncol=m)
    for(n in n_min:nplan){
      for(j in 1:m){
        mat_e_seq[(n-n_min+1),j] <- safeTTest(mat_x_seq[1:n,j], designObj = design_obj)$eValue
      }
    }
    
    ## Compute anytime-valid simultaneous lower confidence bounds for the TDP based on the e-to-p process
    tdp_bounds_p <-seq_TDP_p_s(mat_e_seq=mat_e_seq,alpha=alpha, setSs=list_sets, n_min=11,N=nplan)
    E_vals[["e-to-p"]][[b]] <- tdp_bounds_p$seq_p
    # save the bounds in the output list
    for(i in 1:3){
      TDP_bounds[["e-to-p"]][[i]][(11:nplan),b] <- tdp_bounds_p$seq_TDP[,i]
    }
    
    ## Compute anytime-valid simultaneous lower confidence bounds for the TDP based on the mom process
    tdp_bounds_mom <- seq_TDP_mom_s(mat_e_seq=mat_e_seq,alpha=alpha, setSs=list_sets,  n_min=11,N=nplan)
    E_vals[["mom"]][[b]] <- tdp_bounds_mom$seq_e
    # save the bounds in the output list
    for(i in 1:3){
      TDP_bounds[["mom"]][[i]][(11:nplan),b] <- tdp_bounds_mom$seq_TDP[,i]
    }
    
    ## Compute anytime-valid simultaneous lower confidence bounds for the TDP based on the e-to-p process
    tdp_bounds_arc_mom <- seq_TDP_ARC_mom_s(mat_e_seq=mat_e_seq,alpha=alpha, setSs=list_sets, n_min=11, N=nplan)
    E_vals[["mom_ARD"]][[b]] <- tdp_bounds_arc_mom$seq_e
    # save the bounds in the output list
    for(i in 1:3){
      TDP_bounds[["mom_ARD"]][[i]][(11:nplan),b] <- tdp_bounds_arc_mom$seq_TDP[,i]
    }
    
    ## Compute simultaneous lower confidence bounds for the TDP using all resolution inference
    tdp_bounds_hom <- seq_TDP_hom(mat_x_seq, alpha=alpha, setSs=list_sets, burn_in=10)
    # save the bounds in the output list
    for(i in 1:3){
      TDP_bounds[["ARI"]][[i]][(11:nplan),b] <- tdp_bounds_hom[,i]
    }
  }
  
  # save the TDP bounds
  file_out <- paste("mu=", effect, "_rho=", rho, "_size_R=", size_R,".RData", sep="")
  save(TDP_bounds, E_vals, file=file_out)
  
  # If run in parallel: add description to each list containing the TDP bounds
  desc <- paste("effect", effect, "rho", rho, "sizeR",size_R, sep="_")
  out <- list(TDP_bounds)
  names(out) <- desc
  return(out)
}

###
# Empirical non-coverage rate
###

## Input
# bounds: vector containing lower confidence bounds for the true TDP
# true_pi: the true TDP

## Output
# The empirical non-coverage rate, i.e., the proportion of bounds that are larger than the true TDP.

# Notes
# The function rounds the bounds (10 digits) to avoid errors due to maschine inaccuracy

non_cov_rate <- function(bounds, true_pi){
  bounds <- round(bounds, digits=10)
  return(sum(bounds>true_pi)/1000)
}


# Functions for the analysis of pre-processed fMRI data -------------------
# For this, we assume that the first level analysis has already been done.
# The results should be saved as a four-dimensional array including the effect sizes per voxel and per subject.
# The effect sizes per voxel are given in a three dimensional array, so that the location of the voxel is defined by dim_1, dim_2 and dim_3 of the corresponding array entry
# The fourth dimension of the array are the subjects

###
# Compute mom e-process per voxel
###

## Input
# data_eff: four dimensional array containing the effect sizes per voxel and per subject
# coord: vector of x,y,z coordinate of the considered voxel
# filt: stopping times (in number of already observed subjects) at which the e-process should be computed
# design_obj: the design object for that voxel, can be NULL
# delta_min: if design_obj is NULL, then the minimum relevant effect size for this voxel
# beta: acceptable type 2 error. 1- beta defines the power

## Output
# loc: location of the voxel, i.e., the coordinates given by coord.
# e_proc: the e-process corresponding to this voxel

e_process_per_voxel <- function(coord, data_eff, filt, design_obj=NULL, delta_min=NULL, beta=NULL, alt="greater"){
  coord <- unlist(coord)
  x <- coord[1]
  y <- coord[2]
  z <- coord[3]
  data_eff <- data_eff[x,y,z,]
  gc()
  ts <- length(filt)
  if(sum(data_eff==0)==length(data_eff)){
    # if the effect size is equal to zero for all observations, return e-process = 0
    # This can happens if the considered coordinate is not within the brain
    # That locations that are outside the brain are considered happens because at this step, no brain mask has been used.
    return(list(loc=coord,
                e_proc=rep(0, ts)))
  }
  if(sum(data_eff==data_eff[1])==length(data_eff)){
      # This means the data at this voxel has variance 0 and hence the SafeTTest cannot be computed
      # Return NA in this case
      return(list(loc=coord,
                  e_proc=rep(NA, ts)))
  }
  # use tryCatch function here, to ensure that function can run in parallel in a loop even if errors occur.
  # Indicate errors with NA
  tryCatch({
    # Initialize vector to store e-process
    e_proc <- numeric(ts)
    
    # If design_obj for SafeTTest is not given: define it
    if(is.null(design_obj)){
      design_obj <- designSafeT(deltaMin = d_min, beta=0.2, alternative = alt, pb=F)
    }
    
    # at each stopping time in filt: compute the e-process value
    for(t in 1:ts){
      e_proc[t] <- safeTTest(data_eff[c(1:filt[t])], designObj = design_obj)$eValue
    }
    return(list(loc=coord,
                e_proc=e_proc))
  },
  error=function(e){
    return(list(loc=coord,
                e_proc=rep(NA, ts)))
  })
}

###
# Compute the anytime-valid simultaneous TDP bounds based on the mom e-process in pre-defined ROIs
###

## Input
# e_proc_array: four dimensional array containing the e-processes per voxel.
#               The first three dimensions are the spatial dimensions of the fMRI brain, the last dimensions are the considered stopping times of the e-process
# atlas: three dimensional array with same dimensions as fMRI brain, contains for each voxel a number indicating the brain region (ROI) the voxel belongs to; 0 if it belongs to no ROI
# alpha: confidence level
# ROI: vector indicating for which brain regions in the atlas TDP bounds should be computed. If NULL, bounds are computed for all brain regions given in the atlas
#       numbers in ROI have to coincide with the numbers in atlas
# exclude_0: True/False, indicates whether voxels belonging to no brain region should be removed from set of all considered voxels.
# ROI_only: if NULL: all available voxels are used (without 0 if wanted). Otherwise, a vector specifying which ROIs are used (in this case: bounds are only simultaneously valid in these ROIs)

## Output
# seq_TDP: matrix with TDP bounds without ARD guarantee as entries for each considered stopping time (rows) and each ROI (columns)
# seq_TDP_ARD: same as seq_TDP but with ARD guarantee
# e_proc_mat: matrix containing in the columns the location, ROI and observed e-processes per voxel.
# mat_ell_rej: matrix containing the anytime-valid simultaneous upper confidence bounds for the number of false discoveries (non-increasing)

tdp_bounds_fmri <- function(e_proc_array, atlas, alpha=0.2, ROI=NULL, exclude_0=T, ROI_only=NULL){
  dims <- dim(atlas)
  if(!identical(dims, dim(e_proc_array)[1:3])){stop("dimension of atlas and observations have to be identical")}
  if(length(dim(e_proc_array))!=4){stop("e_proc_array has to be four dimensional")}
  
  ## Initialize a matrix used in the computation of the e-processes with columns for coordinates, ROIs and stopping times
  ts <- dim(e_proc_array)[4]
  
  # Build all coordinate triples at once
  coords <- as.matrix(expand.grid(1:dims[1], 1:dims[2], 1:dims[3]))
  
  # Extract atlas values in the same order, i.e., flatten atlas to two dimension
  atlas_vals <- as.vector(atlas)
  
  # Extract e_proc_array values, flattened in the right order (flatten into two dimensions)
  e_vals <- matrix(aperm(e_proc_array, c(4,1,2,3)), nrow = ts) |> t()
  
  # Combine into final matrix
  e_proc_mat <- cbind(coords, atlas_vals, e_vals)
  
  colnames(e_proc_mat) <- c("x", "y", "z","ROI", paste("n", c(1:ts),sep=""))
  
  ## remove the voxels for which the e-process could not be computed, i.e, e_process is equal to zero everywhere or NA
  # these voxels are not considered to be part of the brain.
  id.na <- which(is.na(e_proc_mat[,5]))
  id.zero <- which(apply(e_proc_mat[,c(5:(ts+4))], 1, function(r) all(r == 0)))
  e_proc_mat <- e_proc_mat[-c(id.na, id.zero),]
  
  ## if voxels belonging to no ROI (atlas=0) are to be removed
  # this might differ from the voxels not considered to be part of the brain above, since it could happen that we find activation in the white matter which is not considered in an atlas 
  if(exclude_0){
    id.zero <- which(e_proc_mat[,4]==0)
    e_proc_mat <- e_proc_mat[-c(id.zero),]
  }
  rm(id.na)
  rm(id.zero)
  
  if(!is.null(ROI_only)){
    id.considered <- which(e_proc_mat[,4]%in% ROI_only)
    e_proc_mat <- e_proc_mat[id.considered,]
  }
  
  ## compute the TDP for the ROIs
  # if no ROI is specified, use all available ones
  if(is.null(ROI)){
    ROI <- sort(unique(e_proc_mat[,4]))
  }
  
  ## initialize matrix for sequential TDP bounds
  mat_TDP_mom <- mat_ell_rej <- mat_TDP_ARD <- matrix(nrow=ts, ncol=length(ROI)) #one column is one ROI
  
  ## Compute the anytime-valid simultaneous lower confidence bounds for the TDP:
  # first stopping time: bounds based on mom process and mom_ARD process are identical 
  for(s in 1:length(ROI)){
    # find the rows in e_proc_mat that correspond to the considered ROI
    id_ROI <- which(e_proc_mat[,4]==ROI[s])
    # 1) Compute upper confidence bounds for the number of false discoveries
    mat_ell_rej[1,s] <- m0_bound_average(vec_e=e_proc_mat[,5], setS=id_ROI, alpha=alpha)
    # 2) Compute the lower confidence bounds for the TDP
    mat_TDP_mom[1,s]<- mat_TDP_ARD[1,s]<- 1-mat_ell_rej[1,s]/length(id_ROI)
  }

  # second observation time points, mom and mom_ARD might differ
  if(ts>1){
    for(n in 2:ts){
      for(s in 1:length(ROI)){
        # find the rows in e_proc_mat that correspond to the considered ROI
        id_ROI <- which(e_proc_mat[,4]==ROI[s])
        # 1) Compute upper confidence bounds for the number of false discoveries using mom process
        mat_TDP_mom[n,s]<-  1-(m0_bound_average(vec_e=e_proc_mat[,(4+n)], setS=id_ROI, alpha=alpha))/length(id_ROI)
        # 2) Compute upper confidence bounds for the number of false discoveries using mom_ARD process
        mat_ell_rej[n,s] <- m0_bound_average(vec_e=e_proc_mat[,(4+n)], setS=id_ROI, alpha=alpha, ell_rej= mat_ell_rej[(n-1),s])
        mat_TDP_ARD[n,s]<- 1-mat_ell_rej[n,s]/length(id_ROI)
      }
    }
  }
  
  return(list(seq_TDP=mat_TDP_mom,
              seq_TDP_ARD = mat_TDP_ARD,
              e_proc_mat = e_proc_mat,
              mat_ell_rej= mat_ell_rej))
  
}

