################## 
# Simulation Study
##################
require(safestats)
require(hommel)

###
# Compute anytime-valid simultaneous upper confidence bound for the number of 
# false discoveries based on e-processes at time n as for a fixed discovery set.
###

## Input:
# vec_e:    Vector of the sorted values of the e-processes at time n corresponding 
#           to all considered elementary hypotheses.
# setS:     Indices of the (sorted) elementary hypotheses that are discoveries.
# alpha:    Confidence level
# ell_rej:  c_alpha at time (n-1) for the discovery set setS. 
#           If ell_rej=NULL: ell_rej is set to be the number of discoveries. 

## Output:
# c_alpha:  Upper (1-alpha) confidence bound for the number of false discoveries 
#           at time n for the discovery set setS.

m0_bound_average <- function(vec_e, setS, alpha, ell_rej=NULL){
  S_e <- vec_e[setS]
  # Define tilde(h)(n-1) (see Lemma 3.2)
  if(!is.null(ell_rej)){
    h <- ell_rej
  }else{
    h <- length(S_e)
  }
  c_alpha <- 0
  
  # Compute the cumulative sums over the ordered e-process values 
  # corresponding to discovieries
  cumsum_S_e <- cumsum(S_e[1:h])
  
  ## If all elementary hypotheses are discoveries:
  # Iteratively decrease h=tilde(h)(n-1),...,1 until the hypothesis of size h
  # with the smallest corresponding e-process value is not rejected
  if(length(setS)==length(vec_e)){
    while(c_alpha==0 & h >=1){
      if(cumsum_S_e[h]<(h/alpha)){
        c_alpha <- h
      }else{
        h <- h-1
      }
    }
  }else{
    ## If some elementary hypotheses are not discoveries:
    # Compute discovery set specific threshold as in Lemma 3.2 
    BarS_e <- vec_e[-setS]
    max_Cumsum_e <- max((1:length(BarS_e))/alpha-cumsum(BarS_e))
    
    # Iteratively decrease h=tilde(h)(n-1),...,1 until the hypothesis of size h
    # intersecting only discoveries with the smallest corresponding e-process value 
    # is not rejected by closed testing
    while(c_alpha== 0 & h>=1){
      # Check: If H_K_H is not rejected by local level-alpha test it 
      #        is also not rejected by closed testing.
      if(cumsum_S_e[h]<(h/alpha)){
        c_alpha <- h
      }else{
        if((cumsum_S_e[h]-h/alpha)<max_Cumsum_e){ 
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


# (Anytime-valid) simultaneous lower confidence bounds for the TDP ------------------------------------

###
# Carefree anytime-valid simultaneous lower confidence bounds based on an e-process 
###

## Input:
# mat_e_seq:  Matrix containing the e-processes corresponding to each of the m 
#             elementary hypotheses. 
# alpha:      Confidence level of the bounds.
# setSs:      List of discovery sets to be considered at each time point 
#             (fixed over time).
#             Each list element contains a vector of indices indicating the 
#             hypotheses in the discovery set.
# n_min:      Minimal sample size for which the bounds are to be computed.
# N:          Maximal sample size for which the bounds are to be computed.

## Output:
# seq_TDP:    Matrix containing the TDP bounds at each considered time point (rows) 
#             for each discovery set (columns).
# seq_e:      Matrix containing the e-processes, equal to mat_e_seq.

## Notes:
# The matrix mat_e_seq has N-n_min+1 rows and m columns.
# Lemma 3.2 is utilized to compute the bounds.
# The bounds are computed at time n=n_min,n_min+1,...,N.

seq_TDP_e <- function(mat_e_seq,alpha, setSs,  n_min,N ){
  m <- ncol(mat_e_seq)
  n_sets <- length(setSs)
  length_sets <- sapply(setSs, length)
  
  # Initialize matrix to save the sequential TDP bounds and 
  # the upper bounds for the number of false rejections (ell_rej)
  mat_TDP <- mat_ell_rej <- matrix(nrow=(N-n_min+1), ncol=n_sets)
  
  ## Compute bounds for every observation time point, starting with n=n_min
  for(n in n_min:N){
    ##  At each time point:
    # sort the e-values and remember their original indices
    id_ordered <- order(mat_e_seq[(n-n_min+1),])
    vec_e_sort <- mat_e_seq[(n-n_min+1),id_ordered]
    
    # update the rejection sets to be in terms of the rejected hypotheses
    # (id_ordered tells me the original position. Therefore, I need to know at which position in ordered trial the old indices are)
    u_setSs <- vector("list", length=n_sets)
    for(s in 1:n_sets){
      u_setSs[[s]] <- which(id_ordered%in%setSs[[s]])
    }
      
    # Compute the (non-decreasing) TDP bounds for each considered discovery set 
     if(n==n_min){
       # At the first time point, no previous rejections are known
      for(s in 1:n_sets){
        mat_ell_rej[1,s] <- m0_bound_average(vec_e=vec_e_sort, 
                                             setS=u_setSs[[s]], 
                                             alpha=alpha)
        mat_TDP[1,s] <- 1-mat_ell_rej[1,s]/length_sets[s]
      }
    }else{
      # After the first time point, 
      # Utilize knowledge about bounds at previous observation time points
      for(s in 1:n_sets){
        mat_ell_rej[(n-n_min+1),s] <- m0_bound_average(vec_e=vec_e_sort,
                                                       setS=u_setSs[[s]], 
                                                       alpha=alpha, 
                                                       ell_rej =mat_ell_rej[(n-n_min),s])
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
# mat_x_seq:  Matrix containing the observed values, 
#             with ncol=m, nrow=n_obs, which are the imput for a one-sample t-test
# alpha:      Confidence level of the bounds.
# setSs:      List of discovery sets to be considered at each time point 
#             (fixed over time).
#             Each list element contains a vector of indices indicating the 
#             hypotheses in the discovery set.
# n_min:      Minimal sample size for which the bounds are to be computed.
# N:          Maximal sample size for which the bounds are to be computed.

## Output
# seq_TDP:    Matrix containing the TDP bounds at each considered time point (rows) 
#             for each discovery set (columns).
# seq_p:      Matrix containing the p-values at each time point.

## Notes
# These bounds are not anytime-valid.
# Computation of the bounds utilizes the computational shortcut by Goeman et al (2019), 
# implemented in the "hommel" package

seq_TDP_hom <- function(mat_x_seq, alpha, setSs, n_min){ 
  n_obs <- nrow(mat_x_seq)
  m <- ncol(mat_x_seq) 
  n_sets <- length(setSs)
  
  # Initialize matrix to save p-values at each time point
  mat_p_seq <- matrix(nrow=(n_obs-(n_min-1)), ncol=m)
  
  # Initialize matrix to save sequential TDP bounds
  mat_TDP <- matrix(nrow=(n_obs-(n_min-1)), ncol=n_sets)
  
  # Compute bounds for every observation time point, starting with n=n_min
  for(n in n_min:n_obs){
    # At each time point:
    # 1) Compute the p-values corresponding to the m elementary null hypotheses
    for(j in 1:m){
      mat_p_seq[(n-(n_min-1)),j] <- t.test(x= mat_x_seq[c(1:n),j],
                                           alternative="greater")$p.value
    }
    # 2) Compute the TDP bounds for each considered discovery set
    hom <- hommel(mat_p_seq[(n-(n_min-1)),], simes=T)
    for(s in 1:n_sets){
      mat_TDP[(n-(n_min-1)),s] <- tdp(hom, ix=setSs[[s]], alpha=alpha)
    }
  }
  return(list(seq_TDP=mat_TDP,
              seq_p = mat_p_seq))
}

# Functions used in the simulation study -----------------------------------------

###
# Simulating one scenario
###

## Input
# rho:        Level of dependency between the elementary hypotheses.
# size_R:     Size of the rejection sets.
# effect:     Mean under the alternative.
# design_obj: Design object needed to compute e-processes. 
#             Does not depent on the observed data.
# seed:       Set a seed to make study reproducible
# nplan:      Maximal sample size for which the bounds are computed.
# B:          Number of Monte-Carlo iterations.
# alpha:      Confidence level of the bounds.
# n_min:      Minimal sample size for which the bounds are computed.
# save_e_proc:Should the computed e-processes (and p-values) be saved seperately?

## Output
# List containing the (anytime-valid) simultaneous lower confidence bounds based 
# on the mom process and ARI
# The list is structured as follows: [[Procedure]][[true TDP]][[TDP_bounds]]
# TDP_bounds are saved in a matrix with N rows and B columns

## Notes:
# If save_e_proc=True: the e-processes for each elementary hypothesis are saved 
# in the working directory.
# The file is named "mu=[effect]_rho=[rho]_size_R=[size_R].RData"
# Only rejection sets with a true TDP of 0.1, 0.5 and 0.9 are computed, 
# Should this be adapted: change m0, list_sets and for-loops over the 
# discovery sets (indexed by i) in the function

sim_tdp_bounds_fun <- function(rho,size_R, effect, design_obj, seed, 
                               nplan=100, B=1000, alpha=0.2, n_min=11,
                               save_e_proc=F){
  set.seed(seed)
  # Initialize parameters for simulation of the data
  m <- 1000
  pi0 <- 0.5 # half of the m elementary hypotheses are true.
  V <- matrix(rho,ncol=m, nrow=m)
  diag(V) <- 1
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
  TDP_bounds <- lapply(c(1:2), make_list_matrix, nplan=nplan, B=B)
  P_E_vals<- lapply(rep(B,2), vector, mode="list")
  names(TDP_bounds)<- names(P_E_vals) <- c("mom", "ARI")
  
  ## Simulate the data and compute the e-processes and TDP bounds
  # Do this for B independent Monte Carlo repetitions
  for(b in 1:B){
    set.seed(seed*b)
    # Randomly draw data from multivariate normal distribution
    mat_x_seq <- mvrnorm(n=nplan, mu=c(rep(0,m*pi0), rep(effect,m*(1-pi0))), Sigma=V) 
    
    # Compute the mom e-processes corresponding to the elementary hypotheses mu_j<=0.
    mat_e_seq <- matrix(nrow=(nplan-n_min+1), ncol=m)
    for(n in n_min:nplan){
      for(j in 1:m){
        mat_e_seq[(n-n_min+1),j] <- safeTTest(mat_x_seq[1:n,j], 
                                              designObj = design_obj)$eValue
      }
    }
    
    # Compute carefree anytime-valid simultaneous lower confidence bounds 
    # for the TDP based on the mom e-process
    tdp_bounds_mom <- seq_TDP_e(mat_e_seq=mat_e_seq,alpha=alpha, 
                                setSs=list_sets, n_min=n_min, N=nplan)
    P_E_vals[["mom"]][[b]] <- tdp_bounds_mom$seq_e
    
    # Save the bounds in the output list
    for(i in 1:3){
      TDP_bounds[["mom"]][[i]][(n_min:nplan),b] <- tdp_bounds_mom$seq_TDP[,i]
    }
    
    # Compute simultaneous lower confidence bounds for the TDP using 
    # all resolution inference
    tdp_bounds_hom <- seq_TDP_hom(mat_x_seq, alpha=alpha, 
                                  setSs=list_sets, n_min=n_min)
    P_E_vals[["ARI"]][[b]] <- tdp_bounds_hom$seq_p
    # Save the bounds in the output list
    for(i in 1:3){
      TDP_bounds[["ARI"]][[i]][(n_min:nplan),b] <- tdp_bounds_hom$seq_TDP[,i]
    }
  }
  # Save the e-processes (p-values) in extra file
  if(save_e_proc){
    file_out <- paste("P_E_values_mu=", effect, "_rho=", rho, "_size_R=", size_R,".RData", sep="")
    save(P_E_vals, file=file_out)
  }
  
  # Add description to each list containing the TDP bounds
  desc <- paste("effect", effect, "rho", rho, "sizeR",size_R, sep="_")
  out <- list(TDP_bounds)
  names(out) <- desc
  return(out)
}

###
# Compute the empirical non-coverage rate
###

## Input
# bounds: vector containing lower confidence bounds for the true TDP
# true_pi: the true TDP

## Output
# The empirical non-coverage rate, i.e., the proportion of lower confidence bounds 
# for the TDP that are larger than the true TDP.

# Notes
# The function uses rounding (10 digits) to avoid errors due to machine inaccuracy

non_cov_rate <- function(bounds, true_pi){
  bounds <- round(bounds, digits=10)
  return(sum(bounds>true_pi)/1000)
}

################## 
# Case Study
##################
require(safestats)
require(hommel)
require(fmri)
require(oro.nifti)

# Functions for the subject-level analysis -------------------
###
# Get estimated effect-size per voxel per subject
###

## Input:
# subject_vec:  Vector containing the names of the nii.gz files per subject, 
#               saved in the folder: data.
# stim:         Convoluted HRF model and condition onset time series,
#               Output of fmri.stimulus (package: fmri)
# conf_vec:     Vector containing the names of the file containing the confounding 
#               time series per subject, saved in the folder: data
# conf_para:    Vector containing the names of the confounding parameters to be used
#               If null: use none
# save_subject: Should the subject specific results be automatically saved ? 

## Ouput
# effects:      A four dimensional array containing the effect-sizes per voxel 
#               for all n subjects.
#               The first three dimensions correspond to the location of the voxels
#               The fourth dimension correspond to the subjects

## Notes
# If the nii.gz files per subject or the confounding time series per subject
#  are not saved in a folder called data or if this folder is a subfolder in the 
#  current working directory: change the file.path() command defining the variable 
#  nii_name and cef_path, respectively.
# If save_subject is true, the subject specific effect size maps are saved in the
#  current working directory
# The effect of interest is the contrast between the first and second condition,
#  i.e., the contrast between the semantic condition and the control condition.

sub_effect_size <- function(subject_vec, stim, conf_vec =NULL,conf_para=NULL, save_subject =F){
  n <- length(subject_vec)
  design_mat <- design_mat_i <- fmri.design(stimulus= stim,
                                            order=2)
  
  
  # Initialize 4 dimensional array to save the effect sizes of all subjects
  effects_per_sub <- array(NA, dim=c(97,115,97,n)) 
  
  ## Estimate the effect-size per voxel for each subject
  for(i in 1:n){
    # Initialize subject-specific design matrix if confounding effects are included
    # in the design matrix
    if(!is.null(conf_para)){
      cef_path <- file.path("data", conf_vec[i])
      cef <- read.delim(cef_path)
      id.cef <- which(names(cef)%in%conf_para)
      design_mat_i <- cbind(design_mat, cef[,id.cef]) 
    }
    # Load the preprocessed BOLD signal
    nii_name <- file.path("data",subject_vec[i])
    dss <- readNIfTI(nii_name)
    gc()
    # Convert nifti object to fmridata object
    dss <- niftiImage2fmri(niftiobj=dss, level = 0.75, mask=NULL, setmask = F, indx = NULL,
                           indy = NULL, indz = NULL, avoidnegs = FALSE)
    gc()
    # Estimate the effect sizes per subject
    out_sub <- fmri.lm(ds=dss, 
                       z=design_mat_i, 
                       mask = NULL,
                       actype = c("smooth"),
                       contrast = c(1,-1), 
                       verbose = FALSE)
    gc()
    # If wanted: save subject specific results
    effects_per_sub[,,,i]<- effect <-out_sub$cbeta/(sqrt(out_sub$var))
    if(save_subject){
      file_out <- paste(subject_vec[i], "_effect_per_voxel.RData", sep="")
      save(effect,file=file_out)
    }
    rm(effect)
    gc()
  }
  return(effects= effects_per_sub)
}




# Functions for the group-level analysis ----------------------------------
###
# Compute mom e-process for one voxel based on results of subject-level analysis
###

## Input
# data_eff:   Four dimensional array containing the effect sizes per voxel and 
#             per subject
# coord:      Vector of x,y,z coordinate of the considered voxel
# filt:       Stopping times at which the e-process should be computed
# design_obj: The design object for that voxel required to compute the e-process, 
#             can be NULL.
# delta_min:  If design_obj is NULL, then the minimum relevant effect size.
# beta:       If design_obj is NULL, the acceptable type 2 error

## Output
# loc:        Location of the voxel, i.e., the coordinates given by coord.
# e_proc:     The e-process corresponding to this voxel

## Notes
# This function is written so that it can be used in a foreach loop

e_process_per_voxel <- function(coord, data_eff, filt, design_obj=NULL, 
                                delta_min=NULL, beta=NULL, alt="greater"){
  coord <- unlist(coord)
  x <- coord[1]
  y <- coord[2]
  z <- coord[3]
  data_eff <- data_eff[x,y,z,]
  gc()
  ts <- length(filt)
  
  # Check if the effect size is equal to zero for all observations.
  # This occurs if the considered coordinate is not within the brain.
  # Return e-process = 0 at all stopping times.
  if(sum(data_eff==0)==length(data_eff)){
    return(list(loc=coord,
                e_proc=rep(0, ts)))
  }
  
  # Check which estimated effect-size has variance 0. 
  # At this voxel, the SafeTTest (i.e., the e-process) cannot be computed
  # Return NA 
  if(sum(data_eff==data_eff[1])==length(data_eff)){
      return(list(loc=coord,
                  e_proc=rep(NA, ts)))
  }
  
  # Use tryCatch function to ensure that function can run in parallel in a loop 
  # even if errors occur.
  # Indicate errors with NA.
  tryCatch({
    # Initialize vector to store e-process
    e_proc <- numeric(ts)
    
    # If design_obj for SafeTTest is not given: define it
    if(is.null(design_obj)){
      design_obj <- designSafeT(deltaMin = d_min, beta=beta, alternative = alt, pb=F)
    }
    
    # at each stopping time: compute the e-process value
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
# Compute anytime-valid simultaneous upper confidence bound for the number of 
# false discoveries based on e-processes at time n as for a fixed discovery set.
###

## Input:
# vec_e:    Vector of the sorted values of the e-processes at time n corresponding 
#           to all considered elementary hypotheses.
# setS:     Indices of the (sorted) elementary hypotheses that are discoveries.
# alpha:    Confidence level
# ell_rej:  c_alpha at time (n-1) for the discovery set setS. 
#           If ell_rej=NULL: ell_rej is set to be the number of discoveries. 

## Output:
# c_alpha:  Upper (1-alpha) confidence bound for the number of false discoveries 
#           at time n for the discovery set setS.

m0_bound_average <- function(vec_e, setS, alpha, ell_rej=NULL){
  S_e <- vec_e[setS]
  # Define tilde(h)(n-1) (see Lemma 3.2)
  if(!is.null(ell_rej)){
    h <- ell_rej
  }else{
    h <- length(S_e)
  }
  c_alpha <- 0
  
  # Compute the cumulative sums over the ordered e-process values 
  # corresponding to discovieries
  cumsum_S_e <- cumsum(S_e[1:h])
  
  ## If all elementary hypotheses are discoveries:
  # Iteratively decrease h=tilde(h)(n-1),...,1 until the hypothesis of size h
  # with the smallest corresponding e-process value is not rejected
  if(length(setS)==length(vec_e)){
    while(c_alpha==0 & h >=1){
      if(cumsum_S_e[h]<(h/alpha)){
        c_alpha <- h
      }else{
        h <- h-1
      }
    }
  }else{
    ## If some elementary hypotheses are not discoveries:
    # Compute discovery set specific threshold as in Lemma 3.2 
    BarS_e <- vec_e[-setS]
    max_Cumsum_e <- max((1:length(BarS_e))/alpha-cumsum(BarS_e))
    
    # Iteratively decrease h=tilde(h)(n-1),...,1 until the hypothesis of size h
    # intersecting only discoveries with the smallest corresponding e-process value 
    # is not rejected by closed testing
    while(c_alpha== 0 & h>=1){
      # Check: If H_K_H is not rejected by local level-alpha test it 
      #        is also not rejected by closed testing.
      if(cumsum_S_e[h]<(h/alpha)){
        c_alpha <- h
      }else{
        if((cumsum_S_e[h]-h/alpha)<max_Cumsum_e){ 
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

###
# Compute the anytime-valid simultaneous bounds for the number of true discoveries 
# based on the mom e-process in pre-defined regions of interest
###

## Input
# e_proc_array: Four dimensional array containing the e-processes per voxel.
#               The first three dimensions are the spatial dimensions, 
#               the last dimensions are the considered stopping times of the e-process
# atlas:        Three dimensional array with same dimensions as the spatial dimensions of e_proc_array, 
#               assigns each voxel an integer indicating the brain region (ROI) the voxel belongs to.
#               Voxels which are assigned zero are not considered to be part of the brain
# alpha:        Confidence level of the bounds for the number of true discoveries
# ROI:          Vector indicating which brain regions in the atlas are discovery sets. 
#               Integers in ROI have to coincide with the integers in atlas.
#               If ROI=NULL, bounds are computed for all brain regions simultaneously.
# exclude_0:    True/False, indicates whether voxels belonging to no brain region in the atlas should 
#               be removed from set of all considered voxels.

## Output
# seq_TD:       Matrix with TD bounds for each considered stopping time (rows) and each ROI (columns)
# e_proc_mat:   Matrix containing in the columns the location, ROI and observed e-processes per voxel.
# size_per_ROI: Vector containing the size of each considered ROI.

td_bounds_fmri <- function(e_proc_array, atlas, alpha=0.2, ROI=NULL, exclude_0=T){
  # Check that the atlas and the e-process array have the same dimensions
  dims <- dim(atlas)
  if(!identical(dims, dim(e_proc_array)[1:3])){stop("dimension of atlas and observations have to be identical")}
  if(length(dim(e_proc_array))!=4){stop("e_proc_array has to have four dimensional")}
  
  # Initialize a matrix used in the computation of the e-processes 
  # columns for coordinates, ROIs and stopping times
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
  
  # remove the voxels for which the e-process could not be computed, 
  # i.e, e_process is NA
  # These voxels are not considered to be part of the brain.
  id.na <- which(is.na(e_proc_mat[,5]))
  e_proc_mat <- e_proc_mat[-c(id.na),]
  rm(id.na)
  
  # if voxels belonging to no ROI (atlas=0) are to be removed:
  # this might differ from the voxels not considered to be part of the brain above, 
  # since we might find activation in the white matter 
  # which is not considered as part of the brain in an atlas 
  if(exclude_0){
    id.zero <- which(e_proc_mat[,4]==0)
    e_proc_mat <- e_proc_mat[-c(id.zero),]
    rm(id.zero)
  }
  
  # compute the TDP for the ROIs
  # if no ROI is specified, use all ROIs in the atlas
  if(is.null(ROI)){
    ROI <- sort(unique(e_proc_mat[,4]))
  }
  
  # Initialize matrix for sequential TDP bounds
  mat_TDP_mom <- mat_m0_bounds <- matrix(nrow=ts, ncol=length(ROI)) #one column corresponds to one ROI
  size_ROI <- numeric(length(ROI))
  
  # Compute the anytime-valid simultaneous lower confidence bounds for the TDP:
  # first stopping time: 
  for(s in 1:length(ROI)){
    # Sort the e-process values at the first stopping time:
    e_proc_sort <- e_proc_mat[order(e_proc_mat[,5]),]
    # find the rows in e_proc_mat that correspond to the considered ROI
    id_ROI <- which(e_proc_sort[,4]==ROI[s])
    size_ROI[s] <-length(id_ROI)
    # 1) Compute upper confidence bounds for the number of false discoveries
    mat_m0_bounds[1,s] <- m0_bound_average(vec_e=e_proc_sort[,5], 
                                           setS=id_ROI, 
                                           alpha=alpha)
    # 2) Compute the lower confidence bounds for the number of true discoveries
    mat_TDP_mom[1,s]<- length(id_ROI)-mat_m0_bounds[1,s]
    rm(e_proc_sort)
  }

  if(ts>1){
    for(n in 2:ts){
      for(s in 1:length(ROI)){
        # Sort the e-process values at the first stopping time:
        e_proc_sort <- e_proc_mat[order(e_proc_mat[,(4+n)]),]
        # find the rows in e_proc_mat that correspond to the considered ROI
        id_ROI <- which(e_proc_sort[,4]==ROI[s])
        #1) Compute upper confidence bounds for the number of false discoveries using mom process
        mat_m0_bounds[n,s] <- m0_bound_average(vec_e=e_proc_sort[,(4+n)], 
                                               setS=id_ROI, 
                                               alpha=alpha, 
                                               ell_rej= mat_m0_bounds[(n-1),s])
        #2) Compute the lower confidence bounds for the number of true discoveries
        mat_TDP_mom[n,s]<- length(id_ROI)-mat_m0_bounds[n,s]
        rm(e_proc_sort)
      }
    }
  }
  
  return(list(seq_TD=mat_TDP_mom,
              e_proc_mat = e_proc_mat,
              size_per_ROI= size_ROI))
  
}
