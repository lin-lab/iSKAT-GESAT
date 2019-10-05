# Last Edited: October 4, 2019
# suppressWarnings for call to cor(V,Z)
# fix error where include.nonsingleton.z was not initialized if
# is_check_genotype == FALSE and is_dosage == FALSE.

# Last Edited: October 14, 2015
# Add a check for missingness in Z
# if is_check_genotype==FALSE and is_dosage== FALSE,
# then Z cannot have missigneness
# iSKAT/GESAT core functions cannot have missingness in all variables (the functions work assuming no missingness)
# and the checks for missingness occur in main.R
# in main.R, imputation for missing values takes place only IF is_check_genotype = TRUE OR is_dosage==TRUE
# thus if both of these are false, Z must have no missingness
# updated main.R such that if Z has missingness and is_check_genotype==FALSE and is_dosage== FALSE, it throws an error
#
#
# Last Edited: January 29, 2014
# v1.* of package onwards names the common variant method as GESAT and the rare variant method as iSKAT
# iSKAT_MAIN_Check_Z_and_V() is new in v1.* of iSKAT package
# also removes common variants

iSKAT <- function(Z, Y, E, X=NULL, type="davies",
                  lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=FALSE, plotfile=NA, scale.Z=FALSE, weights.Z=NULL, weights.V=NULL,
                  out_type="C", impute.method = "fixed", is_check_genotype=TRUE, is_dosage=FALSE, missing_cutoff=0.15,
                  r.corr=(0:10)/10, weights.beta=c(1,25), MAF_cutoff=0.05, SetID=NULL){

  # if weights.beta is specified and weights.Z is not, weights.Z will be set to beta weights
  # if weights.beta is specified and weights.V is not, weights.V will be set to beta weights
  # if weights.Z is specified, weights.beta is ignored for Z
  # if weights.V is specified, weights.beta is ignored for V
  # if scale.Z is T, weights.beta is ignored for Z and weights.Z is also ignored
  # NB: in v1.* package onwards, iSKAT main function has three more arguments r.corr,  weights.beta, MAF_cutoff than GESAT main function
  # the default on scale.Z is also different.
  # in iSKAT the default is to NOT scale Z and use weights.beta=c(1,25) weights

  #-------------------------------------
  # check outcome type
  if(out_type != "C" && out_type != "D"){
    stop("Invalid out_type!. Please use only \"C\" for continous outcome and \"D\" for dichotomous outcome.")
  }


  #-------------------------------------
  # check dimensions
  if(nrow(E)!= nrow(Y)) stop("Dimensions of E and Y don't match.")
  if(nrow(Z)!= nrow(Y)) stop("Dimensions of Z and Y don't match.")
  if(is.null(X)==FALSE){
    if(nrow(X)!= nrow(Y)) stop("Dimensions of X and Y don't match.")
    if(class(X)!= "matrix") stop("X is not a matrix.")
  }
  if(class(Z)!= "matrix") stop("Z is not a matrix.")
  if(class(E)!= "matrix") stop("E is not a matrix.")
  if(class(Y)!= "matrix") stop("Y is not a matrix.")


  #----------------------------------------- added on Oct 14, 2015 (start)
  if(is_check_genotype==FALSE & is_dosage== FALSE){
    if(sum(is.na(Z))!= 0) stop("Z cannot have any missing values if is_check_genotype = FALSE and is_dosage=FALSE.")
  }
  #----------------------------------------- added on Oct 14, 2015 (end)

  if(is.null(weights.Z)==FALSE){
    if(ncol(Z)!= length(weights.Z)) stop("Dimensions of Z and weights.Z don't match.")
    if(sum(weights.Z<0)!=0) stop("weights.Z have to be non-negative and cannot be missing.")
  }

  if(is.null(weights.V)==FALSE){
    if(ncol(Z)!= length(weights.V)) stop("dimensions of Z and weights.V don't match")
    if(sum(weights.V<0)!=0) stop("weights.V have to be non-negative and cannot be missing")
  }


  #-------------------------------------
  # check other arguments
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored, beta.weights are also ignored for Z! To use weights as weights.Z or beta.weights for Z, set scale.Z=F")


  #-------------------------------------
  # check for missing values in Y, E, X
  if(is.null(X)==FALSE){
    if(sum(is.na(X))!= 0) stop("X cannot have any missing values.")
  }

  if(sum(is.na(E))!= 0) stop("E cannot have any missing values.")
  if(sum(is.na(Y))!= 0) stop("Y cannot have any missing values.")


  #-------------------------------------
  # check that X doesn't have intercept

  if(is.null(X)==FALSE){

    if(ncol(X)==1){
      if(checkpolymorphic(X)==FALSE) stop("X should not include intercept and must have more than one level.")
    }else{
      if(sum(apply(X, 2, checkpolymorphic))!= ncol(X)) stop("X should not include intercept and must have more than one level.")

    }
  }


  #-------------------------------------
  # check that E has more than one levels
  if(ncol(E)==1){
    if(checkpolymorphic(E)==FALSE) stop("E must have more than one level.")
  }else{
    if(sum(apply(E, 2, checkpolymorphic))!= ncol(E)) stop("E must have more than one level.")

  }


  #-------------------------------------
  # check Z and impute
  if(is_dosage ==TRUE){
    impute.method="fixed"
  }

  if(is_check_genotype==TRUE | is_dosage==TRUE){
    Z.out <- iSKAT_MAIN_Check_Z_and_V(Z=Z, SetID=SetID, weights.Z=weights.Z, weights.V=weights.V, impute.method=impute.method,  										missing_cutoff=missing_cutoff, MAF_cutoff=MAF_cutoff)
    Z <- as.matrix(Z.out$Z.test)
    weights.Z <- as.vector(Z.out$weights.Z.test)
    weights.V <- as.vector(Z.out$weights.V.test)
    include.nonsingleton.Z <- as.vector(Z.out$include.nonsingleton.Z)
  }

  if(is_check_genotype==FALSE & is_dosage==FALSE){
    include.nonsingleton.Z <- c(1:ncol(Z))
  }

  if(is.null(Z)==TRUE){

    if(is.null(SetID)){
      msg <- sprintf("The Z matrix has no SNPs." )
    } else {
      msg <- sprintf("In %s, the Z matrix has no SNPs.", SetID )
    }

    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=0, n.GE.test=NA))
  }

  if(ncol(Z)<2){

    if(is.null(SetID)){
      msg <- sprintf("The Z matrix has fewer than 2 SNPs. Do not need penalized procedures." )
    } else {
      msg <- sprintf("In %s, the Z matrix has fewer than 2 SNPs. Do not need penalized procedures.", SetID )
    }

    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=NA))
  }


  #-------------------------------------
  # check and assign weights
  MAF_for_betaweights <- colMeans(as.matrix(Z), na.rm=T)/2
  MAF_for_betaweights <- pmin(MAF_for_betaweights, 1-MAF_for_betaweights)

  if(is.null(weights.beta)==FALSE&&is.null(weights.Z)==TRUE) weights.Z <- Beta.Weights(MAF_for_betaweights, weights.beta)
  if(is.null(weights.beta)==FALSE&&is.null(weights.V)==TRUE) weights.V <- Beta.Weights(MAF_for_betaweights, weights.beta)


  #-------------------------------------
  # check V
  if(is.null(include.nonsingleton.Z)==TRUE){

    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix is empty.")
    } else {
      msg <- sprintf("In %s, the GxE matrix is empty.", SetID)
    }

    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=0))

  }

  Ztemp <- NULL
  Ztemp <- as.matrix(Z[,include.nonsingleton.Z])
  if(is.null(weights.V)==FALSE){
    weights.V.temp <- weights.V[include.nonsingleton.Z]
    Ztemp <- t(t(Ztemp) * (weights.V.temp))
  }


  if(ncol(E)==1){
    V <- as.matrix(drop(E)*Ztemp)
  }else{
    V <- NULL
    for(hhh in 1:ncol(E)){
      V <- as.matrix(cbind(V, drop(E[,hhh])*Ztemp))
    }
  }

  cor_vz <- suppressWarnings(cor(V, Z))

  Idx.V.check <- apply(V, 2, checkpolymorphic)&(apply(abs(cor_vz)>0.9999999999,1,sum,na.rm=T)==0)

  if(sum(Idx.V.check)==0){

    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix is empty. " )
    } else {
      msg <- sprintf("In %s, the GxE matrix is empty.", SetID )
    }
    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=0))
  }


  if(sum(Idx.V.check)==1){

    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix has fewer than 2 SNPs. " )
    } else {
      msg <- sprintf("In %s, the GxE matrix has fewer than 2 SNPs.", SetID )
    }

    warning(msg,call.=FALSE)
    return(list(pvalue=NA, param=NULL, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=1))
  }

  V <- as.matrix(V[,Idx.V.check])


  #-------------------------------------
  # Run iSKAT
  Xtilde <- as.matrix(cbind(X,E))
  if(out_type == "C"){
    iSKAT.out <- iSKAT.linear(Y=Y, Xtilde=Xtilde, Z=Z, V=V, type=type,
                              lower=lower, upper=upper, nintervals=nintervals, plotGCV=plotGCV, plotfile=plotfile, scale.Z=scale.Z, 				weights.Z=weights.Z, weights.V=NULL, r.corr=r.corr)
  }

  if(out_type == "D"){
    iSKAT.out <- iSKAT.logistic(Y=Y, Xtilde=Xtilde, Z=Z, V=V, type=type,
                                lower=lower, upper=upper, nintervals=nintervals, plotGCV=plotGCV, plotfile=plotfile, scale.Z=scale.Z, 				weights.Z=weights.Z, weights.V=NULL, r.corr=r.corr)
  }

  return(list(pvalue=iSKAT.out$pvalue, param=iSKAT.out$param,
              lambda=iSKAT.out$lambda, n.G.test=drop(ncol(Z)), n.GE.test=drop(ncol(V))))
}









#
#	Check the Z, and do imputation

#
iSKAT_MAIN_Check_Z_and_V <- function(Z, SetID, weights.Z=NULL, weights.V=NULL, impute.method,  missing_cutoff, MAF_cutoff){

  # check.Z.error = 0 : no snps removed, but some snps possibly imputed
  # check.Z.error = 1 : all snps removed, returns NULL matrix for Z
  # check.Z.error = 2 : some snps removed, remainder snps may have been imputed

  check.Z.error <- 0
  n <- nrow(Z)
  ##############################################
  # Recode Missing to NA

  IDX_MISS <- union(which(is.na(Z)), which(Z == 9))
  if(length(IDX_MISS) > 0){
    Z[IDX_MISS] <- NA
  }

  ###################################################
  # Check missing rates and exclude any SNPs with missing rate > missing_cutoff
  # Also exclude non-polymorphic SNPs
  ########################################## modification #1
  m <- ncol(Z)
  ID_INCLUDE_SNP <- NULL
  for(i in 1:m){
    missing.ratio <- length(which(is.na(Z[,i])))/n
    sd1 <- sd(Z[,i], na.rm=TRUE)
    maf1 <- min(mean(Z[,i], na.rm=TRUE)/2, 1- mean(Z[,i], na.rm=TRUE)/2)
    if(is.na(sd1)==FALSE){
      if(missing.ratio < missing_cutoff && sd1 > 0 && maf1<MAF_cutoff){
        ID_INCLUDE_SNP <- c(ID_INCLUDE_SNP,i)
      }
    }
  }


  if(length(ID_INCLUDE_SNP) == 0){

    if(is.null(SetID)){
      msg <- sprintf("ALL SNPs have either high missing rates or no-variation or exceed MAF_cutoff. ")
    } else {
      msg <- sprintf("In %s, ALL SNPs have either high missing rates or no-variation or exceed MAF_cutoff. ",SetID )
    }
    warning(msg, call.=FALSE)

    re <- list(Z.test=NULL, weights.Z.test=NULL, weights.V.test=NULL, check.Z.error=1, include.nonsingleton.Z=NULL)


  } else{

    if(m - length(ID_INCLUDE_SNP) > 0 ){

      if(is.null(SetID)){
        msg <- sprintf("%d SNPs with either high missing rates or no-variation or exceed MAF_cutoff are excluded!", m - length(ID_INCLUDE_SNP))
      } else {
        msg <- sprintf("In %s, %d SNPs with either high missing rates or no-variation or exceed MAF_cutoff are excluded!",SetID, m - length(ID_INCLUDE_SNP) )
      }

      ########################################## end of modification #1

      warning(msg, call.=FALSE)
      Z <- as.matrix(Z[,ID_INCLUDE_SNP])
      if(is.null(weights.Z)==FALSE) weights.Z <- weights.Z[ID_INCLUDE_SNP]
      if(is.null(weights.V)==FALSE) weights.V <- weights.V[ID_INCLUDE_SNP]
      check.Z.error <- 2
      IDX_MISS <- which(is.na(Z))
    }

    ########################################## modification #2
    # Check for non-singleton Z's before imputation

    ID_INCLUDE_SNP_GE <- NULL

    for(iii in 1:ncol(Z)){
      nonsingleton <- checkvariation2(Z[,iii])

      if(nonsingleton==TRUE){
        ID_INCLUDE_SNP_GE <- c(ID_INCLUDE_SNP_GE,iii)
      }
    }


    ########################################## end of modification #2

    ##########################################
    # Missing Imputation
    if(length(IDX_MISS) > 0){
      if(is.null(SetID)){
        msg <- sprintf("The missing genotype rate is %f. Imputation is applied.", (length(IDX_MISS))/length(Z) )
      } else {
        msg <- sprintf("In %s, the missing genotype rate is %f. Imputation is applied.", SetID, (length(IDX_MISS))/length(Z) )
      }

      warning(msg,call.=FALSE)
      Z <- Impute_iSKAT(Z,impute.method)
    }
    re <- list(Z.test=Z, weights.Z.test=weights.Z, weights.V.test=weights.V, check.Z.error=check.Z.error,
               include.nonsingleton.Z=ID_INCLUDE_SNP_GE)

  }

  return(re)
}

