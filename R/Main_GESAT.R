# Last Edited: October 4, 2019
# suppressWarnings for call to cor(V,Z)

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
# Four changes to GESAT()
# also, fixed two bugs related to weights.V in GESAT() main function call:
# (1)a one line bug to if(ncol(Z)!= length(weights.V))
# (2) weights.V are handled by directly multiplying it to Z, before GE matrix=V is formed
# previously if weights.V were specified, there wd be problems due to bugs
# (3)also now all colliner columns with Z are removed added "&(apply(abs(cor(V,Z))>0.9999999999,1,sum,na.rm=T)==0)"
# (4)also, changed to "davies" and "liu"
# iSKAT_Main_Check_Z() is also modified to add a check: if(is.na(sd1)==FALSE)
# Last Edited: January 27, 2014
# v0.6 onwards allow for binary outcome

GESAT <- function(Z, Y, E, X=NULL, type="davies",
                  lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=FALSE, plotfile=NA, scale.Z=TRUE, weights.Z=NULL, weights.V=NULL,
                  out_type="C", impute.method = "fixed", is_check_genotype=TRUE, is_dosage=FALSE, missing_cutoff=0.15, SetID=NULL){

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
    if(ncol(Z)!= length(weights.V)) stop("dimensions of V and weights.V don't match")  #change!!!
    if(sum(weights.V<0)!=0) stop("weights.V have to be non-negative and cannot be missing")
  }


  #-------------------------------------
  # check other arguments
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")


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
    Z.out <- iSKAT_MAIN_Check_Z(Z=Z, SetID=SetID, weights.Z=weights.Z, weights.V=weights.V, impute.method=impute.method,  									missing_cutoff=missing_cutoff)
    Z <- as.matrix(Z.out$Z.test)
    weights.Z <- as.vector(Z.out$weights.Z.test)
    weights.V <- as.vector(Z.out$weights.V.test)
  }

  if(is.null(Z)==TRUE){

    if(is.null(SetID)){
      msg <- sprintf("The Z matrix has no SNPs." )
    } else {
      msg <- sprintf("In %s, the Z matrix has no SNPs.", SetID )
    }

    warning(msg,call.=FALSE)
    return(list(pvalue=NA, Is_converge=NA, lambda=NA, n.G.test=0, n.GE.test=NA))
  }

  if(ncol(Z)<2){

    if(is.null(SetID)){
      msg <- sprintf("The Z matrix has fewer than 2 SNPs. Do not need penalized procedures." )
    } else {
      msg <- sprintf("In %s, the Z matrix has fewer than 2 SNPs. Do not need penalized procedures.", SetID )
    }

    warning(msg,call.=FALSE)
    return(list(pvalue=NA, Is_converge=NA, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=NA))
  }



  #-------------------------------------
  # check V
  Ztemp <- Z
  if(is.null(weights.V)==FALSE){
    Ztemp <- t(t(Z) * (weights.V))
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
    return(list(pvalue=NA, Is_converge=NA, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=0))	}


  if(sum(Idx.V.check)==1){

    if(is.null(SetID)){
      msg <- sprintf("The GxE matrix has fewer than 2 SNPs. " )
    } else {
      msg <- sprintf("In %s, the GxE matrix has fewer than 2 SNPs.", SetID )
    }

    warning(msg,call.=FALSE)
    return(list(pvalue=NA, Is_converge=NA, lambda=NA, n.G.test=drop(ncol(Z)), n.GE.test=1))
  }

  V <- as.matrix(V[,Idx.V.check])



  #-------------------------------------
  # Run GESAT
  Xtilde <- as.matrix(cbind(X,E))
  if(out_type == "C"){
    iSKAT.out <- GxEscore.linear.GCV(Y=Y, Xtilde=Xtilde, Z=Z, V=V, type=type,
                                     lower=lower, upper=upper, nintervals=nintervals, plotGCV=plotGCV, plotfile=plotfile, scale.Z=scale.Z, 				weights.Z=weights.Z, weights.V=NULL)
  }

  if(out_type == "D"){
    iSKAT.out <- GxEscore.logistic.GCV(Y=Y, Xtilde=Xtilde, Z=Z, V=V, type=type,
                                       lower=lower, upper=upper, nintervals=nintervals, plotGCV=plotGCV, plotfile=plotfile, scale.Z=scale.Z, 				weights.Z=weights.Z, weights.V=NULL)
  }

  return(list(pvalue=iSKAT.out$pvalue, Is_converge=iSKAT.out$Is_converge,
              lambda=iSKAT.out$lambda, n.G.test=drop(ncol(Z)), n.GE.test=drop(ncol(V))))
}














#
#	Check the Z, and do imputation
#     modified significantly from SKAT_MAIN_Check_Z from V0.78
#
iSKAT_MAIN_Check_Z <- function(Z, SetID, weights.Z=NULL, weights.V=NULL, impute.method,  missing_cutoff){

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
  m <- ncol(Z)
  ID_INCLUDE_SNP <- NULL
  for(i in 1:m){
    missing.ratio <- length(which(is.na(Z[,i])))/n
    sd1 <- sd(Z[,i], na.rm=TRUE)
    if(is.na(sd1)==FALSE){  # change in v1.0 of iSKAT package, Jan 29, 2014
      if(missing.ratio < missing_cutoff && sd1 > 0){
        ID_INCLUDE_SNP <- c(ID_INCLUDE_SNP,i)
      }
    }
  }

  if(length(ID_INCLUDE_SNP) == 0){

    if(is.null(SetID)){
      msg <- sprintf("ALL SNPs have either high missing rates or no-variation. ")
    } else {
      msg <- sprintf("In %s, ALL SNPs have either high missing rates or no-variation. ",SetID )
    }
    warning(msg, call.=FALSE)

    re <- list(Z.test=NULL, weights.Z.test=NULL, weights.V.test=NULL, check.Z.error=1)


  } else{

    if(m - length(ID_INCLUDE_SNP) > 0 ){

      if(is.null(SetID)){
        msg <- sprintf("%d SNPs with either high missing rates or no-variation are excluded!", m - length(ID_INCLUDE_SNP))
      } else {
        msg <- sprintf("In %s, %d SNPs with either high missing rates or no-variation are excluded!",SetID, m - length(ID_INCLUDE_SNP) )
      }

      warning(msg, call.=FALSE)
      Z <- as.matrix(Z[,ID_INCLUDE_SNP])
      if(is.null(weights.Z)==FALSE) weights.Z <- weights.Z[ID_INCLUDE_SNP]
      if(is.null(weights.V)==FALSE) weights.V <- weights.V[ID_INCLUDE_SNP]
      check.Z.error <- 2
      IDX_MISS <- which(is.na(Z))
    }



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
    re <- list(Z.test=Z, weights.Z.test=weights.Z, weights.V.test=weights.V, check.Z.error=check.Z.error)

  }

  return(re)
}




# copied without modifcation from SKAT V0.78, renamed Impute_iSKAT()
# Simple Imputation
# Z : an n x p genotype matrix with n samples and p SNPs
# Missing has to be NA: a missing genotype value.

Impute_iSKAT <-function(Z, impute.method){

  p <- dim(Z)[2]

  if(impute.method =="random"){
    for(i in 1:p){
      IDX <- which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1 <- mean(Z[-IDX,i])/2
        Z[IDX,i] <- rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1 <- mean(Z[-IDX,i])/2
        Z[IDX,i] <- 2 * maf1
      }
    }
  } else {
    stop("Error: Imputation method shoud be either \"fixed\" or \"random\"! ")
  }

  return(Z)
}
