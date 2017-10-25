# Jan 29, 2014
# Version 3 keeps only functions not repeated in iSKAT-linear-v3.R and GxE-scoretest-snpset-v12.R 
# remaining functions are identical to those in v2
# Version 2: Feb 1, 2013
# Added three new functions (shared between iSKAT.logistic() and iSKAT.linear()):
# Get_PValue.Lambda_iSKAT()
# Get_Liu_PVal.MOD.Lambda_iSKAT()
# Get_Liu_Params_Mod_Lambda_iSKAT()
# These three functions are only used for the optimal test iSKAT.logistic() and iSKAT.linear()
# Added one line to the function iSKAT_Optiaml_Each_Q()
# Details of the four changes in v2 (all changes apply to both logistic and linear):
# (i)   one-liner change in iSKAT_Optiaml_Each_Q()
# calls
# (ii)  Get_PValue.Lambda_iSKAT()
# which calls
# (iii) Get_Liu_PVal.MOD.Lambda_iSKAT()
# which calls
# (iv)  Get_Liu_Params_Mod_Lambda_iSKAT()
#
# Version 1: April 23, 2012
# iSKAT-Logistic
# modified from GxE-scoretest-logistic-snpset-v20.R and iSKAT-Linear-v1.R
# Functions that are retained from GxE-scoretest-logistic-snpset-v20.R and iSKAT-Linear-v1.R are identical
# Thus any GESAT related functions are the same in  GxE-scoretest-logistic-snpset-v20.R/iSKAT-Linear-v1.R and iSKAT-Logistic-v1.R
# iSKAT-Linear adds the function iSKAT.linear() and all other functions it needs
# All functions added to iSKAT-Linear-v1.R from GxE-scoretest-snpset-v5.R have names starting with iSKAT

#----------------------------------------------------------------------------------------------------------
# Note: these are the functions for linear (GxE-scoretest-snpset-v5.R, iSKAT-Linear-v1.R) 
# and logistic regression (GxE-scoretest-logistic-snpset-v20.R, iSKAT-Logistic-v1.R) respectively:
# GxEscore.logistic.GCV(), iSKAT.logistic()          : main function     
# ridge.logistic()                                   : fit final null model
# chooseridge.logistic() and ridge.select.logistic() : select ridge parameter
# Burdentests.GE.logistic                            : Burden tests

# GxEscore.linear.GCV(), iSKAT.linear()              : main function
# ridge.linear()                                     : fit final null model     
# chooseridge.linear() and ridge.select.linear()     : select ridge parameter
# Burdentests.GE.linear                              : Burden tests

# The analogous functions in each take in exactly the same arguments in linear and logistic codes
# ridge.logistic() returns slightly different values from ridge.linear()
# other analogous functions return the same values in both linear and logistic codes
#----------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# START: Functions specific to logistic
#
#--------------------------------------------------------------------------------------------------------------
iSKAT.logistic <- function(Y, Xtilde, Z, V, type="davies",
        lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=F, plotfile=NA, scale.Z=T, weights.Z=NULL, weights.V=NULL, r.corr=(0:10)/10){

        # Y (n x 1 matrix):  binary outcome variable
        # Xtilde (n x qtilde matrix): are the variables adjusted for (not penalized)
        # Z (n x p matrix):  genetic covariates that are adjusted for (penalized)
        # V (n x p matrix): GxE terms that we are testing
        # n = no. of samples
	    # lower cannot be zero
        # to use a fixed lambda, set nintervals=1
        # NB: if nintervals=1, upper is used and lower is ignored

        # iSKAT.logistic() Modified from GxEscore.logistic.GCV() v20 of code
        # added 1 argument: r.corr(), type is changed to lower-case
        # parts on fitting the null model are identical and unchanged

        if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
        if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
        if(nrow(V)!= nrow(Y)) stop("dimensions of V and Y don't match")
        if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
        if(lower<=0|upper<=0|lower>upper) stop("lower/upper has to be >0, lower<=upper")
        if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")

        n <- drop(nrow(Y))

#---------------------------------------------------------------------------
# fit ridge regression model under the null
# Note that all variables should always be centered as otherwise scaling by scale() will be incorrect
#---------------------------------------------------------------------------
        if(is.null(weights.Z)==F & scale.Z==F){
                Z <- t(t(Z) * (weights.Z))
        }
        Z <- scale(Z, center=T, scale=scale.Z)
        W <- cbind(rep(1,n), scale(Xtilde), Z)         # all the variables in the null model

        if(nintervals>1){
        	lambdahat <- chooseridge.logistic(Y, Xtilde, Z, lambdastart=lower, lambdaend=upper,
                	                        intervals=nintervals, plot=plotGCV, file=plotfile, center.Z=F, scale.Z=F, weights.Z=NULL)
        }else{
                lambdahat <- upper
	}

        model <- ridge.logistic(Y, Xtilde, Z, lambda = lambdahat, center.Z=F, scale.Z=F, weights.Z=NULL)
        Yhat <- model$Yhat
        sqrtvarhat_vec <- sqrt(Yhat*(1-Yhat))
        lambda <- model$lambda
        #varhat_vec <- Yhat*(1-Yhat)
        #varhat <- diag(varhat_vec)
        #sigmahat <- solve(varhat)


#---------------------------------------------------------------------------
# Score statistic
#---------------------------------------------------------------------------
        if(is.null(weights.V)==F){
                V <- t(t(V) * (weights.V))
        }


#---------------------------------------------------------------------------
# new in iSKAT (everything above this line is unchanged from GESAT v20, except for changing "Davies" to "davies" and "Liu" to "liu")
#---------------------------------------------------------------------------
		res <- Y - Yhat

        #Q <- t(Y-Yhat) %*% V %*% t(V) %*% (Y-Yhat)
	    #Q1 <- t(Y-Yhat) %*% V  #v8c
        #Q <- Q1 %*% t(Q1)      #v8c


        p <- ncol(Z)
        qtilde <- ncol(Xtilde) + 1
        # +1 is for intercept,
        # no need to multiply by 2 for lambda as the penalized package has a 0.5 in penalized likelihood
        if(p==1){
                temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda)))
        }else{
                temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda,p))))
        }

        M1 <- sqrtvarhat_vec * W                      # M1 = sqrt(varhat) %*% W
     	M00 <-sqrtvarhat_vec * V		      # M00 = sqrt(varhat) %*% V
	    transM1 <- t(M1)
	    invM1 <- solve(transM1 %*% M1 +temp, transM1) # invM1 = inverse %*% t(M1) = solve(t(M1) %*% M1 +temp) %*% t(M1)
	    V1 <- M00 -  M1 %*% invM1 %*% M00

	results <- iSKAT_Optimal_Logistic(res, V, X1=NULL, kernel=NULL, weights = NULL, pi_1=NULL, method=type, 
										res.out=NULL, n.Resampling =0, r.all=r.corr, V1)


        return(list(pvalue=results$p.value, param=results$param, lambda=drop(lambdahat)))



}







iSKAT_Optimal_Logistic  = function(res, V, X1=NULL, kernel=NULL, weights = NULL, pi_1=NULL , method = NULL
, res.out=NULL, n.Resampling =0, r.all, V1){


	# Modified from SKAT v0.73 SKAT_Optimal_Logistic()
	# Note that arguments like X1, kernel, weights, pi_1, res.out, n.Resampling are never used
	# they are kept for consistency with the corresponding SKAT function
	# all the Z's in original function are renamed V
	# added 1 argument V1

	# if r.all >=0.999 ,then r.all = 0.999
	IDX<-which(r.all >= 0.999)
	if(length(IDX) > 0){
		r.all[IDX]<-0.999	
	}

	n<-dim(V)[1]
	p.m<-dim(V)[2]
	n.r<-length(r.all)
	

	###########################################
	# Compute Q.r and Q.r.res
	##########################################
	out.Q<-iSKAT_Optimal_Get_Q(V, res, r.all, n.Resampling, res.out)
	Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) 

	##################################################
	# Compute P-values 
	#################################################

	out<-iSKAT_Optimal_Get_Pvalue(Q.all, V1 / sqrt(2), r.all, method)

	param<-list(p.val.each=NULL,q.val.each=NULL)
	param$p.val.each<-out$p.val.each[1,]
	param$q.val.each<-Q.all[1,]
	param$rho<-r.all
	param$minp<-min(param$p.val.each)

	id_temp<-which(param$p.val.each == min(param$p.val.each))
	id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
	if(length(id_temp1) > 0){
		param$rho[id_temp1] = 1
	}

	param$rho_est<-param$rho[id_temp]


	p.value<-out$p.value[1]
	p.value.resampling<-NULL
	if(n.Resampling > 0){
		p.value.resampling<-out$p.value[-1]
		#param$pval.each.resample<-out$p.val.each[-1]
	}

 	re<-list(p.value = p.value, p.value.resampling = p.value.resampling
	, Test.Type = method, Q = NA, param=param )  
  	
	return(re)	

}
#--------------------------------------------------------------------------------------------------------------
# END: Functions specific to logistic
#
#--------------------------------------------------------------------------------------------------------------

