# Version 4: June 18, 2015
# change one line of code so that it can handle small p-values better
# changed for Get_Liu_PVal.MOD.Lambda_iSKAT()
# so that instead of using 1- pchisq(), use pchisq(lower.tail=F)

# Version 3: Jan 29, 2014
# Version 3 removed all the functions repeated in GxE-scoretest-snpset-v11.R
# Remaining functions are identical to those in Version 2.
# Version 3 adds a function: checkvariation2()


# Version 2: Feb 2, 2013
# (1) Added three new functions (shared between iSKAT.logistic() and iSKAT.linear()):
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
# (2) Changes to four iSKAT.linear() and GxEscore.linear.GCV() specific functions
# Changes are made so that functions for linear regression are consistent between GxE-scoretest-snpset-v8.R, iSKAT-Linear-v2.R
# (i) Change to iSKAT.linear() so that varhat can accomodate larger p
# The following three functions are also changed so that they are identical to those in GxE-scoretest-snpset-v8.R
# Version 8: Jan 28, 2013
# (ii) change to GxEscore.linear.GCV() such that varhat can accomodate larger p
# Version 7: Dec 20, 2012
# (iii) change to ridge.select.linear() and ridge.linear()
# such that an intercept is now included but Y is not centered
# Version 6: Dec 17, 2012
# (iv) the only change is to add a try() function in ridge.select.linear()
# such that the GCV will only select a model that converges/ matrix is invertible


# Version 1: April 16, 2012
# iSKAT-Linear
# modified from GxE-scoretest-snpset-v5.R
# Functions that are retained from GxE-scoretest-snpset-v5.R are identical
# Thus any GESAT related functions are the same in both GxE-scoretest-snpset-v5.R and iSKAT-Linear-v1.R
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
# START: Functions specific to linear
#
#--------------------------------------------------------------------------------------------------------------
iSKAT.linear <- function(Y, Xtilde, Z, V, type="davies",
        lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=F, plotfile=NA, scale.Z=T, weights.Z=NULL, weights.V=NULL, r.corr=(0:10)/10){

        # Y (n x 1 matrix):  continuous outcome variable
        # Xtilde (n x qtilde matrix): are the variables adjusted for (not penalized)
        # Do NOT include intercept as part of Xtilde (Y is centered in ridge, so no intercept needed)
        # Z (n x p matrix):  genetic covariates that are adjusted for (penalized)
        # V (n x p matrix): GxE terms that we are testing
        # n = no. of samples
        # lower cannot be zero
        # to use a fixed lambda, set nintervals=1
        # NB: if nintervals=1, upper is used and lower is ignored

        # iSKAT.linear() Modified from GxEscore.linear.GCV() v5 of code
        # added 1 argument: r.corr(), type is changed to lower-case
	    # parts on fitting the null model are identical and unchanged

	    # iSKAT.linear() v2 changed varhat estimate to make it consistent with GxEscore.linear.GCV() v8

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
        if(nintervals>1){
        	if(is.null(weights.Z)==F & scale.Z==F){
                	Z <- t(t(Z) * (weights.Z))
        	}
        	Z <-  scale(Z, center=T, scale=scale.Z) 
        	lambdahat <- chooseridge.linear(Y, Xtilde, Z, lambdastart=lower, lambdaend=upper,
                	                        intervals=nintervals, plot=plotGCV, file=plotfile, center.Z=F, scale.Z=F, weights.Z=NULL)
                ridgemodel <- ridge.linear(Y, Xtilde, Z, lambda = lambdahat, center.Z=F, scale.Z=F, weights.Z=NULL)
                Yhat <- ridgemodel$Yhat
        }else{
                lambdahat <- upper
        	ridgemodel <- ridge.linear(Y, Xtilde, Z, lambda = lambdahat, center.Z=T, scale.Z=scale.Z, weights.Z=weights.Z)
        	Yhat <- ridgemodel$Yhat
	}
#---------------------------------------------------------------------------
# Score statistic
#---------------------------------------------------------------------------

        if(is.null(weights.V)==F){
                V <- t(t(V) * (weights.V))
        }


#---------------------------------------------------------------------------
# new in iSKAT (everything above this line is unchanged from GESAT v5, except for changing "Davies" to "davies" and "Liu" to "liu")
#---------------------------------------------------------------------------
	res <- Y - Yhat
        # Q <- t(Y-Yhat) %*% V %*% t(V) %*% (Y-Yhat)
        # Q1 <- t(res) %*% V
        # Q <- Q1 %*% t(Q1)

        #varhat <- drop(var(res))                                # Change in GxE-scoretest-v8.R, iSKAT-Linear-v2.R
	    df1 <- drop(sum(ridgemodel$W * t(ridgemodel$invW)))      # Change in GxE-scoretest-v8.R, iSKAT-Linear-v2.R
        varhat <- drop(var(res)) * (n-1) / (n - df1)             # Change in GxE-scoretest-v8.R, iSKAT-Linear-v2.R
	    V1 <- V - ridgemodel$W %*% ridgemodel$invW %*% V         # V1=t(diag(n)-H) %*% V, H = t(H) for linear regression only
	    results <- iSKAT_Optimal_Linear(res, V, X1=NULL, kernel=NULL, weights = NULL, s2=varhat, method=type, res.out=NULL,     
	    			n.Resampling =0, r.all=r.corr, V1)


        return(list(pvalue=results$p.value, param=results$param, lambda=drop(lambdahat)))
}



iSKAT_Optimal_Linear = function(res, V, X1=NULL, kernel=NULL, weights = NULL, s2, method=NULL
, res.out=NULL, n.Resampling =0, r.all, V1){

	# Modified from SKAT v0.73 SKAT_Optimal_Linear()
	# Note that arguments like X1, kernel, weights, res.out, n.Resampling are never used
	# they are kept for consistency with the corresponding SKAT function
	# all the Z's in original function are renamed V
	# added 1 argument V1

	# if r.all >=0.999 ,then r.all = 0.999. It is just for computation.
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
	Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) / s2

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
# END: Functions specific to linear
#
#--------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# START: common functions shared by both linear and logistic models
# NB: Functions common to both have same names in both scripts
# NB: Common functions are identical in both GxE-scoretest-logistic-snpset-v19.R and GxE-scoretest-snpset-v5.R
#-------------------------------------------------------------------------------------------------------
iSKAT_Optimal_Param<-function(Z1,r.all){

	# copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
	# Function get parameters of optimal test
	n<-dim(Z1)[1]
	p.m<-dim(Z1)[2]
	r.n<-length(r.all)

	z_mean<-rowMeans(Z1)
	Z_mean<-matrix(rep(z_mean,p.m),ncol=p.m,byrow=FALSE)
	cof1<-(t(z_mean) %*% Z1)[1,] / sum(z_mean^2)

	Z.item1<-Z_mean %*% diag(cof1)	
	Z.item2<-Z1 - Z.item1

	# W3.2 Term : mixture chisq
	W3.2.t<-t(Z.item2) %*% Z.item2
	lambda<-Get_Lambda(W3.2.t)
	
	# W3.3 Term : variance of remaining ...
	W3.3.item<-sum((t(Z.item1) %*% Z.item1) * (t(Z.item2) %*% Z.item2)) * 4
	
	# Mixture Parameters
	MuQ<-sum(lambda)
	VarQ<-sum(lambda^2) *2 + W3.3.item
	KerQ<-sum(lambda^4)/(sum(lambda^2))^2 * 12
	Df<-12/KerQ

	# W3.1 Term : tau1 * chisq_1
	tau<-rep(0,r.n)
	for(i in 1:r.n){
		r.corr<-r.all[i]
		#term1<-p.m*r.corr + cof1^2 * (1-r.corr)
		term1<-p.m^2*r.corr + sum(cof1^2) * (1-r.corr)
		tau[i]<-sum(term1) *  sum(z_mean^2)
	}

	out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau)
	return(out)
}



iSKAT_Optiaml_Each_Q<-function(param.m, Q.all, r.all, lambda.all){

	# copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
	#       Function get SKAT statistics with given rho
	#               Q.all is a matrix with n.q x n.r

	n.r<-length(r.all)
	c1<-rep(0,4)
	n.q<-dim(Q.all)[1]

	pval<-matrix(rep(0,n.r*n.q),ncol=n.r)
	pmin.q<-matrix(rep(0,n.r*n.q),ncol=n.r)
	param.mat<-NULL

	for(i in 1:n.r){
		Q<-Q.all[,i]
		r.corr<-r.all[i]
		lambda.temp<-lambda.all[[i]] 
		c1[1]<-sum(lambda.temp)
		c1[2]<-sum(lambda.temp^2)
		c1[3]<-sum(lambda.temp^3)
		c1[4]<-sum(lambda.temp^4)
		param.temp<-iSKAT_Get_Liu_Params_Mod(c1)

		muQ<-param.temp$muQ
		varQ<-param.temp$sigmaQ^2
		df<-param.temp$l

		# get pvalue
		Q.Norm<-(Q - muQ)/sqrt(varQ) * sqrt(2*df) + df
		pval[,i]<- 1-pchisq(Q.Norm,  df = df)

                #----------------------------
                # one-line addition in iSKAT v2
                pval[,i]<-Get_PValue.Lambda_iSKAT(lambda.temp,Q)$p.value
                #---------------------------- 
		param.mat<-rbind(param.mat,c(muQ,varQ,df))
	}

	pmin<-apply(pval,1,min)
	for(i in 1:n.r){
	
		muQ<-param.mat[i,1]
		varQ<-param.mat[i,2]
		df<-param.mat[i,3]

		q.org<-qchisq(1-pmin,df=df)
		q.q<-(q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ
		pmin.q[,i]<-q.q

	}
	
	out<-list(pmin=pmin,pval=pval,pmin.q=pmin.q)
	return(out)

}


iSKAT_Optimal_Integrate_Func_Davies<-function(x,pmin.q,param.m,r.all){
	 # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT

	n.r<-length(r.all)
	n.x<-length(x)

	temp1<-param.m$tau %x% t(x)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)

	re<-rep(0,length(x))
	for(i in 1:length(x)){
		#a1<<-temp.min[i]
		min1<-temp.min[i]
		if(min1 > sum(param.m$lambda) * 10^4){
			temp<-0
		} else {
			min1.temp<- min1 - param.m$MuQ			
			sd1<-sqrt(param.m$VarQ - param.m$VarRemain)/sqrt(param.m$VarQ)
			min1.st<-min1.temp *sd1 + param.m$MuQ
			
			dav.re<-SKAT_davies(min1.st,param.m$lambda,acc=10^(-6))
			temp<-dav.re$Qq
			if(dav.re$ifault != 0){
				stop("dav.re$ifault is not 0")
			}
		}
		if(temp > 1){
			temp=1
		}
		#lambda.record<<-param.m$lambda
		#print(c(min1,temp,dav.re$ifault,sum(param.m$lambda)))
		re[i]<-(1-temp) * dchisq(x[i],df=1)
	}
	return(re)

}


iSKAT_Optimal_PValue_Davies<-function(pmin.q,param.m,r.all){

 # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT
	re<-try(integrate(iSKAT_Optimal_Integrate_Func_Davies, lower=0, upper=30, subdivisions=500,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15), silent = TRUE)
	if(class(re) == "try-error"){
		re<-iSKAT_Optimal_PValue_Liu(pmin.q,param.m,r.all)
		return(re)
	} 

	pvalue<-1-re[[1]]
	if(pvalue < 0){
		pvalue=0
	}
	return(pvalue)

}


iSKAT_Optimal_Integrate_Func_Liu<-function(x,pmin.q,param.m,r.all){
 # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT	
	#x<-1
	#print(length(x))
	#print(x)
	#X1<<-x
	#x<-X1

	n.r<-length(r.all)
	n.x<-length(x)

	temp1<-param.m$tau %x% t(x)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)

	temp.q<-(temp.min - param.m$MuQ)/sqrt(param.m$VarQ)*sqrt(2*param.m$Df) + param.m$Df
	re<-pchisq(temp.q ,df=param.m$Df) * dchisq(x,df=1)
	return(re)

}


iSKAT_Optimal_PValue_Liu<-function(pmin.q,param.m,r.all){
 # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT

	 re<-integrate(iSKAT_Optimal_Integrate_Func_Liu, lower=0, upper=30, subdivisions=500
	,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15)
	
	pvalue<-1-re[[1]]
	return(pvalue)

}


iSKAT_Optimal_Get_Q<-function(Z1, res, r.all, n.Resampling, res.out, res.moments=NULL){
 # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT

	n.r<-length(r.all)
	p.m<-dim(Z1)[2]

	Q.r<-rep(0,n.r)
	Q.r.res<-NULL
	Q.sim<-NULL	
	
	temp<-t(res) %*% Z1
	for(i in 1:n.r){
		r.corr<-r.all[i]
		Q1<-(1-r.corr) * rowSums(temp^2)
		Q2<-r.corr * p.m^2 * rowMeans(temp)^2
		Q.r[i]<-Q1 + Q2
	}
	Q.r = Q.r /2
  	if(n.Resampling > 0){
	
		temp<-t(res.out) %*% Z1
		Q.r.res<-matrix(rep(0,n.Resampling *n.r),ncol=n.r)
		for(i in 1:n.r){
			r.corr<-r.all[i]
			Q1<-(1-r.corr) * rowSums(temp^2)
			Q2<-r.corr * p.m^2 * rowMeans(temp)^2
			Q.r.res[,i]<-Q1 + Q2
		}
		Q.r.res = Q.r.res/2
  	}

	if(!is.null(res.moments)){

		temp<-t(res.moments) %*% Z1
		n.moments<-dim(res.moments)[2]
		Q.sim<-matrix(rep(0,n.moments *n.r),ncol=n.r)
		for(i in 1:n.r){
			r.corr<-r.all[i]
			Q1<-(1-r.corr) * rowSums(temp^2)
			Q2<-r.corr * p.m^2 * rowMeans(temp)^2
			Q.sim[,i]<-Q1 + Q2
		}
		Q.sim = Q.sim/2

	}

	re<-list(Q.r=Q.r, Q.r.res=Q.r.res , Q.sim=Q.sim)
	return(re)


}


iSKAT_Optimal_Get_Pvalue<-function(Q.all, Z1, r.all, method){
 # copied without modification from SKAT 0.73, except to prefix of all names to iSKAT

	n.r<-length(r.all)
	n.q<-dim(Q.all)[1]
	p.m<-dim(Z1)[2]

	lambda.all<-list()
	for(i in 1:n.r){
		r.corr<-r.all[i]
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z2<- Z1 %*% t(L)
		K1<-t(Z2) %*% Z2

		lambda.all[[i]]<-Get_Lambda(K1)
		 
	}

	# Get Mixture param 
	param.m<-iSKAT_Optimal_Param(Z1,r.all)
	Each_Info<-iSKAT_Optiaml_Each_Q(param.m, Q.all, r.all, lambda.all)
	pmin.q<-Each_Info$pmin.q
	pval<-rep(0,n.q)

	if(method == "davies" || method=="optimal"){

		for(i in 1:n.q){
			pval[i]<-iSKAT_Optimal_PValue_Davies(pmin.q[i,],param.m,r.all)
		}


	} else if(method =="liu" || method =="liu.mod"){
		
		for(i in 1:n.q){
			pval[i]<-iSKAT_Optimal_PValue_Liu(pmin.q[i,],param.m,r.all)
		}

	} else {
		stop("Invalid Method!")
	}
	return(list(p.value=pval,p.val.each=Each_Info$pval))

}


iSKAT_Get_Liu_Params_Mod<-function(c1){
 # copied without modification from SKAT 0.73, except to prefix of name to iSKAT
 
 ## Helper function for getting the parameters for the null approximation
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}



Get_Liu_PVal.MOD.Lambda_iSKAT <- function(Q.all, lambda){
        # new in iSKAT v2, copied without modification from SKAT v0.81, renamed to append iSKAT
        param<-Get_Liu_Params_Mod_Lambda_iSKAT(lambda)

        Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
        Q.Norm1<-Q.Norm * param$sigmaX + param$muX
	# change in version 4 the line below
        p.value<- pchisq(Q.Norm1,  df = param$l, ncp=param$d, lower.tail=FALSE)

        return(p.value)

}

Get_Liu_Params_Mod_Lambda_iSKAT <- function(lambda){
  ## new in iSKAT v2, copied without modification from SKAT v0.81, renamed to append iSKAT
  ## Helper function for getting the parameters for the null approximation

  c1<-rep(0,4)
  for(i in 1:4){
        c1[i]<-sum(lambda^i)
  }

  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2

  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0

  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a

  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}


Get_PValue.Lambda_iSKAT<-function(lambda,Q){
        # new in iSKAT v2, copied with a one-line modification from SKAT v0.81, renamed to append iSKAT 
        #print(lambda)
        n1<-length(Q)

        p.val<-rep(0,n1)
        p.val.liu<-rep(0,n1)
        is_converge<-rep(0,n1)
        p.val.liu<-Get_Liu_PVal.MOD.Lambda_iSKAT(Q, lambda)

        for(i in 1:n1){
                out<-SKAT_davies(Q[i],lambda,acc=10^(-6))

                p.val[i]<-out$Qq
                #p.val.liu[i]<-SKAT_liu(Q[i],lambda)

                is_converge[i]<-1

                # check convergence
                if(length(lambda) == 1){
                        p.val[i]<-p.val.liu[i]
                } else if(out$ifault != 0){
                        is_converge[i]<-0
                        p.val[i]<-p.val.liu[i]  # this one-line modification is different from SKAT v0.81
                }

                # check p-value
                if(p.val[i] > 1 || p.val[i] <= 0 ){
                        is_converge[i]<-0
                        p.val[i]<-p.val.liu[i]
                }
        }

        return(list(p.value=p.val, p.val.liu=p.val.liu, is_converge=is_converge))

}



checkvariation <- function(x){
	# added in iSKAT v1
	# on top of checking that x is polymorphic, checks that there are at least 2 indiv in each level of x
        x <- x[is.na(x)==F]
	return(min(summary(as.factor(x)))>1&length(unique(x))>1)
}

checkvariation2 <- function(x){
	# added in iSKAT v3: to allow for dosage data, function returns TRUE only if there are more than 2 levels 
	# or if only 2 levels, there are no singletons
	# on top of checking that x is polymorphic, checks that there are at least 2 indiv in each level of x
        x <- x[is.na(x)==F]
	return((min(summary(as.factor(x)))>1&length(unique(x))==2)|length(unique(x))>2)
}

#-------------------------------------------------------------------------------------------------------
# END: common functions shared by both linear and logistic models
#
#-------------------------------------------------------------------------------------------------------
