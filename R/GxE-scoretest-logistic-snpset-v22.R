# Jan 29, 2014
# v22: change from "Davies" to "davies" and "Liu" to "liu"
# Jan 27, 2014
# v21 is modified from v20
# but v21 only keeps the functions required for logistic regression
# since most of the functions are shared between linear and logistic
# removes old functions, burden tests

#library(penalized)

#--------------------------------------------------------------------------------------------------------------
# START: Functions specific to logistic
#
#--------------------------------------------------------------------------------------------------------------
GxEscore.logistic.GCV <- function(Y, Xtilde, Z, V, type="davies",
        lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, plotGCV=F, plotfile=NA, scale.Z=T, weights.Z=NULL, weights.V=NULL){

        # Y (n x 1 matrix):  binary outcome variable
        # Xtilde (n x qtilde matrix): are the variables adjusted for (not penalized)
        # Z (n x p matrix):  genetic covariates that are adjusted for (penalized)
        # V (n x p matrix): GxE terms that we are testing
        # n = no. of samples
	    # lower cannot be zero
        # to use a fixed lambda, set nintervals=1
        # NB: if nintervals=1, upper is used and lower is ignored

	    # v20: add three new arguments scale.Z=T, weights.Z=NULL, weights.V=NULL
	    # v20: if nintervals=1, upper is used and lower is ignored
        # v19: gives Liu p-value if Davies doesn't converge
        # modified in v16 to introduce davies method
        # v16 added 3 new arguments  lower, upper and type, otherwise unchanged
        # v16 added some initial checks
        # old version (< v16) is similar to GxEscore.logistic(),
        # except the ridge parameter is selected by GCV
        # v16 also removed p-df LRT
        # v16 has davies method, so isn't similar to  GxEscore.logistic()

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

        #Q <- t(Y-Yhat) %*% V %*% t(V) %*% (Y-Yhat)
		Q1 <- t(Y-Yhat) %*% V  #v8c
        Q <- Q1 %*% t(Q1)      #v8c


        p <- ncol(Z)
        qtilde <- ncol(Xtilde) + 1
        # +1 is for intercept,
        # no need to multiply by 2 for lambda as the penalized package has a 0.5 in penalized likelihood
        if(p==1){
                temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda)))
        }else{
                temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda,p))))
        }

        #Wvarhat <- W * varhat_vec                         #v8c,varhat %*% W = W * varhat_vec
        #inverse <- solve(t(W) %*% Wvarhat + temp)         #v8c
        #inverse <- solve(t(W) %*% varhat %*% W + temp)
		# NB: sqrtvarhat_vec = sqrt(varhat_vec)
        M1 <- sqrtvarhat_vec * W                      # M1 = sqrt(varhat) %*% W
        M0 <- t(sqrtvarhat_vec * V)                   # M0 = V' %*% sqrt(varhat)
        #inverse <- solve(t(M1) %*% M1 +temp)         # t(M1) %*% M1 = t(W) %*% varhat %*% W
        #M2 <- M0 - M0 %*% M1 %*% inverse %*% t(M1)
		transM1 <- t(M1)
		invM1 <- solve(transM1 %*% M1 +temp, transM1) # invM1 = inverse %*% t(M1) = solve(t(M1) %*% M1 +temp) %*% t(M1)
		M2 <- M0 - M0 %*% M1 %*% invM1
        M3 <- M2 %*% t(M2)


#---------------------------------------------------------------------------
# p-value from non-central chi-square approximation
#---------------------------------------------------------------------------
        if(type=="liu"){
        #H <- W  %*% inverse %*% t(Wvarhat)       #v8c
        #H <- W  %*% inverse %*% t(W) %*% varhat

        #R1 <- t(diag(n)-H) %*% (V * varhat_vec)  #v8c
        #R <- R1 %*% t(R1)
        #v8c, this is correct, old version shdnt have t(diag(n)-H) but shd be (diag(n)-H)
        #R <- t(diag(n)-H) %*% varhat %*% V %*% t(V) %*% varhat %*% t(diag(n)-H)

        #A <- R %*% sigmahat            #v8c
        #A2 <- A %*% A                  #v8c
        #kappa1 <-  mtrace(A)           #v8c
        #kappa2 <-  2*mtrace(A2)        #v8c
        #kappa3 <-  8*sum(A * t(A2))    #v8c
        #kappa4 <-  48*sum(A2 * t(A2))  #v8c

        #kappa1 <-  mtrace(R %*% sigmahat)
        #kappa2 <-  2*mtrace(R %*% sigmahat %*% R %*% sigmahat)
        #kappa3 <-  8*mtrace(R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat)
        #kappa4 <-  48*mtrace(R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat %*% R %*% sigmahat)

		M4 <- M3 %*% M3
       	kappa1 <-  mtrace(M3)             # = tr(M3)
        kappa2 <-  2*mtrace(M4)           # = 2 tr(M3 M3) = 2 tr(M4)
        kappa3 <-  8*sum(M3 * t(M4))      # = 8 tr(M3 M3 M3) = 8 tr(M3 M4) = 8 sum(M3 * M4')
        kappa4 <-  48*sum(M4 * t(M4))     # = 48 tr(M3 M3 M3 M3) = 48 tr(M4 M4) = 48 sum(M4 * M4')

        approx <- noncentralapproxdirect(kappa2, kappa3, kappa4)
        Q.Norm <-((Q - kappa1)/sqrt(kappa2))*approx$sigmaX +  approx$muX
        pvalue <- pchisq(Q.Norm, df=approx$df, ncp = approx$ncp,  lower.tail=F)
        Is_converge <- 1

        }else{


#---------------------------------------------------------------------------
# p-value from davies
#---------------------------------------------------------------------------
        # sqrt(varhat) %*% W = W * sqrt(varhat_vec)
        # W' %*% sqrt(varhat) = t(W*sqrt(varhat_vec))
        # M1 <- t(sigmahat * sqrt(varhat_vec)) - W  %*% inverse %*% t(W * sqrt(varhat_vec))
        # M2 <- t(V *varhat_vec) %*% M1
        # M3 <- M2 %*% t(M2)

        daviesout <- Get_PValue_GESAT(M3, Q)
        pvalue <- daviesout$p.value
        Is_converge <- daviesout$is_converge

        # new in v19: do nchisq/Liu if Davies doesn't converge
        if(Is_converge<=0){
            M4 <- M3 %*% M3
        	kappa1 <-  mtrace(M3)             # = tr(M3)
        	kappa2 <-  2*mtrace(M4)           # = 2 tr(M3 M3) = 2 tr(M4)
        	kappa3 <-  8*sum(M3 * t(M4))      # = 8 tr(M3 M3 M3) = 8 tr(M3 M4) = 8 sum(M3 * M4')
        	kappa4 <-  48*sum(M4 * t(M4))     # = 48 tr(M3 M3 M3 M3) = 48 tr(M4 M4) = 48 sum(M4 * M4')
			approx <- noncentralapproxdirect(kappa2, kappa3, kappa4)
            Q.Norm <-((Q - kappa1)/sqrt(kappa2))*approx$sigmaX +  approx$muX
            pvalue <- pchisq(Q.Norm, df=approx$df, ncp = approx$ncp,  lower.tail=F)
        }

        }

        return(list(pvalue=pvalue, Is_converge=Is_converge, lambda=drop(model$lambda)))

}



ridge.logistic <- function(Y, Xtilde, Z, lambda=0, center.Z=T, scale.Z=T, weights.Z=NULL){
        # in v20, functions were renamed
        # ridge.logistic() and ridge.select.logistic() were modified from ridge.logistic.lambda() in v19
        # modified in v20 to add three more arguments center.Z=T, scale.Z=T, weights.Z=NULL
        # ridge.logistic() and ridge.select.logistic() are similar but return different things
        # ridge.select.logistic() returns only lambda, GCV, effective.df
        # ridge.logistic() returns only lambda, Yhat, thetahat
	# note that Z should always be centered as scale() behaves weirdly for scaling if not centered
 
        #============================================================
        # Method 2: Use penalized package - do not center Y
        #============================================================
        # when unpenalized is specified as a matrix, intercept has to be specified explicitly
        # ridge.logistic.lambda() includes intercept, so do not add intercept as part of Xtilde
        # uses penalized package (load it externally)
        # Xtilde is a n*qtilde matrix,
        # where the qtilde covariates do not have penalty imposed
        # Z is a n*p matrix where a penalty is imposed on the p covariates
        # Y is the n*1 matrix of binrary outcomes
        # returns a (qtilde+p)*1 matrix for thetahat, n*1 matrix for Yhat

        #if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
        #if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")

        n <- nrow(Y)
        if(is.null(weights.Z)==F & scale.Z==F){
                Z <- t(t(Z) * (weights.Z))
        }

        Z <- scale(Z, center=center.Z, scale=scale.Z)
        test2 <- penalized(response=Y, penalized=Z,  unpenalized=cbind(rep(1,n), scale(Xtilde)), lambda1=0, lambda2=lambda, model="logistic")
        if (test2@converged==F){ # new in v8b
                thetahat <- NA
                Yhat <- NA
        }else{
                thetahat <- coefficients(test2, "all")
                Yhat <- fitted(test2)
	}
	return(list(lambda=lambda, thetahat=thetahat, Yhat = Yhat))
}



ridge.select.logistic <- function(Y, Xtilde, Z, lambda=0, center.Z=T, scale.Z=T, weights.Z=NULL){
	# in v20, functions were renamed
	# ridge.logistic() and ridge.select.logistic() were modified from ridge.logistic.lambda() in v19	
	# modified in v20 to add three more arguments center.Z=T, scale.Z=T, weights.Z=NULL
	# ridge.logistic() and ridge.select.logistic() are similar but return different things
	# ridge.select.logistic() returns only lambda, GCV, effective.df
	# ridge.logistic() returns only lambda, Yhat, thetahat
        # note that Z should always be centered as scale() behaves weirdly for scaling if not centered

        #============================================================
        # Method 2: Use penalized package - do not center Y
        #============================================================
        # when unpenalized is specified as a matrix, intercept has to be specified explicitly
        # ridge.logistic.lambda() includes intercept, so do not add intercept as part of Xtilde
        # uses penalized package (load it externally)
        # Xtilde is a n*qtilde matrix,
        # where the qtilde covariates do not have penalty imposed
        # Z is a n*p matrix where a penalty is imposed on the p covariates
        # Y is the n*1 matrix of binrary outcomes
        # returns a (qtilde+p)*1 matrix for thetahat, n*1 matrix for Yhat

        #if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
        #if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
        
	n <- nrow(Y) 
	if(is.null(weights.Z)==F & scale.Z==F){
		Z <- t(t(Z) * (weights.Z))
	}

	Z <- scale(Z, center=center.Z, scale=scale.Z)
        test2 <- penalized(response=Y, penalized=Z,  unpenalized=cbind(rep(1,n), scale(Xtilde)), lambda1=0, lambda2=lambda, model="logistic")
        if (test2@converged==F){ # new in v8b
                thetahat <- NA
                Yhat <- NA
                GCV <- NA
                effective.df <- NA
        }else{
        	thetahat <- coefficients(test2, "all")
        	Yhat <- fitted(test2)
		sqrtvarhat_vec <- sqrt(Yhat*(1-Yhat))

		if(sum(sqrtvarhat_vec<=0)>0){
			GCV <- NA
                	effective.df <- NA
		}else{
        		W <- cbind(rep(1,n), scale(Xtilde), Z)         # all the variables in the null model
        		#varhat_vec <- Yhat*(1-Yhat)                           #v8c
        		#varhat <- diag(varhat_vec)                            #v8c
        		p <- ncol(Z)
        		qtilde <- ncol(Xtilde) + 1
        		# +1 is for intercept, no need to multiply by 2 for lambda as the penalized package has a 0.5 in penalized likelihood
        		if(p==1){
                		temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda)))
       		 	}else{
                		temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda,p))))
        		}

        		#Wvarhat <- W * varhat_vec
        		#v8c, varhat %*% W =  W * varhat_vec since varhat is diagonal
        		#inverse <- try(solve(t(W)  %*% Wvarhat + temp))          #v8c
        		#inverse <- try(solve(t(W) %*% varhat %*% W + temp))     # try fn new in v8B

			M1 <- sqrtvarhat_vec * W             	        # M1 = sqrt(varhat) %*% W
			MM1<- t(M1) %*% M1			# MM1 = M1' M1 = W' sqrt(varhat) sqrt(varhat) W = W' varhat W
			inverse <- try(solve(MM1+temp))

        		error.inverse <- 0                                       # new in v8b
        		if (is(inverse, "try-error")) error.inverse <- 1         # new in v8b
        		if(error.inverse==1){
                		GCV <- NA
                		effective.df <- NA
        		}else{

                	#H <- W  %*% inverse %*% t(Wvarhat)              	                #v8c
                	#H <- W  %*% inverse %*% t(W) %*% varhat
                	#effective.df <- mtrace(H)-1                      	                # subtract 1 for the intercept
			        equivH <- MM1 %*% inverse		 		        # tr(H) = tr(equivH)
			        effective.df <- mtrace(equivH)-1               		            # subtract 1 for the intercept
                	#GCV <- (t(Y-Yhat) %*% sigmahat %*% (Y-Yhat))/(n*(1-effective.df/n)^2)  #new GCV in v8
			GCV <- (sum(((Y-Yhat)/sqrtvarhat_vec)^2))/(n*(1-effective.df/n)^2)			
        		}

		}	
}
        return(list(lambda=lambda, GCV = GCV, effective.df=effective.df))
	    #return(list(lambda=lambda, thetahat=thetahat, Yhat = Yhat, GCV = GCV, effective.df=effective.df))
}


chooseridge.logistic <- function(Y, Xtilde, Z, lambdastart=1e-20, lambdaend=sqrt(nrow(Y))/log(nrow(Y)),
                         intervals=5, plot=F, file=NA, center.Z=T, scale.Z=T, weights.Z=NULL)
{
	    # renamed chooseridge.logistic in v20 (used to be chooseridge)
	    # modified in v20 to add center.Z=T, scale.Z=T, weights.Z=NULL as arguments
        # modified in v16 to add lambdastart and lambdaend as arguments, otherwise unchanged
        # need lambdastart, lambdaend >=0, lambdastart < lambdaend

        lambda <- c(exp(seq(log(lambdastart),log(lambdaend),length=intervals)))
        output <- c()
        for(ii in 1:length(lambda)){
                temp <- ridge.select.logistic(Y, Xtilde, Z, lambda[ii], center.Z, scale.Z, weights.Z)$GCV
                output <- c(output, temp)
                rm(temp)
        }

        lambdafinal <- lambda[which.min(output)]
          if(plot==T&is.na(file)==T){
                plot(lambda, output, xlab="lambda", ylab="GCV")
                abline(v=lambdafinal, col="red")
          }

        if(plot==T&is.na(file)==F){
                pdf(file)
                plot(lambda, output, xlab="lambda", ylab="GCV")
                abline(v=lambdafinal, col="red")
                dev.off()
        }

        return(lambdafinal)
}

#--------------------------------------------------------------------------------------------------------------
# END: Functions specific to logistic
#-------------------------------------------------------------------------------------------------------------






