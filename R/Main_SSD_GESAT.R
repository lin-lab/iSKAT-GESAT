# Jan 29, 2014
# in v1.0 of iSKAT package there are two sets of SSD functions, one for GESAT(), one for iSKAT()

# modified from SKAT V0.78 SKAT.SSD.OneSet()
GESAT.SSD.OneSet <- function(SSD.INFO, SetID, ...){
	
	id1 <- which(SSD.INFO$SetInfo$SetID == SetID)
	if(length(id1) == 0){
		MSG <- sprintf("Error: cannot find set id [%s] from SSD!", SetID)
		stop(MSG)
	}	
	Set_Index <- SSD.INFO$SetInfo$SetIndex[id1]

	Z <- Get_Genotypes_SSD(SSD.INFO, Set_Index)
	re <- GESAT(Z, ...)
	
	return(re)
}


# modified from SKAT V0.78 SKAT.SSD.OneSet_SetIndex()
GESAT.SSD.OneSet_SetIndex <- function(SSD.INFO, SetIndex, ...){
	
	id1 <- which(SSD.INFO$SetInfo$SetIndex == SetIndex)
	if(length(id1) == 0){
		MSG <- sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
		stop(MSG)
	}	
	SetID <- SSD.INFO$SetInfo$SetID[id1]


	Z <- Get_Genotypes_SSD(SSD.INFO, SetIndex)
	re <- GESAT(Z, ...)
	return(re)
}


# modified from SKAT V0.78 SKAT.SSD.All()
# how GESAT.SSD.ALL handles errors differently from SKAT.SSD.All()
GESAT.SSD.All <- function(SSD.INFO, ...){
	
	N.Set <- SSD.INFO$nSets
	OUT.Pvalue <- rep(NA, N.Set)
	OUT.n.G.test <- rep(NA, N.Set)
	OUT.n.GE.test <- rep(NA, N.Set)
	OUT.Error <- rep(-1, N.Set)
	OUT.lambda <- rep(NA, N.Set)
	OUT.converge <- rep(NA, N.Set)


	for(i in 1:N.Set){

		Is.Error <- TRUE
		try1 <- try(Get_Genotypes_SSD(SSD.INFO, i),silent = TRUE)

		if(class(try1) != "try-error"){
			Z <- try1
			Is.Error <- FALSE
		} else {
			err.msg <- geterrmessage()
			msg <- sprintf("Error to get genotypes of %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
			warning(msg,call.=FALSE)
		}
	

		if(!Is.Error){
			Is.Error <- TRUE
			try2 <- try(GESAT(Z, ...), silent = TRUE)
			
			if(class(try2) != "try-error"){
				re <- try2
				Is.Error <- FALSE
			} else {

				err.msg <- geterrmessage()
				msg <- sprintf("Error to run GESAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
				warning(msg,call.=FALSE)
			}
		}
		

		if(!Is.Error){

			OUT.Pvalue[i] <- re$pvalue
			OUT.converge[i] <- re$Is_converge
			OUT.lambda[i] <- re$lambda
			OUT.n.G.test[i] <- re$n.G.test
			OUT.n.GE.test[i] <- re$n.GE.test
			OUT.Error[i] <- 0

		}else{

			OUT.Pvalue[i] <- NA
			OUT.converge[i] <- NA
			OUT.lambda[i] <- NA
			OUT.n.G.test[i] <- NA
			OUT.n.GE.test[i] <- NA
			OUT.Error[i] <- 1

		}
	}

	
	out.tbl <- data.frame(SetID=SSD.INFO$SetInfo$SetID, pvalue=OUT.Pvalue, Is_converge=OUT.converge, lambda=OUT.lambda, 
					n.G.test=OUT.n.G.test, n.GE.test=OUT.n.GE.test, Error=OUT.Error)



	return(out.tbl)	
}

