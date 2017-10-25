 \name{iSKAT.SSD.All}
 \alias{iSKAT.SSD.All}
 \title{interaction Sequence/SNP-set Kernel Association Test}
 \description{
	Iteratively test for interactions between a set of SNPS/genes and Environment for  SNP sets in SSD file. 
 }
 \usage{

	iSKAT.SSD.All(SSD.INFO, \dots)
 }
\arguments{
      \item{SSD.INFO}{an SSD_INFO object returned from Open_SSD.   }
      \item{\dots}{further arguments to be passed to ``iSKAT''. }
}
\value{
	\item{SetID}{SetID.}
	\item{pvalue}{the p-value of iSKAT.}	
	\item{rho_est}{the estimated \eqn{\rho} parameter.}  
	\item{lambda}{the chosen tuning parameter.}
	\item{n.G.test}{the number of SNPs in the genotype matrix used in the test. It can be different from ncol(Z) when some markers are monomorphic or have higher missing rates than the missing_cutoff. }  
	\item{n.GE.test}{the number of columns in the ZxE interaction matrix used in the test. It can be different from ncol(Z)*ncol(E) when some markers are monomorphic or have higher missing rates than the missing_cutoff. It can also be different from n.G.test*ncol(E) as the resulting ZxE variable might have no variability or might be perfectly collinear with columns of Z. Note that singletons are adjusted for in the main effects but are not tested for interactions. Likewise, columns of ZxE perfectly collinear with columns of Z are not tested for interactions. } 
  \item{Error}{1= Error, 0=No Error.}
  
}
\details{
Returns a data frame, where each row gives pvalue, rho_est, lambda, n.G.test, n.GE.test, Error for that particular SETID. 
Please see iSKAT for details.                     

}


\author{Xinyi (Cindy) Lin}

