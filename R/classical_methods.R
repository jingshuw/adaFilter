#' Compute partial conjunction p-values

#' @inheritParams adaFilter

#' @param p.vec a vector of p-values of length \code{n}
#' @param method the method to construct partial conjunction p-values. 
#'        There are currently 3 options:  "Bonferroni",  "Fisher" and  
#'        "Simes". The default is "Bonferroni".
#'
#' @return a vector of partial conjunction p-values 
#'
#' @export



computePC <- function(p.vec,
					  r,
					  method = c("Bonferroni", "Fisher", "Simes")) {

	method <- match.arg(method, c("Bonferroni", "Fisher", "Simes"))

	n <- sum(!is.na(p.vec))

	p.selected <- sort(p.vec, decreasing = T)[1:(n - r + 1)]

	if (method == "Bonferroni") {
		return((n-r+1) * p.selected[n-r+1])
	} else if (method == "Fisher") {
		df <- 2 * (n - r + 1)
		return(pchisq( -2*sum(log(p.selected)), df, lower.tail=FALSE))
	} else if (method == "Simes") {
		return(min((n -r + 1)/((n - r + 1) : 1) * p.selected))
	}
}

#' Classical procedures for partial conjunction hypotheses
#' @inheritParams adaFilter
#' @inheritParams computePC
#'
#' @return a list of objects
#' \describe{
#' \item{decision}{a length M vector of 1 and 0s indicating whether the corresponding PC hypothesis is rejected or not}
#' \item{PC.pvalues}{the length M vector of partial conjunction p-values}
#' }
#'
#' @export

ClassicalMTP <- function(p.matrix,
						 alpha,
						 r,
						 method = c("Bonferroni", "Fisher", "Simes"),
						 type.I.err = c("FWER", "FDR", "PFER")) {
  	method <- match.arg(method, c("Bonferroni", "Fisher", "Simes"))
	type.I.err <- match.arg(type.I.err, c("FWER", "FDR", "PFER"))
	
	n.noNA <- rowSums(!is.na(p.matrix))
	sh.idx <- which((n.noNA - r + 1) > 0)
	
	bhpc.pvalues <- apply(p.matrix[sh.idx,], 1, computePC, r, method)
	
	if (type.I.err == "FWER")
	  adjust.bhpc.p <- p.adjust(bhpc.pvalues, "bonferroni")
	else if (type.I.err == "PFER")
	  adjust.bhpc.p <- bhpc.pvalues * length(bhpc.pvalues)
	else
	  adjust.bhpc.p <- p.adjust(bhpc.pvalues, "BH")
	
	decision = rep(FALSE, nrow(p.matrix))
	decision[sh.idx] = adjust.bhpc.p <= alpha
	
	BHPC.p = rep(NA, nrow(p.matrix))
	BHPC.p[sh.idx] = bhpc.pvalues
	
	return(list(decision = decision,
	            PC.pvalues = BHPC.p))
	
}

