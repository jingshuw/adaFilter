#' Compute partial conjunction p-values

#' @inheritParams AdaFilteringPC

#' @param p.vec a vector of p-values of length n
#' @param method the method to construct partial conjunction p-values. 
#'        There are currently 3 options using either Bonferroni or Fisher or 
#'        Simes derived partial conjunction p-value
#'
#' @return the BHPC p-value
#'
#' @export



ComputeBHPC <- function(p.vec,
						r,
						method = c("Bonferroni", "Fisher", "Simes")) {

	method <- match.arg(method, c("Bonferroni", "Fisher", "Simes"))

	n <- length(p.vec)

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
#' @inheritParams AdaFilteringPC
#' @inheritParams ComputeBHPC
#'
#' @return a list of objects
#' \describe{
#' \item{decision}{a length M vector of 1 and 0s indicating whether the corresponding PC hypothesis is rejected or not}
#' \item{BHPC.p}{the length M vector of partial conjunction p-values}
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

	bhpc.pvalues <- apply(p.matrix, 1, ComputeBHPC, r, method)

	if (type.I.err == "FWER")
		adjust.bhpc.p <- p.adjust(bhpc.pvalues, "bonferroni")
	else if (type.I.err == "PFER")
		adjust.bhpc.p <- bhpc.pvalues * length(bhpc.pvalues)
	else
		adjust.bhpc.p <- p.adjust(bhpc.pvalues, "BH")

	return(list(decision = adjust.bhpc.p <= alpha,
	            BHPC.p = bhpc.pvalues))

}
