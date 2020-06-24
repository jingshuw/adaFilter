#' Screening for partial conjunction hypotheses controlling
#' using the adaptive filtering procedures
#'
#' @description Provides the adaptive filtering multiple testing procedures
#' for partial conjunction (PC) hypotheses. The procedures can more efficiently reject non-null
#' partial conjunction hypotheses than previously methods
#' especially when the number of
#' candidate hypotheses is large while most of them are complete nulls
#'
#' @docType package
#' @name adaFilter_package
NULL

#' The main function for adaptive filtering procedures
#'
#' @param p.matrix an \code{M * n} matrix of p-values. \code{M} is the
#' number of PC hypotheses and \code{n} is the number of studies. Missing values are allowed.
#' @param r the required replicability level, i.e. the smallest number of nonnull individual hypotheses
#'          for a nonnull PC hypotheses
#' @param type.I.err one of "FDR", "FWER" or "PFER", with the default type I error "FDR". For FWER and PFER control, the adaFilter Bonferroni procedure is performed;
#' for FDR control, the adaFilter BH procedure is performed.
#' @param alpha control level of the chosen simultaneous error. Default is \code{0.05}.
#' @param fast logical value, default is FALSE. If \code{fast} is TRUE, then only results for hypotheses whose filter p-values are less than \code{alpha} are computed, and for other hypotheses the results are NA. This can speed up computation when \code{M} is super large.
#'
#' @return a data frame of size \code{M * 5} where the \code{5} columns are:
#' \describe{
#' \item{decision}{a length M vector of 1 and 0s indicating whether the corresponding PC hypothesis is rejected or not. For a hypothesis, its \code{decision} element is 1 if its \code{adjusted.p} element is not greater than \code{alpha}.}
#' \item{adjusted.p}{the adaFilter adjusted p-values for the selection p-values. If \code{type.I.error} is "FDR", then return the adaFilter BH adjusted p-values, otherwise, return the adaFilter Bonferroni adjusted p-values.}
#' \item{selection.p}{the length M vector of selection p-values}
#' \item{filter.p}{the length M vector of filtering p-values}
#' \item{adj.number}{The adaFilter adjustment number. For more details, see reference.}
#' }
#'
#' @references {
#' J. Wang, L. Gui, W. J. Su, C. Sabatti and A. B. Owen (2020). Detecting Multiple Replicating 
#' Signals using Adaptive Filtering Procedures
#' }
#'
#' @examples
#' ## Simulate a p-value matrix of size 1000 * 2 with 1000 partial conjunction hypotheses and 2 studies
#' data <- GenPMat()
#'
#' ## set r = 2, which means that the partial conjunction hypothesis is nonnull 
#' ## only when in both studies the corresponding hypothesis is nonnull
#' ## Control for FDR at level 0.1
#' result <- adaFilter(data$pvalue.mat, 2, "FDR", 0.1)
#' ## print the true positives
#' print(which(result$decision == 1 & data$truth.pc == 1))
#' ## plot the histogram of the selection p-values (which should be conservative)
#' hist(result$selection.p, breaks = 30,
#' main = "Histogram of the selection p-values")
#'
#' @import MASS
#' @export

adaFilter <- function(p.matrix,
					  r,
					  type.I.err = c("FDR", "FWER", "PFER"),
					  alpha = 0.05,
					  fast = F) {

	type.I.err <- match.arg(type.I.err, c("FWER", "FDR", "PFER"))

	n <- ncol(p.matrix)

	n.noNA <- rowSums(!is.na(p.matrix))
	if (fast) {
		sh.idx <- which((n.noNA >= r) & 
						rowSums(p.matrix <= alpha/(n.noNA - r + 1),na.rm = T) >= (r - 1))
	} else
		sh.idx <- which(n.noNA >= r)

	if(length(sh.idx) == 0){
		if (fast) {
			warning("All p-values are too large! No hypotheses can be rejected. Return NA")
		} else
			warning(paste("No hypotheses have at least", r, "non-missing values. return NA"))
		return(NA)
	} else{
		p.mat.sh <- p.matrix[sh.idx, ]

		M <- sum(n.noNA >= r)
		N <- length(sh.idx)
		n.noNA <- n.noNA[sh.idx]


		p.mat.sh <- t(apply(p.mat.sh, 1, sort, method = "quick",
							na.last = T))

		selection.p <- p.mat.sh[, r] * (n.noNA - r + 1)
		filter.p <- p.mat.sh[, r - 1] * (n.noNA - r + 1) 
		filter.p[(n.noNA < r) & (n.noNA >= r-1)] <- NA


		#### compute AdaFilter adjusted number
		sorted.s <- sort(selection.p, method = "quick", index.return = T)
		names(sorted.s$x) <- paste0("S", 1:N)
		names(filter.p) <- paste0("F", 1:N)
		temp <- sort(c(filter.p, sorted.s$x), method = "quick")
		temp[1:(2 * N)] <- 1:(2* N)
		adj.number <- temp[names(sorted.s$x)] - 1:N
		


		if (type.I.err != "FDR") {
			adjusted.p <- adj.number * sorted.s$x
			adjusted.p <- cummin(adjusted.p[N:1])[N:1]
			if (type.I.err == "FWER")
				adjusted.p <- pmin(1, adjusted.p)
	
		} else {
			adjusted.p <- adj.number / (1:N) * sorted.s$x
			adjusted.p <- pmin(1, cummin(adjusted.p[N:1])[N:1])
		}

		## reorder adj.number and adjusted.p
		adj.number[sorted.s$ix] <- adj.number
		adjusted.p[sorted.s$ix] <- adjusted.p
		
		decision <- adjusted.p <= alpha


		results <- data.frame(matrix(NA, nrow(p.matrix), 5))
		colnames(results) <- c("decision", "adjusted.p", "selection.p",
							   "filter.p", "adj.number")
		results[sh.idx, ] <- cbind(decision, adjusted.p, selection.p, filter.p, adj.number)

		rownames(results) <- rownames(p.matrix)


	return(results)
	}

}



