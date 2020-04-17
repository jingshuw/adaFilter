#' Screening for partial conjunction hypotheses controlling
#' using the adaptive filtering procedures
#'
#' @description Provides the adaptive filtering multiple testing procedures
#' for partial conjunction hypotheses. The procedures can more efficiently reject non-null
#' partial conjunction hypotheses than previously methods
#' especially when the number of
#' candidate hypotheses is large while most of them are complete nulls
#'
#' @docType package
#' @name adaFilter_package
NULL

#' The main function for adaptive filtering procedures
#'
#' @param p.matrix an M*n matrix of p-values. M is the
#' number of PC hypotheses and n is the number of studies
#' @param alpha control level of the chosen simultaneous error
#' @param r the smallest number of nonnull individual hypotheses
#'          for a nonnull partial conjunction hypotheses
#' @param type.I.err for FWER and PFER control,
#' the adaptive filtering Bonferroni procedure is used;
#' for FDR control, the adaptive filtering BH procedure is used
#'
#' @return a list of objects
#' \describe{
#' \item{decision}{a length M vector of 1 and 0s indicating whether the corresponding PC hypothesis is rejected or not}
#' \item{select.p}{the length M vector of selection p-values}
#' \item{filter.p}{the length M vector of filtering p-values}
#' }
#'
#' @references {
#' J. Wang, C. Sabatti and A. B. Owen (2016). Adaptive Filtering
#' Multiple Testing Procedures for Partial Conjunction
#' Hypotheses.
#' }
#'
#' @examples
#' ## Simulate a p-value matrix of size 1000 * 2 with 1000 partial conjunction hypotheses and 2 studies
#' data <- GenPMat()
#'
#' ## set r = 2, which means that the partial conjunction hypothesis is nonnull 
#' ## only when in both studies the corresponding hypothesis is nonnull
#' ## Control for FDR at level 0.1
#' result <- AdaFilteringPC(data$pvalue.mat, 0.1, 2, "FDR")
#' ## print the true positives
#' print(which(result$decision == 1 & data$truth.pc == 1))
#' ## plot the histogram of the selection p-values (which should be conservative)
#' hist(result$select.p, breaks = 30,
#' main = "Histogram of the Selection P-values")
#'
#' @import MASS
#' @export

AdaFilteringPC <- function(p.matrix,
						   alpha,
						   r,
						   type.I.err = c("FWER", "FDR", "PFER")) {

	type.I.err <- match.arg(type.I.err, c("FWER", "FDR", "PFER"))

	n <- ncol(p.matrix)
	N <- nrow(p.matrix)

  n.noNA <- rowSums(!is.na(p.matrix))
  sh.idx <- which((n.noNA - r + 1) > 0 & rowSums(p.matrix <= alpha/(n.noNA - r + 1),na.rm = T) >= (r - 1))
  
  if(length(sh.idx) == 0){
    decision = rep(FALSE, N)
    combined.p.final = matrix(NA, nrow = 2, ncol = N)
  } else{

  p.mat.sh <- p.matrix[sh.idx, ]


  N.ori <- N
  N <- length(sh.idx)
  n.noNA <- n.noNA[sh.idx]


  p.mat.sh <- t(apply(p.mat.sh, 1, sort, method = "quick",
                      na.last = T))
  combined.p <- t(p.mat.sh[, c(r, r-1)] * (n.noNA - r + 1))



	sorted.f <- sort(combined.p[2, ], method = "quick")


	if (type.I.err != "FDR") {

		m.prime <- FindEqualIdx(sorted.f, alpha/(1:N), N)

		m <- ifelse(sorted.f[m.prime] <= alpha/(m.prime - 1),
					m.prime, m.prime - 1)
		
		decision <- rep(FALSE,N.ori)
		decision[sh.idx] <- combined.p[1, ] <= alpha/m
		M <- m
	} else {
		sorted.s <- sort(combined.p[1, ], method = "quick")
		combined.mat <- matrix(0, 2 * N, 3)
		i <- 1
		j <- 1
		k <- 1
		while(i <= N && j <= N) {
			if (sorted.f[i] < sorted.s[j]) {
				combined.mat[k, 1] <- sorted.f[i]
				combined.mat[k, 2] <- i
				combined.mat[k, 3] <- j - 1
				i <- i + 1
			} else if (sorted.f[i] == sorted.s[j]) {
				combined.mat[k, 1] <- sorted.f[i]
				combined.mat[k, 2] <- i
				combined.mat[k, 3] <- j
				i <- i + 1
				j <- j + 1
			} else {
				combined.mat[k, 1] <- sorted.s[j]
				combined.mat[k, 2] <- i - 1
				combined.mat[k, 3] <- j
				j <- j + 1
			}
			k <- k + 1
		}

		if (i== N + 1 && j <= N) {
			combined.mat[k:(k + N - j), 1] <- sorted.s[j:N]
			combined.mat[k:(k + N - j), 2] <- rep(N, N - j + 1)
			combined.mat[k:(k + N - j), 3] <- j:N
			k <- k + N -j + 1
		} else if (j == N + 1 && i <= N) {
			combined.mat[k:(k + N - i), 1] <- sorted.f[i:N]
			combined.mat[k:(k + N - i), 3] <- rep(N, N - i + 1)
			combined.mat[k:(k + N - i), 2] <- i:N
			k <- k + N -i + 1
		}
		combined.mat <- combined.mat[1:(k - 1), ]

		thres <- combined.mat[, 3]/combined.mat[, 2] * alpha
		lower.judge <- combined.mat[, 1] < thres
		combined <- cbind(c(combined.mat[-1 , 1], 1), thres)
		upper <- t(apply(combined, 1, function(v) c(min(v), which.min(v))))

		candidates <- cbind(combined.mat[lower.judge, 1],
							upper[lower.judge, , drop = F])
		if (nrow(candidates) == 0)
			gamma0 <- 0
		else
			gamma0 <- mean(candidates[nrow(candidates), 1:2])
    
		decision <- rep(FALSE,N.ori)
		decision[sh.idx] <- combined.p[1, ] <= gamma0

    M <- sum(combined.p[2, ] <= gamma0)
    
    combined.p.final = matrix(NA, nrow = 2, ncol = N.ori)
    combined.p.final[,sh.idx] = combined.p
	  }
  }
  
	return(list(decision = decision,
				select.p = combined.p[1, ],
				filter.p = combined.p[2, ]))
}

#' v1 is non-decreasing, v2 is non-increasing, find the intersection point
#' @keywords internal
FindEqualIdx <- function(v1, v2, end, start = 1) {
	if (end - start < 2)
		return(end)
	midpoint <- ceiling((end + start) / 2)
	if (v1[midpoint] > v2[midpoint]) {
		return(FindEqualIdx(v1, v2, midpoint, start))
	} else if (v1[midpoint] < v2[midpoint]) {
		return(FindEqualIdx(v1, v2, end, midpoint))
	} else
		return(midpoint)
}



