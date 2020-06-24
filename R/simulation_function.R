#' Generate p-value matrix for simulation
#'
#' @param M number of partial conjunction hypotheses 
#' @param n number of studies
#' @param r the required minimum number of nonnull studies 
#'          for a nonnull partial conjunction hypothesis
#' @param all.zero.frac the fraction of complete nulls
#' @param alternative.frac the total fraction nonnull partial conjunction
#' hypotheses
#' @param mu a vector of possible values of effect size of nonnull individual 
#' hypotheses
#' @param rho correlation within block of screening tests
#' @param block.size the block size of the block dependency structure
#' @param one.sided either construct one-sided p-value or not
#'
#' @return a list of objects
#' \describe{
#' \item{truth.pc}{a length M vector of 1 and 0s indicating whether the 
#' corresponding hypothesis is a true nonnull partial conjunction hypothesis 
#' or not}
#' \item{pvalue.mat}{the individual p-value matrix of size M * n}
#' \item{zmat}{The corresponding z-value matrix of size M * n}
#' }
#' 
#' @references {
#' J. Wang, L. Gui, W. J. Su, C. Sabatti and A. B. Owen (2020). Detecting Multiple Replicating Signals using Adaptive Filtering Procedures
#' }
#'
#' @export


GenPMat <- function(M = 1000,
					n = 2,
					r = 2,
					all.zero.frac = 0.5,
					alternative.frac = 0.1,
					mu = c(3, 4, 6),
					rho = 0,
					block.size = 10,
                    one.sided = F) {
	combs <- as.matrix(expand.grid(data.frame(matrix(rep(c(0, 1), 
														 n), nrow = 2))))
	category <- rowSums(combs)
	n.cases <- sum(category >= r)
	weights <- rep((1 - alternative.frac - all.zero.frac)/(2^n - n.cases - 1),
				   nrow(combs))
	weights[category == 0] <- all.zero.frac
	weights[category >= r] <- rep(alternative.frac/n.cases, n.cases)

	mean.matrix <- combs[sample(1:nrow(combs), M, replace = T, prob = weights), ]



	truth.pc <- rowSums(mean.matrix) >= r


    ## This is the signal matrix
    mean.matrix <- apply(mean.matrix, 1, function(v) {

                             v[v == 1] <- sample(mu, sum(v == 1), replace = T)
                             return(v)
                   })
    mean.matrix <- t(mean.matrix)




    noise.matrix <- matrix(rnorm(n * M), ncol = n)


    ## generate dependency across screening tests
    if (rho != 0) {
        Sigma <- matrix(rep(rho, block.size * block.size), nrow = block.size)
        diag(Sigma) <- 1
        svd.Sigma <- svd(Sigma)
        rot <- t(t(svd.Sigma$u) * sqrt(svd.Sigma$d))

        for (i in 1:(M/block.size)) {
            noise.matrix[(i - 1) * block.size + 1:block.size, ] <- rot %*%
                noise.matrix[(i - 1) * block.size + 1:block.size, ]
        }
    }

	zmat <- noise.matrix + mean.matrix
	signs <- 2 * (zmat >= 0) - 1
	raw.pvalues <- pnorm(abs(zmat), lower.tail = F)

    if (one.sided) {
        pvalue.matrix <- signs * raw.pvalues + (1 - signs)/2
	} else {
        pvalue.matrix <- 2 * raw.pvalues 
    }

    return(list(truth.pc = truth.pc, 
				pvalue.mat = pvalue.matrix,
                zmat = zmat))
}



