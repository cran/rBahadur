#' Binary random variates with Diagonal Plus Low Rank (dplr) correlations
#'
#' Generate second Bahadur order multivariate Bernoulli random variates with 
#' Diagonal Plus Low Rank (dplr) correlation structures.
#'
#' @importFrom stats runif
#'
#' @param n number of observations
#' @param mu vector of means
#' @param U outer product component matrix
#'
#' @details This generates multivariate Bernoulli (MVB) random vectors with mean 
#' vector 'mu' and correlation matrix \eqn{C = D + U U^T} where \eqn{D} is a diagonal
#'  matrix with values dictated by 'U'. 'mu' must take values in the open unit interval 
#'  and 'U' must induce a valid second Bahadur order probability distribution. That is, 
#'  there must exist an MVB probability distribution with first moments 'mu' and 
#'  standardized central second moments \eqn{C} such that all higher order central 
#'  moments are zero.
#'
#' @return An \eqn{n}-by-\eqn{m} matrix of binary random variates, where \eqn{m} is 
#' the length of 'mu'.
#' @export
#'
#' @examples
#' set.seed(1)
#' h2_0 = .5; m = 200; n = 1000; r =.5; min_MAF=.1
#'
#' ## draw standardized diploid allele substitution effects
#' beta <- scale(rnorm(m))*sqrt(h2_0 / m)
#'
#' ## draw allele frequencies
#' AF <- runif(m, min_MAF, 1 - min_MAF)
#'
#' ## compute unstandardized effects
#' beta_unscaled <- beta/sqrt(2*AF*(1-AF))
#'
#' ## generate corresponding haploid quantities
#' beta_hap <- rep(beta, each=2)
#' AF_hap <- rep(AF, each=2)
#'
#' ## compute equilibrium outer product covariance component
#' U <- am_covariance_structure(beta, AF, r)
#'
#' ## draw multivariate Bernoulli haplotypes
#' H <- rb_dplr(n, AF_hap, U)
#'
#' ## convert to diploid genotypes
#' G <- H[,seq(1,ncol(H),2)] + H[,seq(2,ncol(H),2)]
#'
#' ## empirical allele frequencies vs target frequencies
#' emp_afs <- colMeans(G)/2
#' plot(AF, emp_afs)
#'
#' ## construct phenotype
#' heritable_y <-  G%*%beta_unscaled
#' nonheritable_y <-  rnorm(n, 0, sqrt(1-h2_0))
#' y <- heritable_y + nonheritable_y
#'
#' ## empirical h2 vs expected equilibrium h2
#' (emp_h2 <- var(heritable_y)/var(y))
#' h2_eq(r, h2_0)

rb_dplr <- function(n, mu, U) {

  M <- length(mu)

  k <- matrix(NaN, nrow=n, ncol=M)
  rand_U <- matrix(runif(M*n), nrow=n, ncol=M)

  # initial step
  p <- rep(mu[1],n)
  k[ ,1] <- as.numeric(rand_U[,1] <= p)

  tmp_bool <- (k[ ,1]==0)
  p <- tmp_bool*(1-p) + (!tmp_bool)*p
  Bk0 <- tmp_bool*(1-mu[1]) + (!tmp_bool)*mu[1]
  Bk1 <- tmp_bool*(-1) + (!tmp_bool)*1

  x <- Bk1*U[1]/p
  c <- 1

  # recursive steps
  for (m in 2:(M-1)) {
    p <- mu[m] + x * U[m]
    if (any(p < 0 | p > 1)) {
      stop('Infeasible probabilities. Are you sure specified parameters correspond to 
           a valid Bahadur order-2 MVB distribution?')
    }
    k[ ,m] <- (rand_U[ ,m] <= p)

    tmp_bool <- (k[ ,m]==0)
    p <- tmp_bool*(1-p) + (!tmp_bool)*p
    Bk0 <- tmp_bool*(1-mu[m]) + (!tmp_bool)*mu[m]
    Bk1 <- tmp_bool*(-1) + (!tmp_bool)*1

    x <- (x*Bk0 + c*Bk1*U[m])/p
    c <- (Bk0/p)*c
  }
  p <- mu[M] + x*U[M]
  if (any(p < 0 | p > 1)) {
    stop(mu[M],' ',p,'Infeasible probabilities. Are you sure specified parameters
         correspond to a valid Bahadur order-2 MVB distribution?')
  }
  k[ ,M] <-(rand_U[ ,M] <= p)

  return(k)
}


