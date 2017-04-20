#' @title Simulation of doubled haploid breeding values
#'
#' @description A function for the fast simulation of double haploid breeding values.
#'
#' @name dhbv
#'
#' @inheritParams embv
#'
#' @param n_rep Integer. Number of replications.
#'
#' @details This function simulates the breeding values of `n_rep` double haploid lines generated
#' by doubling of gametes meiotically derived from `individual`.
#'
#' @return A double vector.
#'
#' @examples
#' ## Simulate some data
#' set.seed(123L)
#' n_chr <- 10L
#' len <- rep(x = 300.0, times = n_chr)
#' n_loci <- rep(x = 1000L, times = n_chr)
#' alleles <- c(0L, 1L)
#' positions <- lapply(1L:n_chr, function(i) sort(runif(n_loci[i], min = 0.0, max = len[i])))
#' effects <- lapply(n_loci, rnorm)
#' ind <- replicate(n = 2L, lapply(n_loci, sample, x = alleles, replace = TRUE), simplify = FALSE)
#' n_rep <- 1e3L
#' bv <- embvr::dhbv(ind, positions, effects, n_rep, seed = 456L)
#' sd(bv)  ## estimate of square root of segregation variance
#'
#' ## Check the volatility of the estimate by simple bootstrapping.
#' sd(replicate(1e3L, sd(sample(bv, n_rep, TRUE))))
NULL
