#' @title Estimation of expected maximum haploid breeding values (EMBV)
#'
#' @description This is (the R wrapper of) a C++ routine for the fast estimation of
#' expected maximum haploid breeding values.
#'
#' @name embv
#'
#' @param individual The genotype of the individual for which the EMBV should be computed.
#' The structure is a nested list, where the first level refers to the two gametes and the
#' second level to chromatids (integer vectors) within gametes.
#'
#' @param positions The genetic positions of loci. A list, where each element is a
#' (strictly) increasingly sorted double vector.
#'
#' @param effects The effects of loci. A list, where each element is a double vector.
#'
#' @param n_gam Integer. The number of gametes (possibly) contributed by the individual.
#'
#' @param se_level Double. Standard error threshold value.
#'
#' @param min_rep Integer. Minimum number of replications.
#'
#' @param max_rep Integer. Maximum number of replications.
#'
#' @param m Integer. Interference parameter (\code{m=0} means no interference).
#'
#' @param p Double. Proportion of chiasmata coming from non-interference process.
#'
#' @param seed Integer. Seed for conducting simulation.
#'
#' @details There will be at least \code{min_rep} and at maximum \code{max_rep} replicates.
#' If there are already more than \code{min_rep} replicates and the empirical standard
#' error falls below \code{se_level}, the iteration stops.
#'
#' @seealso \code{\link[Meiosis]{crossover}} for a description of the crossover process.
#' @return A list with the following elements:
#'   \itemize{
#'     \item embv: estimate of the EMBV
#'     \item se: standard error
#'     \item n: number of replicates
#'     \item time (ms): execution time in milliseconds
#'   }
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
#' embvr::embv(ind, positions, effects, n_gam = 20L, se_level = 1.0, min_rep = 10L,
#'             max_rep = 1000L, m = 0L, p = 0.0, seed = 456L)
NULL
