## library('embvr')
## context('embv')

## test_that('Test passed', {

##   ref <- structure(list(embv = 90.3773111768218,
##                         se = 0.999016102897148,
##                         n = 321L,
##                         `time (ms)` = 111.264,
##                         seed = 456L),
##                    .Names = c("embv", "se", "n", "time (ms)", "seed"))

##   set.seed(123L)
##   n_chr <- 10L
##   len <- rep(x = 300.0, times = n_chr)
##   n_loci <- rep(x = 1000L, times = n_chr)
##   alleles <- c(0L, 1L)
##   positions <- lapply(1L:n_chr, function(i) sort(runif(n_loci[i], min = 0.0, max = len[i])))
##   effects <- lapply(n_loci, rnorm)
##   ind <- replicate(n = 2L, lapply(n_loci, sample, x = alleles, replace = TRUE), simplify = FALSE)
##   val <- embvr::embv(ind, positions, effects, n_gam = 20L, se_level = 1.0, min_rep = 10L,
##                     max_rep = 1000L, m = 0L, p = 0.0, seed = 456L)

##   expect_equal(ref$embv, val$embv, tolerance = 1e-3)
##   expect_equal(ref$se, val$se, tolerance = 1e-3)
##   expect_equal(ref$n, val$n)
## })
