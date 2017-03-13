// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Meiosis)]]
#include <Rcpp.h>
#include <random>
#include <vector>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <chrono>
#include <Meiosis.h>


//' @rdname bcgv
// [[Rcpp::export]]
Rcpp::List bcgv(std::vector<std::vector<std::vector<double>>> individual,
                std::vector<std::vector<double>> positions,
                std::vector<std::vector<double>> locus_effects,
                const unsigned int n_gam,
                const double se_level,
                const unsigned int min_rep,
                const unsigned int max_rep,
                const unsigned int m,
                const double p,
                const int seed
                )
{
  auto t1 = std::chrono::high_resolution_clock::now();

  if (min_rep > max_rep) Rcpp::stop("'min_rep' must be smaller than or equal to 'max_rep'");
  if (p < 0.0 || p > 1.0) Rcpp::stop("'p' must be in [0, 1].");

  engine.seed(seed);
  const auto& xof = Meiosis::crossover<std::vector<double>>;

  // Prepare xodat parameters.
  const bool obligate_chiasma = true;
  const std::size_t n_chr = positions.size();
  std::vector<double> L(n_chr), Lstar(n_chr);
  const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
  for (std::size_t i{0}; i != n_chr; ++i) {
    const auto& pos = positions[i];
    L[i] = std::ceil(pos.back());
    Lstar[i] = Meiosis::calc_Lstar(L[i], m, p, epsilon);
  }

  // Transform genotypes by multiplication with locus effects and compute sums.
  std::vector<std::vector<double>> sums(individual.size(), std::vector<double>(n_chr));
  for (std::size_t g{0}; g != individual.size(); ++g) {
    auto& gamete = individual[g];
    for (std::size_t i{0}; i != n_chr; ++i) {
      auto& chromatid = gamete[i];
      auto& leff = locus_effects[i];
      std::transform(chromatid.begin(), chromatid.end(),
                     leff.begin(), chromatid.begin(), std::multiplies<double>());
      sums[g][i] = std::accumulate(chromatid.begin(), chromatid.end(), 0.0);
    }
  }

  // Online statistics.
  std::size_t n = 0;
  double mean = 0.0;
  double M2 = 0.0;
  double se = std::numeric_limits<double>::max();
  double delta;
  double delta2;

  while (n < min_rep || (se >= se_level && n < max_rep)) {
    double max_value = std::numeric_limits<double>::lowest();

    for (std::size_t g{0}; g != n_gam; ++g) {
      double value = 0.0;
      for (std::size_t i{0}; i != n_chr; ++i) {
        const auto& xlocations = xof(L[i], m, p, obligate_chiasma, Lstar[i], engine2);
        value += meisum(individual[0][i], individual[1][i], xlocations,
                        positions[i], sums[0][i], sums[1][i], engine2);
      }
      if (value > max_value) max_value = value;
    }

    // Update online statistics.
    ++n;
    delta = max_value - mean;
    mean += delta / n;
    delta2 = max_value - mean;
    M2 += delta * delta2;
    se = std::sqrt(M2 / (n * (n - 1)));
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  return Rcpp::List::create(Rcpp::Named("bcgv") = mean,
                            Rcpp::Named("se") = se,
                            Rcpp::Named("n") = n,
                            Rcpp::Named("time (ms)") = duration / 1000.0);
}


double meisum(const t_alleles& patalle,
              const t_alleles& matalle,
              const t_locations& xlocations,
              const t_locations& pos,
              const double patsum,
              const double matsum,
              std::mt19937& engine)
{
  std::uniform_int_distribution<int> dist(0, 1);
  bool cur_allele = (dist(engine) != 0);

  if (xlocations.size() == 0) {
    if (cur_allele == 0) {
      return matsum;
    }
    else {
      return patsum;
    }
  }

  std::size_t n_loci = pos.size();
  t_locations xloc(xlocations);
  xloc.push_back(*(pos.end() - 1) + 1e-6);
  double sum = 0.0;

  std::size_t ix1 = 0;
  std::size_t ix2;
  for (std::size_t j = 0; j < xloc.size(); ++j){
    const auto it = std::upper_bound(pos.begin() + ix1, pos.end(), xloc[j]);
    ix2 = it - pos.begin();
    if (cur_allele){
      for (std::size_t i = ix1; i != ix2; ++i) {
        sum += patalle[i];
      }
    } else{
      for (std::size_t i = ix1; i != ix2; ++i) {
        sum += matalle[i];
      }
    }
    cur_allele = !cur_allele;
    ix1 = ix2;
  }
  return sum;
}

