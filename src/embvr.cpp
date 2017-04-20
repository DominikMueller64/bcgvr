// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Meiosis)]]
#include <Rcpp.h>
#include <random>
#include <vector>
#include <tuple>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <chrono>
#include <Meiosis.h>

#include "embvr.h"

std::random_device rdev;
std::mt19937 engine;

t_xoparam get_xoparam(const t_positions& positions)
{
  const bool obligate_chiasma = true;
  const std::size_t n_chr = positions.size();
  vec L(n_chr), Lstar(n_chr);
  for (std::size_t i = 0; i < n_chr; ++i) {
    L[i] = std::ceil(positions[i].back());
    Lstar[i] = L[i];
  }

  return std::make_tuple(n_chr, L, Lstar, obligate_chiasma);
}

std::vector<vec> transform_geno(t_individual& individual,
                                const t_effects& effects,
                                const std::size_t n_chr)

{
  // Transform genotypes by multiplication with locus effects and compute sums.
  std::vector<vec> sums(individual.size(), vec(n_chr));
  for (std::size_t g = 0; g < individual.size(); ++g) {
    auto& gamete = individual[g];
    for (std::size_t i = 0; i < n_chr; ++i) {
      auto& chromatid = gamete[i];
      auto& leff = effects[i];
      std::transform(chromatid.begin(), chromatid.end(),
                     leff.begin(), chromatid.begin(), std::multiplies<double>());
      sums[g][i] = std::accumulate(chromatid.begin(), chromatid.end(), 0.0);
    }
  }
  return sums;
}

// Simulation of meiosis, effects are directly accumulated.
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

// Evaluation of embv.
t_tuple embv_routine(t_individual individual,
                     t_positions positions,
                     t_effects effects,
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

  engine.seed(seed);
  const auto& xof = Meiosis::crossover<std::vector<double>>;

  const auto xoparam = get_xoparam(positions);
  const auto n_chr = std::get<0>(xoparam);
  const auto& L = std::get<1>(xoparam);
  const auto& Lstar = std::get<2>(xoparam);
  const auto obligate_chiasma = std::get<3>(xoparam);

  const auto& sums = transform_geno(individual, effects, n_chr);

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
        const auto& xlocations = xof(L[i], m, p, obligate_chiasma, Lstar[i], engine);
        value += meisum(individual[0][i], individual[1][i], xlocations,
                        positions[i], sums[0][i], sums[1][i], engine);
      }
      if (value > max_value) max_value = value;
    }

    // Update online statistics.
    ++n;
    delta = max_value - mean;
    mean += delta / n;
    delta2 = max_value - mean;
    M2 += delta * delta2;
    se = std::sqrt(M2 / (n * (n - 1))); // nan for n = 1.
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

  return std::make_tuple(mean, se, n, duration / 1000.0, seed);
}

// Wrapper for interfacing `embv` with R.
//' @rdname embv
// [[Rcpp::export]]
Rcpp::List embv(t_individual individual,
                t_positions positions,
                t_effects effects,
                const int n_gam,
                const double se_level,
                const int min_rep,
                const int max_rep,
                const int m = 0,
                const double p = 1.0,
                Rcpp::Nullable<Rcpp::IntegerVector> seed = R_NilValue)
{
  if (individual.size() != 2) Rcpp::stop("'individual' must have two gametes (list of two).");
  if (n_gam < 1) Rcpp::stop("'n_gam' must be >= 1.");
  if (se_level <= 0.0) Rcpp::stop("'se_level' must be > 0.0.");
  if (min_rep < 1) Rcpp::stop("'min_rep' must be > 0.");
  if (max_rep < 1) Rcpp::stop("'max_rep' must be > 0.");
  if (min_rep > max_rep) Rcpp::stop("'min_rep' must be smaller than or equal to 'max_rep'.");

  if (m < 0) Rcpp::stop("'m' must be >= 0.");
  if (p < 0.0 || p > 1.0) Rcpp::stop("'p' must be in [0, 1].");

  int seed_;
  if (seed.isNotNull()) {
    seed_ = Rcpp::as<std::vector<int>>(seed)[0];
  } else {
    seed_ = std::random_device{}();
  }

  const auto& ret = embv_routine(individual, positions, effects, n_gam, se_level,
                                 min_rep, max_rep, m, p, seed_);

  return Rcpp::List::create(Rcpp::Named("embv") = std::get<0>(ret),
                            Rcpp::Named("se") = std::get<1>(ret),
                            Rcpp::Named("n") = std::get<2>(ret),
                            Rcpp::Named("time (ms)") = std::get<3>(ret),
                            Rcpp::Named("seed") = std::get<4>(ret));
}



// Double haploid breedinv values
//' @rdname dhbv
// [[Rcpp::export]]
Rcpp::NumericVector dhbv(t_individual individual,
                         t_positions positions,
                         t_effects effects,
                         const unsigned int n_rep,
                         const unsigned int m = 0,
                         const double p = 1.0,
                         Rcpp::Nullable<Rcpp::IntegerVector> seed = R_NilValue
                         )
{
  if (individual.size() != 2) Rcpp::stop("'individual' must have two gametes (list of two).");
  if (n_rep < 2) Rcpp::stop("'n_rep' must be >= 2.");
  if (m < 0) Rcpp::stop("'m' must be >= 0.");
  if (p < 0.0 || p > 1.0) Rcpp::stop("'p' must be in [0, 1].");

  int seed_;
  if (seed.isNotNull()) {
    seed_ = Rcpp::as<std::vector<int>>(seed)[0];
  } else {
    seed_ = std::random_device{}();
  }

  engine.seed(seed_);
  const auto& xof = Meiosis::crossover<std::vector<double>>;

  const auto xoparam = get_xoparam(positions);
  const auto n_chr = std::get<0>(xoparam);
  const auto& L = std::get<1>(xoparam);
  const auto& Lstar = std::get<2>(xoparam);
  const auto obligate_chiasma = std::get<3>(xoparam);

  const auto& sums = transform_geno(individual, effects, n_chr);


  Rcpp::NumericVector res(n_rep);

  // // Online statistics.
  // double M2 = 0.0;
  // double mean = 0.0;
  // double delta, delta2;

  for (std::size_t n = 0; n < n_rep; ++n) {
    double value = 0.0;
    for (std::size_t i = 0; i < n_chr; ++i) {
      const auto& xlocations = xof(L[i], m, p, obligate_chiasma, Lstar[i], engine);
      value += meisum(individual[0][i], individual[1][i], xlocations,
                      positions[i], sums[0][i], sums[1][i], engine);
    }
    res[n] = 2.0 * value;
    // // Update online statistics.
    // delta = value - mean;
    // // Rcpp::Rcout << "delta: " << delta <<std::endl;
    // mean += delta / (n + 1);
    // delta2 = value - mean;
    // M2 += delta * delta2;
  }

  // return std::sqrt(4.0 * M2 / n_rep); // nan for n = 1.
  return res;
}

