#include <Rcpp.h>
#include <random>
#include <vector>
#include <tuple>

typedef std::vector<double> t_alleles;
typedef std::vector<double> t_locations;
typedef std::tuple<double, double, int, double, int> t_tuple;

t_tuple bcgv_routine(std::vector<std::vector<std::vector<double>>> individual,
                     std::vector<std::vector<double>> positions,
                     std::vector<std::vector<double>> effects,
                     const unsigned int n_gam,
                     const double se_level,
                     const unsigned int min_rep,
                     const unsigned int max_rep,
                     const unsigned int m,
                     const double p,
                     const int seed);

double meisum(const t_alleles& patalle,
              const t_alleles& matalle,
              const t_locations& xlocations,
              const t_locations& pos,
              const double patsum,
              const double matsum,
              std::mt19937& engine);
