#include <Rcpp.h>
#include <random>
#include <vector>

std::random_device rdev;
std::mt19937 engine;

typedef std::vector<double> t_alleles;
typedef std::vector<double> t_locations;


Rcpp::List bcgv(std::vector<std::vector<std::vector<double>>> individual,
                std::vector<std::vector<double>> positions,
                std::vector<std::vector<double>> locus_effects,
                const int n_gam,
                const double se_level,
                const int min_rep,
                const int max_rep,
                const int m = 0,
                const double p = 0.0,
                const int seed = std::random_device{}();
                );


double meisum(const t_alleles& patalle,
              const t_alleles& matalle,
              const t_locations& xlocations,
              const t_locations& pos,
              const double patsum,
              const double matsum,
              std::mt19937& engine);
