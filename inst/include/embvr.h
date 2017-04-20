#include <Rcpp.h>
#include <random>
#include <vector>
#include <tuple>

typedef std::vector<double> vec;
typedef std::vector<double> t_alleles;
typedef std::vector<double> t_locations;
typedef std::tuple<double, double, int, double, int> t_tuple;
typedef std::tuple<std::size_t, std::vector<double>, std::vector<double>, bool > t_xoparam;
typedef std::vector<std::vector<t_alleles>> t_individual;
typedef std::vector<std::vector<double>> t_effects;
typedef std::vector<std::vector<double>> t_positions;
