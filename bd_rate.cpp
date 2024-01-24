#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <cxxopts.hpp>

#include "ap.h"
#include "integration.h"
#include "interpolation.h"

cxxopts::ParseResult parse_args(int argc, char *argv[]) {
  cxxopts::Options options("bd_rate", "BD-Rate calculater");

  // clang-format off
  options.add_options()
      ("anchor_rate", "Anchor rates", cxxopts::value<std::vector<double> >())
      ("anchor_metric", "Anchor metrics", cxxopts::value<std::vector<double> >())
      ("test_rate", "Test rates", cxxopts::value<std::vector<double> >())
      ("test_metric", "Test metrics", cxxopts::value<std::vector<double> >())
      ("min_overlap", "Minimum overlap", cxxopts::value<double>()->default_value("0.5"))
      ("method", "Interpolation method", cxxopts::value<std::string>()->default_value("akima"))
      ("h,help", "Print help")
      ;
  // clang-format on

  cxxopts::ParseResult result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  auto method = result["method"].as<std::string>();
  if (method != "akima" && method != "cubic") {
    std::cerr << "Unknown method: " << method << std::endl;
    exit(1);
  }

  return result;
}

void check_overlap(const std::vector<double> &x_A, const std::vector<double> &x_B,
                   const float min_overlap) {
  const auto total_x_min = std::min(*std::min_element(x_A.begin(), x_A.end()),
                                    *std::min_element(x_B.begin(), x_B.end()));
  const auto total_x_max = std::max(*std::max_element(x_A.begin(), x_A.end()),
                                    *std::max_element(x_B.begin(), x_B.end()));

  const auto overlap_x_min = std::max(*std::min_element(x_A.begin(), x_A.end()),
                                      *std::min_element(x_B.begin(), x_B.end()));
  const auto overlap_x_max = std::min(*std::max_element(x_A.begin(), x_A.end()),
                                      *std::max_element(x_B.begin(), x_B.end()));

  const auto overlap = std::max(overlap_x_max - overlap_x_min, 0.) / (total_x_max - total_x_min);

  if (overlap == 0) {
    std::cerr << "Curves do not overlap. BD cannot be calculated." << std::endl;
    exit(1);
  } else if (overlap < min_overlap) {
    std::cerr << "Insufficient curve overlap: " << overlap << ". Minimum overlap: " << min_overlap
              << ". " << std::endl;
    exit(1);
  }
}

void prepare_input(const std::vector<double> &rate, const std::vector<double> &metric,
                   alglib::real_1d_array &x, alglib::real_1d_array &y) {
  assert(rate.size() == metric.size());

  // sort (rates, metrics) in ascending order of rates
  std::vector<size_t> indices(rate.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    indices[i] = i;
  }
  std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) { return rate[i] < rate[j]; });

  x.setlength(indices.size());
  y.setlength(indices.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    x[i] = metric[indices[i]];
    y[i] = std::log10(rate[indices[i]]);
  }
}

int main(int argc, char *argv[]) {
  auto args = parse_args(argc, argv);

  auto anchor_rate = args["anchor_rate"].as<std::vector<double>>();
  auto anchor_metric = args["anchor_metric"].as<std::vector<double>>();
  auto test_rate = args["test_rate"].as<std::vector<double>>();
  auto test_metric = args["test_metric"].as<std::vector<double>>();
  auto min_overlap = args["min_overlap"].as<double>();
  auto method = args["method"].as<std::string>();

  // prepare input
  alglib::real_1d_array x_A, y_A;
  prepare_input(anchor_rate, anchor_metric, x_A, y_A);
  alglib::real_1d_array x_B, y_B;
  prepare_input(test_rate, test_metric, x_B, y_B);

  // check overlap
  const auto total_x_min = std::min(x_A[0], x_B[0]);
  const auto total_x_max = std::max(x_A[x_A.length() - 1], x_B[x_B.length() - 1]);

  const auto overlap_x_min = std::max(x_A[0], x_B[0]);
  const auto overlap_x_max = std::min(x_A[x_A.length() - 1], x_B[x_B.length() - 1]);

  const auto overlap = std::max(overlap_x_max - overlap_x_min, 0.) / (total_x_max - total_x_min);
  if (overlap < min_overlap) {
    std::cerr << "Insufficient curve overlap: " << overlap << ". Minimum overlap: " << min_overlap
              << ". " << std::endl;
    exit(1);
  }

  // build spline
  alglib::spline1dinterpolant spline_A, spline_B;
  if (method == "akima") {
    alglib::spline1dbuildakima(x_A, y_A, spline_A);
    alglib::spline1dbuildakima(x_B, y_B, spline_B);
  } else if (method == "cubic") {
    alglib::spline1dbuildcubic(x_A, y_A, spline_A);
    alglib::spline1dbuildcubic(x_B, y_B, spline_B);
  }

  // calc integration
  double inter_A = alglib::spline1dintegrate(spline_A, overlap_x_max) -
                   alglib::spline1dintegrate(spline_A, overlap_x_min);
  double inter_B = alglib::spline1dintegrate(spline_B, overlap_x_max) -
                   alglib::spline1dintegrate(spline_B, overlap_x_min);

  double bd_diff = (inter_B - inter_A) / (overlap_x_max - overlap_x_min);
  double bd_rate = (exp10(bd_diff) - 1) * 100;

  std::cout << bd_rate << "\n";

  return 0;
}
