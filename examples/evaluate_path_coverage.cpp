#include "cuda_kernels.h"
#include "obs.hpp"
#include "triangle_struct.hpp"
#include "utils.hpp"
#include "viewpoint_struct.hpp"

#include <iostream>
#include <cmath>
#include <string>

int main() {
    // test solutions
    // std::cout << "k_1_f_1 coverage=" << compute_coverage_path("test/k_1_0_f_1_0_.csv") << std::endl;
    // std::cout << "k_1000_f_1 coverage=" << compute_coverage_path("test/k_1000_0_f_1_0_.csv") << std::endl;
    std::vector<std::string> files = {
        "2m_local",
        "2m_global",
        "4m_local",
        "4m_global",
        "8m_local",
        "8m_global",
        "16m_local",
        "16m_global"
    };

    for (std::string file : files) {
        std::vector<bool> coverage;
        std::cout << std::to_string(compute_coverage_path( file + "/k_100_0_f_1_0_.csv", coverage)) << std::endl;
    }


    // std::vector<std::vector<float>> path_coverage;
    // std::vector<float> weights = {-3.0f, -2.5f, -2.0f, -1.5f, -1.0f, -0.5f, 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f};
    // for (size_t k = 0; k < weights.size(); k++) {
    //     for (size_t f = 0; f < weights.size(); f++) {
    //         float kw = std::pow(10, weights[k]);
    //         float fw = std::pow(10, weights[f]);
    //         std::string file = "k_" + getnum(kw) + "_f_" + getnum(fw) + "_.csv";
    //         std::cout << "Computing coverage for " << file << " ";
    //         std::vector<float> coverage = {compute_coverage_path("pareto_front_2m/" + file)};
    //         std::cout << " coverage=" << std::to_string(coverage[0]) << std::endl;
    //         path_coverage.push_back(coverage);
    //     }
    // }
    // saveCSV("../knot_ocp/pareto_front/pareto_front_2m_local/cov.csv", path_coverage);
}