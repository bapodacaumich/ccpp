#include "cuda_kernels.h"
#include "obs.hpp"
#include "triangle_struct.hpp"
#include "utils.hpp"
#include "viewpoint_struct.hpp"

#include <iostream>
#include <cmath>
#include <string>

int main(int argc, char** argv) {
    // test solutions
    // std::cout << "k_1_f_1 coverage=" << compute_coverage_path("test/k_1_0_f_1_0_.csv") << std::endl;
    // std::cout << "k_1000_f_1 coverage=" << compute_coverage_path("test/k_1000_0_f_1_0_.csv") << std::endl;
    // std::vector<std::string> files = {
    //     "2m_local",
    //     "2m_global",
    //     "4m_local",
    //     "4m_global",
    //     "8m_local",
    //     "8m_global",
    //     "16m_local",
    //     "16m_global"
    // };


    // for (std::string file : files) {
    //     std::vector<bool> coverage;
    //     std::cout << std::to_string(compute_coverage_path( file + "/k_100_0_f_1_0_.csv", coverage)) << std::endl;
    // }


    std::vector<std::vector<float>> path_coverage;
    // std::vector<float> weights = {-3.0f, -2.5f, -2.0f, -1.5f, -1.0f, -0.5f, 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f};
    std::string folder = argv[1];
    folder += "/";

    std::vector<std::vector<float>> weights;
    loadCSV("../knot_ocp/pareto_front/" + folder + "wcomb.csv", weights, 2, ' ');
    std::cout << weights.size() << std::endl;
    std::cout << weights[0].size() << std::endl;
    for (size_t i =0; i < weights.size(); i++) {
        std::cout << weights[i][0] << " " << weights[i][1] << std::endl;
    }

    // std::vector<std::vector<bool>> coverages;
    for (size_t i = 0; i < weights.size(); i++) {
        std::vector<bool> coverage_per_face;
        float kw = static_cast<float>(weights[i][0]);
        float fw = static_cast<float>(weights[i][1]);
        std::string file = folder + "k_" + getnum(kw) + "_f_" + getnum(fw) + ".csv";
        std::cout << "Computing coverage for " << file << " ";
        std::vector<float> coverage = {compute_coverage_path(file, coverage_per_face)};
        std::cout << " coverage=" << std::to_string(coverage[0]) << std::endl;
        // coverages.push_back(coverage_per_face);
        path_coverage.push_back(coverage);
    }

    // std::vector<size_t> common_uncovered;
    // for (size_t i = 0; i < coverages[0].size(); i++) {
    //     bool is_common = true;
    //     for (size_t j = 0; j < coverages.size(); j++) {
    //         if (coverages[j][i]) {
    //             is_common = false;
    //             break;
    //         }
    //     }
    //     if (is_common) {
    //         common_uncovered.push_back(i);
    //     }
    // }
    // std::cout << "common uncovered idxs: ";
    // for (size_t i = 0; i < common_uncovered.size(); i++) {
    //     std::cout << common_uncovered[i] << "\n";
    // }
    saveCSV("../knot_ocp/pareto_front/" + folder + "cov.csv", path_coverage);
}