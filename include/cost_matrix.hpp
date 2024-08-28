#ifndef COST_MATRIX_HPP
#define COST_MATRIX_HPP

#include "vec3_struct.hpp"

#include <vector>
#include <stddef.h>

class CostMatrix {
    public:
        CostMatrix(size_t rrtz_iterations);
        void generateCostMatrix();
        void generatePathMatrix();
        void loadViewpoints(std::string filename);
        void saveSimpleCostMatrix(std::string filename);
        void saveCostMatrix(std::string filename);
        void savePathMatrix(std::string filename);

    private:
        size_t rrtz_iterations; // number of iterations to run rrtz
        size_t n_vp; // number of viewpoints

        // coverage viewpoints
        std::vector<vec3> viewpoints;

        // (n, n, n) matrix of costs from viewpoint i to j to k.
        std::vector<std::vector<std::vector<float>>> cost_matrix; 
        std::vector<std::vector<float>> simple_cost_matrix;

        // (n, n, pathlength) matrix (symmetric) paths from viewpoint i to j.
        std::vector<std::vector<std::vector<vec3>>> path_matrix;

};

#endif // COST_MATRIX_HPP