#ifndef COST_MATRIX_HPP
#define COST_MATRIX_HPP

#include "vec3_struct.hpp"

#include <vector>

class CostMatrix {
    public:
        CostMatrix();
        CostMatrix(size_t rrtz_iterations);
        void generateCostMatrix();
        void generatePathMatrix();
        void loadViewpoints(std::string filename, vec3 start);
        void loadSimpleCostMatrix(std::string filename);
        void saveSimpleCostMatrix(std::string filename);
        void loadCostMatrix(std::string filename);
        void saveCostMatrix(std::string filename);
        void loadPathMatrix(std::string filename);
        void savePathMatrix(std::string filename);
        float getCost(size_t i, size_t j, size_t k);
        float getSimpleCost(size_t i, size_t j);
        std::vector<vec3> getPath(size_t i, size_t j);
        size_t getNVP() { return this->n_vp; }

        size_t rrtz_iterations; // number of iterations to run rrtz
        size_t n_vp; // number of viewpoints

        // coverage viewpoints
        std::vector<vec3> viewpoints;
        std::vector<vec3> viewpoint_dirs;

        // (n, n, n) matrix of costs from viewpoint j to k given i->j->k.
        std::vector<std::vector<std::vector<float>>> cost_matrix; 

        // (n, n) matrix of costs from viewpoint i to j.
        std::vector<std::vector<float>> simple_cost_matrix;

        // (n, n, pathlength) matrix (symmetric) paths from viewpoint i to j.
        std::vector<std::vector<std::vector<vec3>>> path_matrix;

};

#endif // COST_MATRIX_HPP