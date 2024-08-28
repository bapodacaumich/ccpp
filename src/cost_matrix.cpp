#include "cost_matrix.hpp"
#include "limit_struct.hpp"
#include "obs.hpp"
#include "rrtz.hpp"
#include "utils.hpp"

#include <algorithm>
#include <string>
#include <stddef.h>
#include <vector>


CostMatrix::CostMatrix(size_t rrtz_iterations) {
    this->rrtz_iterations = rrtz_iterations;
}

void CostMatrix::generateCostMatrix() {
    // generate the cost matrix from the path matrix (i to j to k)
    // clear and resize cost matrix
    this->cost_matrix.clear();
    this->cost_matrix.resize(this->n_vp, std::vector<std::vector<float>>(this->n_vp, std::vector<float>(this->n_vp, INFINITY)));
    for (size_t i = 0; i < this->n_vp; i++) {
        for (size_t j = 0; j < this->n_vp; j++) {
            for (size_t k = 0; k < this->n_vp; k++) {
                if (i != j && j != k && i != k) {
                    // if any of the viewpoints are the same, leave cost at infinity
                    // otherwise, calculate the cost as the sum of the path lengths
                    if (this->simple_cost_matrix[i][j] != INFINITY && this->simple_cost_matrix[j][k] != INFINITY) {
                        vec3 last_dir = 
                            this->path_matrix[i][j].at(path_matrix[i][j].size()-1) -
                            this->path_matrix[i][j].at(path_matrix[i][j].size()-2);
                        vec3 next_dir = 
                            this->path_matrix[j][k].at(1) - 
                            this->path_matrix[j][k].at(0);
                        this->cost_matrix[i][j][k] = 
                            this->simple_cost_matrix[i][j] + 
                            this->simple_cost_matrix[j][k] + 
                            heading_change(last_dir, next_dir);
                    }
                }
            }
        }
    }
}

void CostMatrix::generatePathMatrix() {
    // generate the path matrix from the viewpoints

    // empty matrices
    this->path_matrix.clear();
    this->path_matrix.resize(this->n_vp, std::vector<std::vector<vec3>>(this->n_vp));
    this->simple_cost_matrix.clear();
    this->simple_cost_matrix.resize(this->n_vp, std::vector<float>(this->n_vp, INFINITY));


    // set wide limits
    Limit limits = { -10.0f, 10.0f, -10.0f, 30.0f, -10.0f, 10.0f };

    // load obstacles
    std::vector<OBS> obsVec;
    loadConvexStationOBS(obsVec);

    for (size_t vp_idx0 = 0; vp_idx0 < this->n_vp; vp_idx0++) {
        for (size_t vp_idx1 = 0; vp_idx1 < this->n_vp; vp_idx1++) {
            // if the viewpoints are the same, add empty path and infinite cost to matrix
            if (vp_idx0 == vp_idx1) { continue; }

            if (vp_idx0 > vp_idx1) { 
                // get reverse path for flipped viewpoints (We already planned this path in reverse)
                std::vector<vec3> path = this->path_matrix[vp_idx1][vp_idx0];
                std::reverse(path.begin(), path.end());
                this->path_matrix[vp_idx0][vp_idx1] = path;
                this->simple_cost_matrix[vp_idx0][vp_idx1] = this->simple_cost_matrix[vp_idx1][vp_idx0];
            } else {
                // setup and run rrtz
                std::vector<vec3> path;
                vec3 start = viewpoints[vp_idx0];
                vec3 goal = viewpoints[vp_idx1];
                RRTZ rrtz = RRTZ(start, goal, obsVec, limits, this->rrtz_iterations);
                if (!rrtz.run(path)) {
                    // if rrtz runs, add path and path cost to matrix
                    this->path_matrix[vp_idx0][vp_idx1] = path;
                    this->simple_cost_matrix[vp_idx0][vp_idx1] = rrtz.getBestCost();
                } else {
                    // if rrtz fails, add empty path and infinite cost to matrix
                    std::cout << "Failed to find path from " << this->viewpoints[vp_idx0].toString() << " to " << this->viewpoints[vp_idx1].toString() << std::endl;
                }
            }
        }
    }
}

void CostMatrix::loadViewpoints(std::string filename) {
    // load the viewpoints from a file
    std::vector<std::vector<float>> viewpoint_data;
    loadCSV(filename, viewpoint_data, 3, ',');
    this->n_vp = viewpoint_data.size();

    this->viewpoints.clear();
    for (size_t i = 0; i < viewpoint_data.size(); i++) {
        vec3 vp = vec3(viewpoint_data[i][0], viewpoint_data[i][1], viewpoint_data[i][2]);
        this->viewpoints.push_back(vp);
    }
}

void CostMatrix::saveSimpleCostMatrix(std::string filename) {
    // save the simple cost matrix to a file
    std::vector<std::vector<float>> simple_cost_matrix_data; // flattened simple cost matrix
    for (size_t i = 0; i < this->n_vp; i++) {
        for (size_t j = 0; j < this->n_vp; j++) {
            std::vector<float> row;
            row.push_back(i);
            row.push_back(j);
            row.push_back(this->simple_cost_matrix[i][j]);
            simple_cost_matrix_data.push_back(row);
        }
    }
}

void CostMatrix::saveCostMatrix(std::string filename) {
    // save the cost matrix to a file
    std::vector<std::vector<float>> cost_matrix_data; // flattened cost matrix
    for (size_t i = 0; i < this->n_vp; i++) {
        for (size_t j = 0; j < this->n_vp; j++) {
            for (size_t k = 0; k < this->n_vp; k++) {
                std::vector<float> row;
                row.push_back(i);
                row.push_back(j);
                row.push_back(k);
                row.push_back(this->cost_matrix[i][j][k]);
                cost_matrix_data.push_back(row);
            }
        }
    }
}

void CostMatrix::savePathMatrix(std::string filename) {
    // save the path matrix to a file
    std::vector<std::vector<float>> path_matrix_data; // flattened path matrix
    for (size_t i = 0; i < this->n_vp; i++) {
        for (size_t j = 0; j < this->n_vp; j++) {
            std::vector<float> row;
            row.push_back(i);
            row.push_back(j);
            for (size_t k = 0; k < this->path_matrix[i][j].size(); k++) {
                row.push_back(this->path_matrix[i][j][k].x);
                row.push_back(this->path_matrix[i][j][k].y);
                row.push_back(this->path_matrix[i][j][k].z);
                path_matrix_data.push_back(row);
            }
        }
    }
}