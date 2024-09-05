#include "cost_matrix.hpp"
#include "limit_struct.hpp"
#include "obs.hpp"
#include "rrtz.hpp"
#include "utils.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <stddef.h>
#include <vector>


CostMatrix::CostMatrix() {
    this->rrtz_iterations = 2000;
}

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
                if (i != j && j != k && i != k && j != 0 && k != 0) {
                    // if any of the viewpoints are the same, leave cost at infinity
                    // if j or k are the start viewpoint, leave cost at infinity
                    // otherwise, calculate the cost as the sum of the path lengths
                    if (this->simple_cost_matrix[i][j] != INFINITY && this->simple_cost_matrix[j][k] != INFINITY) {
                        vec3 last_dir = 
                            this->path_matrix[i][j].at(path_matrix[i][j].size()-1) -
                            this->path_matrix[i][j].at(path_matrix[i][j].size()-2);
                        vec3 next_dir = 
                            this->path_matrix[j][k].at(1) - 
                            this->path_matrix[j][k].at(0);
                        this->cost_matrix[i][j][k] = 
                            // this->simple_cost_matrix[i][j] +
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
    Limit limits = { -5.0f, 10.0f, -5.0f, 15.0f, -5.0f, 10.0f };

    // load obstacles
    std::vector<OBS> obsVec;
    loadConvexStationOBS(obsVec);


    std::ostringstream message;
    auto start = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    for (size_t vp_idx0 = 0; vp_idx0 < this->n_vp; vp_idx0++) {
        for (size_t vp_idx1 = 0; vp_idx1 < this->n_vp; vp_idx1++) {
            // auto prev_now = now;
            now = std::chrono::high_resolution_clock::now();
            // Calculate the duration
            std::chrono::duration<double> duration = now - start;
            double time_taken = duration.count();

            // Output the duration in seconds
            size_t leftinrow = this->n_vp - vp_idx1;
            size_t leftincol = this->n_vp - vp_idx0;
            size_t itemsleft = leftinrow + leftincol * leftincol / 2;
            size_t totalitems = (this->n_vp + 1) * (this->n_vp + 1) / 2;

            double avg_loop_period = time_taken / (totalitems - itemsleft);
            double seconds_remaining = avg_loop_period * (itemsleft);
            double minutes_remaining = seconds_remaining / 60.0;
            seconds_remaining = std::fmod(seconds_remaining, 60.0);
            message << " Time remaining: " << int(minutes_remaining) << "m " << int(seconds_remaining) << "s";
            double progress = 1.0 - static_cast<double>(itemsleft) / static_cast<double>((this->n_vp + 1) * (this->n_vp + 1) / 2);
            displayProgressBar(progress, 50, message);
            message.str("");
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
                vec3 start = this->viewpoints[vp_idx0].pose;
                vec3 goal = this->viewpoints[vp_idx1].pose;
                RRTZ rrtz = RRTZ(start, goal, obsVec, limits, this->rrtz_iterations);
                if (rrtz.run(path)) {
                    // if rrtz runs, add path and path cost to matrix
                    this->path_matrix[vp_idx0][vp_idx1] = path;
                    this->simple_cost_matrix[vp_idx0][vp_idx1] = rrtz.getBestCost();
                // } else {
                //     // if rrtz fails, add empty path and infinite cost to matrix
                //     std::cout << "Failed to find path from " << this->viewpoints[vp_idx0].toString() << " to " << this->viewpoints[vp_idx1].toString() << std::endl;
                }
            }
        }
    }
}

void CostMatrix::loadViewpoints(std::string filename, Viewpoint start) {
    // load the viewpoints from a file
    std::vector<std::vector<float>> viewpoint_data;
    loadCSV(filename, viewpoint_data, 7, ',');
    this->n_vp = viewpoint_data.size() + 1;

    this->viewpoints.clear();

    // add start viewpoint
    this->viewpoints.push_back(start);
    for (size_t i = 0; i < viewpoint_data.size(); i++) {
        viewpoints.push_back(Viewpoint(
            vec3(viewpoint_data[i][0], viewpoint_data[i][1], viewpoint_data[i][2]),
            vec3(viewpoint_data[i][3], viewpoint_data[i][4], viewpoint_data[i][5]),
            viewpoint_data[i][6]
        ));
        std::cout << "Viewpoint " << i << " mm: " << viewpoints.back().module_idx << std::endl;
    }
    std::cout << "n_vp = " << this->n_vp << " should equal " << this->viewpoints.size() << std::endl;
}

void CostMatrix::loadSimpleCostMatrix(std::string filename) {
    // load the simple cost matrix from a file
    std::vector<std::vector<float>> simple_cost_matrix_data;
    loadCSV(filename, simple_cost_matrix_data, 3, ',');
    this->simple_cost_matrix.clear();
    this->simple_cost_matrix = std::vector<std::vector<float>>(this->n_vp, std::vector<float>(this->n_vp, INFINITY));

    for (size_t i = 0; i < simple_cost_matrix_data.size(); i++) {
        size_t idx0 = simple_cost_matrix_data[i][0];
        size_t idx1 = simple_cost_matrix_data[i][1];
        float cost = simple_cost_matrix_data[i][2];
        this->simple_cost_matrix[idx0][idx1] = cost;
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
    saveCSV(filename, simple_cost_matrix_data);
}

void CostMatrix::loadCostMatrix(std::string filename) {
    // load cost matrix from file
    std::vector<std::vector<float>> cost_matrix_data;
    loadCSV(filename, cost_matrix_data, 4, ',');
    this->cost_matrix.clear();
    this->cost_matrix = std::vector<std::vector<std::vector<float>>>(this->n_vp, std::vector<std::vector<float>>(this->n_vp, std::vector<float>(this->n_vp, INFINITY)));

    for (size_t i = 0; i < cost_matrix_data.size(); i++) {
        size_t idx0 = cost_matrix_data[i][0];
        size_t idx1 = cost_matrix_data[i][1];
        size_t idx2 = cost_matrix_data[i][2];
        float cost = cost_matrix_data[i][3];
        cost_matrix[idx0][idx1][idx2] = cost;
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
    saveCSV(filename, cost_matrix_data);
}

void CostMatrix::loadPathMatrix(std::string filename) {
    // load path matrix from file
    std::vector<std::vector<float>> path_matrix_data;
    loadCSV(filename, path_matrix_data, 3, ',', true); // raw so number of columns is not checked
    this->path_matrix.clear();
    this->path_matrix = std::vector<std::vector<std::vector<vec3>>>(this->n_vp, std::vector<std::vector<vec3>>(this->n_vp));

    for (size_t i = 0; i < path_matrix_data.size(); i++) {
        size_t idx0 = path_matrix_data[i][0];
        size_t idx1 = path_matrix_data[i][1];

        if ( (path_matrix_data[i].size()-2) % 3 != 0) {
            std::cout << idx0 << " " << idx1 << (path_matrix_data[i].size()-2) % 3 << std::endl;
        }

        std::vector<vec3> path;
        for (size_t j = 2; j < path_matrix_data[i].size(); j+=3) {
            vec3 point = vec3(path_matrix_data[i][j], path_matrix_data[i][j+1], path_matrix_data[i][j+2]);
            path.push_back(point);
            // this->path_matrix[idx0][idx1].push_back(point);
        }
        this->path_matrix[idx0][idx1] = path;
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
            }
            path_matrix_data.push_back(row);
        }
    }
    saveCSV(filename, path_matrix_data);
}

float CostMatrix::getCost(size_t i, size_t j, size_t k) {
    // get the cost from viewpoint i to j to k
    if (i >= this->n_vp || j >= this->n_vp || k >= this->n_vp || // if indexing outside of the cost matrix
        i == j || j == k || i == k || // will be infinity anyways
        j == 0 || k == 0) { // if j or k are the start viewpoint
        std::cout << "Invalid Coordinates: " << i << ", " << j << ", " << k << std::endl;
        return INFINITY;
    }

    return this->cost_matrix[i][j][k];
}

float CostMatrix::getSimpleCost(size_t i, size_t j) {
    // get the simple cost from viewpoint i to j
    if (i >= this->n_vp || j >= this->n_vp ||
        i == j || j == 0) { // if i and j point at the same viewpoint or j is the start viewpoint
        std::cout << "Invalid Coordinates: " << i << ", " << j << std::endl;
        return INFINITY;
    }

    return this->simple_cost_matrix[i][j];
}

std::vector<vec3> CostMatrix::getPath(size_t i, size_t j) {
    // get the path from viewpoint i to j
    if (i >= this->n_vp || j >= this->n_vp) {
        std::cout << "Invalid Coordinates: " << i << ", " << j << std::endl;
        return std::vector<vec3>();
    }

    return this->path_matrix[i][j];
}