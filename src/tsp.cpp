#include "cost_matrix.hpp"
#include "tsp.hpp"

#include <algorithm>
#include <iostream>
#include <string>

TSP::TSP() {
    this->cm = CostMatrix();
    this->path = std::vector<size_t>(1, 0); // start with 'start' viewpoint
}

TSP::TSP(CostMatrix cm) {
    this->cm = cm;
    this->path = std::vector<size_t>(1, 0); // start with 'start' viewpoint
    this->n_vp = cm.getNVP();
}

void TSP::loadCM(int vgd, vec3 start) {
    std::string viewpoint_file = "../data/coverage_viewpoint_sets/coverage_" + std::to_string(vgd) + "m_vp_set.csv";
    std::string cost_matrix_file = "../data/tsp/" + std::to_string(vgd) + "m_cost_matrix.csv";
    std::string path_matrix_file = "../data/tsp/" + std::to_string(vgd) + "m_path_matrix.csv";
    this->cm.loadViewpoints(viewpoint_file, start);
    this->cm.loadCostMatrix(cost_matrix_file);
    this->cm.loadPathMatrix(path_matrix_file);
    this->n_vp = cm.getNVP();
}

void TSP::greedyInit() {
    // Greedy algorithm
    // https://en.wikipedia.org/wiki/Greedy_algorithm
    // https://en.wikipedia.org/wiki/Travelling_salesman_problem
    for (size_t i = 1; i < this->n_vp; i++) {
    // for (size_t i = 1; i < 5; i++) {
        this->nodes.push_back(i);
    }
    while (this->nodes.size() > 0) {
    // for (size_t j = 0; j < 5; j++) {
    //     if (this->nodes.size() == 0) {
    //         break;
    //     }
        std::cout << "Nodes left: " << this->nodes.size() << std::endl;
        for (size_t i = 0; i < this->nodes.size(); i++) {
            std::cout << this->nodes[i] << " ";
        }
        std::cout << std::endl;
        // get best node to insert to path
        float best_cost;
        size_t node_idx;
        size_t idx_to_insert = this->nearest(best_cost, node_idx);

        std::cout << "Node to insert: " << this->nodes[node_idx] << " at index: " << idx_to_insert << " with cost: " << best_cost << std::endl;

        // get iterators to perform insertion
        auto it_erase = this->nodes.begin() + node_idx;
        auto it_insert = this->path.begin() + idx_to_insert;

        // insert node into path
        this->path.insert(it_insert, this->nodes[node_idx]);

        // erase node from nodes
        this->nodes.erase(it_erase);
        std::cout << "Greedy path: ";
        for (size_t i = 0; i < this->path.size(); i++) {
            std::cout << this->path[i] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Greedy path: ";
    for (size_t i = 0; i < this->path.size(); i++) {
        std::cout << this->path[i] << " ";
    }
    std::cout << std::endl << "Cost: " << this->pathCost() << std::endl;
}

void TSP::twoOpt() {
    // 2-opt algorithm
    // https://en.wikipedia.org/wiki/2-opt
    // https://en.wikipedia.org/wiki/Travelling_salesman_problem
    float last_cost = std::numeric_limits<float>::max();
    float new_cost = this->pathCost();
    while (new_cost < last_cost) {
        for (size_t idx0 = 1; idx0 < this->path.size() - 2; idx0++) {
            for (size_t idx1 = 3; idx1 < this->path.size(); idx1++) {
                // get iterators to perform swap
                auto it0 = this->path.begin() + idx0;
                auto it1 = this->path.begin() + idx1;

                // swap
                std::reverse(it0, it1);

                // calculate new cost
                last_cost = new_cost;
                new_cost = this->pathCost();

                // if new cost is worse, swap back
                if (new_cost >= last_cost) {
                    std::reverse(it0, it1);
                }
            }
        }
    }

    std::cout << "Two Opt Path: ";
    for (size_t i = 0; i < this->path.size(); i++) {
        std::cout << this->path[i] << " ";
    }
    std::cout << std::endl << "Cost: " << new_cost << std::endl;
}

void TSP::getPath(std::vector<std::vector<float>>& path) {
    // get path as vector of vec3
    for (size_t i = 0; i < this->path.size(); i++) {
        std::vector<float> vp;
        vp.push_back(this->cm.viewpoints[this->path[i]].x);
        vp.push_back(this->cm.viewpoints[this->path[i]].y);
        vp.push_back(this->cm.viewpoints[this->path[i]].z);
        vp.push_back(this->cm.viewpoint_dirs[this->path[i]].x);
        vp.push_back(this->cm.viewpoint_dirs[this->path[i]].y);
        vp.push_back(this->cm.viewpoint_dirs[this->path[i]].z);
        path.push_back(vp);
    }
}

float TSP::insertionCost(size_t idx, size_t insert_idx) {
    // cost to insert idx at insert_idx position in path
    // Before:
    // ( insertion_idx - 1 ) --- ( insertion_idx )
    // After
    // ( insertion_idx - 1 ) --- ( idx ) --- ( insertion_idx )

    // can't insert at first position (should never happen, but protect)
    if (insert_idx == 0) {
        return INFINITY;
    }

    // out of bounds
    if (insert_idx > this->path.size()) {
        return INFINITY;
    }

    // insert at end of path (easy case)
    if (insert_idx == this->path.size()) {
        if (this->path.size() == 1) {
            return this->cm.getSimpleCost(this->path[0], idx); // args should be 0 and idx
        } else {
            return this->cm.getCost(this->path[insert_idx-2], this->path[insert_idx-1], idx);
        }
    }

    // insert somewhere in the middle of the path

    // we will look at the difference between pre and post costs of path subsection
    float before_cost = 0.0f;
    float after_cost = 0.0f;

    // inserting one after the start (affeects middle leg cost)
    if (insert_idx == 1) {
        // there is only one node before insertion point
        // std::cout << "One after start, checking vp idxs: " << this->path[insert_idx-1] << " " << this->path[insert_idx] << " " << idx << std::endl;
        before_cost += this->cm.getSimpleCost(this->path[insert_idx-1], this->path[insert_idx]);
        after_cost += this->cm.getSimpleCost(this->path[insert_idx-1], idx) + this->cm.getCost(this->path[insert_idx-1], idx, this->path[insert_idx]);
    } else {
        // there are two consecutive nodes before insertion point
        // std::cout << "more than one after start, checking vp idxs: " << this->path[insert_idx-2] << " " << this->path[insert_idx-1] << " " << this->path[insert_idx] << " " << idx << std::endl;
        before_cost += this->cm.getCost(this->path[insert_idx-2], this->path[insert_idx-1], this->path[insert_idx]);
        after_cost += this->cm.getCost(this->path[insert_idx-2], this->path[insert_idx-1], idx) + this->cm.getCost(this->path[insert_idx-1], idx, this->path[insert_idx]);
    }

    // segment after insertion point ( insertion_idx ) --- ( insertion_idx + 1 )
    if (insert_idx < this->path.size() - 2) {
        // std::cout << "Not inserting right before end: " << this->path[insert_idx-1] << " " << this->path[insert_idx] << " " << this->path[insert_idx+1] << std::endl;
        before_cost += this->cm.getCost(this->path[insert_idx-1], this->path[insert_idx], this->path[insert_idx+1]);
        after_cost += this->cm.getCost(idx, this->path[insert_idx], this->path[insert_idx+1]);
    }

    return after_cost - before_cost;
}
size_t TSP::nearestNeighbor(size_t idx, float& best_cost) {
    // takes vp_idx and returns the best place to insert it into path (for best_cost cost)
    // initialize best values
    best_cost = std::numeric_limits<float>::max();
    size_t best_idx = -1; // largest ulong (probably out of range)

    // find best place to insert idx into path lowest cost
    for (size_t i = 1; i < this->path.size() + 1; i++) {
        float cost = this->insertionCost(idx, i);
        if (cost < best_cost) {
            best_cost = cost;
            best_idx = i;
        }
    }
    return best_idx;
}

size_t TSP::nearest(float& best_cost, size_t& node_idx) {
    // returns the best position in path to insert node_idx -- will cost best_cost
    best_cost = std::numeric_limits<float>::max();
    size_t best_idx = -1; // place in this->path to insert viewpoint
    node_idx = -1; // node index in this->nodes to insert into path

    // edge conditions : path size is 1
    // find the best vp idx in nodes to insert into path with lowest cost
    for (size_t i = 0; i < this->nodes.size(); i++) {
        // find best place to insert idx into path lowest cost
        float insertion_cost = 0.0f; // cost after inserting 'nodes[i]' into path at closest/cheapest location
        size_t best_insertion_idx = this->nearestNeighbor(this->nodes[i], insertion_cost);
        std::cout << "Insertion cost for node " << this->nodes[i] << " at index " << best_insertion_idx << " is " << insertion_cost << std::endl;

        if (insertion_cost < best_cost) {
            best_cost = insertion_cost;
            best_idx = best_insertion_idx;
            node_idx = i;
        }
    }
    return best_idx;
}

float TSP::pathCost() {
    // if no viewpoints or one viewpoint, there is no traversal so cost is 0
    if (this->path.size() == 0 || this->path.size() == 1) {
        return 0.0f;
    }

    if (this->path.size() == 2) {
        return this->cm.getSimpleCost(this->path[0], this->path[1]);
    }

    // accumulate cost (get simple cost for first two viewpoints)
    float cost = this->cm.getSimpleCost(this->path[0], this->path[1]);

    // set up traversal vars
    size_t prev_idx = this->path[0]; // this should always be 0
    size_t this_idx = this->path[1];
    size_t next_idx;
    for (size_t i = 2; i < this->path.size(); i++) {
        next_idx = this->path[i];

        // accumulate cost
        cost += this->cm.getCost(prev_idx, this_idx, next_idx);

        // set up traversal vars for next iteration
        prev_idx = this_idx;
        this_idx = next_idx;
    }
    return cost;
}

float TSP::pathCost(std::vector<size_t>& path) {
    // if no viewpoints or one viewpoint, there is no traversal so cost is 0
    if (path.size() == 0 || path.size() == 1) {
        return 0.0f;
    }

    if (path.size() == 2) {
        return this->cm.getSimpleCost(path[0], path[1]);
    }

    // accumulate cost (get simple cost for first two viewpoints)
    float cost = this->cm.getSimpleCost(path[0], path[1]);

    // set up traversal vars
    size_t prev_idx = path[0]; // this should always be 0
    size_t this_idx = path[1];
    size_t next_idx;
    for (size_t i = 2; i < path.size(); i++) {
        next_idx = path[i];

        // accumulate cost
        cost += this->cm.getCost(prev_idx, this_idx, next_idx);

        // set up traversal vars for next iteration
        prev_idx = this_idx;
        this_idx = next_idx;
    }
    return cost;
}