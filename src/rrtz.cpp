#include "limit_struct.hpp"
#include "node3d_struct.hpp"
#include "obs.hpp"
#include "plane_struct.hpp"
#include "rrtz.hpp"
#include "utils.hpp"
#include "vec3_struct.hpp"

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <limits>

RRTZ::RRTZ() {
    this->start = vec3();
    this->goal = vec3();
    this->max_nodes = 0;
    this->min_iter = 0;
    this->max_iter = 0;
    this->tree_nodes = std::vector<Node3D>();
    this->obs = std::vector<OBS>();
}

RRTZ::RRTZ(
    vec3 start,
    vec3 goal,
    std::vector<OBS> obs,
    Limit limits,
    size_t max_nodes,
    size_t min_iter,
    size_t max_iter,
    bool debug
    ) {

    // save problem config
    this->start = start;
    this->goal = goal;
    this->max_nodes = max_nodes;
    this->min_iter = min_iter; // not using TODO: implement convergence criteria
    this->max_iter = max_iter; // not using 

    // initialize tree
    this->tree_nodes = std::vector<Node3D>();
    this->tree_nodes.push_back(Node3D(start, start, 0, 0, 0)); // root node

    // get env
    this->obs = obs;
    this->limits = limits;

    // initialize best solution
    this->best_cost = std::numeric_limits<float>::infinity();
    this->best_idx = -1;

    // pass in debug setting
    this->debug = debug;
}

bool RRTZ::run(std::vector<vec3>& path) {
    // initialize best solution
    while (!this->terminate()) {
        // random number between 0 and 1
        float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
        float cost = std::numeric_limits<float>::infinity();

        if (this->debug) {std::cout << "\rIteration=" << this->tree_nodes.size();}

        if (r < 0.5) { 
            // extend to goal
            Plane plane = this->sample_plane(this->goal);
            if (this->extend(plane, cost)) {
                // success to goal --> save index of this leaf node
                if (cost < this->best_cost) {
                    this->best_cost = cost;
                    if (this->debug) {std::cout << " Improved Cost=" << this->best_cost << std::endl;}
                    this->best_idx = this->tree_nodes.size() - 1;
                }
            }
        } else {
            // extend to random plane
            bool success = false;
            while (!success) {
                Plane plane = this->sample_plane();
                success = this->extend(plane, cost);
            }
        }
    }

    // reconstruct path
    if (this->best_cost != std::numeric_limits<float>::infinity()) {
        this->reconstruct_path(this->best_idx, path);
        return true;
    }
    return false;
}

bool RRTZ::terminate() {
    // TODO: implement convergence criteria condition
    // currently: end at max_iterations
    return this->tree_nodes.size() >= this->max_nodes;
}

void RRTZ::reconstruct_path(size_t parent_idx, std::vector<vec3>& path) {
    path.push_back(this->tree_nodes[parent_idx].end);
    while (parent_idx != 0) {
        path.insert(path.begin(),this->tree_nodes[parent_idx].origin);
        parent_idx = this->tree_nodes[parent_idx].parent_idx;
    }
}
Plane RRTZ::sample_plane() {
    vec3 pose = vec3(
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * (this->limits.xmax - this->limits.xmin) + this->limits.xmin,
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * (this->limits.ymax - this->limits.ymin) + this->limits.ymin,
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * (this->limits.zmax - this->limits.zmin) + this->limits.zmin
    );
    vec3 normal = vec3(
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX),
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX),
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX)
    );
    return Plane(normal, pose);
}

Plane RRTZ::sample_plane(vec3 point) {
    vec3 normal = vec3(
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX),
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX),
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX)
    );
    return Plane(normal, point);
}

bool RRTZ::extend(Plane p, float& cost) {
    // return true if the tree is successfully extended
    // populate cost with cost to extended-to node

    // generate a node from start to the plane point
    if (!this->collision(this->start, p.point)) {
        Node3D start_node = Node3D(this->start, p.point, this->tree_nodes.size(), 0, 0);
        this->add_node_to_tree(start_node, cost);
        return true;
    }

    // try connecting plane with other nodes
    std::vector<Node3D*> near_nodes = std::vector<Node3D*>();
    std::vector<vec3> intersection_points = std::vector<vec3>();

    // find lowest cost candidate node
    Node3D *best_node_ptr = nullptr;
    vec3 *best_int_point_ptr = nullptr;
    float best_cost = std::numeric_limits<float>::infinity();

    if (this->get_near_nodes(p, near_nodes, intersection_points)) {
        for (size_t i = 0; i < near_nodes.size(); i++) {
            if (!this->collision(intersection_points[i], p.point)) {
                float node_cost = (*near_nodes[i]).cost + heading_change(*near_nodes[i], p.point-intersection_points[i]);
                if (node_cost < best_cost) {
                    best_cost = cost;
                    best_node_ptr = near_nodes[i];
                    best_int_point_ptr = &intersection_points[i];
                }
            }
        }
    }

    // if a collisionfree candidate node was found, add it to the tree
    if (best_node_ptr != nullptr) {
        Node3D new_node = Node3D(*best_int_point_ptr, p.point, this->tree_nodes.size(), (*best_node_ptr).idx);
        this->add_node_to_tree(new_node, cost);
        return true;
    }
    return false;
}

void RRTZ::print_node(Node3D node) {
    std::cout << " origin=" << node.origin.toString() << " end=" 
    << node.end.toString() << " vector=" << node.vector.toString() 
    << " cost=" << node.cost << " idx=" << node.idx << " parent idx="
    << node.parent_idx << std::endl;
}

bool RRTZ::get_near_nodes(Plane plane, std::vector<Node3D*>& near_nodes, std::vector<vec3>& int_points) {
    for (size_t i = 0; i < this->tree_nodes.size(); i++) {
        vec3 intPoint;
        if (ray_int_plane(this->tree_nodes[i], plane, 1e-9f, intPoint)) {
            near_nodes.push_back(&(this->tree_nodes[i]));
            int_points.push_back(intPoint);
        }
    }
    return !near_nodes.empty();
}

void RRTZ::add_node_to_tree(Node3D& node, float& cost) {
    node.cost = this->tree_nodes[node.parent_idx].cost + heading_change(this->tree_nodes[node.parent_idx], node);
    cost = node.cost;
    this->tree_nodes.push_back(node);
}

bool RRTZ::check_bounds(vec3 point) {
    return point.x >= this->limits.xmin && point.x <= this->limits.xmax &&
           point.y >= this->limits.ymin && point.y <= this->limits.ymax &&
           point.z >= this->limits.zmin && point.z <= this->limits.zmax;
}

bool RRTZ::collision(vec3 origin, vec3 end) {
    for (size_t i = 0; i < this->obs.size(); i++) {
        if (this->obs[i].collision(origin, end)) {
            return true;
        }
    }
    return false;
}