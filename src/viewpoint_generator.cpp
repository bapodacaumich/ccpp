#include "cone_camera.hpp"
#include "cuda_kernels.h"
#include "obs.hpp"
#include "utils.hpp"
#include "viewpoint_generator.hpp"
#include "viewpoint_struct.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

ViewpointGenerator::ViewpointGenerator() {
    // default constructor -- useless to do this
    this->structure = std::vector<OBS>();
    this->vgd = 2.0f;
    this->cam = ConeCamera();
    this->structure = std::vector<OBS>();
    this->unfiltered_viewpoints = std::vector<Viewpoint>();
    this->coverage_viewpoints = std::vector<Viewpoint>();

    std::cout << "Warning: default constructor for ViewpointGenerator assumes no inspectable object and is useless" << std::endl;
}
ViewpointGenerator::ViewpointGenerator(
    std::vector<OBS> structure
    ) {

    // initialize members
    this->structure = structure;
    this->vgd = 2.0f;
    this->cam = ConeCamera();
    this->unfiltered_viewpoints = std::vector<Viewpoint>();
    this->coverage_viewpoints = std::vector<Viewpoint>();

    // initialize viewpoint generator
    this->initialize();
}
ViewpointGenerator::ViewpointGenerator(
    std::vector<OBS> structure,
    float vgd
    ) {

    // initialize members
    this->structure = structure;
    this->vgd = vgd;
    this->cam = ConeCamera();
    this->unfiltered_viewpoints = std::vector<Viewpoint>();
    this->coverage_viewpoints = std::vector<Viewpoint>();

    // initialize viewpoint generator
    this->initialize();
}

ViewpointGenerator::ViewpointGenerator(
    std::vector<OBS> structure,
    float vgd,
    ConeCamera cam
    ) {

    // initialize members
    this->structure = structure;
    this->vgd = vgd;
    this->cam = cam;
    this->unfiltered_viewpoints = std::vector<Viewpoint>();
    this->coverage_viewpoints = std::vector<Viewpoint>();

    // initialize viewpoint generator
    this->initialize();
}

void ViewpointGenerator::initialize() {

    // count number of mesh faces
    std::cout << "Initializing ViewpointGenerator..." << std::endl;
    this->countMeshFaces();

    // then populate unfiltered_viewpoints
    std::cout << "Populating Viewpoints..." << std::endl;
    if (!this->populateViewpoints()) { return; }; // return if no viewpoints were added

    // once we know how many viable viewpoints there are, we can initialize the coverage map
    std::cout << "Allocating Coverage Map..." << std::endl;
    this->coverage_map = new bool*[this->unfiltered_viewpoints.size()];
    for (size_t i = 0; i < this->unfiltered_viewpoints.size(); i++) {
        this->coverage_map[i] = new bool[this->num_mesh_faces];
        for (size_t j = 0; j < this->num_mesh_faces; j++) {
            this->coverage_map[i][j] = false;
        }
    }

    // populate unfiltered viewpoint x mesh face coverage matrix
    std::cout << "Populating Coverage Map..." << std::endl;
    this->populateCoverage();
}

std::vector<Viewpoint> ViewpointGenerator::getCoverageViewpoints() {
    // run greedy algorithm to select viewpoints and put in coverage_viewpoints
    this->greedy();
    return this->coverage_viewpoints;
}

void ViewpointGenerator::greedy() {
    /*
    * greedy algorithm to select viewpoints from unfiltered_viewpoints
    * and put into filtered_viewpoints based on coverage
    */
    // initialize marginal gain and vector of ordered pointers for each viewpoint
    int* marginal_gain = new int[this->unfiltered_viewpoints.size()];
    std::vector<int*> sorted_marginal_gain;
    for (size_t i = 0; i < this->unfiltered_viewpoints.size(); i++) {
        marginal_gain[i] = std::numeric_limits<int>::max();
        sorted_marginal_gain.push_back(&marginal_gain[i]);
    }

    // initialize coverage map for filtered viewpoints
    bool* filtered_coverage = new bool[this->num_mesh_faces];
    for (size_t i = 0; i < this->num_mesh_faces; i++) {
        filtered_coverage[i] = false; // initialize all elements to false
    }

    // iterate over all viewpoints (upper bound for number of viewpoints we need to check)
    std::cout << std::endl;
    // for (size_t i = 0; i < this->unfiltered_viewpoints.size(); i++) {
    for (size_t i = 0; i < 5; i++) {
        // update marginal gains and decending pointers
        this->updateMarginalGainPtrs(filtered_coverage, sorted_marginal_gain);

        std::cout << marginal_gain[0] << std::endl;
        std::cout << marginal_gain[1] << std::endl;
        std::cout << marginal_gain[2] << std::endl;
        std::cout << marginal_gain[3] << std::endl;

        // get viewpoint index of maximal marginal gain and add to unfiltered viewpoints
        size_t idx_to_add = sorted_marginal_gain[0] - marginal_gain;
        coverage_viewpoints.push_back(this->unfiltered_viewpoints[idx_to_add]);

        // remove first (maximal) element of sorted_marginal_gain because 
        //  we added the viewpoint it is associated with to coverage viewpoints
        if (!sorted_marginal_gain.empty()) {
            sorted_marginal_gain.erase(sorted_marginal_gain.begin());
        } else {
            std::cout << "Error: sorted_marginal_gain is empty (should not happen)" << std::endl;
        }

        // check if we have covered all faces
        if (allTrue(filtered_coverage, this->num_mesh_faces)) { break; } // from utils
        else {
            size_t num_true = 0;
            numTrue(filtered_coverage, this->num_mesh_faces, num_true); // from utils
            std::cout << "Faces Covered: " << num_true << " / " << this->num_mesh_faces << std::endl;
        }
    }

    delete marginal_gain;
    delete filtered_coverage;
}

void ViewpointGenerator::updateMarginalGainPtrs(
    const bool* filtered_coverage,         // coverage map for filtered viewpoints (already covered)
    std::vector<int*> sorted_marginal_gain  // array of pointers to marginal gain sorted by decending marginal gain (not const so will be updated)
    ) {

    // initialize the pointer that indicates the last element that needs to be sorted in the sorted marginal gain array
    auto sort_end = sorted_marginal_gain.begin();

    // update marginal gain as we search for the maximal element
    for (size_t sm_idx = 0; sm_idx < sorted_marginal_gain.size(); sm_idx++) {
        int gain = 0;
        for (size_t face_idx = 0; face_idx < this->num_mesh_faces; face_idx++) {
            // if the viewpoint covers the face and the face is not already covered, increment gain
            if (!filtered_coverage[sm_idx] && this->coverage_map[sm_idx][face_idx]) { gain++; }
        }
        // update the marginal gain array (marginal_gain)
        *sorted_marginal_gain[sm_idx] = gain;

        // increment sort_end
        sort_end++;

        // check if this is the largest element by comparing it to the next element (except if on last element)
        if (sm_idx != sorted_marginal_gain.size() - 1 && gain > *sorted_marginal_gain[sm_idx+1]) {
            break;
        }

    }

    // sort the subset of the marginal gain array that still needs to be sorted
    std::sort(
        sorted_marginal_gain.begin(), 
        sort_end,
        [](int* a, int* b) -> bool { 
            return *a > *b;
        }
    );
}

void ViewpointGenerator::populateCoverage() {

    // populate the filtered or unfiltered coverage map
    for (size_t vp_idx=0; vp_idx < this->unfiltered_viewpoints.size(); vp_idx++) {
        std::cout << "\rPopulating Coverage: " << vp_idx << "/" << this->unfiltered_viewpoints.size() << " viewpoints";
        for (size_t face_idx=0; face_idx < this->num_mesh_faces; face_idx++) {

            // check if viewpoint covers face
            std::vector<vec3*> points_ptr = {
                &(this->all_faces[face_idx]->a),
                &(this->all_faces[face_idx]->b),
                &(this->all_faces[face_idx]->c),
                };

            // if viewpoint sees all points of the triangle, then it covers the triangle
            if (this->collision(unfiltered_viewpoints[vp_idx].pose, points_ptr)) {
                this->coverage_map[vp_idx][face_idx] = true;
            } else {
                this->coverage_map[vp_idx][face_idx] = false;
            }
        }
    }
}

void ViewpointGenerator::countMeshFaces() {
    // count the mesh faces and populate all_faces (vector of pointers to each face)
    this->num_mesh_faces = 0;
    for (size_t obs_idx = 0; obs_idx < this->structure.size(); obs_idx++) {
        this->num_mesh_faces += this->structure[obs_idx].faces.size();
        for (size_t face_idx = 0; face_idx < this->structure[obs_idx].faces.size(); face_idx++) {
            this->all_faces.push_back(&(this->structure[obs_idx].faces[face_idx]));
        }
    }
}

bool ViewpointGenerator::populateViewpoints() {
    // size_t num_successful_viewpoints = 0;
    // size_t num_failed_viewpoints = 0;
    // for (size_t obs_idx = 0; obs_idx < this->structure.size(); obs_idx++) {
    //     for (size_t face_idx = 0; face_idx < this->structure[obs_idx].faces.size(); face_idx++) {
    //         vec3 centroid = this->structure[obs_idx].faces[face_idx].getCentroid();
    //         vec3 normal = this->structure[obs_idx].faces[face_idx].n;
    //         Viewpoint vp = Viewpoint(
    //             centroid + normal * this->vgd, // TODO: dynamic vgd for obstacle avoidance
    //             normal * -1.0f
    //             );
    //         if (!this->collision(vp)) { // make sure viewpoint is not inside an obstacle
    //             std::vector<vec3*> points_ptr = {
    //                 &(this->structure[obs_idx].faces[face_idx].a),
    //                 &(this->structure[obs_idx].faces[face_idx].b),
    //                 &(this->structure[obs_idx].faces[face_idx].c),
    //                 };
    //             if (this->collision(vp.pose, points_ptr)) { // make sure viewpoint can see triangle
    //                 num_successful_viewpoints++;
    //                 this->unfiltered_viewpoints.push_back(vp);
    //                 std::cout << "Successfully added viewpoint at " << vp.pose.toString() << std::endl;
    //             } else {
    //                 num_failed_viewpoints++;
    //                 std::cout << "Viewpoint at " << vp.pose.toString() << " cannot see own triangle" << std::endl;
    //             }
    //         } else { // generate new viewpoint somehow --> using intersection?
    //             num_failed_viewpoints++;
    //             std::cout << "Collision detected for viewpoint at " << vp.pose.toString() << std::endl;
    //         }
    //     }
    // }
    // std::cout << "Successfully added " << num_successful_viewpoints << "/" << num_failed_viewpoints + num_successful_viewpoints << " viewpoints" << std::endl;
    // return num_successful_viewpoints > 0;

    size_t num_successful_viewpoints = 0;
    size_t num_failed_viewpoints = 0;
    std::vector<Viewpoint> sampled_viewpoints = std::vector<Viewpoint>();
    for (size_t obs_idx = 0; obs_idx < this->structure.size(); obs_idx++) {
        for (size_t face_idx = 0; face_idx < this->structure[obs_idx].faces.size(); face_idx++) {
            vec3 centroid = this->structure[obs_idx].faces[face_idx].getCentroid();
            vec3 normal = this->structure[obs_idx].faces[face_idx].n;
            Viewpoint vp = Viewpoint(
                centroid + normal * this->vgd, // TODO: dynamic vgd for obstacle avoidance
                normal * -1.0f
                );
            sampled_viewpoints.push_back(vp);
        }
    }

    bool** collisions = new bool*[sampled_viewpoints.size()];
    for (size_t i = 0; i < sampled_viewpoints.size(); i++) {
        collisions[i] = new bool[this->num_mesh_faces];
        for (size_t j = 0; j < this->num_mesh_faces; j++) {
            collisions[i][j] = false;
        }
    }

    // cuda_kernel_intersect_triangles(
    //     sampled_viewpoints,
    //     this->all_faces,
    //     collisions
    // );

    std::cout << "Collision Matrix:" << std::endl;
    for (size_t i = 0; i < sampled_viewpoints.size(); i++) {
        std::cout << sampled_viewpoints[i].pose.toString() << ": ";
        for (size_t j = 0; j < this->num_mesh_faces; j++) {
            std::cout << collisions[i][j] << " ";
        }
        std::cout << std::endl;
    }

    for (size_t i = 0; i < sampled_viewpoints.size(); i++) {
        delete[] collisions[i];
    }
    delete[] collisions;


    std::cout << "Successfully added " << num_successful_viewpoints << "/" << num_failed_viewpoints + num_successful_viewpoints << " viewpoints" << std::endl;
    return num_successful_viewpoints > 0;
}

bool ViewpointGenerator::collision(Viewpoint vp) {
    for (size_t obs_idx = 0; obs_idx < this->structure.size(); obs_idx++) {
        if (this->structure[obs_idx].collision(vp.pose)) { return true; }
    }
    return false;
}

bool ViewpointGenerator::collision(vec3 pose, const std::vector<vec3*>& points) {
    // must displace centroid in direction of normal by small amount to avoid self-collision
    for (size_t obs_idx = 0; obs_idx < this->structure.size(); obs_idx++) {
        for (size_t i=0; i < points.size(); i++) {
            if (this->structure[obs_idx].collision(pose, *points[i] + (pose - *points[i]) * 1e-9)) { return true; }
        }
    }
    return false;
}

ViewpointGenerator::~ViewpointGenerator() {
    for (size_t i = 0; i < this->unfiltered_viewpoints.size(); i++) {
        delete[] this->coverage_map[i];
    }
    delete[] this->coverage_map;
}