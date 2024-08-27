#include "cone_camera.hpp"
#include "cuda_kernels.h"
#include "obs.hpp"
#include "utils.hpp"
#include "viewpoint_coverage_gain_struct.hpp"
#include "viewpoint_generator.hpp"
#include "viewpoint_struct.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

ViewpointGenerator::ViewpointGenerator() {
    // default constructor -- useless to do this
    this->vgd = 2.0f;
    this->cam = ConeCamera();
    this->structure = std::vector<OBS>();
    this->unfiltered_viewpoints = std::vector<Viewpoint>();
    this->coverage_viewpoints = std::vector<Viewpoint>();
    this->coverage_map = std::vector<std::vector<bool>>();

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
    this->coverage_map = std::vector<std::vector<bool>>();

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
    this->coverage_map = std::vector<std::vector<bool>>();

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
    this->coverage_map = std::vector<std::vector<bool>>();

    // initialize viewpoint generator
    this->initialize();
}

void ViewpointGenerator::initialize() {

    // count number of mesh faces
    std::cout << "Initializing Viewpoint Generator..." << std::endl;
    this->countMeshFaces();

    // then populate unfiltered_viewpoints
    std::cout << "Populating Viewpoints..." << std::endl;
    if (!this->populateViewpoints()) { return; }; // return if no viewpoints were added

    // compute incidence angle between each viewpoint and each face
    cuda_kernel_inc_angle(this->unfiltered_viewpoints, this->all_faces, this->inc_angle_map);
}

void ViewpointGenerator::printIncidenceAngles() {
    // print viewpoints and faces with corresponding indices, then print incidence angles organized wrt indices
    for (size_t i=0; i < this->unfiltered_viewpoints.size(); i++) {
        std::cout << "Viewpoint " << i << ": " << this->unfiltered_viewpoints[i].pose.toString() << std::endl;
    }

    for (size_t i=0; i < this->num_mesh_faces; i++) {
        std::cout << "Face " << i << ": " << this->all_faces[i]->toString() << std::endl;
    }

    // print inc_angles
    std::cout << "Incidence Angles:" << std::endl;
    for (size_t vp_idx = 0; vp_idx < this->unfiltered_viewpoints.size(); vp_idx++) {
        std::cout << "Viewpoint " << vp_idx << ": "; 
        for (size_t face_idx = 0; face_idx < this->num_mesh_faces; face_idx++) {
            std::cout << this->inc_angle_map[vp_idx][face_idx] << ", ";
        }
        std::cout << std::endl;
    }
}

void ViewpointGenerator::printCoverageMap() {
    std::cout << "Coverage Map:" << std::endl;
    for (size_t i = 0; i < this->coverage_map.size(); i++) {
        std::cout << "Viewpoint " << i << ": ";
        for (size_t j = 0; j < this->coverage_map[i].size(); j++) {
            std::cout << this->coverage_map[i][j] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void ViewpointGenerator::saveCoverageMap(std::string& filename) {
    // save coverage map to file
    std::vector<std::vector<float>> data;
    for (size_t i = 0; i < this->coverage_map.size(); i++) {
        std::vector<float> row;
        for (size_t j = 0; j < this->coverage_map[i].size(); j++) {
            row.push_back(this->coverage_map[i][j] ? 1.0f : 0.0f);
        }
        data.push_back(row);
    }
    saveCSV("../data/coverage_maps/" + filename, data);
}

void ViewpointGenerator::loadCoverageMap(std::string& filename) {
    // assume right coverage_map - viewpoints - faces are loaded
    // load coverage map from file
    std::vector<std::vector<float>> data;
    loadCSV("../data/coverage_maps/" + filename, data, this->num_mesh_faces);
    for (size_t i = 0; i < data.size(); i++) {
        std::vector<bool> row;
        for (size_t j = 0; j < data[i].size(); j++) {
            row.push_back(data[i][j] == 1.0f);
        }
        this->coverage_map.push_back(row);
    }
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

    this->setUpCoverageGain();

    this->filtered_coverage = std::vector<bool>(this->num_mesh_faces, false);

    // iterate over number of viewpoints (most time we can add viewpoints to this->coverage_viewpoints)
    std::cout << "Number of unfiltered viewpoints=" << this->unfiltered_viewpoints.size() << std::endl;
    for (size_t i = 0; i < this->unfiltered_viewpoints.size(); i++) {
        // each iteration sort the viewpoint-gain objects
        this->sortUpdateMarginalGain();

        // get viewpoint with maximal gain
        this->coverage_viewpoints.push_back(this->vpcg_unfiltered.begin()->vp);
        this->vpcg_filtered.push_back(*(this->vpcg_unfiltered.begin()));

        // remove first element from vp_coverage_gains
        this->vpcg_unfiltered.erase(this->vpcg_unfiltered.begin());

        // update filtered coverage map
        this->updateCoverage();

        // check if we have covered all faces
        if (allTrue(this->filtered_coverage) || allZeroGain(this->vpcg_unfiltered)) { break; }
    }
    this->sortUpdateMarginalGain();
    size_t num_covered = 0;
    for (size_t i = 0; i < this->filtered_coverage.size(); i++) {
        if (this->filtered_coverage[i]) { num_covered++; }
    }
    std::cout << "Number of faces covered=" << num_covered << "/" << this->num_mesh_faces << " or " << filtered_coverage.size() << std::endl;
}

void ViewpointGenerator::sortUpdateMarginalGain() {

    // update marginal gain as we search for the maximal element
    for (size_t i = 0; i < this->vpcg_unfiltered.size(); i++) {
        int gain = 0;
        for (size_t face_idx = 0; face_idx < this->num_mesh_faces; face_idx++) {
            // if the viewpoint covers the face and the face is not already covered, increment gain
            if (!(this->filtered_coverage[face_idx]) && this->coverage_map[vpcg_unfiltered[i].vp_map_idx][face_idx]) { 
                // define gain function here
                // 1 is added to avoid division by zero and normalize best incidence angle to gain = 1
                // gain += 1/(inc_angle_map[vpcg_unfiltered[i].vp_map_idx][face_idx] + 1);
                gain++;
            }
        }
        // update the marginal gain array (marginal_gain)
        vpcg_unfiltered[i].gain = gain;

        // check if this is the largest element by comparing it to the next element (except if on last element)
        if (i < vpcg_unfiltered.size() - 1 && gain >= vpcg_unfiltered[i+1].gain) {
            break;
        }
    }

    // sort the subset of the marginal gain array that still needs to be sorted
    std::sort(
        vpcg_unfiltered.begin(), 
        vpcg_unfiltered.end(),
        [](VP_Coverage_Gain a, VP_Coverage_Gain b) -> bool { 
            return a.gain > b.gain;
        }
    );
}

void ViewpointGenerator::updateCoverage() {
    // update the filtered coverage map based on viewpoints added to this->coverage_viewpoints
    for (size_t i = 0; i < this->vpcg_filtered.size(); i++) {
        for (size_t j = 0; j < this->num_mesh_faces; j++) {
            if (this->coverage_map[vpcg_filtered[i].vp_map_idx][j]) {
                this->filtered_coverage[j] = true;
            }
        }
    }
}

void ViewpointGenerator::populateCoverage() {
    std::cout << "Populating Coverage Map..." << std::endl;
    // populate the filtered or unfiltered coverage map
    getCoverage(this->unfiltered_viewpoints, this->all_faces, this->coverage_map);
}


void ViewpointGenerator::setUpCoverageGain() {
    this->vpcg_filtered = std::vector<VP_Coverage_Gain>();
    this->vpcg_unfiltered = std::vector<VP_Coverage_Gain>();
    for (size_t i = 0; i < this->coverage_map.size(); i++) {
        VP_Coverage_Gain vpcg;
        vpcg.vp = this->unfiltered_viewpoints[i];
        vpcg.gain = std::numeric_limits<int>::max();
        vpcg.coverage = this->coverage_map[i];
        vpcg.vp_map_idx = i;
        // incidence angle in radians:
        getIncidenceAngle(this->unfiltered_viewpoints[i].viewdir, *(this->all_faces[i]), vpcg.inc_angle);
        this->vpcg_unfiltered.push_back(vpcg);
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

    size_t num_successful_viewpoints = 0;
    size_t num_failed_viewpoints = 0;
    std::vector<Viewpoint> sampled_viewpoints = std::vector<Viewpoint>();
    std::vector<size_t> sampled_face_indices = std::vector<size_t>();
    for (size_t face_idx = 0; face_idx < this->num_mesh_faces; face_idx++) {
        vec3 centroid = this->all_faces[face_idx]->getCentroid();
        vec3 normal = this->all_faces[face_idx]->n;
        Viewpoint vp = Viewpoint(
            centroid + normal * this->vgd, // TODO: dynamic vgd for obstacle avoidance
            normal * -1.0f
            );
        sampled_viewpoints.push_back(vp);
        sampled_face_indices.push_back(face_idx);

    }

    bool* collisions = new bool[sampled_viewpoints.size()];
    for (size_t i = 0; i < sampled_viewpoints.size(); i++) {
        collisions[i] = false;
    }

    // // debug int points
    // vec3 **int_points = new vec3*[sampled_viewpoints.size()];
    // for (size_t i = 0; i < sampled_viewpoints.size(); i++) {
    //     int_points[i] = new vec3[this->num_mesh_faces * 3];
    //     for (size_t j = 0; j < this->num_mesh_faces * 3; j++) {
    //         int_points[i][j].set(0.0f, 0.0f, 0.0f);
    //     }
    // }


    // check each viewpoint-face combination for ray casting collisions
    std::vector<bool> in_collision;
    cuda_kernel_collision_points(
        sampled_viewpoints,
        this->all_faces,
        vec3(-30.0f, -30.0f, -30.0f),
        in_collision
    );

    for (size_t i = 0; i < sampled_viewpoints.size(); i++) {
        if (!in_collision[i]) {
            // std::cout << "Viewpoint " << sampled_viewpoints[i].pose.toString() << " is not in collision" << std::endl;
            this->unfiltered_viewpoints.push_back(sampled_viewpoints[i]);
            num_successful_viewpoints++;
        } else {
            // std::cout << "Viewpoint " << i << " : "<< sampled_viewpoints[i].pose.toString() << " is in collision" << std::endl;
            num_failed_viewpoints++;
        }
    }

    // debug
    // std::cout << "Collision Matrix:" << std::endl;
    // for (size_t i = 0; i < sampled_viewpoints.size(); i++) {
    //     std::cout << "Viewpoint " << sampled_viewpoints[i].pose.toString();
    //     std::cout << " paired with " << this->all_faces[sampled_face_indices[i]]->toString() << ": ";
    //     if (collisions[i]) {
    //         std::cout << "collisions: ";
    //         for (size_t j = 0; j < this->num_mesh_faces * 3; j++) {
    //             if (int_points[i][j].x != std::numeric_limits<float>::infinity()) {
    //                 std::cout << this->all_faces[size_t(j/3)]->toString() << ", ";
    //                 std::cout << ", Intersection=" << int_points[i][j].toString() << ", ";
    //             }
    //         }
    //         std::cout << std::endl;
    //     } else {
    //         std::cout << "no collision" << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    // for (size_t i = 0; i < sampled_viewpoints.size(); i++) {
    //     delete[] int_points[i];
    // }
    // delete[] int_points;

    // delete[] collisions;


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