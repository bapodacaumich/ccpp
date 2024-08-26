#include "viewpoint_generator.hpp"
#include "utils.hpp"
#include "vec3_struct.hpp"
#include "triangle_struct.hpp"
#include "cuda_kernels.h"

#include <iostream>
#include <limits>
#include <vector>

int main(int argc, char** argv) {
    // argc is the number of arguments passed to the program
    // argv is an array of strings containing the arguments

    // check for correct number of arguments
    if (argc != 1 && argc != 2) {
        std::cout << "Usage: ./rrtz or ./rrtz max_nodes" << std::endl;
        return 1;
    }

    size_t max_nodes = 500;
    if (argc == 2) {
        std::cout << "Running with max nodes=" << argv[1] << std::endl;
        max_nodes = std::stoi(argv[1]);
    }

    std::vector<OBS> obsVec;
    loadStationOBS(obsVec);
    // std::vector<OBS> obsVec;
    // // load in cube data and convert to triangles
    // std::vector<std::vector<std::vector<float>>> cubeData;
    // loadCube(cubeData);
    // std::vector<Triangle> triCubeFaces;
    // vecToTri(cubeData, triCubeFaces);

    // // create obstacle
    // OBS obs = OBS(triCubeFaces);
    // obsVec.push_back(obs);

    float vgd = 5.0f;
    ViewpointGenerator vg(obsVec, vgd);
    std::string save_file = "station.csv";

    // vg.populateCoverage();
    // vg.saveCoverageMap(save_file);

    // save coverage map
    vg.loadCoverageMap(save_file);

    std::cout << "Running Greedy.." << std::endl;
    std::vector<Viewpoint> coverage_viewpoints = vg.getCoverageViewpoints();
    std::vector<std::vector<float>> viewpoint_data_save;
    for (size_t i = 0; i < coverage_viewpoints.size(); i++) {
        std::vector<float> vp_data;
        vp_data.push_back(coverage_viewpoints[i].pose.x);
        vp_data.push_back(coverage_viewpoints[i].pose.y);
        vp_data.push_back(coverage_viewpoints[i].pose.z);
        vp_data.push_back(coverage_viewpoints[i].viewdir.x);
        vp_data.push_back(coverage_viewpoints[i].viewdir.y);
        vp_data.push_back(coverage_viewpoints[i].viewdir.z);
        viewpoint_data_save.push_back(vp_data);
    }
    saveCSV("../data/coverage_viewpoint_sets/station_viewpoints_coverage.csv", viewpoint_data_save);
    for (size_t i = 0; i < coverage_viewpoints.size(); i++) {
        std::cout << "Viewpoint " << i << ": pose=" << coverage_viewpoints[i].pose.toString() << " viewdir=" << coverage_viewpoints[i].viewdir.toString() << std::endl;
    }

    return 0;
}