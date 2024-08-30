#include "cost_matrix.hpp"
#include "cuda_kernels.h"
#include "triangle_struct.hpp"
#include "tsp.hpp"
#include "utils.hpp"
#include "vec3_struct.hpp"
#include "viewpoint_generator.hpp"

#include <iostream>
#include <limits>
#include <vector>

int main(int argc, char** argv) {
    // argc is the number of arguments passed to the program
    // argv is an array of strings containing the arguments

    // check for correct number of arguments
    if (argc != 1 && argc != 2) {
        std::cout << "Usage: ./run_ccpp or ./run_ccpp vgd" << std::endl;
        return 1;
    }

    float vgd = 4.0f;
    if (argc == 2) {
        std::cout << "Running with VGD=" << argv[1] << std::endl;
        vgd = std::stof(argv[1]);
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

    ViewpointGenerator vg(obsVec, vgd);
    std::string save_file = "station_remeshed_coverage_" + std::to_string(static_cast<int>(vgd)) + "m.csv";

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
    // saveCSV("../data/coverage_viewpoint_sets/coverage_" + std::to_string(static_cast<int>(vgd)) + "m_vp_set.csv", viewpoint_data_save);

    // save bool vector of final coverage of mesh faces
    std::vector<bool> final_coverage;
    vg.getFilteredCoverage(final_coverage);
    vg.missedCoverage();
    std::vector<std::vector<float>> final_coverage_data;
    for (size_t i = 0; i < final_coverage.size(); i++) {
        std::vector<float> coverage_data;
        coverage_data.push_back(final_coverage[i]);
        final_coverage_data.push_back(coverage_data);
    }
    // saveCSV("../data/coverage_viewpoint_sets/coverage_" + std::to_string(static_cast<int>(vgd)) + "m_coverage.csv", final_coverage_data);
    std::cout << "NUM VIEWPOINTS: " << coverage_viewpoints.size() << std::endl;

    // std::cout << "\n-----------------Coverage Viewpoints-----------------\n";
    // for (size_t i = 0; i < coverage_viewpoints.size(); i++) {
    //     std::cout << "Viewpoint " << i << ": pose=" << coverage_viewpoints[i].pose.toString() << " viewdir=" << coverage_viewpoints[i].viewdir.toString() << std::endl;
    // }

    // // Compute cost matrix for travelling salesman problem for all viewpoints
    // vec3 start = vec3(1.8f, 4.7f, 2.7f);
    // size_t rrtz_iter = 2000;
    // CostMatrix cm(rrtz_iter);
    // std::cout << "loading viewpoints" << std::endl;
    // cm.loadViewpoints(
    //     "../data/coverage_viewpoint_sets/coverage_" + std::to_string(static_cast<int>(vgd)) + "m_vp_set.csv",
    //     start
    // );
    // // std::cout << "generating paths" << std::endl;
    // // cm.generatePathMatrix();
    // // std::cout << "saving path matrix" << std::endl;
    // // cm.savePathMatrix("../data/tsp/" + std::to_string(static_cast<int>(vgd)) + "m_path_matrix.csv");
    // std::cout << "loading path matrix" << std::endl;
    // cm.loadPathMatrix("../data/tsp/" + std::to_string(static_cast<int>(vgd)) + "m_path_matrix.csv");
    // // std::cout << "saving simple cost matrix" << std::endl;
    // // cm.saveSimpleCostMatrix("../data/tsp/" + std::to_string(static_cast<int>(vgd)) + "m_simple_cost_matrix.csv");
    // std::cout << "loading simple cost matrix" << std::endl;
    // cm.loadSimpleCostMatrix("../data/tsp/" + std::to_string(static_cast<int>(vgd)) + "m_simple_cost_matrix.csv");
    // // std::cout << "generating cost matrix" << std::endl;
    // // cm.generateCostMatrix();
    // // std::cout << "saving cost matrix" << std::endl;
    // // cm.saveCostMatrix("../data/tsp/" + std::to_string(static_cast<int>(vgd)) + "m_cost_matrix.csv");
    // std::cout << "loading cost matrix" << std::endl;
    // cm.loadCostMatrix("../data/tsp/" + std::to_string(static_cast<int>(vgd)) + "m_cost_matrix.csv");

    // // TSP tsp(cm);
    // std::cout << "creating TSP object" << std::endl;
    // TSP tsp(cm);
    // tsp.greedyInit();
    // tsp.twoOpt();
    
    // // get path
    // std::vector<std::vector<float>> path;
    // tsp.getPath(path);

    // // view path
    // for (size_t i = 0; i < path.size(); i++) {
    //     std::cout << "Viewpoint " << i << ": ";
    //     for (size_t j = 0; j < path[i].size(); j++) {
    //         std::cout << path[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // // save path
    // saveCSV("../data/ordered_viewpoints/" + std::to_string(static_cast<int>(vgd)) + "m.csv", path);


    return 0;
}