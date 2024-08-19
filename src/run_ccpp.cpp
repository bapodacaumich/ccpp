#include "viewpoint_generator.hpp"
#include "utils.hpp"
#include "vec3_struct.hpp"

#include <iostream>
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
    float vgd = 2.0f;
    ViewpointGenerator vg(obsVec, vgd);
    std::cout << "Running Greedy.." << std::endl;
    std::vector<Viewpoint> coverage_viewpoints = vg.getCoverageViewpoints();
    for (size_t i = 0; i < coverage_viewpoints.size(); i++) {
        std::cout << "Viewpoint " << i << ": pose=" << coverage_viewpoints[i].pose.toString() << " viewdir=" << coverage_viewpoints[i].viewdir.toString() << std::endl;
    }

    return 0;
}