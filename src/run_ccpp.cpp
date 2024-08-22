#include "viewpoint_generator.hpp"
#include "utils.hpp"
#include "vec3_struct.hpp"
#include "triangle_struct.hpp"
#include "cuda_kernels.h"

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

    // TESTING CUDA KERNELS:

    std::vector<Viewpoint> viewpoints;
    std::vector<Triangle*> triangles;

    Triangle tri0 = Triangle(vec3(0.0f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f));
    Triangle tri1 = Triangle(vec3(0.0f, 0.0f, 1.0f), vec3(1.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 1.0f));
    Triangle tri2 = Triangle(vec3(0.0f, 0.0f, 2.0f), vec3(1.0f, 0.0f, 2.0f), vec3(0.0f, 1.0f, 2.0f));

    triangles.push_back(&tri0);
    triangles.push_back(&tri1);
    triangles.push_back(&tri2);

    Viewpoint vp0 = Viewpoint(vec3(0.25f, 0.25f, 3.0f), vec3(0.0f, 0.0f, -1.0f));
    Viewpoint vp1 = Viewpoint(vec3(0.25f, 0.25f, 1.5f), vec3(0.0f, 0.0f, -1.0f));
    Viewpoint vp2 = Viewpoint(vec3(0.25f, 0.25f, 0.5f), vec3(0.0f, 0.0f, -1.0f));

    viewpoints.push_back(vp0);
    viewpoints.push_back(vp1);
    viewpoints.push_back(vp2);

    vec3 ***int_points = new vec3**[triangles.size() * triangles.size() * 3];
    bool ***collisions = new bool**[triangles.size() * triangles.size() * 3];
    for (size_t i = 0; i < triangles.size(); i++) {
        int_points[i] = new vec3*[triangles.size()];
        collisions[i] = new bool*[triangles.size()];
        for (size_t j = 0; j < triangles.size(); j++) {
            int_points[i][j] = new vec3[3];
            collisions[i][j] = new bool[3];
            for (size_t k = 0; k < 3; k++) {
                collisions[i][j][k] = false;
                int_points[i][j][k].set(0.0f, 0.0f, 0.0f);
            }
        }
    }
    cuda_kernel_intersect_triangles(vp0, triangles, collisions, int_points);

    // for (size_t i = 0; i < viewpoints.size(); i++) {
    //     for (size_t j = 0; j < triangles.size(); j++) {
    //         for (size_t k = 0; k < 3; k++) {
    //             std::cout << "Viewpoint =" << viewpoints[i].pose.toString() << " | Vertex=";

    //             if (k == 0) {
    //                 std::cout << triangles[j]->a.toString();
    //             } else if (k == 1) {
    //                 std::cout << triangles[j]->b.toString();
    //             } else {
    //                 std::cout << triangles[j]->c.toString();
    //             }

    //             std::cout << " | Collision=" << collisions[i][j][k];
    //             std::cout << " | IntPoint=" << int_points[i][j][k].toString() << std::endl;
    //         }
    //     }
    // }

    for (size_t i = 0; i < viewpoints.size(); i++) {
        for (size_t j = 0; j < triangles.size(); j++) {
            delete[] int_points[i][j];
            delete[] collisions[i][j];
        }
        delete[] int_points[i];
        delete[] collisions[i];
    }

    delete[] int_points;
    delete[] collisions;



    // // std::vector<OBS> obsVec;
    // // loadStationOBS(obsVec);
    // std::vector<OBS> obsVec;
    // // load in cube data and convert to triangles
    // std::vector<std::vector<std::vector<float>>> cubeData;
    // loadCube(cubeData);
    // std::vector<Triangle> triCubeFaces;
    // vecToTri(cubeData, triCubeFaces);

    // // create obstacle
    // OBS obs = OBS(triCubeFaces);
    // obsVec.push_back(obs);

    // float vgd = 5.0f;
    // ViewpointGenerator vg(obsVec, vgd);



    // std::cout << "Running Greedy.." << std::endl;
    // std::vector<Viewpoint> coverage_viewpoints = vg.getCoverageViewpoints();
    // for (size_t i = 0; i < coverage_viewpoints.size(); i++) {
    //     std::cout << "Viewpoint " << i << ": pose=" << coverage_viewpoints[i].pose.toString() << " viewdir=" << coverage_viewpoints[i].viewdir.toString() << std::endl;
    // }

    return 0;
}