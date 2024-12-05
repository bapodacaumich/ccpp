#include "cuda_kernels.h"

#include "utils.hpp"
#include "vec3_struct.hpp"
#include "viewpoint_struct.hpp"

#include <iostream>
#include <string>
#include <vector>

/* 
    * This file is used to evaluate the path saturation of the different pathfinding algorithms.
    * path saturation is the measure of 'how many frames' each face gets given a camera specification such as framerate
    * path saturation also measures the 'quality' of each frame by evaluating the incidence angle between the camera and the face
*/

int main(int argc, char** argv) {
    // build path data based on framerate

    if (argc < 2) {
        std::cout << "Usage: ./evaluate_path_saturation <fps>" << std::endl;
        return 1;
    }

    float fps = std::stof(argv[1]);
    float dt = 1.0f / fps;

    std::vector<std::string> folders = {
        "cw_opt_packaged_10",
        "cw_opt_packaged_50",
        "cw_opt_packaged_var",
        "pf_final_ko",
        "pf_final_so"
    };

    std::vector<std::string> savefolders = {
        "ivt_10",
        "ivt_50",
        "ivt_var",
        "ocp_ko",
        "ocp_so"
    };

    std::vector<std::string> conditions = {
        "2m_global",
        "2m_local",
        "4m_global",
        "4m_local",
        "8m_global",
        "8m_local",
        "16m_global",
        "16m_local"
    };

    std::string dir = "../knot_ocp/";

    std::string savedir = "../visualization_python/saturation/";

    // iterate through each folder and each condition to evaluate path saturation
    for (size_t f = 0; f < folders.size(); f++) {
        std::string folder = folders[f];
        std::string savefolder = savefolders[f];
        for (std::string condition : conditions) {
            std::cout << "Evaluating path saturation for " << folder << " " << condition << std::endl;

            // first load path
            std::string pathfile = dir + folder + "/" + condition + ".csv";
            std::vector<std::vector<float>> path_data;
            loadCSV(pathfile, path_data, 7);

            // TODO: interpolate path to desired framerate
            std::vector<std::vector<float>> path_data_fps;
            float t = 0;
            for (size_t after_idx = 1; after_idx < path_data.size(); after_idx++) {
                size_t prev_idx = after_idx - 1;
                while (t < path_data[after_idx][6]) {
                    float alpha = (t - path_data[prev_idx][6]) / (path_data[after_idx][6] - path_data[prev_idx][6]);
                    vec3 viewdir = slerp(
                        vec3(path_data[prev_idx][3], path_data[prev_idx][4], path_data[prev_idx][5]),
                        vec3(path_data[after_idx][3], path_data[after_idx][4], path_data[after_idx][5]),
                        alpha
                    );
                    std::vector<float> path_step = {
                        lerp(path_data[prev_idx][0], path_data[after_idx][0], alpha),
                        lerp(path_data[prev_idx][1], path_data[after_idx][1], alpha),
                        lerp(path_data[prev_idx][2], path_data[after_idx][2], alpha),
                        viewdir.x,
                        viewdir.y,
                        viewdir.z,
                        t
                    };
                    path_data_fps.push_back(path_step);
                    t += dt;
                }
            }

            // compute saturation
            std::vector<std::vector<float>> saturation_map;
            compute_saturation_path(path_data_fps, saturation_map);

            // save saturation
            std::string savefile = savedir + savefolder + "/" + condition + "_sat.csv";
            saveCSV(savefile, saturation_map); // saturation = {count, avg_incidence_angle, min_incidence_angle}
            std::cout << std::endl;
        }
    }

}