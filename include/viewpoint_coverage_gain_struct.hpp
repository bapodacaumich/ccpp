#ifndef VIEWPOINT_COVERAGE_GAIN_STRUCT_HPP
#define VIEWPOINT_COVERAGE_GAIN_STRUCT_HPP

#include "viewpoint_struct.hpp"
#include <vector>

struct VP_Coverage_Gain {
    Viewpoint vp;
    int gain;
    std::vector<bool> coverage;
    size_t vp_map_idx;
    float inc_angle;
};

#endif // VIEWPOINT_COVERAGE_GAIN_STRUCT_HPP