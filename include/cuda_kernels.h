#ifndef CUDA_KERNELS_H
#define CUDA_KERNELS_H

#include "vec3_struct.hpp"
#include "viewpoint_struct.hpp"
#include "triangle_struct.hpp"

#include <vector>

// given a viewpoint: determine visible (can see all three vertices) triangles -- true means not visible (collision during raycast)
// populate int_points with intersection points between rays cast to each vertex of triangle in question and other triangles
extern "C" void cuda_kernel_coverage(
    const Viewpoint& viewpoint,
    const std::vector<Triangle*>& triangles,
    std::vector<bool>& collisions, // empty vector to be filled with true/false for each triangle
    vec3** int_points
);

// given a list of viewpoints mapped to triangles, determine if each viewpoint can see its respective triangle
extern "C" void cuda_kernel_many(
    const std::vector<Viewpoint>& viewpoints,
    const std::vector<size_t>& triangle_indices,
    const std::vector<Triangle*>& triangles,
    bool* collisions,
    vec3** int_points
);

extern "C" void cuda_kernel_inc_angle(
    const std::vector<Viewpoint>& viewpoint,
    const std::vector<Triangle*>& triangles,
    std::vector<float>& inc_angles // populate with incidence angles for each triangle
);

#endif // CUDA_KERNELS_H