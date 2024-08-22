#ifndef CUDA_KERNELS_H
#define CUDA_KERNELS_H

#include "vec3_struct.hpp"
#include "viewpoint_struct.hpp"
#include "triangle_struct.hpp"

#include <vector>

extern "C" void cuda_kernel_intersect_triangles(
    const Viewpoint& viewpoint,
    const std::vector<Triangle*>& triangles,
    bool*** collisions,
    vec3*** int_points
);

#endif // CUDA_KERNELS_H