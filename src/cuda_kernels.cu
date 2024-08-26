#include "cuda_kernels.h"
#include "vec3_struct.hpp"
#include "triangle_struct.hpp"
#include "viewpoint_struct.hpp"

#include <vector>
#include <limits>
#include <cmath>

// one origin, mapped to many end points
__global__ void ray_int_plane(
    bool *result, // flattened 3d
    vec3 *int_points, // flattened 3d
    const vec3 origin,  // vp (vec3)
    const vec3 *ends,    // n_vp (1dim)
    const Triangle *tri,// n_tri (1dim)
    size_t n_tri
    ) {

    // epsilon for floating point comparison
    float eps = 1e-6f;

    // get indices
    size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_tri
    size_t tri_idx = blockIdx.y * blockDim.y + threadIdx.y; // n_tri
    size_t tri_pt_idx = blockIdx.z * blockDim.z + threadIdx.z; // 3
    size_t ray_idx = vp_idx * 3 + tri_pt_idx; // n_ray
    size_t res_idx = tri_pt_idx * n_tri * n_tri + tri_idx * n_tri + vp_idx;

    if (vp_idx > n_tri - 1 || tri_idx > n_tri - 1 || tri_pt_idx > 2) { return; }


    // instantiate ray
    vec3 end = ends[ray_idx];
    vec3 vec = end - origin;

    // look for any intersections between the ray and triangle
    vec3 e1 = tri[tri_idx].b - tri[tri_idx].a;
    vec3 e2 = tri[tri_idx].c - tri[tri_idx].a;
    vec3 h = vec.cross(e2);
    float a = e1.dot(h);

    // if ray is parallel to triangle
    if (a > -eps && a < eps) {
        result[res_idx] = false;
        return;
    }

    float f = 1 / a;
    vec3 s = origin - tri[tri_idx].a;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) {
        result[res_idx] = false;
        return;
    }
    vec3 q = s.cross(e1);
    float v = f * vec.dot(q);
    if (v < 0.0f || u + v > 1.0f) {
        result[res_idx] = false;
        return;
    }

    // find intersection point
    float t = f * e2.dot(q);
    vec3 intPoint = origin + vec * t;
    int_points[res_idx] = intPoint;

    // check if intersection point is between origin and end
    vec3 vec_dir = vec/vec.norm();
    if ((intPoint-origin).dot(vec_dir) < vec.norm() - eps && (intPoint-origin).dot(vec_dir) > 0) {
        result[res_idx] = true;
        return;
    }

    result[res_idx] = false;
    return;
}

// dims: viewpoints (x dim) x faces (y dim) x 3 (tri dim)
// many origins, each mapped to an end point
__global__ void ray_int_plane_many(
    bool *result, // flattened 3d
    vec3 *int_points, // flattened 3d
    const vec3 *starts,  // n_vp (vec3)
    const vec3 *ends,    // n_vp (1dim)
    size_t n_vp,
    const Triangle *tri,// n_tri (1dim)
    size_t n_tri
    ) {

    // epsilon for floating point comparison
    float eps = 1e-6f;

    // get indices
    size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_tri
    size_t tri_idx = blockIdx.y * blockDim.y + threadIdx.y; // n_tri
    size_t tri_pt_idx = blockIdx.z * blockDim.z + threadIdx.z; // 3
    size_t ray_idx = vp_idx * 3 + tri_pt_idx; // n_ray
    size_t res_idx = tri_pt_idx * n_tri * n_vp + tri_idx * n_vp + vp_idx;

    if (vp_idx > n_vp - 1 || tri_idx > n_tri - 1 || tri_pt_idx > 2) { return; }


    // instantiate ray
    vec3 origin = starts[ray_idx];
    vec3 end = ends[ray_idx];
    vec3 vec = end - origin;

    // look for any intersections between the ray and triangle
    vec3 e1 = tri[tri_idx].b - tri[tri_idx].a;
    vec3 e2 = tri[tri_idx].c - tri[tri_idx].a;
    vec3 h = vec.cross(e2);
    float a = e1.dot(h);

    // if ray is parallel to triangle
    if (a > -eps && a < eps) {
        result[res_idx] = false;
        return;
    }

    float f = 1 / a;
    vec3 s = origin - tri[tri_idx].a;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) {
        result[res_idx] = false;
        return;
    }
    vec3 q = s.cross(e1);
    float v = f * vec.dot(q);
    if (v < 0.0f || u + v > 1.0f) {
        result[res_idx] = false;
        return;
    }

    // find intersection point
    float t = f * e2.dot(q);
    vec3 intPoint = origin + vec * t;
    int_points[res_idx] = intPoint;

    // check if intersection point is between origin and end
    vec3 vec_dir = vec/vec.norm();
    if ((intPoint-origin).dot(vec_dir) < vec.norm() - eps && (intPoint-origin).dot(vec_dir) > 0) {
        result[res_idx] = true;
        return;
    }

    result[res_idx] = false;
    return;
}

__global__ void collision_or(bool* vp_collision, const bool* ray_tri_collision, size_t n_vp, size_t n_tri) {
    // for each viewpoint-triangle correspondance, check if rays to each vertex collide with any other triangle. if so, write in true
    // get viewpoint index
    size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_vp

    if (vp_idx > n_vp - 1) { return; } // n_tri = n_vp

    vp_collision[vp_idx] = false;
    for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
        for (size_t tri_pt_idx = 0; tri_pt_idx < 3; tri_pt_idx++) {
            size_t res_idx = tri_pt_idx * n_tri * n_vp + tri_idx * n_vp + vp_idx;
            if (ray_tri_collision[res_idx]) {
                vp_collision[vp_idx] = true;
            }
        }
    }
    return;
}

__global__ void inc_angle(
    float *angles, // flattened 2d
    const vec3 *poses, // n_vp (vec3)
    const vec3 *centroids, // n_tri (vec3)
    const vec3 *normals, // n_tri (vec3)
    size_t n_vp,
    size_t n_tri
    ) {

    // get indices
    size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_vp
    size_t tri_idx = blockIdx.y * blockDim.y + threadIdx.y; // n_tri
    size_t res_idx = tri_idx * n_vp + vp_idx; // n_tri * n_vp

    if (vp_idx > n_vp - 1 || tri_idx > n_tri - 1) { return; }

    // calculate angle
    vec3 vec = centroids[tri_idx] - poses[vp_idx];
    vec3 norm = normals[tri_idx];
    float angle = acos(fabsf(vec.dot(norm))/(vec.norm()*norm.norm()));
    angles[res_idx] = angle;
    return;
}

extern "C" void cuda_kernel_coverage(
    const Viewpoint& vp, 
    const std::vector<Triangle*>& faces,
    std::vector<bool>& collisions,
    vec3** int_points
    ) {

    // put viewpoints into arrays
    size_t n_tri = faces.size();
    size_t n_ray = n_tri * 3;

    vec3 *ends = new vec3[n_ray];
    Triangle *tri = new Triangle[n_tri];
    bool *result_arr = new bool[n_tri];
    bool *intersection_arr = new bool[n_ray * n_tri];
    vec3 *result_int_points = new vec3[n_ray * n_tri];

    // thread, block size
    size_t thread_x = 16;
    size_t thread_y = 16;
    size_t thread_z = 3;

    // put faces into array
    for (size_t i = 0; i < n_tri; i++) {
        tri[i] = *faces[i];
    }

    // put viewpoints into array
    for (size_t tri_pt_idx = 0; tri_pt_idx < 3; tri_pt_idx++) {
        for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
            size_t ray_idx = tri_idx * 3 + tri_pt_idx; // n_ray
            if (tri_pt_idx == 0) {
                ends[ray_idx] = tri[tri_idx].a;
            } else if (tri_pt_idx == 1) {
                ends[ray_idx] = tri[tri_idx].b;
            } else if (tri_pt_idx == 2) {
                ends[ray_idx] = tri[tri_idx].c;
            }
        }
    }

    // allocate gpu memory
    vec3 *d_ends;
    Triangle *d_tri;

    bool *d_intersections; // collisions per ray per triangle
    bool *d_result; // collisions per triangle
    vec3 *d_int_points; // intersection points

    cudaMalloc(&d_ends, n_ray * sizeof(vec3));
    cudaMalloc(&d_tri, n_tri * sizeof(Triangle));
    cudaMalloc(&d_result, n_tri * sizeof(bool));
    cudaMalloc(&d_intersections, n_ray * n_tri * sizeof(bool));
    cudaMalloc(&d_int_points, n_ray * n_tri * sizeof(vec3));

    // copy data to gpu
    cudaMemcpy(d_ends, ends, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tri, tri, n_tri * sizeof(Triangle), cudaMemcpyHostToDevice);

    // set thread and block size
    dim3 threadsPerBlock(thread_x, thread_y, thread_z);

    // 2D blocks because 3d blocks can account for the 3rd dim by themselves
    dim3 numBlocks(int((n_tri + threadsPerBlock.x - 1) / threadsPerBlock.x), (n_tri + threadsPerBlock.y - 1) / threadsPerBlock.y, 1);
    ray_int_plane<<<numBlocks, threadsPerBlock>>>(d_intersections, d_int_points, vp.pose, d_ends, d_tri, n_tri);

    cudaMemcpy(intersection_arr, d_intersections, n_ray * n_tri * sizeof(bool), cudaMemcpyDeviceToHost);

    // same as above but without the 3rd dimension
    threadsPerBlock.x = 1024;
    threadsPerBlock.y = 1;
    threadsPerBlock.z = 1;
    numBlocks.x = (n_tri + threadsPerBlock.x - 1) / threadsPerBlock.x;
    collision_or<<<numBlocks, threadsPerBlock>>>(d_result, d_intersections, n_tri, n_tri);

    cudaMemcpy(result_arr, d_result, n_tri * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(result_int_points, d_int_points, n_ray * n_tri * sizeof(vec3), cudaMemcpyDeviceToHost);

    cudaFree(d_ends);
    cudaFree(d_tri);
    cudaFree(d_result);
    cudaFree(d_intersections);
    cudaFree(d_int_points);

    for (size_t vp_idx=0; vp_idx < n_tri; vp_idx++) {
        collisions.push_back(result_arr[vp_idx]);
    }

    if (int_points != nullptr) {
        for (size_t vp_idx = 0; vp_idx < n_tri; vp_idx++) {
            for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
                for (size_t tri_pt_idx = 0; tri_pt_idx < 3; tri_pt_idx++) {
                    size_t res_idx = tri_pt_idx * n_tri * n_tri + tri_idx * n_tri + vp_idx;
                    size_t ray_idx = tri_idx * 3 + tri_pt_idx;
                    if (intersection_arr[res_idx]) {
                        int_points[vp_idx][ray_idx] = result_int_points[res_idx];
                    } else {
                        int_points[vp_idx][ray_idx] = vec3(
                            std::numeric_limits<float>::infinity(), 
                            std::numeric_limits<float>::infinity(), 
                            std::numeric_limits<float>::infinity()
                        );
                    }
                    // int_points[vp_idx][ray_idx] = result_int_points[res_idx];
                }
            }
        }
    }

    delete[] ends;
    delete[] tri;
    delete[] result_arr;
    delete[] intersection_arr;
    delete[] result_int_points;
}

extern "C" void cuda_kernel_many(
    const std::vector<Viewpoint>& viewpoints,
    const std::vector<size_t>& triangle_indices,
    const std::vector<Triangle*>& faces,
    bool* collisions,
    vec3** int_points
    ) {

    // put viewpoints into arrays
    size_t n_tri = faces.size();
    size_t n_vp = viewpoints.size();
    size_t n_ray = n_vp * 3;

    vec3 *starts = new vec3[n_ray];
    vec3 *ends = new vec3[n_ray];
    Triangle *tri = new Triangle[n_tri];
    bool *result_arr = new bool[n_vp];
    bool *intersection_arr = new bool[n_ray * n_tri];
    vec3 *result_int_points = new vec3[n_ray * n_tri];

    // thread, block size
    size_t thread_x = 16;
    size_t thread_y = 16;
    size_t thread_z = 3;

    // put faces into array
    for (size_t i = 0; i < n_tri; i++) {
        tri[i] = *faces[i];
    }

    // put viewpoints into array
    for (size_t tri_pt_idx = 0; tri_pt_idx < 3; tri_pt_idx++) {
        for (size_t vp_idx = 0; vp_idx < n_vp; vp_idx++) {
            size_t ray_idx = vp_idx * 3 + tri_pt_idx; // n_ray
            size_t tri_idx = triangle_indices[vp_idx];
            starts[ray_idx] = viewpoints[vp_idx].pose;
            if (tri_pt_idx == 0) {
                ends[ray_idx] = tri[tri_idx].a;
            } else if (tri_pt_idx == 1) {
                ends[ray_idx] = tri[tri_idx].b;
            } else if (tri_pt_idx == 2) {
                ends[ray_idx] = tri[tri_idx].c;
            }
        }
    }

    // allocate gpu memory
    vec3 *d_starts;
    vec3 *d_ends;
    Triangle *d_tri;

    bool *d_intersections; // collisions per ray per triangle
    bool *d_result; // collisions per triangle
    vec3 *d_int_points; // intersection points

    cudaMalloc(&d_starts, n_ray * sizeof(vec3));
    cudaMalloc(&d_ends, n_ray * sizeof(vec3));
    cudaMalloc(&d_tri, n_tri * sizeof(Triangle));
    cudaMalloc(&d_intersections, n_ray * n_tri * sizeof(bool));
    cudaMalloc(&d_result, n_vp * sizeof(bool));
    cudaMalloc(&d_int_points, n_ray * n_tri * sizeof(vec3));

    // copy data to gpu
    cudaMemcpy(d_starts, starts, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ends, ends, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tri, tri, n_tri * sizeof(Triangle), cudaMemcpyHostToDevice);

    // set thread and block size
    dim3 threadsPerBlock(thread_x, thread_y, thread_z);

    // 2D blocks because 3d blocks can account for the 3rd dim by themselves
    dim3 numBlocks(int((n_vp + threadsPerBlock.x - 1) / threadsPerBlock.x), (n_tri + threadsPerBlock.y - 1) / threadsPerBlock.y, 1);
    ray_int_plane_many<<<numBlocks, threadsPerBlock>>>(
        d_intersections, 
        d_int_points, 
        d_starts, 
        d_ends, 
        n_vp,
        d_tri, 
        n_tri
    );

    cudaDeviceSynchronize();

    cudaMemcpy(intersection_arr, d_intersections, n_ray * n_tri * sizeof(bool), cudaMemcpyDeviceToHost);

    // same as above but without the 3rd dimension
    threadsPerBlock.x = 1024;
    threadsPerBlock.y = 1;
    threadsPerBlock.z = 1;

    // reusing numBlocks
    numBlocks.x = (n_tri + threadsPerBlock.x - 1) / threadsPerBlock.x;
    numBlocks.y = 1; // numBlocks.z is already 1
    collision_or<<<numBlocks, threadsPerBlock>>>(d_result, d_intersections, n_vp, n_tri);

    cudaDeviceSynchronize();

    cudaMemcpy(result_arr, d_result, n_vp * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(result_int_points, d_int_points, n_ray * n_tri * sizeof(vec3), cudaMemcpyDeviceToHost);

    cudaFree(d_starts);
    cudaFree(d_ends);
    cudaFree(d_tri);
    cudaFree(d_result);
    cudaFree(d_intersections);
    cudaFree(d_int_points);

    for (size_t vp_idx=0; vp_idx < n_vp; vp_idx++) {
        collisions[vp_idx] = result_arr[vp_idx];
    }

    if (int_points != nullptr) {
        for (size_t vp_idx = 0; vp_idx < n_vp; vp_idx++) {
            for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
                for (size_t tri_pt_idx = 0; tri_pt_idx < 3; tri_pt_idx++) {
                    size_t res_idx = tri_pt_idx * n_tri * n_vp + tri_idx * n_vp + vp_idx;
                    size_t ray_idx = tri_idx * 3 + tri_pt_idx;
                    if (intersection_arr[res_idx]) {
                        int_points[vp_idx][ray_idx] = result_int_points[res_idx];
                    } else {
                        int_points[vp_idx][ray_idx] = vec3(
                            std::numeric_limits<float>::infinity(), 
                            std::numeric_limits<float>::infinity(), 
                            std::numeric_limits<float>::infinity()
                        );
                    }
                }
            }
        }
    }

    delete[] starts;
    delete[] ends;
    delete[] tri;
    delete[] result_arr;
    delete[] intersection_arr;
    delete[] result_int_points;
}

extern "C" void cuda_kernel_inc_angle(
    std::vector<Viewpoint>& viewpoints,
    std::vector<Triangle*>& faces,
    std::vector<std::vector<float>>& inc_angles // n_vp x n_tri
    ) {

    // put viewpoints into arrays
    size_t n_tri = faces.size();
    size_t n_vp = viewpoints.size();

    vec3 *poses = new vec3[n_vp];
    vec3 *centroids = new vec3[n_tri];
    vec3 *normals = new vec3[n_tri];
    float *angles = new float[n_vp * n_tri];

    // put faces into array
    for (size_t vp_idx = 0; vp_idx < n_vp; vp_idx++) {
        poses[vp_idx] = viewpoints[vp_idx].pose;
    }

    for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
        centroids[tri_idx] = faces[tri_idx]->getCentroid();
        normals[tri_idx] = faces[tri_idx]->n;
    }

    // initialize pointers for gpu memory
    vec3 *d_poses;
    vec3 *d_centroids;
    vec3 *d_normals;
    float *d_angles;

    // allocate memory on gpu
    cudaMalloc(&d_poses, n_vp * sizeof(vec3));
    cudaMalloc(&d_centroids, n_tri * sizeof(vec3));
    cudaMalloc(&d_normals, n_tri * sizeof(vec3));
    cudaMalloc(&d_angles, n_vp * n_tri * sizeof(float));

    // copy data to gpu
    cudaMemcpy(d_poses, poses, n_vp * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_centroids, centroids, n_tri * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_normals, normals, n_tri * sizeof(vec3), cudaMemcpyHostToDevice);

    // set thread and block size
    size_t thread_x = 32;
    size_t thread_y = 32;
    size_t thread_z = 1;
    dim3 threadsPerBlock(thread_x, thread_y, thread_z);
    dim3 numBlocks(
        int((n_vp + threadsPerBlock.x - 1) / threadsPerBlock.x), 
        int((n_tri + threadsPerBlock.y - 1) / threadsPerBlock.y),
        1
    );

    // calculate incidence angles
    inc_angle<<<numBlocks, threadsPerBlock>>>(
        d_angles, 
        d_poses, 
        d_centroids, 
        d_normals, 
        n_vp, 
        n_tri
    );

    // copy data from gpu
    cudaMemcpy(angles, d_angles, n_vp * n_tri * sizeof(float), cudaMemcpyDeviceToHost);

    // free memory from gpu
    cudaFree(d_poses);
    cudaFree(d_centroids);
    cudaFree(d_normals);
    cudaFree(d_angles);

    // put data into inc_angles
    for (size_t vp_idx = 0; vp_idx < n_vp; vp_idx++) {
        std::vector<float> vp_angles;
        for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
            vp_angles.push_back(angles[tri_idx * n_vp + vp_idx]);
        }
        inc_angles.push_back(vp_angles);
    }

    delete[] poses;
    delete[] centroids;
    delete[] normals;
    delete[] angles;
}