#include "cuda_kernels.h"
#include "vec3_struct.hpp"
#include "triangle_struct.hpp"
#include "viewpoint_struct.hpp"

#include <vector>

// dims: viewpoints (x dim) x faces (y dim) x 3 (tri dim)
__global__ void ray_int_plane(
    bool *result, // flattened 3d
    vec3 *int_points, // flattened 3d
    const vec3 origin,  // vp (vec3)
    const vec3 *ends,    // n_vp (1dim)
    const Triangle *tri,// n_tri (1dim)
    size_t n_tri
    ) {

    // epsilon for floating point comparison
    float eps = 1e-7f;

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

    // int_points[res_idx] = ends[ray_idx];

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

__global__ void collision_or(bool* collision, bool* intersections, size_t n_tri) {
    // get indices
    size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_ray
    size_t tri_idx = blockIdx.y * blockDim.y + threadIdx.y; // n_tri
    // size_t tri_pt_idx = blockIdx.z * blockDim.z + threadIdx.z; // 3
    size_t v0_idx = 0 * n_tri * n_tri + tri_idx * n_tri + vp_idx;
    size_t v1_idx = 1 * n_tri * n_tri + tri_idx * n_tri + vp_idx;
    size_t v2_idx = 2 * n_tri * n_tri + tri_idx * n_tri + vp_idx;

    if (v0_idx >= n_tri * n_tri) { return; }

    collision[v0_idx] = intersections[v0_idx] || intersections[v1_idx] || intersections[v2_idx];

    return;
}


extern "C" void cuda_kernel_intersect_triangles(
    const Viewpoint& vp, 
    const std::vector<Triangle*>& faces,
    bool*** collisions,
    vec3*** int_points
    ) {

    // put viewpoints into arrays
    size_t n_tri = faces.size();
    size_t n_ray = n_tri * 3;

    vec3 *ends = new vec3[n_ray];
    Triangle *tri = new Triangle[n_tri];
    std::cout << "n_tri = " << n_tri << std::endl;
    bool *result_arr = new bool[n_ray * n_tri];
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
    // bool *d_result; // collisions per triangle
    vec3 *d_int_points; // intersection points

    cudaMalloc(&d_ends, n_ray * sizeof(vec3));
    cudaMalloc(&d_tri, n_tri * sizeof(Triangle));
    cudaMalloc(&d_intersections, n_ray * n_tri * sizeof(bool));
    // cudaMalloc(&d_result, n_tri * n_tri * sizeof(bool));
    cudaMalloc(&d_int_points, n_ray * n_tri * sizeof(vec3));

    // copy data to gpu
    cudaMemcpy(d_ends, ends, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tri, tri, n_tri * sizeof(Triangle), cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(thread_x, thread_y, thread_z);
    // std::cout << "numBlocks = " << numBlocks.x << " " << numBlocks.y << " " << numBlocks.z << std::endl;
    std::cout << "threadsPerBlock = " << threadsPerBlock.x << " " << threadsPerBlock.y << " " << threadsPerBlock.z << std::endl;
    // 2D blocks because 3d blocks can account for the 3rd dim by themselves
    dim3 numBlocks(int((n_tri + threadsPerBlock.x - 1) / threadsPerBlock.x), (n_tri + threadsPerBlock.y - 1) / threadsPerBlock.y, 1);
    ray_int_plane<<<numBlocks, threadsPerBlock>>>(d_intersections, d_int_points, vp.pose, d_ends, d_tri, n_tri);

    // // same as above but without the 3rd dimension
    // threadsPerBlock.x = 32;
    // threadsPerBlock.y = 32;
    // threadsPerBlock.z = 1;
    // collision_or<<<numBlocks, threadsPerBlock>>>(d_result, d_intersections, n_tri);

    cudaMemcpy(result_arr, d_intersections, n_ray * n_tri * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(result_int_points, d_int_points, n_ray * n_tri * sizeof(vec3), cudaMemcpyDeviceToHost);
    // cudaMemcpy(saved_indices, d_saved_indices, thread_x * thread_y * sizeof(size_t), cudaMemcpyDeviceToHost);
    // cudaMemcpy(thread_indices, d_thread_indices, thread_x * thread_y * sizeof(vec3), cudaMemcpyDeviceToHost);
    // cudaMemcpy(result_arr, d_result, thread_x * thread_y * thread_z * sizeof(bool), cudaMemcpyDeviceToHost);
    // cudaMemcpy(result_int_points, d_int_points, thread_x * thread_y * thread_z * sizeof(vec3), cudaMemcpyDeviceToHost);

    cudaFree(d_tri);
    cudaFree(d_ends);
    cudaFree(d_intersections);
    cudaFree(d_int_points);
    // cudaFree(d_result);
    // cudaFree(d_saved_indices);
    // cudaFree(d_thread_indices);
    std::cout << "Ray start=" << vp.pose.toString() << std::endl;

    for (size_t tri_pt_idx = 0; tri_pt_idx < 3; tri_pt_idx++) {
        for (size_t vp_idx = 0; vp_idx < n_tri; vp_idx++) {
            for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
                // size_t saved_idx = k * thread_x * thread_y + j * thread_x + i;
                size_t ray_idx = vp_idx * 3 + tri_pt_idx;
                size_t res_idx = tri_pt_idx * n_tri * n_tri + tri_idx * n_tri + vp_idx;
                std::cout << "Ray end=" << ends[ray_idx].toString() << " | Tri z=" << tri[tri_idx].a.z;
                std::cout << " >>>> Collision=" << result_arr[res_idx];
                std::cout << " | IntPoint=" << result_int_points[res_idx].toString() << std::endl;

                //     std::cout << " oob idx=";
                //     std::cout << result_int_points[res_idx].x << std::endl;
                // }
                // std::cout << res_idx << std::endl;
                // std::cout << "kernel ran?=" << result_arr[res_idx] << " | cuda idxs =";
                // std::cout << result_int_points[res_idx].toString() << " | ";
                // size_t res_idx_cuda = result_int_points[res_idx].x * n_vp * n_tri + result_int_points[res_idx].y * n_vp + result_int_points[res_idx].z;
                // std::cout << "cuda idx =" << res_idx_cuda << std::endl;
                // collisions[i][j][k] = result_arr[res_idx];
                // int_points[i][j][k] = result_int_points[res_idx];
            }
            // collisions[i][j] = result_arr[j * n_vp + i];
        }
    }

}