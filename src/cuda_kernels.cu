#include "cuda_kernels.h"
#include "vec3_struct.hpp"
#include "triangle_struct.hpp"
#include "viewpoint_struct.hpp"

#include <vector>
#include <limits>
#include <cmath>

__device__ void pinhole_camera(
    bool& out_of_frame,  //n_points
    vec3 pose,
    vec3 viewdir,
    vec3 point,
    float hfov=M_PI/4.0f, // rad
    float vfov=M_PI/4.0f // rad
    ) {
    // calculate the angle between the view direction and the vector from the viewpoint to the intersection point
    // size_t point_idx = blockIdx.y * blockDim.y + threadIdx.y; // n_points

    // if (point_idx > n_points - 1) { return; }

    // set up variables
    // vec3 vec = point[point_idx] - pose;
    vec3 vec = point - pose;
    // vec3 viewdir = viewdir;
    float norm_dot = vec.dot(viewdir);

    // check if point is behind camera
    if (norm_dot <= 0) {
        out_of_frame = true;
        // out_of_frame[point_idx] = false;
        return;
    }

    // project point onto view plane
    float d = 1.0f; // distance from viewpoint to view plane
    float w = 2 * d * tanf(hfov/2);
    float h = 2 * d * tanf(vfov/2);
    vec3 point_proj = (vec/norm_dot - viewdir) * d;
    vec3 v_hat = vec3(
        viewdir.z * cosf( atan2f(viewdir.y, viewdir.x) ),
        viewdir.y * sinf( atan2f(viewdir.y, viewdir.x) ),
        sqrtf(viewdir.x * viewdir.x + viewdir.y * viewdir.y) // 0 to 1
    );
    vec3 u_hat = v_hat.cross(viewdir);

    // check if point is within field of view
    if (abs(point_proj.dot(u_hat)) < w/2 && abs(point_proj.dot(v_hat)) < h/2) {
        out_of_frame = false;
        // visible[point_idx] = true;
        return;
    }

    out_of_frame = true;
    // visible[point_idx] = false;
    return;
}

// __global__ void pinhole_camera(
//     bool *visible,  // n_vp x n_points
//     vec3 *pose,
//     vec3 *viewdirs,
//     size_t n_vp,
//     vec3 *point,
//     size_t n_points,
//     float hfov=M_PI/4.0f, // rad
//     float vfov=M_PI/4.0f // rad
//     ) {
//     // calculate the angle between the view direction and the vector from the viewpoint to the intersection point
//     size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_vp
//     size_t point_idx = blockIdx.y * blockDim.y + threadIdx.y; // n_points
//     size_t res_idx = point_idx * n_vp + vp_idx; // n_points * n_vp

//     if (vp_idx > n_vp - 1 || point_idx > n_points - 1) { return; }

//     // set up variables
//     vec3 vec = point[point_idx] - pose[vp_idx];
//     vec3 viewdir = viewdirs[vp_idx];
//     float norm_dot = vec.dot(viewdir);

//     // check if point is behind camera
//     if (norm_dot <= 0) {
//         visible[res_idx] = false;
//         return;
//     }

//     // project point onto view plane
//     float d = 1.0f; // distance from viewpoint to view plane
//     float w = 2 * d * tanf(hfov/2);
//     float h = 2 * d * tanf(vfov/2);
//     vec3 point_proj = (vec/norm_dot - viewdir) * d;
//     vec3 v_hat = vec3(
//         viewdir.z * cosf( atan2f(viewdir.y, viewdir.x) ),
//         viewdir.y * sinf( atan2f(viewdir.y, viewdir.x) ),
//         sqrtf(viewdir.x * viewdir.x + viewdir.y * viewdir.y) // 0 to 1
//     );
//     vec3 u_hat = v_hat.cross(viewdir);

//     // check if point is within field of view
//     if (abs(point_proj.dot(u_hat)) < w/2 && abs(point_proj.dot(v_hat)) < h/2) {
//         visible[res_idx] = true;
//         return;
//     }

//     visible[res_idx] = false;
//     return;
// }

// one viewpoint, mapped to many end points
__global__ void ray_int_tri(
    bool *result, // flattened 3d
    vec3 *int_points, // flattened 3d
    const vec3 origin,  // vp (vec3)
    const vec3 viewdir, // vp (vec3)
    const vec3 *ends,    // n_vp (1dim)
    const Triangle *tri, // n_tri (1dim)
    size_t n_tri
    ) {
    // true for not visiblw

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

    // check if triangle is facing away from camera
    float norm_dot = tri[tri_idx].n.dot(origin - end);
    // if (norm_dot > 0) {
    //     result[res_idx] = true;
    //     return;
    // }

    // look for any intersections between the ray and triangle
    vec3 e1 = tri[tri_idx].b - tri[tri_idx].a;
    vec3 e2 = tri[tri_idx].c - tri[tri_idx].a;
    vec3 h = vec.cross(e2);
    float a = e1.dot(h);

    // if ray is parallel to triangle
    if (a > -eps && a < eps) {
        pinhole_camera(result[res_idx], origin, viewdir, end);
        // result[res_idx] = false;
        return;
    }

    float f = 1 / a;
    vec3 s = origin - tri[tri_idx].a;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) {
        pinhole_camera(result[res_idx], origin, viewdir, end);
        // result[res_idx] = false;
        return;
    }
    vec3 q = s.cross(e1);
    float v = f * vec.dot(q);
    if (v < 0.0f || u + v > 1.0f) {
        pinhole_camera(result[res_idx], origin, viewdir, end);
        // result[res_idx] = false;
        return;
    }

    // find intersection point
    float t = f * e2.dot(q);
    if (t < eps) {
        pinhole_camera(result[res_idx], origin, viewdir, end);
        // result[res_idx] = false;
        return;
    }
    vec3 intPoint = origin + vec * t;
    int_points[res_idx] = intPoint;

    // check if intersection point is between origin and end
    vec3 vec_dir = vec/vec.norm();
    if ((intPoint-origin).dot(vec_dir) < vec.norm() - eps && 
        (intPoint-origin).dot(vec_dir) > 0
        ) {
        result[res_idx] = true;
        return;
    }

    pinhole_camera(result[res_idx], origin, viewdir, end);
    // result[res_idx] = false;
    return;
}

// dims: viewpoints (x dim) x faces (y dim)
// many origins, each mapped to an end point
__global__ void ray_int_tri_many_2d(
    bool *result, // flattened 2d
    vec3 *int_points, // flattened 2d
    const vec3 *starts,  // n_vp (vec3)
    // const vec3 *viewdirs, // n_vp (vec3)
    const vec3 *ends,    // n_vp (vec3)
    size_t n_vp,
    const Triangle *tri,// n_tri (1dim)
    size_t n_tri
    ) {

    // epsilon for floating point comparison
    float eps = 1e-6f;

    // get indices
    size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_tri
    size_t tri_idx = blockIdx.y * blockDim.y + threadIdx.y; // n_tri
    size_t res_idx = tri_idx * n_vp + vp_idx; // n_tri * n_vp

    if (vp_idx > n_vp - 1 || tri_idx > n_tri - 1) { return; }

    // instantiate ray
    vec3 origin = starts[vp_idx];
    vec3 end = ends[vp_idx];
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
    if (t < eps) {
        result[res_idx] = false;
        return;
    }

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
__global__ void ray_int_tri_many(
    bool *result, // flattened 3d
    vec3 *int_points, // flattened 3d
    const vec3 *starts,  // n_ray (vec3)
    const vec3 *viewdirs, // n_vp (vec3)
    const vec3 *ends,    // n_ray (1dim)
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
    vec3 viewdir = viewdirs[ray_idx];
    vec3 end = ends[ray_idx];
    vec3 vec = end - origin;

    // look for any intersections between the ray and triangle
    vec3 e1 = tri[tri_idx].b - tri[tri_idx].a;
    vec3 e2 = tri[tri_idx].c - tri[tri_idx].a;
    vec3 h = vec.cross(e2);
    float a = e1.dot(h);

    // if ray is parallel to triangle
    if (a > -eps && a < eps) {
        pinhole_camera(result[res_idx], origin, viewdir, end);
        // result[res_idx] = false;
        return;
    }

    float f = 1 / a;
    vec3 s = origin - tri[tri_idx].a;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) {
        pinhole_camera(result[res_idx], origin, viewdir, end);
        // result[res_idx] = false;
        return;
    }
    vec3 q = s.cross(e1);
    float v = f * vec.dot(q);
    if (v < 0.0f || u + v > 1.0f) {
        pinhole_camera(result[res_idx], origin, viewdir, end);
        // result[res_idx] = false;
        return;
    }

    // find intersection point
    float t = f * e2.dot(q);

    if (t < eps) {
        pinhole_camera(result[res_idx], origin, viewdir, end);
        // result[res_idx] = false;
        return;
    }

    vec3 intPoint = origin + vec * t;
    int_points[res_idx] = intPoint;

    // check if intersection point is between origin and end
    vec3 vec_dir = vec/vec.norm();
    if ((intPoint-origin).dot(vec_dir) < vec.norm() - eps && (intPoint-origin).dot(vec_dir) > 0) {
        result[res_idx] = true;
        return;
    }

    pinhole_camera(result[res_idx], origin, viewdir, end);
    // result[res_idx] = false;
    return;
}

__global__ void ray_int_plane(
    bool *result, // n_rays
    vec3 *int_points, // n_rays
    vec3 *ray_starts,  // n_rays (vec3)
    vec3 *ray_ends,   // n_rays (vec3)
    vec3 plane_point,
    vec3 plane_normal,
    size_t n_rays
    ) {
    // get ray index 
    size_t ray_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_rays

    // check if out of bounds on rays
    if (ray_idx > n_rays - 1) { return; }

    // get ray vecs to test
    vec3 origin_to_point = plane_point - ray_starts[ray_idx];
    vec3 end_to_point = plane_point - ray_ends[ray_idx];

    // get projections onto plane normal
    float origin_proj = origin_to_point.dot(plane_normal);
    float end_proj = end_to_point.dot(plane_normal);
    float abs_origin_proj = fabsf(origin_proj);
    float abs_end_proj = fabsf(end_proj);

    // test if ray intersects plane (end and origin are on opposite sides of plane)
    if (abs_origin_proj < 1e-6f || abs_end_proj < 1e-6f) {
        result[ray_idx] = false;
        int_points[ray_idx] = vec3(
            INFINITY,
            INFINITY,
            INFINITY
        );
        return;
    }

    // must be opposite signs if they are on opposite sides of the plane
    if (origin_proj * end_proj > 0) {
        result[ray_idx] = false;
        int_points[ray_idx] = vec3(
            INFINITY,
            INFINITY,
            INFINITY
        );
        return;
    }
    float fac = abs_origin_proj / (abs_origin_proj + abs_end_proj);
    int_points[ray_idx] = ray_starts[ray_idx] * (1 - fac) + ray_ends[ray_idx] * fac;
    result[ray_idx] = true;
}

__global__ void collision_odd(bool* vp_collision, const bool* ray_tri_collision, size_t n_vp, size_t n_tri) {
    // for each viewpoint-triangle correspondance, check if rays to each vertex collide with any other triangle. if so, write in true
    // get viewpoint index
    size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_vp

    if (vp_idx > n_vp - 1) { return; } // n_tri = n_vp

    size_t count = 0;
    vp_collision[vp_idx] = false;
    for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
        size_t res_idx = tri_idx * n_vp + vp_idx;
        if (ray_tri_collision[res_idx]) {
            count++;
        }
    }
    if ((count % 2) == 1) {
        vp_collision[vp_idx] = true;
    } else {
        vp_collision[vp_idx] = false;
    }
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

__global__ void collision_or_2d(bool* vp_collision, const bool* ray_tri_collision, size_t n_vp, size_t n_tri) {
    // for each viewpoint-triangle correspondance, check if rays to each vertex collide with any other triangle. if so, write in true
    // get viewpoint index
    size_t vp_idx = blockIdx.x * blockDim.x + threadIdx.x; // n_vp

    if (vp_idx > n_vp - 1) { return; } // n_tri = n_vp

    vp_collision[vp_idx] = false;
    for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
        size_t res_idx = tri_idx * n_vp + vp_idx;
        if (ray_tri_collision[res_idx]) {
            vp_collision[vp_idx] = true;
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
    // vec3 vec = centroids[tri_idx] - poses[vp_idx];
    vec3 vec = poses[vp_idx] - centroids[tri_idx];
    vec3 norm = normals[tri_idx];
    float angle = acos(vec.dot(norm)/(vec.norm()*norm.norm()));
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
    ray_int_tri<<<numBlocks, threadsPerBlock>>>(
        d_intersections, 
        d_int_points, 
        vp.pose, 
        vp.viewdir,
        d_ends, 
        d_tri, 
        n_tri
    );

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

extern "C" void cuda_kernel_ray_int_plane(
    const std::vector<vec3>& ray_starts,
    const std::vector<vec3>& ray_ends,
    const vec3& plane_point,
    const vec3& plane_normal,
    bool* collisions,
    vec3* int_points
) {
    size_t n_rays = ray_starts.size();

    vec3 *starts = new vec3[n_rays];
    vec3 *ends = new vec3[n_rays];
    bool *result_arr = new bool[n_rays];
    vec3 *result_int_points = new vec3[n_rays];

    // thread, block size
    size_t thread_x = 1024;

    // put viewpoints into array
    for (size_t i = 0; i < n_rays; i++) {
        starts[i] = ray_starts[i];
        ends[i] = ray_ends[i];
    }

    vec3 *d_starts;
    vec3 *d_ends;
    bool *d_result;
    vec3 *d_int_points;

    cudaMalloc(&d_starts, n_rays * sizeof(vec3));
    cudaMalloc(&d_ends, n_rays * sizeof(vec3));
    cudaMalloc(&d_result, n_rays * sizeof(bool));
    cudaMalloc(&d_int_points, n_rays * sizeof(vec3));

    cudaMemcpy(d_starts, starts, n_rays * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ends, ends, n_rays * sizeof(vec3), cudaMemcpyHostToDevice);

    // set up thread and block size
    dim3 threadsPerBlock(thread_x, 1, 1);
    dim3 numBlocks((n_rays + threadsPerBlock.x - 1) / threadsPerBlock.x, 1, 1);

    // run kernel
    ray_int_plane<<<numBlocks, threadsPerBlock>>>(
        d_result,
        d_int_points,
        d_starts,
        d_ends,
        plane_point,
        plane_normal,
        n_rays
    );

    cudaDeviceSynchronize();

    // copy results back to host
    cudaMemcpy(result_arr, d_result, n_rays * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(result_int_points, d_int_points, n_rays * sizeof(vec3), cudaMemcpyDeviceToHost);

    // free gpu memory
    cudaFree(d_starts);
    cudaFree(d_ends);
    cudaFree(d_result);
    cudaFree(d_int_points);

    // copy memory into output
    for (size_t i = 0; i < n_rays; i++) {
        collisions[i] = result_arr[i];
        int_points[i] = result_int_points[i];
    }

    delete[] starts;
    delete[] ends;
    delete[] result_arr;
    delete[] result_int_points;
}

extern "C" void cuda_kernel_many_ray(
    const std::vector<vec3>& start_ray,
    const std::vector<vec3>& end_ray,
    const std::vector<Triangle*>& faces,
    bool* collisions,
    vec3** int_points
    ) {

    // put viewpoints into arrays
    size_t n_tri = faces.size();
    size_t n_ray = start_ray.size();

    vec3 *starts = new vec3[n_ray];
    vec3 *ends = new vec3[n_ray];
    Triangle *tri = new Triangle[n_tri];
    bool *result_arr = new bool[n_ray];
    bool *intersection_arr = new bool[n_ray * n_tri];
    vec3 *result_int_points = new vec3[n_ray * n_tri];

    // thread, block size
    size_t thread_x = 32;
    size_t thread_y = 32;
    size_t thread_z = 1;

    // put faces into array
    for (size_t i = 0; i < n_tri; i++) {
        tri[i] = *faces[i];
    }

    // put viewpoints into array
    for (size_t ray_idx = 0; ray_idx < n_ray; ray_idx++) {
        starts[ray_idx] = start_ray[ray_idx];
        ends[ray_idx] = end_ray[ray_idx];
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
    cudaMalloc(&d_result, n_ray * sizeof(bool));
    cudaMalloc(&d_int_points, n_ray * n_tri * sizeof(vec3));

    // copy data to gpu
    cudaMemcpy(d_starts, starts, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ends, ends, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tri, tri, n_tri * sizeof(Triangle), cudaMemcpyHostToDevice);

    // set thread and block size
    dim3 threadsPerBlock(thread_x, thread_y, thread_z);

    // 2D blocks because 3d blocks can account for the 3rd dim by themselves
    dim3 numBlocks(int((n_ray + threadsPerBlock.x - 1) / threadsPerBlock.x), (n_tri + threadsPerBlock.y - 1) / threadsPerBlock.y, 1);
    ray_int_tri_many_2d<<<numBlocks, threadsPerBlock>>>(
        d_intersections, 
        d_int_points,
        d_starts, 
        d_ends, 
        n_ray,
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
    collision_or_2d<<<numBlocks, threadsPerBlock>>>(d_result, d_intersections, n_ray, n_tri);

    cudaDeviceSynchronize();

    cudaMemcpy(result_arr, d_result, n_ray * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(result_int_points, d_int_points, n_ray * n_tri * sizeof(vec3), cudaMemcpyDeviceToHost);

    cudaFree(d_starts);
    cudaFree(d_ends);
    cudaFree(d_tri);
    cudaFree(d_result);
    cudaFree(d_intersections);
    cudaFree(d_int_points);

    for (size_t vp_idx=0; vp_idx < n_ray; vp_idx++) {
        collisions[vp_idx] = result_arr[vp_idx];
    }

    if (int_points != nullptr) {
        for (size_t ray_idx = 0; ray_idx < n_ray; ray_idx++){
            for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
                size_t res_idx = tri_idx * n_ray + ray_idx;
                int_points[ray_idx][tri_idx] = result_int_points[res_idx];
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
    vec3 *viewdirs = new vec3[n_ray];
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
            viewdirs[ray_idx] = viewpoints[vp_idx].viewdir;
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
    vec3 *d_viewdirs;
    vec3 *d_ends;
    Triangle *d_tri;

    bool *d_intersections; // collisions per ray per triangle
    bool *d_result; // collisions per triangle
    vec3 *d_int_points; // intersection points

    cudaMalloc(&d_starts, n_ray * sizeof(vec3));
    cudaMalloc(&d_viewdirs, n_ray * sizeof(vec3));
    cudaMalloc(&d_ends, n_ray * sizeof(vec3));
    cudaMalloc(&d_tri, n_tri * sizeof(Triangle));
    cudaMalloc(&d_intersections, n_ray * n_tri * sizeof(bool));
    cudaMalloc(&d_result, n_vp * sizeof(bool));
    cudaMalloc(&d_int_points, n_ray * n_tri * sizeof(vec3));

    // copy data to gpu
    cudaMemcpy(d_starts, starts, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_viewdirs, viewdirs, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ends, ends, n_ray * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tri, tri, n_tri * sizeof(Triangle), cudaMemcpyHostToDevice);

    // set thread and block size
    dim3 threadsPerBlock(thread_x, thread_y, thread_z);

    // 2D blocks because 3d blocks can account for the 3rd dim by themselves
    dim3 numBlocks(int((n_vp + threadsPerBlock.x - 1) / threadsPerBlock.x), (n_tri + threadsPerBlock.y - 1) / threadsPerBlock.y, 1);
    ray_int_tri_many<<<numBlocks, threadsPerBlock>>>(
        d_intersections, 
        d_int_points, 
        d_starts, 
        d_viewdirs,
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
    cudaFree(d_viewdirs);
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
    delete[] viewdirs;
    delete[] ends;
    delete[] tri;
    delete[] result_arr;
    delete[] intersection_arr;
    delete[] result_int_points;
}

extern "C" void cuda_kernel_collision_points(
    const std::vector<Viewpoint>& viewpoints,
    const std::vector<Triangle*>& faces,
    const vec3 free_space_point,
    std::vector<bool>& in_collision // number of collisions
    ) {

    // put viewpoints into arrays
    size_t n_tri = faces.size();
    size_t n_vp = viewpoints.size();

    vec3 *starts = new vec3[n_vp];
    // vec3 *viewdirs = new vec3[n_vp];
    vec3 *ends = new vec3[n_vp];
    Triangle *tri = new Triangle[n_tri];
    bool *result_arr = new bool[n_vp];
    bool *intersection_arr = new bool[n_vp * n_tri];
    vec3 *result_int_points = new vec3[n_vp * n_tri];

    // thread, block size
    size_t thread_x = 32;
    size_t thread_y = 32;
    size_t thread_z = 1;

    // put faces into array
    for (size_t i = 0; i < n_tri; i++) {
        tri[i] = *faces[i];
    }

    // put viewpoints into array
    for (size_t vp_idx = 0; vp_idx < n_vp; vp_idx++) {
        starts[vp_idx] = viewpoints[vp_idx].pose;
        // viewdirs[vp_idx] = viewpoints[vp_idx].viewdir;
        ends[vp_idx] = free_space_point;
    }

    // allocate gpu memory
    vec3 *d_starts;
    // vec3 *d_viewdirs;
    vec3 *d_ends;
    Triangle *d_tri;

    bool *d_intersections; // collisions per ray per triangle
    bool *d_result; // collisions per triangle
    vec3 *d_int_points; // intersection points

    cudaMalloc(&d_starts, n_vp * sizeof(vec3));
    // cudaMalloc(&d_viewdirs, n_vp * sizeof(vec3));
    cudaMalloc(&d_ends, n_vp * sizeof(vec3));
    cudaMalloc(&d_tri, n_tri * sizeof(Triangle));
    cudaMalloc(&d_intersections, n_vp * n_tri * sizeof(bool));
    cudaMalloc(&d_result, n_vp * sizeof(bool));
    cudaMalloc(&d_int_points, n_vp * n_tri * sizeof(vec3));

    // copy data to gpu
    cudaMemcpy(d_starts, starts, n_vp * sizeof(vec3), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_viewdirs, viewdirs, n_vp * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ends, ends, n_vp * sizeof(vec3), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tri, tri, n_tri * sizeof(Triangle), cudaMemcpyHostToDevice);

    // set thread and block size
    dim3 threadsPerBlock(thread_x, thread_y, thread_z);

    // 2D blocks because 3d blocks can account for the 3rd dim by themselves
    dim3 numBlocks(int((n_vp + threadsPerBlock.x - 1) / threadsPerBlock.x), (n_tri + threadsPerBlock.y - 1) / threadsPerBlock.y, 1);
    ray_int_tri_many_2d<<<numBlocks, threadsPerBlock>>>(
        d_intersections, 
        d_int_points, 
        d_starts, 
        // d_viewdirs,
        d_ends,
        n_vp,
        d_tri, 
        n_tri
    );

    cudaDeviceSynchronize();

    cudaMemcpy(intersection_arr, d_intersections, n_vp * n_tri * sizeof(bool), cudaMemcpyDeviceToHost);
    // cudaMemcpy(result_int_points, d_int_points, n_vp * n_tri * sizeof(vec3), cudaMemcpyDeviceToHost);

    // same as above but without the 3rd dimension
    threadsPerBlock.x = 1024;
    threadsPerBlock.y = 1;
    threadsPerBlock.z = 1;

    // reusing numBlocks
    numBlocks.x = (n_tri + threadsPerBlock.x - 1) / threadsPerBlock.x;
    numBlocks.y = 1; // numBlocks.z is already 1
    collision_odd<<<numBlocks, threadsPerBlock>>>(
        d_result, 
        d_intersections, 
        n_vp, 
        n_tri
    );

    cudaDeviceSynchronize();

    cudaMemcpy(result_arr, d_result, n_vp * sizeof(bool), cudaMemcpyDeviceToHost);

    cudaFree(d_starts);
    // cudaFree(d_viewdirs);
    cudaFree(d_ends);
    cudaFree(d_tri);
    cudaFree(d_result);
    cudaFree(d_intersections);
    cudaFree(d_int_points);

    for (size_t vp_idx=0; vp_idx < n_vp; vp_idx++) {
        in_collision.push_back(result_arr[vp_idx]);
    }

    // if (int_points != nullptr) {
    //     for (size_t vp_idx = 0; vp_idx < n_vp; vp_idx++) {
    //         for (size_t tri_idx = 0; tri_idx < n_tri; tri_idx++) {
    //             for (size_t tri_pt_idx = 0; tri_pt_idx < 3; tri_pt_idx++) {
    //                 size_t res_idx = tri_pt_idx * n_tri * n_vp + tri_idx * n_vp + vp_idx;
    //                 size_t ray_idx = tri_idx * 3 + tri_pt_idx;
    //                 if (intersection_arr[res_idx]) {
    //                     int_points[vp_idx][ray_idx] = result_int_points[res_idx];
    //                 } else {
    //                     int_points[vp_idx][ray_idx] = vec3(
    //                         std::numeric_limits<float>::infinity(), 
    //                         std::numeric_limits<float>::infinity(), 
    //                         std::numeric_limits<float>::infinity()
    //                     );
    //                 }
    //             }
    //         }
    //     }
    // }

    delete[] starts;
    // delete[] viewdirs;
    delete[] ends;
    delete[] tri;
    delete[] result_arr;
    delete[] intersection_arr;
    delete[] result_int_points;
}

extern "C" void cuda_kernel_inc_angle(
    const std::vector<Viewpoint>& viewpoints,
    const std::vector<Triangle*>& faces,
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