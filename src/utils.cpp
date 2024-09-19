#include "cuda_kernels.h"
#include "limit_struct.hpp"
#include "node3d_struct.hpp"
#include "obs.hpp"
#include "rrtz.hpp"
#include "triangle_coverage_struct.hpp"
#include "triangle_struct.hpp"
#include "utils.hpp"
#include "vec3_struct.hpp"
#include "viewpoint_struct.hpp"
#include "viewpoint_coverage_gain_struct.hpp"

#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

bool ray_int_plane(Node3D node, Plane plane, float eps, vec3& intPoint) {
    vec3 origin_to_point = plane.point - node.origin;
    // vec3 end_to_point = node.end - node.origin;
    vec3 end_to_point = node.end - plane.point;
    float origin_dot = origin_to_point.dot(plane.normal);
    float end_dot = end_to_point.dot(plane.normal);
    float abs_origin_dot = fabsf(origin_dot);
    float abs_end_dot = fabsf(end_dot);
    if (abs_origin_dot > eps && abs_end_dot > eps && (origin_dot > 0 ) ^ (end_dot > 0)) {
        float fac = abs_origin_dot / (abs_origin_dot + abs_end_dot);
        intPoint = node.end * fac + node.origin * (1 - fac);
        return true;
    }
    return false;
}

bool ray_int_triangle(vec3 origin, vec3 vector, vec3 end, Triangle tri, vec3& intPoint, float eps) {
    // look for any intersections between the ray and triangle (before end-point)
    vec3 e1 = tri.b - tri.a;
    vec3 e2 = tri.c - tri.a;
    vec3 h = vector.cross(e2);
    float a = e1.dot(h);

    // if ray is parallel to triangle
    if (fabsf(a) < eps) { return false; }

    float f = 1 / a;
    vec3 s = origin - tri.a;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) { return false; }
    vec3 q = s.cross(e1);
    float v = f * vector.dot(q);
    if (v < 0.0f || u + v > 1.0f) { return false; }

    // now find intersection point
    float t = f * e2.dot(q);
    intPoint = origin + vector * t; // will save garbage answer to intpoint if t <= eps

    // check if intersection point is between start and endpoint
    vec3 origin_end = end - origin;
    vec3 origin_int = intPoint - origin;
    float origin_end_dot = origin_end.dot(origin_int/origin_int.norm());
    float origin_int_norm = origin_int.norm();


    return origin_int_norm > 0 && origin_int_norm < origin_end_dot;
}

bool ray_int_triangle(vec3 origin, vec3 vector, Triangle tri, vec3& intPoint, float eps) {
    // look for any intersections between the ray and triangle (no end point)
    vec3 e1 = tri.b - tri.a;
    vec3 e2 = tri.c - tri.a;
    vec3 h = vector.cross(e2);
    float a = e1.dot(h);

    // if ray is parallel to triangle
    if (fabsf(a) < eps) { return false; }

    float f = 1 / a;
    vec3 s = origin - tri.a;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) { return false; }
    vec3 q = s.cross(e1);
    float v = f * vector.dot(q);
    if (v < 0.0f || u + v > 1.0f) { return false; }

    // now find intersection point
    float t = f * e2.dot(q);
    intPoint = origin + vector * t; // will save garbage answer to intpoint if t <= eps
    // std::cout << " inFunction: " << intPoint.x << " " << intPoint.y << " " << intPoint.z << std::endl;
    return t > eps;
}

bool ray_int_triangle(Node3D node, Triangle tri, vec3& intPoint, float eps) {
    vec3 e1 = tri.b - tri.a;
    vec3 e2 = tri.c - tri.a;
    vec3 h = node.vector.cross(e2);
    float a = e1.dot(h);

    // if ray is parallel to triangle
    if (fabsf(a) < eps) { return false; }

    float f = 1 / a;
    vec3 s = node.origin - tri.a;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f) { return false; }
    vec3 q = s.cross(e1);
    float v = f * node.vector.dot(q);
    if (v < 0.0f || u + v > 1.0f) { return false; }

    // now find intersection point
    float t = f * e2.dot(q);
    intPoint = node.origin + node.vector * t; // will save garbage answer to intpoint if t < eps
    // std::cout << " inFunction: " << intPoint.x << " " << intPoint.y << " " << intPoint.z << std::endl;
    return t > eps;
}

float heading_change(Node3D node, Node3D next_node) {
    if (node.vector.norm() < 1e-9f || next_node.vector.norm() < 1e-9f) { return 0.0f; }
    float heading_change = acosf(node.vector.dot(next_node.vector) / (node.vector.norm() * next_node.vector.norm() + 1e-9f));
    return heading_change;
}

float heading_change(Node3D node, vec3 vector) {
    return acosf(node.vector.dot(vector) / (node.vector.norm() * vector.norm() + 1e-9f));
}

float heading_change(vec3 v0, vec3 v1) {
    return acosf(v0.dot(v1) / (v0.norm() * v1.norm() + 1e-9f));
}

bool loadCSV(const std::string& filename, std::vector<std::vector<float>>& data, int rowlen, char delimiter, bool raw){
    /*
    load a csv file into a vector of vectors
    args:
    - filename: std::string, path to csv file
    - data: std::vector<std::vector<float>>, output data
    */
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::vector<float> row;
        std::stringstream ss(line);
        std::string cell;
        int num_cells = 0;
        while (std::getline(ss, cell, delimiter) && num_cells < rowlen)
        {
            try {
                row.push_back(std::stof(cell)); // Convert string to float and add to row
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number format: " << cell << std::endl;
                return false;
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range: " << cell << std::endl;
                return false;
            }
            if (!raw) {
                ++num_cells;
            }
            // ++num_cells;
        }
        data.push_back(row);
    }

    file.close();
    return true;
}

void saveCSV(const std::string& filename, const std::vector<std::vector<float>>& data) {
    /*
    * Save 2d std::vector float to a csv file
    * @param filename: std::string, path to save file
    * @param data: std::vector<std::vector<float>>, data to save
    */


    // Open the file in output mode
    std::ofstream file(filename);
    
    // Check if the file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    // Iterate over the rows
    size_t n_rows = data.size();
    for (size_t i = 0; i < n_rows; ++i) {
        // Iterate over the columns
        size_t n_cols = data[i].size();
        for (size_t j = 0; j < n_cols; ++j) {
            file << std::fixed << std::setprecision(6) << data[i][j]; // Write value
            if (j < n_cols - 1) {
                file << ","; // Separate values with a comma
            }
        }
        file << "\n"; // End of row
    }

    // Close the file
    file.close();
}

void loadCube(std::vector<std::vector<std::vector<float>>>& data, float xs, float xf) {
    /*
    * instantiate cube triangles into data
    * @param data: std::vector<std::vector<std::vector<float>>>, output data
    */
    data = {
        {{xf, 1,-1}, {xf,-1,-1}, {xs, 1,-1}, { 0, 0,-1}},  {{xf,-1,-1},  {xs, 1,-1},  {xs,-1,-1}, {0, 0, -1}}, // base
        {{xs, 1, 1}, {xs,-1, 1}, {xs, 1,-1}, {-1, 0, 0}},  {{xs,-1, 1},  {xs, 1,-1},  {xs,-1,-1}, {-1, 0, 0}}, // front
        {{xf, 1,-1}, {xs, 1,-1}, {xs, 1, 1}, { 0, 1, 0}},  {{xf, 1,-1},  {xs, 1, 1},  {xf, 1, 1}, { 0, 1, 0}}, // left
        {{xf, 1, 1}, {xf,-1, 1}, {xs, 1, 1}, { 0, 0, 1}},  {{xf,-1, 1},  {xs, 1, 1},  {xs,-1, 1}, { 0, 0, 1}}, // top
        {{xf,-1, 1}, {xs,-1, 1}, {xf,-1,-1}, { 0,-1, 0}},  {{xf,-1,-1},  {xs,-1,-1},  {xs,-1, 1}, { 0,-1, 0}}, // right
        {{xf, 1,-1}, {xf,-1,-1}, {xf,-1, 1}, { 1, 0, 0}},  {{xf, 1,-1},  {xf, 1, 1},  {xf,-1, 1}, { 1, 0, 0}}  // back
    };
}

void convertFlatToTriangle(const std::vector<std::vector<float>>& flatData, std::vector<Triangle>& triangles, size_t module_idx) {
    /* convert flat 2d object to triangles assuming each triangle is flattened into x1,y1,z1,x2,y2,z2,x3,y3,z3
    * each triangle also has normal -- [3 points, 9 numbers] [1 normal, 3 numbers]
    */
    for (size_t i = 0; i < flatData.size(); i++) {
        triangles.push_back(Triangle(
            vec3(flatData[i][0], flatData[i][1],  flatData[i][2]),  // v0
            vec3(flatData[i][3], flatData[i][4],  flatData[i][5]),  // v1
            vec3(flatData[i][6], flatData[i][7],  flatData[i][8]),  // v2
            vec3(flatData[i][9], flatData[i][10], flatData[i][11]), // normal
            module_idx                                              // module index
        ));
    }
}

void loadCubeOBS(std::vector<OBS>& obsVec) {
    /*
    * instantiate station into obstacle objects
    * @param obsVec: std::vector<OBS>, output data
    */
    // load in cube data
    std::vector<std::vector<std::vector<float>>> cubeData;
    loadCube(cubeData); // -1 to 1 cube

    // load cube data into triangles
    std::vector<Triangle> triCubeFaces;
    vecToTri(cubeData, triCubeFaces);

    // load cube data into a single obstacle
    obsVec.push_back(OBS(triCubeFaces));
}

void loadConvexStationOBS(std::vector<OBS>& obsVec, float scale) {
    /*
    * instantiate station into obstacle objects
    * @param obsVec: std::vector<OBS>, output data
    */
    size_t num_obstacles = 15;
    std::string model_dir = "../data/model_convex/";

    for (size_t i=0; i < num_obstacles; ++i) {
        // load triangle mesh data
        std::string filename = model_dir + std::to_string(i) + "_faces_normals.csv"; // 9 length vector
        std::vector<std::vector<float>> tri_data;

        // each row is a triangle (3 points = 9 numbers) and normals (3 numbers)
        loadCSV(filename, tri_data, 12,',');

        // convert flat data to triangle objects
        std::vector<Triangle> tris;
        convertFlatToTriangle(tri_data, tris, i);

        vec3 offset = vec3(2.529f, 4.821, 2.591);

        // change offset to reflect center of gravity
        offset -= vec3(2.92199515f, 5.14701097f, 2.63653781f);

        for (size_t j = 0; j < tris.size(); j++) {
            tris[j].n = tris[j].n / tris[j].n.norm();
            tris[j].a += offset;
            tris[j].b += offset;
            tris[j].c += offset;
            tris[j].a *= scale;
            tris[j].b *= scale;
            tris[j].c *= scale;
        }

        OBS obs = OBS(tris);
        obsVec.push_back(obs);
    }
}

void loadStationOBS(std::vector<OBS>& obsVec, float scale) {
    /*
    * instantiate station into obstacle objects
    * @param obsVec: std::vector<OBS>, output data
    */
    size_t num_obstacles = 10;
    std::string model_dir = "../data/model_remeshed/";

    for (size_t i=0; i < num_obstacles; ++i) {
        // load triangle mesh data
        std::string filename = model_dir + std::to_string(i) + "_faces_normals.csv"; // 9 length vector
        std::vector<std::vector<float>> tri_data;

        // each row is a triangle (3 points = 9 numbers) and normals (3 numbers)
        loadCSV(filename, tri_data, 12, ' ');

        // convert flat data to triangle objects
        std::vector<Triangle> tris;
        convertFlatToTriangle(tri_data, tris, i);

        // change offset to reflect center of gravity
        vec3 offset = vec3(2.92199515f, 5.14701097f, 2.63653781f) * -1;

        for (size_t j = 0; j < tris.size(); j++) {
            tris[j].n = tris[j].n / tris[j].n.norm();
            tris[j].a += offset;
            tris[j].b += offset;
            tris[j].c += offset;
            tris[j].a *= scale;
            tris[j].b *= scale;
            tris[j].c *= scale;
        }

        OBS obs = OBS(tris);
        obsVec.push_back(obs);
    }
}

void vecToTri(const std::vector<std::vector<std::vector<float>>>& data, std::vector<Triangle>& tris) {
    /*
    * Convert a vector of vectors of vectors to a vector of triangles
    * @param data: std::vector<std::vector<std::vector<float>>>, input data
    * @param vec: std::vector<Triangles>, output data
    */
    for (size_t i = 0; i < data.size(); ++i) {
        tris.push_back(Triangle(
            vec3(data[i][0][0], data[i][0][1], data[i][0][2]),
            vec3(data[i][1][0], data[i][1][1], data[i][1][2]),
            vec3(data[i][2][0], data[i][2][1], data[i][2][2]),
            vec3(data[i][3][0], data[i][3][1], data[i][3][2])
        ));
    }
}

bool allTrue(const std::vector<TriangleCoverage>& arr, size_t module_idx) {
    /*
    * Check if all elements.covered in an array are true
    * @param arr: const std::vector<TriangleCoverage>&, input array
    * @return bool, true if all elements are true
    */
    for (size_t i = 0; i < arr.size(); ++i) {
        if (!arr[i].covered && arr[i].module_idx == module_idx) {
            return false;
        }
    }
    return true;
}

bool allTrue(const std::vector<TriangleCoverage>& arr) {
    /*
    * Check if all elements.covered in an array are true
    * @param arr: const std::vector<TriangleCoverage>&, input array
    * @return bool, true if all elements are true
    */
    for (size_t i = 0; i < arr.size(); ++i) {
        if (!arr[i].covered) {
            return false;
        }
    }
    return true;
}

bool allTrue(const std::vector<bool>& arr) {
    /*
    * Check if all elements in an array are true
    * @param arr: const bool*, input array
    * @param len: size_t, length of array
    * @return bool, true if all elements are true
    */
    std::cout << "Checking all true, no moduleidx" << std::endl;
    for (size_t i = 0; i < arr.size(); ++i) {
        if (!arr[i]) {
            return false; 
        }
    }
    return true;
}

void numTrue(const std::vector<bool>& arr, size_t& num_true) {
    /*
    * Count the number of true elements in an array
    * @param arr: const bool*, input array
    * @param len: size_t, length of array
    * @param num_true: size_t&, output number of true elements
    */
    num_true = 0;
    for (size_t i = 0; i < arr.size(); ++i) {
        if (arr[i]) { 
            ++num_true; 
        }
    }
}

bool allZeroGain(const std::vector<VPCoverageGain>& arr) {
    /*
    * Check if all elements.gain in an array are zero
    * @param arr: const std::vector<VPCoverageGain>&, input array
    * @param all_zero: bool&, output true if all elements are zero
    */
    for (size_t i = 0; i < arr.size(); ++i) {
        if (arr[i].gain > 0) {
            return false;
        }
    }
    return true;
}

void getCoverage(const std::vector<Viewpoint>& viewpoints, const std::vector<Triangle*>& triangles, std::vector<std::vector<bool>>& coverage_map) {
    /*
    * Get the coverage map for a set of viewpoints and triangles
    * @param viewpoints: const std::vector<Viewpoint>&, input viewpoints
    * @param triangles: const std::vector<Triangle>&, input triangles
    * @param coverage_map: bool**, output coverage map
    */
    std::ostringstream message;
    auto now = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < viewpoints.size(); ++i) {
        std::vector<bool> coverage;
        cuda_kernel_coverage(viewpoints[i], triangles, coverage, nullptr);
        for (size_t j = 0; j < coverage.size(); j++) {
            coverage[j] = !coverage[j];
        }
        coverage_map.push_back(coverage);
        auto prev_now = now;
        now = std::chrono::high_resolution_clock::now();
        // Calculate the duration
        std::chrono::duration<double> duration = now - prev_now;
        double loop_period = duration.count();

        // Output the duration in seconds
        double seconds_remaining = loop_period * (viewpoints.size() - i);
        double minutes_remaining = seconds_remaining / 60.0;
        seconds_remaining = std::fmod(seconds_remaining, 60.0);
        message << " Time remaining: " << int(minutes_remaining) << "m " << int(seconds_remaining) << "s";
        displayProgressBar(static_cast<double>(i) / viewpoints.size(), 150, message);
        message.str("");
    }
}

void displayProgressBar(double progress, int width, std::ostringstream& message) {
    // Clear the current line
    std::cout << "\r";

    // Calculate the number of '#' characters
    int pos = static_cast<int>(width * progress);
    
    // Draw the progress bar
    std::cout << "[";
    for (int i = 0; i < width; ++i) {
        if (i < pos) 
            std::cout << "#";
        else 
            std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(1) << (progress * 100) << "%";

    std::cout << message.str();
    
    // Flush the output to ensure it updates immediately
    std::cout.flush();
}

void getIncidenceAngle(vec3 viewdir, Triangle tri, float& angle) {
    /*
    * Get the incidence angle between a view direction and a triangle (always positive and between 0 and pi)
    * @param viewdir: vec3, input view direction
    * @param tri: Triangle, input triangle
    * @param angle: float&, output incidence angle
    */
    vec3 normal = tri.n;
    angle = acosf(fabsf(viewdir.dot(normal)) / (viewdir.norm() * normal.norm() + 1e-9f));
}

void pinhole_camera_test(
    bool& visible, 
    vec3 pose, 
    vec3 viewdir, 
    vec3 point,
    float hfov, // rad
    float vfov // rad
    ) {
    // calculate the angle between the view direction and the vector from the viewpoint to the intersection point
    vec3 vec = point - pose;
    float norm_dot = vec.dot(viewdir);

    // check if point is behind camera
    if (norm_dot <= 0) {
        visible = false;
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
        visible = true;
        return;
    }

    visible = false;
    return;
}

void cw_acceleration(vec3& acceleration, vec3 point, vec3 velocity) {
    // compute cost to oppose CW disturbance from start to end at speed.
    // n = number of discretization steps to take -- inclusive

    // this->speed;
    float mu = 3.986e14; // gravitational parameter
    float a = 6.6e3; // semi-major axis
    float n = std::sqrt(mu / (a * a * a)); // orbital rate 
    float g0 = 9.80665; // standard gravity
    float m = 5; // mass of the spacecraft - 5 kg
    float Isp = 80; // specific impulse

    std::vector<vec3> displacement;

    // compute drift acceleration
    acceleration.x = 3 * n * n * point.x + 2 * n * velocity.y;
    acceleration.y = -2 * n * velocity.x;
    acceleration.z = -1 * n * n * point.z;

    // // compute cost
    // float impulse = 0;
    // float dt = (end - start).norm() / speed / (N - 1);
    // for (size_t i = 0; i < N; i++) {
    //     impulse += acceleration[i].norm() * std::sqrt(dt);
    // }

    // float cost = impulse * m / (g0 * Isp);

    // return cost;
}

float cw_cost(vec3 start, vec3 end, float speed, size_t N) {
    // compute cost to oppose CW disturbance from start to end at speed.
    // n = number of discretization steps to take -- inclusive

    float mu = 3.986e14; // gravitational parameter
    float a = 6.6e3; // semi-major axis
    float n = std::sqrt(mu / (a * a * a)); // orbital rate 
    float g0 = 9.80665; // standard gravity
    float m = 5; // mass of the spacecraft - 5 kg
    float Isp = 80; // specific impulse

    std::vector<vec3> displacement;
    vec3 velocity = (end - start) / (end - start).norm() * speed;
    std::vector<vec3> acceleration;

    // calculate the displacement vector
    for (size_t i = 0; i < N; i++) {
        float t = static_cast<float>(i) / static_cast<float>(N);
        vec3 point = start * (1 - t) + end * t;
        displacement.push_back(point);
    }

    // compute drift acceleration
    for (size_t i = 0; i < N; i++) {
        vec3 acc;
        acc.x = 3 * n * n * displacement[i].x + 2 * n * velocity.y;
        acc.y = -2 * n * velocity.x;
        acc.z = -1 * n * n * displacement[i].z;
        acceleration.push_back(acc);
    }

    // compute cost
    float impulse = 0;
    float dt = (end - start).norm() / speed / (N);
    for (size_t i = 0; i < N; i++) {
        impulse += acceleration[i].norm() * std::sqrt(dt);
    }

    float cost = impulse * m / (g0 * Isp);

    return cost;
}

float fuel_cost(vec3 pose, vec3 v0, vec3 v1, float speed, float dt) {
    vec3 cw_acc;
    cw_acceleration(cw_acc, pose, v0);

    vec3 dv = v1 / v1.norm() - v0 / v0.norm() - cw_acc * dt; // change in velocity
    float mf = 5;
    float Isp = 80;
    float m0 = std::exp(dv.norm() / (9.81 * Isp));
    return m0 - mf;
}