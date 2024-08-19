#ifndef UTILS_HPP
#define UTILS_HPP

#include "node3d_struct.hpp"
#include "obs.hpp"
#include "triangle_struct.hpp"
#include "vec3_struct.hpp"

#include <string>
#include <vector>

static vec3* TEMP = new vec3();

bool ray_int_plane(Node3D node, Plane plane, float eps=1e-9f, vec3& intPoint= * TEMP);
bool ray_int_triangle(vec3 origin, vec3 vector, vec3 end, Triangle tri, vec3& intPoint= * TEMP, float eps=1e-9f);
bool ray_int_triangle(vec3 origin, vec3 vector, Triangle tri, vec3& intPoint= * TEMP, float eps=1e-9f);
bool ray_int_triangle(Node3D node, Triangle tri, vec3& intPoint= * TEMP, float eps=1e-9f);
float heading_change(Node3D node, Node3D next_node);
float heading_change(Node3D node, vec3 vector);
bool loadCSV(const std::string& filename, std::vector<std::vector<float>>& data, int rowlen);
void saveCSV(const std::string& filename, const std::vector<std::vector<float>>& data);
void loadCube(std::vector<std::vector<std::vector<float>>>& data, float xs=-1, float xf=1);
void convertFlatToTriangle(const std::vector<std::vector<float>>& flatData, std::vector<Triangle>& tris);
void loadStationOBS(std::vector<OBS>& obsVec);
void vecToTri(const std::vector<std::vector<std::vector<float>>>& data, std::vector<Triangle>& tris);
bool allTrue(const bool* arr, size_t len);
void numTrue(const bool* arr, size_t len, size_t& num_true);

#endif // UTILS_HPP