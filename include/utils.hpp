#ifndef UTILS_HPP
#define UTILS_HPP

#include "node3d_struct.hpp"
#include "obs.hpp"
#include "triangle_struct.hpp"
#include "vec3_struct.hpp"
#include "viewpoint_struct.hpp"
#include "viewpoint_coverage_gain_struct.hpp"

#include <string>
#include <vector>

static vec3* TEMP = new vec3();

bool ray_int_plane(Node3D node, Plane plane, float eps=1e-5f, vec3& intPoint= * TEMP);
bool ray_int_triangle(vec3 origin, vec3 vector, vec3 end, Triangle tri, vec3& intPoint= * TEMP, float eps=1e-5f);
bool ray_int_triangle(vec3 origin, vec3 vector, Triangle tri, vec3& intPoint= * TEMP, float eps=1e-5f);
bool ray_int_triangle(Node3D node, Triangle tri, vec3& intPoint= * TEMP, float eps=1e-5f);
float heading_change(Node3D node, Node3D next_node);
float heading_change(Node3D node, vec3 vector);
float heading_change(vec3 v0, vec3 v1);
bool loadCSV(const std::string& filename, std::vector<std::vector<float>>& data, int rowlen, char delimiter=',', bool raw=false);
void saveCSV(const std::string& filename, const std::vector<std::vector<float>>& data);
void loadCube(std::vector<std::vector<std::vector<float>>>& data, float xs=-1, float xf=1);
void convertFlatToTriangle(const std::vector<std::vector<float>>& flatData, std::vector<Triangle>& tris, size_t module_idx=0);
void loadCubeOBS(std::vector<OBS>& obsVec);
void loadConvexStationOBS(std::vector<OBS>& obsVec);
void loadStationOBS(std::vector<OBS>& obsVec);
void vecToTri(const std::vector<std::vector<std::vector<float>>>& data, std::vector<Triangle>& tris);
bool allTrue(const std::vector<bool>& arr);
void numTrue(const std::vector<bool>& arr, size_t& num_true);
bool allZeroGain(const std::vector<VP_Coverage_Gain>& arr);
void getCoverage(const std::vector<Viewpoint>& viewpoints, const std::vector<Triangle*>& triangles, std::vector<std::vector<bool>>& coverage_map);
void displayProgressBar(double progress, int width, std::ostringstream& message);
void getIncidenceAngle(vec3 viewdir, Triangle tri, float& angle);
void pinhole_camera_test( bool& visible, vec3 pose, vec3 viewdir, vec3 point, float hfov, float vfov);

#endif // UTILS_HPP