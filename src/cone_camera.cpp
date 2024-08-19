#include "cone_camera.hpp"
#include <cmath>

ConeCamera::ConeCamera() {
    this->fov = M_PI/4; // radians - 45 deg
}

ConeCamera::ConeCamera(float fov) {
    this->fov = fov * M_PI / 180; // input in degrees
}

bool ConeCamera::cover_point(Viewpoint vp, vec3 point) {
    vec3 vp_to_point = point - vp.pose;
    float angle = acos(vp_to_point.dot(vp.viewdir) / vp_to_point.norm());
    return angle < this->fov;
}