#ifndef VIEWPOINT_HPP
#define VIEWPOINT_HPP

#include "vec3_struct.hpp"

struct Viewpoint {
    vec3 pose;
    vec3 viewdir; // always unit vec

    Viewpoint() : pose(0.0f, 0.0f, 0.0f), viewdir(1.0f, 0.0f, 0.0f) { this->viewdir.normalize(); }

    Viewpoint(vec3 p) : pose(p), viewdir(1.0f, 0.0f, 0.0f) { viewdir.normalize(); }

    Viewpoint(vec3 p, vec3 v) : pose(p), viewdir(v) { viewdir.normalize(); }
};

#endif