#ifndef CONE_CAMERA_HPP
#define CONE_CAMERA_HPP

#include "vec3_struct.hpp"
#include "viewpoint_struct.hpp"

class ConeCamera {
    public:

        ConeCamera();
        ConeCamera(float fov);
        bool cover_point(Viewpoint vp, vec3 point); // check if a point is covered by the camera
        float get_fov();

    private:

        float fov;
};

#endif // CONE_CAMERA_HPP