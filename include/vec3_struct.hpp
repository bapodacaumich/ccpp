#ifndef VEC3_STRUCT_HPP
#define VEC3_STRUCT_HPP

#include <math.h>
#include <string>

struct vec3 {
    float x, y, z;

    vec3 () : x(0.0f), y(0.0f), z(0.0f) {}

    vec3 (const vec3& v) : x(v.x), y(v.y), z(v.z) {}

    vec3 (float a, float b, float c) : x(a), y(b), z(c) {}

    void set(float a, float b, float c) {
        this->x = a;
        this->y = b;
        this->z = c;
    }
    std::string toString() {
        return "(" + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + ")";
    }
    void normalize() {
        float n = norm();
        x /= n;
        y /= n;
        z /= n;
    }
    inline vec3 operator+(const vec3& v) const {
        return vec3(x + v.x, y + v.y, z + v.z);
    }
    inline vec3 operator-(const vec3& v) const {
        return vec3(x - v.x, y - v.y, z - v.z);
    }
    inline vec3 operator*(const float& f) const {
        return vec3(x * f, y * f, z * f);
    }
    inline vec3 operator/(const float& f) const {
        return vec3(x / f, y / f, z / f);
    }
    inline vec3 operator+=(const vec3& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    inline vec3 operator-=(const vec3& v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    inline float dot(const vec3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }
    inline vec3 cross(const vec3& v) const {
        return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    inline float norm() const {
        return sqrtf(x * x + y * y + z * z);
    }
};

#endif // VEC3_STRUCT_HPP