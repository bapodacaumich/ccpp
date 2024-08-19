#ifndef TRIANGLE_STRUCT_HPP
#define TRIANGLE_STRUCT_HPP

#include "vec3_struct.hpp"

#include <iostream>
#include <string>

struct Triangle {
    vec3 a, b, c, n;
    Triangle () : a(), b(), c(), n() {}

    Triangle (const Triangle& t) : a(t.a), b(t.b), c(t.c), n(t.n) {}

    Triangle (vec3 a, vec3 b, vec3 c) : a(a), b(b), c(c){
        this->n = (b - a).cross(c - a);
        (this->n).normalize();
    }

    Triangle (vec3 a, vec3 b, vec3 c, vec3 n): a(a), b(b), c(c), n(n) {
        (this->n).normalize();
    }

    vec3 getCentroid() {
        return (this->a + this->b + this->c) / 3.0f;
    }

    std::string toString() {
        return "Triangle (v0, v1, v2, n): " + 
                this->a.toString() + ", " + 
                this->b.toString() + ", " + 
                this->c.toString() + ", " + 
                this->n.toString();
    }
};

#endif // TRIANGLE_STRUCT_HPP