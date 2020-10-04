#pragma once
#include "geometry.h"


class Material {
    Vec3f albedo;
    Vec3f diffuse_color;
    float specular_exponent;

public:
    Vec3f get_albedo() const;
    Vec3f get_diffuse_color() const;
    float get_specular_exponent() const;
    Material(const Vec3f &a, const Vec3f &color, const float &spec);
    Material();
};
