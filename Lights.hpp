#pragma once
#include "geometry.h"

class Light {
    Vec3f position;
    float intensity;

public:
    Light(const Vec3f &p, const float &i);
    Vec3f get_position() const;
    float get_intensity() const;
};