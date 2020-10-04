#include "Lights.hpp"

Light::Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}

Vec3f Light::get_position() const {
    return position;
}

float Light::get_intensity() const {
    return intensity;
}