#include "Material.hpp"

Material::Material(const Vec3f &a, const Vec3f &color, const float &spec) : albedo(a), diffuse_color(color), specular_exponent(spec) {}

Material::Material(): albedo(1,0,0), diffuse_color(), specular_exponent(){}

Vec3f Material::get_albedo() const {
    return albedo;
}
Vec3f Material::get_diffuse_color() const {
    return diffuse_color;
}
float Material::get_specular_exponent() const {
    return specular_exponent;
}