#pragma once
#include "geometry.h"
#include "Material.hpp"
#include <limits>


#pragma ide diagnostic ignored "openmp-use-default-none"


class Object {
public:
    virtual bool ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const = 0;
    virtual Material get_material(Vec3f hit) const = 0;
};

class Sphere: public Object {
    Vec3f center;
    float radius;
    Material material;
    std::vector<Vec3f>envmap;

public:
    Sphere(const Vec3f &c, const float &r, const Material &m);
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const override;
    Material get_material(Vec3f hit) const override;
};

class Panorama: public Object {
    int envmap_width, envmap_height;
    std::vector<Vec3f> envmap;
    Vec3f center;
    float radius;
    Material material;

public:
    void open_img(const char * img);
    Panorama(const Vec3f &c, const float &r, const Material &m, const char * img);
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const override;
    Material get_material(Vec3f hit) const override;
};

bool plane_intersect(const Vec3f &orig, const Vec3f &dir, const Vec3f &norm, const float &d, Vec3f &p);

class Plane: public Object {
    float y1_lim;
    float x1_lim;
    float x2_lim;
    float z1_lim;
    float z2_lim;
    Material first;
    Material second;

public:
    Plane(const float &y1, const float &x1, const float &x2, const float &z1, const float &z2,
          const Material &first, const Material &second);
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const override;
    Material get_material(Vec3f hit) const override;
};

class PlaneCircle: public Object {
    Vec3f center;
    float r;
    Material material;
    int envmap_width, envmap_height;
    std::vector<Vec3f> envmap;

public:
    void open_img(const char * img);
    PlaneCircle(const Vec3f &c, const float &r, const Material &m, const char * img);
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const override;
    Material get_material(Vec3f hit) const override;
};

class Cylinder: public Object {
    Vec3f center;
    float radius;
    Vec3f vec;
    float height;
    Material material;

public:
    Cylinder(const Vec3f &c, const float &r, const Vec3f &v, const float &h, const Material &m);
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const override;
    Material get_material(Vec3f hit) const override;
};

class Cube: public Object {
    Vec3f center;
    Vec3f vec1;
    Vec3f vec2;
    float a;
    float b;
    float h;
    Material material;

public:
    Cube(const Vec3f &c, const Vec3f &v1, const Vec3f &v2, const float &a, const float &b,\
    const float &h, const Material &m);
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const override;
    Material get_material(Vec3f hit) const override;
};
