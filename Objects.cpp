#include "Objects.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

Sphere::Sphere(const Vec3f &c, const float &r, const Material &m): center(c), radius(r), material(m) {
    envmap.push_back(Vec3f(1., 1., 1.));
    envmap.push_back(Vec3f(1., 1., 1.));
}

bool Sphere::ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const {
    Vec3f L = center - orig;
    float tca = L * dir; //the center projection on a ray (>=0, if exists
    // real projection, <0, if imagine; |tca| = distance to projection)
    float d2 = L * L - tca * tca; //distance from the center to the point of the projection
    if (d2 > radius * radius) return false; //if distance from the center to the point of the projection
    // then ray doesn't intersect
    float thc = sqrtf(radius * radius - d2); //distance from the point of the projection to the sphere
    float t = tca - thc;
    float t1 = tca + thc;
    if (t < 0) {
        t = t1;
    }
    hit = orig + dir * t;
    N = (hit - center).normalize();
    return t >= 0;
}

Material Sphere::get_material(Vec3f hit) const {
    return material;
}

void Panorama::open_img(const char * img) {
    int n = -1;
    unsigned char *pixmap = stbi_load(img, &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    #pragma omp parallel for
    for (int j=envmap_height-1; j>=0; j--) {
    #pragma omp parallel for
        for (int i=0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);
}

Panorama::Panorama(const Vec3f &c, const float &r, const Material &m, const char * img): center(c), radius(r), material(m){
    open_img(img);
}

bool Panorama::ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const {
    Vec3f L = center - orig;
    float tca = L * dir; //the center projection on a ray (>=0, if exists
    // real projection, <0, if imagine; |tca| = distance to projection)
    float d2 = L * L - tca * tca; //distance from the center to the point of the projection
    if (d2 > radius * radius) return false; //if distance from the center to the point of the projection
    // then ray doesn't intersect
    float thc = sqrtf(radius * radius - d2); //distance from the point of the projection to the sphere
    float t = tca - thc;
    float t1 = tca + thc;
    if (t < 0) {
        t = t1;
    }
    hit = orig + dir * t;
    N = (center - hit).normalize();
    return t >= 0;
}


Material Panorama::get_material(Vec3f hit) const {
    float phi;
    float teta;
    Vec3f hit_c(hit - center);
    if (hit_c.x == 0 && hit_c.z > 0) phi = 0;
    else if (hit_c.x > 0) phi = M_PI / 2. - atan(hit_c.z / hit_c.x);
    else if (hit_c.x == 0 && hit_c.z < 0) phi = M_PI;
    else if (hit_c.x < 0) phi = 3 * M_PI / 2. + atan(hit_c.z / -hit_c.x);
    if (hit_c.z == 0 && hit_c.x == 0 && hit_c.y > 0) teta = M_PI / 2.;
    else if (fabs(hit_c.z) < 1e-2 && fabs(hit_c.x) < 1e-2 && hit_c.y < 0) teta = -M_PI / 2.;
    else teta = atan(hit_c.y / sqrtf(hit_c.z * hit_c.z + hit_c.x * hit_c.x));
    int i = floor(phi / (2 * M_PI) * (envmap_width - 1));
    int j = floor((M_PI / 2. - teta) / (M_PI) * (envmap_height - 1));
    return Material(material.get_albedo(), envmap[i+j*envmap_width], material.get_specular_exponent());
}

bool plane_intersect(const Vec3f &orig, const Vec3f &dir, const Vec3f &norm, const float &d, Vec3f &p) {
    float c = dir * norm;
    if (c == 0) {
        return false;
    }
    float alpha = (d - norm * orig) / c;
    if (alpha < 0) {
        return false;
    }
    else {
        p = orig + dir * alpha;
        return true;
    }
}

Plane::Plane(const float &y1, const float &x1, const float &x2, const float &z1, const float &z2,
             const Material &first, const Material &second):
             y1_lim(y1),
             x1_lim(x1),
             x2_lim(x2),
             z1_lim(z1),
             z2_lim(z2),
             first(first),
             second(second) {}

bool Plane::ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const {
    if (fabs(dir.y) > 1e-3)  {
        float dist = -(orig.y - y1_lim) / dir.y;
        hit = orig + dir * dist;
        if (dist > 0 && hit.x > x1_lim && hit.x < x2_lim && hit.z > z1_lim && hit.z < z2_lim) {
            N = Vec3f(0., 1., 0.);
            return true;
        }
        return false;
    }
    return false;
}

Material Plane::get_material(Vec3f hit) const {
    Material material = int(1.*hit.x+1000) & 1 ? first : second;
    return material;
}

void PlaneCircle::open_img(const char * img) {
    int n = -1;
    unsigned char *pixmap = stbi_load(img, &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    #pragma omp parallel for
    for (int j=envmap_height-1; j>=0; j--) {
    #pragma omp parallel for
        for (int i=0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);
}

PlaneCircle::PlaneCircle(const Vec3f &c, const float &r, const Material &m, const char * img): center(c), r(r), material(m){
    open_img(img);
}

bool PlaneCircle::ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const {
    if (fabs(dir.y) > 1e-3) {
        float d = -(orig.y - center.y) / dir.y;
        hit = orig + dir * d;
        Vec3f pt_c = hit - center;
        if (d > 0 && (pow(pt_c.x, 2) + pow(pt_c.z, 2)) <= pow(r, 2)) {
            N = Vec3f(0., 1., 0.);
            return true;
        }
        return false;
    }
    return false;
}

Material PlaneCircle::get_material(Vec3f hit) const {
    float min_len = std::min(envmap_height, envmap_width);
    Vec3f pt_c = hit - center;
    int i = floor((pt_c.x / r) * (min_len / 2.) + (envmap_width / 2.));
    int j = floor((pt_c.z / r) * (min_len / 2.) + (envmap_height / 2.));
    return Material(material.get_albedo(), envmap[i+j*envmap_width], material.get_specular_exponent());
}

Cylinder::Cylinder(const Vec3f &c, const float &r, const Vec3f &v, const float &h, const Material &m) : \
    center(c), radius(r), vec(v), height(h), material(m) {}

bool Cylinder::ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const {
    Vec3f c1 = (center - vec * height * 0.5);
    Vec3f c2 = (center + vec * height * 0.5);

    float c = dir * vec;
    if (c == 0) {
        Vec3f pt;
        if (plane_intersect(orig - dir * radius, dir, dir, center * dir, pt)) {
            float r = ((pt - center) - vec * ((pt - center) * vec)).norm();
            if (abs((pt - center) * vec) <= height / 2. && r <= radius) {
                float tmp = 0.;
                if (abs(tmp = (center - orig) * dir <= radius)) hit = orig + dir * (sqrt(radius * radius - r * r) + tmp);
                else hit = orig + dir * (tmp - sqrtf(radius * radius - r * r));
                N = (hit - center).normalize();
                return true;
            } else return false;
        } else return false;
    }
    Vec3f p1, p2;
    float h = 0;
    bool pln1 = plane_intersect(orig, dir, vec, c1 * vec, p1), pln2 = plane_intersect(orig, dir, vec, c2 * vec, p2);
    if(!pln1 && !pln2) return false;
    if(pln1 && !pln2) {
        p2 = p1;
        p1 = orig;
        c2 = c1;
    }
    if(!pln1 && pln2) {
        p1 = orig;
    }
    if (pln1 && pln2) {
        if (vec * dir < 0) {
            Vec3f tmp = p2;
            p2 = p1;
            p1 = tmp;
            c2 = c1;
        }
    }
    Vec3f vec1 = vec;
    if (vec * dir < 0) vec1 = -vec;
    h = (p2 - p1) * vec1;
    c1 = c2 - (vec1 * h);
    if ((p1 - c1).norm() <= radius) {
        if (h - height <= 1e-4) {
            hit = p1;
            N = -vec1;
        } else if ((p2 - c2).norm() < sqrtf(radius * radius)) {
            hit = p2;
            N = vec1;
        } else {

        }
        return true;
    }
    else {
        float t = 0;
        Vec3f p2pr = p2 - vec1 * h;
        Vec3f tmp_vec = (p2pr - p1).normalize();
        Vec3f L = c1 - p1;
        float tca = L * tmp_vec;
        float d2 = L * L - tca * tca;
        if (d2 > radius * radius) return false;
        float thc = sqrtf(radius * radius - d2);
        t = tca - thc;
        float t1 = tca + thc;
        if (t < 0) t = t1;
        if (t < 0 || t > (p2pr - p1).norm()) return false;
        hit = p1 + dir * (t / (tmp_vec * dir));
        N = (p1 + tmp_vec * t - c1).normalize();
        return true;
    }
}

Material Cylinder::get_material(Vec3f hit) const {
    return material;
}

Cube::Cube(const Vec3f &c, const Vec3f &v1, const Vec3f &v2, const float &a, const float &b,\
    const float &h, const Material &m): center(c), vec1(v1), a(a), b(b), h(h), vec2(v2), material(m){}

bool Cube::ray_intersect(const Vec3f &orig, const Vec3f &dir, Vec3f &hit, Vec3f &N) const {
    Vec3f vec12 = -vec1;
    Vec3f vec22 = -vec2;
    Vec3f vec3 = cross(vec1, vec2).normalize();
    Vec3f vec32 = - vec3;

    Vec3f L = orig - center;
    bool in = fabs(L * vec1) < h && fabs(L * vec2) < a && fabs(L * vec3) < b;
    bool intersect = false;
    Vec3f p;
    float dist = std::numeric_limits<float>::max();
    Vec3f vecs[6] = {vec1, vec12, vec2, vec22, vec3, vec32};
    float l = 0;
    for (int i = 0; i < 6; i++) {
        switch(i / 2) {
            case 0: l = h; break;
            case 1: l = a; break;
            case 2: l = b; break;
        }
        Vec3f c_proec = center + vecs[i] * l;
        //if ray intersects one of planes and intersection point belongs to corresponding parallelepiped face return true
        if (plane_intersect(orig, dir, vecs[i], c_proec * vecs[i], p) && \
            fabs((p - c_proec) * vec1) <= h && fabs((p - c_proec) * vec2) <= a && \
            fabs((p - c_proec) * vec3) <= b) {
            intersect = true;
            if ((p - orig).norm() < dist) {
                dist = (p - orig).norm();
                if (in) N = -vecs[i];
                else N = vecs[i];
                hit = p;
            }
        }
    }
    return intersect;
}

Material Cube::get_material(Vec3f hit) const {
    return material;
}