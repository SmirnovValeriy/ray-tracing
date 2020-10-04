#pragma ide diagnostic ignored "openmp-use-default-none"

#include <iostream>
#include <cstdio>
#include <limits>
#include <fstream>
#include <vector>
#include <cstring>
#include "geometry.h"
#include <omp.h>
#include "Lights.hpp"
#include "Material.hpp"
#include "Objects.hpp"
#include "Time.hpp"
#include "stb_image_write.h"
#include "stb_image.h"

//macros
#define BCKGRND_COLOR Vec3f(0.9, 0.9, 0.9)
#define X1_LIMIT -100
#define X2_LIMIT 100
#define Y1_LIMIT -4
#define Z1_LIMIT -100
#define Z2_LIMIT 100
#define PAN_RAD 8

using namespace std;


//parameters of the scenes
int width = 720;
int height = 720;
double fov = 3. * M_PI / 4.;
int scene = 2;


//reflect function
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N*2.f*(I*N);
}

//function that checks if the ray with given origin and direction intersects some objects of the scene
bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const vector<Object*> &objects, \
Vec3f &hit, Vec3f &N, Material &material) {
    //check if intersects spheres
    float dist = numeric_limits<float>::max();
    for(size_t i=0; i < objects.size(); i++) {
        float dist_i;
        Vec3f hit0;
        Vec3f N0;
        if (objects[i]->ray_intersect(orig, dir, hit0, N0) && (dist_i = (hit0-orig).norm()) < dist) {
            dist = dist_i;
            hit = hit0;
            N = N0;
            material = objects[i]->get_material(hit0);
        }
    }
    if (scene == 1) return dist < 1000;
    return true;
}

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const vector<Object*> &objects,\
const vector<Light> &lights, size_t depth=0) {
    Vec3f point, N;
    Material material;

    //check if our ray intersects some objects
    //depth of reflections = 3
    if (depth > 3 || !scene_intersect(orig, dir, objects, point, N, material)) {
        return BCKGRND_COLOR; // background color
    }

    //get reflection of other objects
    Vec3f reflect_dir = reflect(dir, N).normalize();
    Vec3f reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    Vec3f reflect_color = cast_ray(reflect_orig, reflect_dir, objects, lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;

    //getting lighted areas and shadows
    for (size_t i = 0; i < lights.size(); i++) {
        Vec3f light_dir = (lights[i].get_position() - point).normalize();
        float light_distance = (lights[i].get_position() - point).norm();

        Vec3f shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
        Vec3f shadow_pt, shadow_N;
        Material tmp;
        if (scene_intersect(shadow_orig, light_dir, objects, shadow_pt, shadow_N, tmp) &&
            (shadow_pt - shadow_orig).norm() < light_distance)
            continue;
        diffuse_light_intensity += lights[i].get_intensity() * max(0.f, light_dir * N);
        specular_light_intensity +=
                powf(std::max(0.f, -reflect(-light_dir, N) * dir),
                     material.get_specular_exponent()) * lights[i].get_intensity();
    }
    return material.get_diffuse_color() * diffuse_light_intensity * material.get_albedo()[0] +
           Vec3f(1., 1., 1.) * specular_light_intensity * material.get_albedo()[1] + reflect_color * material.get_albedo()[2];
}

//rendering the scene function
void getImg(const vector<Object*> &objects, const vector<Light> &lights, const char img_file[], int threads) {

    vector<Vec3f> framebuffer(width * height);

    //matrix for getting rid of stepping
    float matrix[4 * 2] = {
            0.5, -0.25,
            -0.25, 0.5,
            1.25, 0.5,
            0.5, 1.25
    };

    //parallelization
    //getting rid of stepping (only for scene #1)
    if (scene == 1) {
        #pragma omp parallel for num_threads(threads)
        for (size_t j = 0; j < height; j++) {
            #pragma omp parallel for num_threads(threads)
            for (size_t i = 0; i < width; i++) {
                framebuffer[i + j * width] = Vec3f(0., 0., 0.);
                for (size_t sample = 0; sample < 4; sample++) {
                    float x = (2 * (i + matrix[2 * sample]) / (float) width - 1) * tan(fov / 2.) * width /
                            (float) height;
                    float y = -(2 * (j + matrix[2 * sample + 1]) / (float) height - 1) * tan(fov / 2.);
                    Vec3f dir = Vec3f(x, y, -1).normalize(); //unit vector passing this pixel
                    framebuffer[i + j * width] = framebuffer[i + j * width] + \
            cast_ray(Vec3f(0, 0, 0), dir, objects, lights);
                }
                framebuffer[i + j * width] = framebuffer[i + j * width] * 0.2;
            }
        }
    }
    // without getting rid of stepping (for scene #2)
    else {
        #pragma omp parallel for num_threads(threads)
        for (size_t j = 0; j < height; j++) {
        #pragma omp parallel for num_threads(threads)
            for (size_t i = 0; i < width; i++) {
                float x = (2 * (i + 0.5) / (float) width - 1) * tan(fov / 2.) * width / (float) height;
                float y = -(2 * (j + 0.5) / (float) height - 1) * tan(fov / 2.);
                Vec3f dir = Vec3f(x, y, -1).normalize(); //unit vector passing this pixel
                framebuffer[i + j * width] = cast_ray(Vec3f(0, 0, 0), dir, objects, lights);
            }
        }
    }


    //using stb library to save the scene in .jpg format
    std::vector<unsigned char> pixmap(width * height * 3);
    #pragma omp parallel for
    for (size_t i = 0; i < height * width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; j++) {
            pixmap[i * 3 + j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    stbi_write_jpg(img_file, width, height, 3, pixmap.data(), 100);
}

int main(int argc, char *argv[]) {
    char * img_file; //filename for image to write in
    strcpy(img_file = (char *)malloc(12), "scene.jpg");
    int threads = 8; //default number of threads
    scene = 1; //default scene

    Timer Timer; //execution time control

    for(int i = 0; i < argc; i++) {
        if(!strcmp(argv[i], "-out")) {
            strcpy(img_file = (char *)malloc(strlen(argv[i + 1])), argv[i + 1]);
            // if argument -out is given get filename
        }
        if(!strcmp(argv[i], "-threads")) {
            threads = std::max(threads, atoi(argv[i + 1]));
            //if -threads argument is given get threads number and check if it is greater than default
        }
        if(!strcmp(argv[i], "-scene")) {
            scene = atoi(argv[i + 1]);
            //if argument -scene is given get scene number
        }
    }

    //push objects to a scene
    vector<Object*> objects;

    //push lights to a scene
    vector<Light> lights;

    //set of materials
    Material plane(Vec3f(0.9,  0.1, 0.2), Vec3f(0.45, 0.15, 0.15), 5.);
    Material gray_asphalt(Vec3f(0.9, 0.9, 0.4), Vec3f(0.8, 0.8, 0.8) * 0.3, 10000.);
    Material gray(Vec3f(0.9, 0.9, 0.4), Vec3f(0.8, 0.8, 0.8), 10000.);
    Material wet_asphalt(Vec3f(1., 0.2, 0.), Vec3f(80., 80., 80.) * (1. / 255.) * 0.3, 5.);
    Material yellow_taxi(Vec3f(0.9,  0.9, 0.05), Vec3f(240., 145., 7.) * (1. / 255.), 20.);
    Material red_phone(Vec3f(0.6,  0.1, 0.1), Vec3f(165., 39., 27.) * (1./ 255.), 10.);
    Material big_ben(Vec3f(0.9, 0.1, 0.01), Vec3f(100., 92., 78.)*(1. / 255.), 5.);
    Material green_taxi(Vec3f(0.9,  0.9, 0.05), Vec3f(18., 61., 43.) * (1. / 255.), 20.);
    Material blue(Vec3f(0.9, 0.9, 0.4), Vec3f(0.0, 0.2, 0.8), 10000.);
    Material red(Vec3f(0.9, 0.9, 0.4), Vec3f(0.9, 0.2, 0.1), 10000.);

    if (scene == 1) {
        //London scene
        height = 540;
        width = 720;
        fov = M_PI / 2.;
        Sphere sphere1(Vec3f(-3.0, Y1_LIMIT + 1.4, -14), 1.4, yellow_taxi);
        Sphere sphere2(Vec3f(-0.8, Y1_LIMIT + 1.4, -20), 1.4, green_taxi);
        objects.push_back(&sphere1);
        objects.push_back(&sphere2);

        Vec3f center3(7.0, Y1_LIMIT + 6.5, -25);
        Vec3f center4(7.0, Y1_LIMIT + 14.5, -25);
        Cylinder cylinder1(center3, 2.4, Vec3f(0, 1, 0), 13, big_ben);
        Cylinder cylinder2(center4, 2., Vec3f(0, 1, 0), 3, big_ben);
        objects.push_back(&cylinder1);
        objects.push_back(&cylinder2);

        Vec3f center1(-8.0, Y1_LIMIT + 2.5, -23);
        Vec3f dir1 = Vec3f(-2, 0, 1.8).normalize();
        Vec3f dir2 = Vec3f(0., 1., 0.);
        Cube cube1(center1, dir2, dir1, 4., 1.4, 2.5, red_phone);
        objects.push_back(&cube1);
        Vec3f center2 = center1 + dir1 * 5.4;
        center2.y = Y1_LIMIT + 1.4;
        Cube cube2(center2, dir2, dir1, 1.4, 1.4, 1.4, red_phone);
        objects.push_back(&cube2);

        Plane plane1(Y1_LIMIT, X1_LIMIT, X2_LIMIT, Z1_LIMIT, Z2_LIMIT, wet_asphalt, gray_asphalt);
        objects.push_back(&plane1);

        lights.push_back(Light(Vec3f( 0.5, 6.,  0.5), 1.8));
        lights.push_back(Light(Vec3f( -5., 10.,  -0.5), 1.3));

        Timer.start(); //start timer
        //render the scene
        getImg(objects, lights, img_file, threads);
    }
    else if (scene == 2) {
        //NASA scene
        height = 960;
        width = 960;
        fov = 3. * M_PI / 4.;
        //create panorama object with texture from "./spacestation.jpg"
        Panorama panorama(Vec3f(0., 0., -1.),
                          PAN_RAD,
                          Material(),
                          "./spacestation.jpg");
        //create circle plane object with texture from "./nasa_centre.jpg"
        PlaneCircle circle(Vec3f(0., -1.0, -1.6), 3, plane, "./nasa_center.jpg");
        Sphere sphere1(Vec3f(0., -0.39, -3.2), 0.6, blue);
        Sphere sphere2(Vec3f(-1.4, -0.39, -2.4), 0.6, gray);
        Sphere sphere3(Vec3f(1.4, -0.39, -2.4), 0.6, red);
        objects.push_back(&sphere1);
        objects.push_back(&sphere2);
        objects.push_back(&sphere3);
        objects.push_back(&panorama);
        objects.push_back(&circle);
        lights.push_back(Light(Vec3f( 4, 5.,  -3.5), 1.));
        lights.push_back(Light(Vec3f( 0., 0.,  0.), 0.7));

        Timer.start(); //start timer
        //render the scene
        getImg(objects, lights, img_file, threads);
    }
    else return 0; //finish if scene number is not equal to 1 or 2

    Timer.stop(); //stop the timer
    Timer.print_time(); //print the rendering time
    return 0;
}
