/*
CSCI 420
Assignment 3 Raytracer

Name: Jiny Yang
*/

#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <math.h>

#include "opencv2/core/core.hpp"
#include "opencv2/core/utility.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"

#include "pic.h"
#include "Ray.h"

using point = vec3;


#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10
#define PI 3.14159265
#define SUPERSAMPLE_SIZE 20
#define RECURSION 0

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode= MODE_JPEG;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480
//#define WIDTH 320
//#define HEIGHT 240

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

typedef struct _Bary
{
	double coord[3];
} Bary;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

bool intersectSphere(Ray r, Sphere sphere, double& t);
bool intersectTriangle(Ray r, Triangle triangle, double& t);
bool intersectTriangle(Ray r, Triangle tri, double& t, double dist);
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
vec3 illuminationTriangle(Triangle tri, point p, int t_index, int depth);

const point CAMERA_POS = point(0, 0, 0);

inline double area(point a, point b, point c) { return 0.5 * ((b - a).crossProduct(c - a)).length();};

Bary barycentric(point p, point p1, point p2, point p3) {
	double total_area = area(p1, p2, p3);
	Bary b;
	b.coord[0] = area(p, p2, p3) / total_area;
	b.coord[1] = area(p1, p, p3) / total_area;
	b.coord[2] = area(p1, p2, p) / total_area;
	return b;
}

vec3 illuminationSphere(Sphere s, Ray r, point p, int s_index, int depth) {
	if (depth >= RECURSION) {
		return vec3(0, 0, 0);
	}
	//send shadow ray L from light source
	vec3 color;
	vec3 center = vec3(s.position[0], s.position[1], s.position[2]);
	// N = unit surface normal
	vec3 N = (p - center) / s.radius;
	//negate if ray originates inside of sphere!
	//If the distance from the point to the center of the sphere is less than the radius of the sphere,
	if ((r.orig - vec3(s.position)).length() < s.radius) {
		N.x = -N.x;
		N.y = -N.y;
		N.z = -N.z;
	}
	N.normalize();
	double finalColor[3] = { 0, 0, 0 };
	double localPhongColor[3] = { 0, 0, 0 };
	double reflectedColor[3] = { 0, 0, 0 };

	for (int i = 0; i < num_lights; i++) {
		Light light = lights[i];
		double shadow = 0;
		vec3 L, V;
		//L = unit vector to light
		L = vec3(light.position) - p;
		L.normalize();
		//generate shadow ray L
		//go through loop to see if it intersects with other objects
		Ray shadowRay = Ray(p, L);
		for (int i = 0; i < num_spheres; i++) {
			if (i == s_index)
				continue;
			double t;
			if (intersectSphere(shadowRay, spheres[i], t)) {
				shadow = 1;
				break;
			}		
		}
		for (int i = 0; i < num_triangles && shadow < 1; i++) {
			double t;
			if (intersectTriangle(shadowRay, triangles[i], t)) {
				shadow = 1;
				break;
			}
		}
		if (shadow == 1) {
			continue;
		}
		//V = unit vector to camera
		vec3 cop = CAMERA_POS;
		V = cop - p;
		V.normalize();
		//R = unit reflected vector
		//2(n · l) n − l
		vec3 R = (2 * (N.dotProduct(L)*N)) - L;
		R.normalize();
		//I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ α)
		localPhongColor[0] += light.color[0] * s.color_diffuse[0] * std::max(L.dotProduct(N), 0.) +
			s.color_specular[0] * pow(std::max(R.dotProduct(V), 0.), s.shininess) ;
		localPhongColor[1] += light.color[1] * s.color_diffuse[1] * std::max(L.dotProduct(N), 0.) +
			s.color_specular[1] * pow(std::max(R.dotProduct(V), 0.), s.shininess);
		localPhongColor[2] += light.color[2] * s.color_diffuse[2] * std::max(L.dotProduct(N), 0.) +
			s.color_specular[2] * pow(std::max(R.dotProduct(V), 0.), s.shininess);
		//make reflected ray and get color from the reflected object
		double t = DBL_MAX;
		double t_temp;
		point intPoint;
		vec3 reflected = vec3(0,0,0);
		bool intersect = false;
		Ray reflectedRay = Ray(p, R);
		//do intersect tests on the objects in scene
		if (s.color_specular[0] != 0 || s.color_specular[1] != 0 && s.color_specular[0] != 0) {
			for (int i = 0; i < num_spheres; i++) {
				Sphere s = spheres[i];
				if (intersectSphere(reflectedRay, s, t_temp) && t_temp < t) {
					t = t_temp;
					intPoint = point(reflectedRay.r(t));
					reflected = illuminationSphere(s, reflectedRay, intPoint, i, depth + 1);
					intersect = true;
				}
			}
			for (int t_index = 0; t_index < num_triangles; t_index++) {
				Triangle tri = triangles[t_index];
				if (intersectTriangle(reflectedRay, tri, t_temp) && t_temp < t) {
					t = t_temp;
					intPoint = point(reflectedRay.r(t));
					reflected = illuminationTriangle(tri, intPoint, t_index, depth + 1);
					intersect = true;
				}
			}
			if (intersect) {
				reflectedColor[0] += reflected.x;
				reflectedColor[1] += reflected.y;
				reflectedColor[2] += reflected.z;
			}
		}
		localPhongColor[0] += ambient_light[0];
		localPhongColor[1] += ambient_light[1];
		localPhongColor[2] += ambient_light[2];
	}
	finalColor[0] = (1 - s.color_specular[0]) * localPhongColor[0] + s.color_specular[0] * reflectedColor[0];
	finalColor[1] = (1 - s.color_specular[1]) * localPhongColor[1] + s.color_specular[1] * reflectedColor[1];
	finalColor[2] = (1 - s.color_specular[2]) * localPhongColor[2] + s.color_specular[2] * reflectedColor[2];

	if (finalColor[0] > 1) {
		finalColor[0] = 1;
	}
	if (finalColor[1] > 1) {
		finalColor[1] = 1;
	}
	if (finalColor[2] > 1) {
		finalColor[2] = 1;
	}

	return vec3(finalColor[0], finalColor[1], finalColor[2]);
}


vec3 illuminationSphere(Sphere s, Ray r, point p, int s_index) {
	//send shadow ray L from light source
	vec3 color;
	vec3 center = vec3(s.position[0], s.position[1], s.position[2]);
	// N = unit surface normal
	vec3 N = (p - center) / s.radius;
	//negate if ray originates inside of sphere!
	//If the distance from the point to the center of the sphere is less than the radius of the sphere,
	if ((r.orig - vec3(s.position)).length() < s.radius) {
		N.x = -N.x;
		N.y = -N.y;
		N.z = -N.z;
	}
	N.normalize();
	double finalColor[3] = { 0, 0, 0 };
	double localPhongColor[3] = { 0, 0, 0 };
	double reflectedColor[3] = { 0, 0, 0 };

	for (int i = 0; i < num_lights; i++) {
		Light light = lights[i];
		double shadow = 0;
		vec3 L, V;
		//L = unit vector to light
		L = vec3(light.position) - p;
		L.normalize();
		//generate shadow ray L
		//go through loop to see if it intersects with other objects
		Ray shadowRay = Ray(p, L);
		for (int i = 0; i < num_spheres; i++) {
			if (i == s_index)
				continue;
			double t;
			if (intersectSphere(shadowRay, spheres[i], t)) {
				shadow = 1;
				break;
			}
		}
		for (int i = 0; i < num_triangles && shadow < 1; i++) {
			double t;
			if (intersectTriangle(shadowRay, triangles[i], t)) {
				shadow = 1;
				break;
			}
		}
		if (shadow == 1) {
			continue;
		}
		//V = unit vector to camera
		vec3 cop = CAMERA_POS;
		V = cop - p;
		V.normalize();
		//R = unit reflected vector
		//2(n · l) n − l
		vec3 R = (2 * (N.dotProduct(L)*N)) - L;
		R.normalize();
		//I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ α)
		localPhongColor[0] += light.color[0] * s.color_diffuse[0] * std::max(L.dotProduct(N), 0.) +
			s.color_specular[0] * pow(std::max(R.dotProduct(V), 0.), s.shininess) ;
		localPhongColor[1] += light.color[1] * s.color_diffuse[1] * std::max(L.dotProduct(N), 0.) +
			s.color_specular[1] * pow(std::max(R.dotProduct(V), 0.), s.shininess);
		localPhongColor[2] += light.color[2] * s.color_diffuse[2] * std::max(L.dotProduct(N), 0.) +
			s.color_specular[2] * pow(std::max(R.dotProduct(V), 0.), s.shininess);
		//make reflected ray and get color from the reflected object
		localPhongColor[0] += ambient_light[0];
		localPhongColor[1] += ambient_light[1];
		localPhongColor[2] += ambient_light[2];
	}
	finalColor[0] = localPhongColor[0];
	finalColor[1] = localPhongColor[1];
	finalColor[2] = localPhongColor[2];

	if (finalColor[0] > 1) {
		finalColor[0] = 1;
	}
	if (finalColor[1] > 1) {
		finalColor[1] = 1;
	}
	if (finalColor[2] > 1) {
		finalColor[2] = 1;
	}

	return vec3(finalColor[0], finalColor[1], finalColor[2]);
}

vec3 illuminationTriangle(Triangle tri, point p, int t_index, int depth) {
	if (depth >= RECURSION) {
		return vec3(0, 0, 0);
	}
	vec3 color;
	double finalColor[3] = { 0, 0, 0 };
	double localPhongColor[3] = { 0, 0, 0 };
	double reflectedColor[3] = { 0, 0, 0 };
	double color_diffuse[3];
	double color_specular[3];
	//get baricentric params
	Bary bary = barycentric(p, tri.v[0].position, tri.v[1].position, tri.v[2].position);
	color_diffuse[0] = bary.coord[0] * tri.v[0].color_diffuse[0] + bary.coord[1] * tri.v[1].color_diffuse[0] + bary.coord[2] * tri.v[2].color_diffuse[0];
	color_diffuse[1] = bary.coord[0] * tri.v[0].color_diffuse[1] + bary.coord[1] * tri.v[1].color_diffuse[1] + bary.coord[2] * tri.v[2].color_diffuse[1];
	color_diffuse[2] = bary.coord[0] * tri.v[0].color_diffuse[2] + bary.coord[1] * tri.v[1].color_diffuse[2] + bary.coord[2] * tri.v[2].color_diffuse[2];
	color_specular[0] = bary.coord[0] * tri.v[0].color_specular[0] + bary.coord[1] * tri.v[1].color_specular[0] + bary.coord[2] * tri.v[2].color_specular[0];
	color_specular[1] = bary.coord[0] * tri.v[0].color_specular[1] + bary.coord[1] * tri.v[1].color_specular[1] + bary.coord[2] * tri.v[2].color_specular[1];
	color_specular[2] = bary.coord[0] * tri.v[0].color_specular[2] + bary.coord[1] * tri.v[1].color_specular[2] + bary.coord[2] * tri.v[2].color_specular[2];

	for (int i = 0; i < num_lights; i++) {
		Light light = lights[i];
		double shadow = 0;
		vec3 L, V, N, R;
		vec3 a = tri.v[0].position;
		vec3 b = tri.v[1].position;
		vec3 c = tri.v[2].position;
		//Make a normal to the plane containing triangle
		N = (b - a).crossProduct(c - a);
		N.x = bary.coord[0] * tri.v[0].normal[0] + bary.coord[1] * tri.v[1].normal[0] + bary.coord[2] * tri.v[2].normal[0];
		N.y = bary.coord[0] * tri.v[0].normal[1] + bary.coord[1] * tri.v[1].normal[1] + bary.coord[2] * tri.v[2].normal[1];
		N.z = bary.coord[0] * tri.v[0].normal[2] + bary.coord[1] * tri.v[1].normal[2] + bary.coord[2] * tri.v[2].normal[2];
		N.normalize();
		//L = unit vector to light
		L = vec3(light.position) - p;
		double L_length = L.length();
		L.normalize();
		//go through loop to see if it intersects with other objects
		//if does not intersect, color = (0,0,0)
		//else add color
		//generate shadow ray L
		Ray shadowRay = Ray(p, L);
		for (int i = 0; i < num_spheres; i++) {
			double t;
			if (intersectSphere(shadowRay, spheres[i], t)) {
				shadow = 1;
				break;
			}
		}
		for (int i = 0; i < num_triangles && shadow < 1; i++) {
			double t;
			if (i == t_index){
				continue; 
			} 
			if (intersectTriangle(shadowRay, triangles[i], t, L_length)) {
				shadow = 1;
				break;
			}
		}
		if (shadow == 1) {
			continue;
		}
		//V = unit vector to camera
		vec3 cop = CAMERA_POS;
		V = cop - p;
		V.normalize();
		//R = unit reflected vector
		//2(n · l) n − l
		R = (2 * (N.dotProduct(L)*N)) - L;
		R.normalize();
		//I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ α)

		double shininess = bary.coord[0] * tri.v[0].shininess + bary.coord[1] * tri.v[1].shininess + bary.coord[2] * tri.v[2].shininess;
		localPhongColor[0] += light.color[0] * color_diffuse[0] * std::max(L.dotProduct(N), 0.) +
			color_specular[0] * pow(std::max(R.dotProduct(V), 0.), shininess);
		localPhongColor[1] += light.color[1] * color_diffuse[1] * std::max(L.dotProduct(N), 0.) +
			color_specular[1] * pow(std::max(R.dotProduct(V), 0.), shininess);
		localPhongColor[2] += light.color[2] * color_diffuse[2] * std::max(L.dotProduct(N), 0.) +
			color_specular[2] * pow(std::max(R.dotProduct(V), 0.), shininess);
		//make reflected ray and get color from the reflected object
		double t = DBL_MAX;
		double t_temp;
		point intPoint;
		vec3 reflected = vec3(0,0,0);
		bool intersect = false;
		Ray reflectedRay = Ray(p, R);
		//do intersect tests on the objects in scene
		if (color_specular[0] != 0 || color_specular[1] != 0 && color_specular[2] != 0) {
			for (int i = 0; i < num_spheres; i++) {
				Sphere s = spheres[i];
				if (intersectSphere(reflectedRay, s, t_temp) && t_temp < t) {
					t = t_temp;
					intPoint = point(reflectedRay.r(t));
					reflected = illuminationSphere(s, reflectedRay, intPoint, i, depth + 1);
					intersect = true;
				}
			}
			for (int i = 0; i < num_triangles; i++) {
				Triangle tri = triangles[i];
				if (intersectTriangle(reflectedRay, tri, t_temp) && t_temp < t) {
					t = t_temp;
					intPoint = point(reflectedRay.r(t));
					reflected = illuminationTriangle(triangles[i], intPoint, i, depth + 1);
					intersect = true;
				}
			}
			if (intersect) {
				reflectedColor[0] += reflected.x;
				reflectedColor[1] += reflected.y;
				reflectedColor[2] += reflected.z;
			}
		}
	}
	localPhongColor[0] += ambient_light[0];
	localPhongColor[1] += ambient_light[1];
	localPhongColor[2] += ambient_light[2];

	finalColor[0] = (1 - color_specular[0]) * localPhongColor[0] + color_specular[0] * reflectedColor[0];
	finalColor[1] = (1 - color_specular[1]) * localPhongColor[1] + color_specular[1] * reflectedColor[1];
	finalColor[2] = (1 - color_specular[2]) * localPhongColor[2] + color_specular[2] * reflectedColor[2];

	if (finalColor[0] > 1) {
		finalColor[0] = 1;
	}
	if (finalColor[1] > 1) {
		finalColor[1] = 1;

	}
	if (finalColor[2] > 1) {
		finalColor[2] = 1;
	}
	return vec3(finalColor[0], finalColor[1], finalColor[2]);
}
vec3 illuminationTriangle(Triangle tri, point p, int t_index) {

	vec3 color;
	double finalColor[3] = { 0, 0, 0 };
	double localPhongColor[3] = { 0, 0, 0 };
	double color_diffuse[3];
	double color_specular[3];
	//get baricentric params
	Bary bary = barycentric(p, tri.v[0].position, tri.v[1].position, tri.v[2].position);
	color_diffuse[0] = bary.coord[0] * tri.v[0].color_diffuse[0] + bary.coord[1] * tri.v[1].color_diffuse[0] + bary.coord[2] * tri.v[2].color_diffuse[0];
	color_diffuse[1] = bary.coord[0] * tri.v[0].color_diffuse[1] + bary.coord[1] * tri.v[1].color_diffuse[1] + bary.coord[2] * tri.v[2].color_diffuse[1];
	color_diffuse[2] = bary.coord[0] * tri.v[0].color_diffuse[2] + bary.coord[1] * tri.v[1].color_diffuse[2] + bary.coord[2] * tri.v[2].color_diffuse[2];
	color_specular[0] = bary.coord[0] * tri.v[0].color_specular[0] + bary.coord[1] * tri.v[1].color_specular[0] + bary.coord[2] * tri.v[2].color_specular[0];
	color_specular[1] = bary.coord[0] * tri.v[0].color_specular[1] + bary.coord[1] * tri.v[1].color_specular[1] + bary.coord[2] * tri.v[2].color_specular[1];
	color_specular[2] = bary.coord[0] * tri.v[0].color_specular[2] + bary.coord[1] * tri.v[1].color_specular[2] + bary.coord[2] * tri.v[2].color_specular[2];

	for (int i = 0; i < num_lights; i++) {
		Light light = lights[i];
		double shadow = 0;
		vec3 L, V, N, R;
		vec3 a = tri.v[0].position;
		vec3 b = tri.v[1].position;
		vec3 c = tri.v[2].position;
		//Make a normal to the plane containing triangle
		N = (b - a).crossProduct(c - a);
		N.x = bary.coord[0] * tri.v[0].normal[0] + bary.coord[1] * tri.v[1].normal[0] + bary.coord[2] * tri.v[2].normal[0];
		N.y = bary.coord[0] * tri.v[0].normal[1] + bary.coord[1] * tri.v[1].normal[1] + bary.coord[2] * tri.v[2].normal[1];
		N.z = bary.coord[0] * tri.v[0].normal[2] + bary.coord[1] * tri.v[1].normal[2] + bary.coord[2] * tri.v[2].normal[2];
		N.normalize();
		//L = unit vector to light
		L = vec3(light.position) - p;
		double L_length = L.length();
		L.normalize();
		//go through loop to see if it intersects with other objects
		//if does not intersect, color = (0,0,0)
		//else add color
		//generate shadow ray L
		Ray shadowRay = Ray(p, L);
		for (int i = 0; i < num_spheres; i++) {
			double t;
			if (intersectSphere(shadowRay, spheres[i], t)) {
				shadow = 1;
				break;
			}
		}
		for (int i = 0; i < num_triangles && shadow < 1; i++) {
			double t;
			if (i == t_index) {
				continue;
			}
			if (intersectTriangle(shadowRay, triangles[i], t, L_length)) {
				shadow = 1;
				break;
			}
		}
		if (shadow == 1) {
			continue;
		}
		//V = unit vector to camera
		vec3 cop = CAMERA_POS;
		V = cop - p;
		V.normalize();
		//R = unit reflected vector
		//2(n · l) n − l
		R = (2 * (N.dotProduct(L)*N)) - L;
		R.normalize();
		//I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ α)
		double shininess = bary.coord[0] * tri.v[0].shininess + bary.coord[1] * tri.v[1].shininess + bary.coord[2] * tri.v[2].shininess;
		localPhongColor[0] += light.color[0] * color_diffuse[0] * std::max(L.dotProduct(N), 0.) +
			color_specular[0] * pow(std::max(R.dotProduct(V), 0.), shininess);
		localPhongColor[1] += light.color[1] * color_diffuse[1] * std::max(L.dotProduct(N), 0.) +
			color_specular[1] * pow(std::max(R.dotProduct(V), 0.), shininess);
		localPhongColor[2] += light.color[2] * color_diffuse[2] * std::max(L.dotProduct(N), 0.) +
			color_specular[2] * pow(std::max(R.dotProduct(V), 0.), shininess);
	}
	localPhongColor[0] += ambient_light[0];
	localPhongColor[1] += ambient_light[1];
	localPhongColor[2] += ambient_light[2];

	finalColor[0] = localPhongColor[0];
	finalColor[1] = localPhongColor[1];
	finalColor[2] = localPhongColor[2];

	if (finalColor[0] > 1) {
		finalColor[0] = 1;
	}
	if (finalColor[1] > 1) {
		finalColor[1] = 1;

	}
	if (finalColor[2] > 1) {
		finalColor[2] = 1;
	}
	return vec3(finalColor[0], finalColor[1], finalColor[2]);
}

//return the color of surface the ray hit
vec3 trace(point camera_pos, Ray ray, bool& intersect) {
	double t = DBL_MAX;
	double t_temp;
	point intPoint;
	intersect = false;
	vec3 finalColor;
	//do intersect tests on the objects in scene
	for (int i = 0; i < num_spheres; i++) {
		Sphere s = spheres[i];
		if (intersectSphere(ray, s, t_temp) && t_temp < t){
			t = t_temp;
			//compute position of ray intersection with the t
			intPoint = point(ray.r(t));
			if(!RECURSION)
				finalColor = illuminationSphere(s, ray, intPoint, i);
			else
				finalColor = illuminationSphere(s, ray, intPoint, i, 0);
			intersect = true;
		}
	}
	for (int t_index = 0; t_index < num_triangles; t_index++) {
		Triangle tri = triangles[t_index];
		if (intersectTriangle(ray, tri, t_temp) && t_temp < t) {
			t = t_temp;
			//compute position of ray intersection with the t
			intPoint = point(ray.r(t));
			double color_specular[3];
			if(!RECURSION)
				finalColor = illuminationTriangle(tri, intPoint, t_index);
			else
				finalColor = illuminationTriangle(tri, intPoint, t_index, 0);
			intersect = true;
		}	
	}
	if (t == DBL_MAX) {
		//background color
		return vec3(1, 1, 1);
	}
	if (finalColor.x > 1)
		finalColor.x = 1;
	if (finalColor.y > 1)
		finalColor.y = 1;
	if (finalColor.z > 1)
		finalColor.z = 1;
	return finalColor;
}

//ray-sphere intersection
bool intersectSphere(Ray r, Sphere sphere, double& t)
{
	double t0, t1;
	point center = point(sphere.position);
	vec3 oc = center - r.orig;
	double radius2 = pow(sphere.radius,2);
	double t_oa = oc.dotProduct(r.dir);
	double d2 = oc.dotProduct(oc) - pow(t_oa, 2);
	if (d2 > radius2) {
		return false;
	}
	double t_pa= sqrt(radius2 - d2);
	//Check if t0, t1 > 0
	t0 = t_oa - t_pa;
	t1 = t_oa + t_pa;
	t = std::min(t0,t1);
	if (t0 < 0 && t1 < 0) {
		return false;
	}
	else if (t0 < 0) {
		t = t1;
	}
	else if (t1 < 0) {
		t = t0;
	}
	return true;
}


//ray-triangle intersection
bool intersectTriangle(Ray r, Triangle tri, double& t)
{
	vec3 a = tri.v[0].position;
	vec3 b = tri.v[1].position;
	vec3 c = tri.v[2].position;
	//Make a normal to the plane containing triangle
	vec3 N = (b - a).crossProduct(c - a);
	//Derive t from plane and ray eq
	t = (a - r.orig).dotProduct(N) / r.dir.dotProduct(N);
	//do a intersection test with plane
	//The Ray And The Triangle Are Parallel
	if (fabs(N.dotProduct(r.dir)) < 0.0000001) {
		return false;
	}
	//The Triangle Is "Behind" The Ray
	if (t < 0) {
		return false;
	}
	//get point where ray intersect the plane
	point p = r.r(t);

	//(b−a)×(x−a)⋅n > 0
	//(c−b)×(x−b)⋅n > 0
	//(a−c)×(x−c)⋅n > 0
	if ((b - a).crossProduct(p - a).dotProduct(N) > 0 &&
		(c - b).crossProduct(p - b).dotProduct(N) > 0 &&
		(a - c).crossProduct(p - c).dotProduct(N) > 0)
		return true;
	return false;
}


//ray-triangle intersection
bool intersectTriangle(Ray r, Triangle tri, double& t, double dist)
{
	vec3 a = tri.v[0].position;
	vec3 b = tri.v[1].position;
	vec3 c = tri.v[2].position;
	//Make a normal to the plane containing triangle
	vec3 N = (b - a).crossProduct(c - a);
	//Derive t from plane and ray eq
	t = (a - r.orig).dotProduct(N) / r.dir.dotProduct(N);
	//do a intersection test with plane
	//The Ray And The Triangle Are Parallel
	if (fabs(N.dotProduct(r.dir)) < 0.0000001) {
		return false;
	}
	//The Triangle Is "Behind" The Ray
	if (t < 0) {
		return false;
	}
	if (t > dist) {
		return false;
	}
	//get point where ray intersect the plane
	point p = r.r(t);
	//test if inside triangle 
	//(b−a)×(x−a)⋅n > 0
	//(c−b)×(x−b)⋅n > 0
	//(a−c)×(x−c)⋅n > 0
	if ((b - a).crossProduct(p - a).dotProduct(N) > 0 &&
		(c - b).crossProduct(p - b).dotProduct(N) > 0 &&
		(a - c).crossProduct(p - c).dotProduct(N) > 0)
		return true;
	return false;
}


//generates ray within pixel in a random position
Ray generateRay(double i, double j) {
	// pixel position in x, y planes
	double x_pos = (i + rand() / (1. + RAND_MAX )) / (double)(WIDTH - 1);
	double y_pos = (j + rand() / (1. + RAND_MAX )) / (double)(HEIGHT - 1);

	double aspect = (double)WIDTH / HEIGHT;
	double scale = tan(PI * 0.5 * fov / 180.);
	
	vec3 l = vec3(2 * aspect *scale, 0, 0);//length vector of the in the viewport
	vec3 h = vec3(0, 2.*scale, 0); 	//height vector of the viewport
	
	point center = vec3(0, 0, -1); 	//the target position, center of image plane
	vec3 bottom_left = center - h / 2 - l / 2;
	vec3 raydir = bottom_left + l * x_pos + h * y_pos;
	raydir.normalize();//normalize raydir
	point cop = CAMERA_POS; //the eye position
	Ray r = Ray(cop,raydir);
	return r;
}

void draw_scene()
{
	Sphere s = spheres[0];

  unsigned int x,y;
  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++){
		Ray rays[SUPERSAMPLE_SIZE];
		// for each pixel in a screen make a number of rays
		vec3 aliased_pix = vec3(0,0,0);
		for (int r_idx = 0; r_idx < SUPERSAMPLE_SIZE; r_idx++)
		{
			rays[r_idx] = generateRay(x, y) ;
			// Trace the projection of the ray created.
			// use illumination, intersection code inside of trace func.
			bool intersect;
			vec3 ray_pixel = trace(CAMERA_POS, rays[r_idx], intersect);
			aliased_pix.x += ray_pixel.x;
			aliased_pix.y += ray_pixel.y;
			aliased_pix.z += ray_pixel.z;
		}
		aliased_pix = aliased_pix/ SUPERSAMPLE_SIZE;
		plot_pixel(x, y, (int)(aliased_pix.x * 255.f), (int)(aliased_pix.y * 255.f), (int)(aliased_pix.z *  255.f));
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

/* Write a jpg image from buffer*/
void save_jpg()
{
	if (filename == NULL)
		return;

	// Allocate a picture buffer // 
	cv::Mat3b bufferBGR = cv::Mat::zeros(HEIGHT, WIDTH, CV_8UC3); //rows, cols, 3-channel 8-bit.
	printf("File to save to: %s\n", filename);

	// unsigned char buffer[HEIGHT][WIDTH][3];
	for (int r = 0; r < HEIGHT; r++) {
		for (int c = 0; c < WIDTH; c++) {
			for (int chan = 0; chan < 3; chan++) {
				unsigned char red = buffer[r][c][0];
				unsigned char green = buffer[r][c][1];
				unsigned char blue = buffer[r][c][2];
				bufferBGR.at<cv::Vec3b>(r,c) = cv::Vec3b(blue, green, red);
			}
		}
	}
	if (cv::imwrite(filename, bufferBGR)) {
		printf("File saved Successfully\n");
	}
	else {
		printf("Error in Saving\n");
	}
}

void parse_check(char *expected,char *found)
{
  if(stricmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);

}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(stricmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(stricmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(stricmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

static int once = 0;

void idle()
{
  //hack to make it only draw once
  if(!once)
  {
	draw_scene();
	if(mode == MODE_JPEG)
		save_jpg();
	once = 1;
  }
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;
  
  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

