#include "vec3.h"
#include "math.h"

void normalize() {
	
}

vec3 vec3::normalize()
{
	double length = sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
	x /= length;
	y /= length;
	z /= length;
	return *this;
}

double vec3::length()
{
	double length = sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
	return length;
}




vec3 vec3::operator/(double t)
{
	return vec3(x/t, y/t, z/t);
}

vec3 vec3::operator+(const vec3 & v)
{
	return vec3(x+v.x, y+v.y, z+ v.z);
}

vec3 vec3::operator-(const vec3 & v)
{
	return vec3(x - v.x, y - v.y, z - v.z);
}

vec3 vec3::crossProduct(vec3 b)
{
	vec3 a = *this;
	vec3 c;
	c.x = a.y * b.z - a.z * b.y;
	c.y = -(a.x * b.z - a.z * b.x);
	c.z = a.x * b.y - a.y * b.x;
	return c;
}

double vec3::dotProduct(vec3 b){
	vec3 a = *this;
	double c;
	c = a.x * b.x + a.y * b.y +a.z * b.z;
	return c;
}