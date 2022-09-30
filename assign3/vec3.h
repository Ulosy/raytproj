
#include <math.h>

class vec3
{
public:
	double x;
	double y;
	double z;
	vec3() {};
	vec3(double x_in, double y_in, double z_in) : x(x_in), y(y_in), z(z_in)
	{};
	vec3(double* arr) : x(arr[0]), y(arr[1]), z(arr[2])
	{};

	//normalize the vector
	vec3 normalize();

	// scalar multiplication
	vec3 operator*( double t) {
		return vec3(x*t, y*t, z*t);
	}

	friend vec3 operator* (double t, const vec3& v ) {
		return vec3(v.x*t, v.y*t, v.z*t);
	}
	
	vec3 operator/(double t);

	//vector-vector addition
	vec3 operator+(const vec3 &);

	vec3 operator-(const vec3 &);
	//cross product
	vec3 crossProduct(vec3 b);

	//dot product
	double dotProduct(vec3 b);

	double length();
};

