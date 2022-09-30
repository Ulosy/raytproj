#include "vec3.h"
#pragma once

using point = vec3;

class Ray
{
public:
	point orig;
	vec3 dir;
	Ray() {};
	Ray(point e, vec3 dir) : orig(e), dir(dir) {

	};
	//does r(t) = p + dt to return a point in space
	point r(double t) {
		return orig + (dir * t);
	}
};


