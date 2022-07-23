/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// DirectionaLight is a child class of Light and implements a light
// with constant value in a given direction. The light has no position
// and can be treated as infinitely far away.

#pragma once
#include "light.h"

class DirectionalandPointLight : public Light {
public:
	
	Vector direction;
	Vertex position;
	Colour intensity;

	DirectionalandPointLight();
	DirectionalandPointLight(Vector dir, Colour col, Vertex pos);

	bool get_direction(Vertex &surface, Vector &dir);
	bool get_position(Vertex& pos);
	bool get_direction_point(Vertex &surface, Vector &dir, float &r2);
	void get_intensity(Vertex& surface, Colour& intensity, float &r2);

};
