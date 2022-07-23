/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include "directional_and_point_light.h"
#include <iostream>

DirectionalandPointLight::DirectionalandPointLight()
{
	Light();
}

DirectionalandPointLight::DirectionalandPointLight(Vector dir, Colour col, Vertex pos)
{
	Light();
	position = pos;
	direction = dir;
	direction.normalise();
	intensity = col;
	
}

bool DirectionalandPointLight::get_direction(Vertex &surface, Vector &dir)
{
	dir = direction;

	return true;
}

bool DirectionalandPointLight::get_position(Vertex &pos)
{

	pos = position;

	return true;
}

bool DirectionalandPointLight::get_direction_point(Vertex &surface, Vector &dir, float &r2)
{
	
	dir.x = surface.x - position.x;
	dir.y = surface.y - position.y;
	dir.z = surface.z - position.z;
	
	//r2 is the radius squared of the point light to the object
	r2 = dir.len_sqr();
	r2 = r2 / 30.0f;;
	dir.normalise();
	
	return true;
}

void DirectionalandPointLight::get_intensity(Vertex &surface, Colour &level, float &r2)
{

	// Calculate intensity of point light with respect to distance
	level.r = intensity.r / (3.14159265359f * r2);
	level.g = intensity.g / (3.14159265359f * r2);
	level.b = intensity.b / (3.14159265359f * r2);
}
