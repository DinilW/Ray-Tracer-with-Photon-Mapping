/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Material is the base class for materials.

#pragma once

#include "vector.h"
#include "colour.h"

class Material {
public:

	float reflection;
	float refraction;
	float refracI;



	virtual void compute_base_colour(Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}
	virtual void compute_light_colour(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}

	virtual void compute_russian(float(&power)[3], char& russian)
	{
		russian = 'a';
	}

	virtual void compute_russian_caustic(float(&power)[3], char& russian)
	{
		russian = 'a';
	}

	virtual void get_ambient(float(&amb)[3])
	{

		amb[0] = 0.0f;
		amb[1] = 0.0f;
		amb[2] = 0.0f;
	}


};
