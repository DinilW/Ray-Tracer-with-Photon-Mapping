/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include "phong.h"

#include <math.h>
#include <algorithm>
#include <iostream>




// A simple Phong based lighting model

void Phong::compute_base_colour(Colour &result)
{
	result.r = ambient.r;
	result.g = ambient.g;
	result.b = ambient.b;
}

void Phong::compute_light_colour(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
{

	float diff;

	Vector tolight;
	Vector toviewer;

	result.r=0.0f;
	result.g=0.0f;
	result.b=0.0f;

	tolight = ldir;
	tolight.negate();

	toviewer = viewer;
	toviewer.negate();

	diff = normal.dot(tolight);
	
	if (diff < 0.0f) // light is behind surface
	{
		return;
	}

	// diffuse

	result.r += diffuse.r * diff;
	result.g += diffuse.g * diff;
	result.b += diffuse.b * diff;

	// the specular component

	Vector r;
	
	normal.reflection(tolight, r);
	r.normalise();

	float h;

	h = r.dot(toviewer);

	if (h > 0.0f)
	{
		float p = (float)pow(h, power);

		result.r += specular.r * p;
		result.g += specular.g * p;
		result.b += specular.b * p;
	}



}

// russian roulette for non refractive objects
void Phong::compute_russian( float (&power)[3], char &russian)
{
	float Pd;
	float Ps;
	float epsilon;
	float maxdiffrgb, maxspecrgb, maxpower;
	maxdiffrgb = std::max(std::max(diffuse.r * power[0], diffuse.g * power[1]), diffuse.b * power[2]);
	maxspecrgb= std::max(std::max(specular.r * power[0], specular.g * power[1]), specular.b * power[2]);
	maxpower = std::max(std::max(power[0], power[1]),power[2]);
	Pd = maxdiffrgb / maxpower;
	Ps = maxspecrgb / maxpower;
	epsilon = ((float)rand() / (RAND_MAX));
	if (epsilon < Pd)
	{
		russian = 'd';
		power[0] = (power[0] * diffuse.r) / Pd;
		power[1] = (power[1] * diffuse.g) / Pd;
		power[2] = (power[2] * diffuse.b) / Pd;
		//std::cout << power[0] << " " << power[1] << " " << power[2] << std::endl;
	}
	else if (Pd < epsilon && epsilon < (Ps + Pd))
	{
		russian = 's';
		power[0] = (power[0] * specular.r) / Ps;
		power[1] = (power[1] * specular.g) / Ps;
		power[2] = (power[2] * specular.b) / Ps;
		//std::cout << power[0] << " " << power[1] << " " << power[2] << std::endl;
	}
	else
	{
		russian = 'a';
	}
	//std::cout << power[0] << " " << power[1] << " " << power[2] << std::endl;
}

//russian roulette for when object has refractive properties
void Phong::compute_russian_caustic(float(&power)[3], char& russian)
{
	float Pd;
	float Ps;
	float Pt;
	float epsilon;
	float max_transmission, maxdiffrgb, maxspecrgb, maxpower;
	max_transmission = std::max(std::max(refraction * power[0], refraction * power[1]), refraction * power[2]);
	maxdiffrgb = std::max(std::max(diffuse.r * power[0], diffuse.g * power[1]), diffuse.b * power[2]);
	maxspecrgb = std::max(std::max(specular.r * power[0], specular.g * power[1]), specular.b * power[2]);
	maxpower = std::max(std::max(power[0], power[1]), power[2]);
	Pd = maxdiffrgb / maxpower;
	Ps = maxspecrgb / maxpower;
	Pt = max_transmission / maxpower;
	epsilon = ((float)rand() / (RAND_MAX));
	if (epsilon < Pd)
	{
		russian = 'd';
		power[0] = (power[0] * diffuse.r) / Pd;
		power[1] = (power[1] * diffuse.g) / Pd;
		power[2] = (power[2] * diffuse.b) / Pd;
		//std::cout << power[0] << " " << power[1] << " " << power[2] << std::endl;
	}
	else if (Pd < epsilon && epsilon < (Ps + Pd))
	{
		russian = 's';
		power[0] = (power[0] * specular.r) / Ps;
		power[1] = (power[1] * specular.g) / Ps;
		power[2] = (power[2] * specular.b) / Ps;
		//std::cout << power[0] << " " << power[1] << " " << power[2] << std::endl;
	}
	else if ((Ps + Pd) < epsilon && epsilon < (Pt + Ps + Pd))
	{
		russian = 't';
		power[0] = (power[0] * refraction) / Pt;
		power[1] = (power[1] * refraction) / Pt;
		power[2] = (power[2] * refraction) / Pt;
	}
	else
	{
		russian = 'a';
	}
}




