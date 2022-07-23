/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include "scene.h"

Scene::Scene()
{
	object_list = 0;
	light_list = 0;
}

void Scene::trace(Ray ray, Object *objects, Hit &hit)
{
	Hit current_hit;

	hit.flag = false;
	hit.t = 0.0f;
	hit.what = 0;

	while (objects != 0)
	{
		Hit hit_current;

		objects->intersection(ray, hit_current);

		if (hit_current.flag == true)
		{
			if (hit.flag == false)
			{
				hit = hit_current;

			} else if (hit_current.t < hit.t)
			{
				hit = hit_current;
			}
		}

		objects = objects->next;
	}
}

void Scene::raytrace(Ray ray, int level, Object *objects, Light *lights, Colour &colour)
{
	
	// a default colour if we hit nothing.
	float red   = 0.2f;
	float green = 0.2f;
	float blue  = 0.2f;

	// check we've not recursed too far.
	level = level - 1;
	if (level < 0)
	{
		return;
	}

	// first step, find the closest primitive

	Hit best_hit;
	trace(ray, objects, best_hit);

	// if we found a primitive then compute the colour we should see
	if (best_hit.flag)
	{

		best_hit.what->material->compute_base_colour(colour);

		Light *light = lights;

		while (light != (Light *)0)
		{
			Vector viewer;
			Vector ldir;

			viewer.x = -best_hit.position.x;
			viewer.y = -best_hit.position.y;
			viewer.z = -best_hit.position.z;
			viewer.normalise();

			bool lit;

			lit = light->get_direction(best_hit.position, ldir);

			if (lit)
			{
				Colour intensity;
				Colour scaling;

				float r2 = 12.5264f;
				light->get_intensity(best_hit.position, scaling, r2);

				best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);

				intensity.scale(scaling);

				colour.add(intensity);
			}

			light = light->next;
		}

		// TODO: compute reflection ray if material supports it.
		if(1)
		{
		}

		// TODO: compute refraction ray if material supports it.
		if(1)
		{
		}
		
	}
	else
	{
		colour.r = 0.0f;
		colour.g = 0.0f;
		colour.b = 0.0f;
	}
}
