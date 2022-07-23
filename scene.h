/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include "colour.h"
#include "ray.h"
#include "object.h"
#include "light.h"
#include "hit.h"

// Scene is a class that is used to build a scene database of objects
// and lights and then trace a ray through it.

class Scene {
public:

	Object *object_list;
	Light *light_list;

	Scene();

	// Trace a ray through the scene and find the closest if any object
	// intersection in front of the ray.
	void trace(Ray ray, Object *objects, Hit &hit);
	
	// Trace a ray through the scene and return its colour. This function
	// is the one that should recurse down the reflection/refraction tree.
	void raytrace(Ray ray, int level, Object *objects, Light *lights, Colour &colour);

};
