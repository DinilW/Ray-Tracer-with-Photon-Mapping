 /***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab two.
 * You should not need to modify this code.
 * It creates a framebuffer, loads an triangle mesh object, calls the drawing function to render the object and then outputs the framebuffer as a ppm file.
 *
 * On linux.bath.ac.uk:
 *
 * Compile the code using g++ -o lab4 main_lab4.cpp framebuffer.cpp polymesh.cpp sphere.cpp phong.cpp directional_and_point_light.cpp -lm
 *
 * Execute the code using lab4
 *
 * This will produce an image file called test.ppm. You can convert this a png file for viewing using
 *
 *
 *
 * 
 */

#include "framebuffer.h"
#include "ray.h"
#include "hit.h"
#include "polymesh.h"
#include "sphere.h"
#include "light.h"
#include "directional_and_point_light.h"
#include "material.h"
#include "phong.h"
#include "src/tree.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm> 
#include <math.h>

using namespace std;
//Point is our photon
struct Point
{
    float pos[3];
    float power[3];
    char type;
    Vector normal;

    Point() {}
    Point(float x, float y, float z) {
        pos[0] = x;
        pos[1] = y;
        pos[2] = z;
    }

    Point(float x, float y, float z, float a, float b, float c, Vector norm) {
        pos[0] = x;
        pos[1] = y;
        pos[2] = z;
        power[0] = a;
        power[1] = b;
        power[2] = c;
        normal = norm;
    }

    Point(float x, float y, float z, float a, float b, float c, Vector norm, char typePhoton) {
        pos[0] = x;
        pos[1] = y;
        pos[2] = z;
        power[0] = a;
        power[1] = b;
        power[2] = c;
        normal = norm;
        type = typePhoton;
    }



    float operator [] (int i) const {
        return pos[i];
    }

    bool operator == (const Point& p) const {
        return pos[0] == p[0] && pos[1] == p[1] && pos[2] == p[2];
    }

    friend std::ostream& operator << (std::ostream& s, const Point& p) {
        return s << '(' << p[0] << ", " << p[1] << ", " << p[2] << ')';
    }
};


typedef KD::Core<3, Point> CORE;

void object_test(Ray ray, Object *objects, Hit &best_hit)
{
  Object *obj = objects;

  best_hit.flag = false;


  while(obj != 0)
  {
    Hit obj_hit;
    obj_hit.flag=false;
	  
    obj->intersection(ray, obj_hit);

    
    if (obj_hit.flag)
    {
      if (obj_hit.t > 0.0f)
      {
        if (best_hit.flag == false)
	{
	  best_hit = obj_hit;
	} else if (obj_hit.t < best_hit.t)
	{
	  best_hit = obj_hit;
	}
      }
    }
    
    obj = obj->next;
  }


  return;
}

//reflection calculation
void reflection(Hit &best_hit, Ray ray, Ray &rray)
{
    Colour rcolour;
    best_hit.normal.reflection(ray.direction, rray.direction);
    rray.position.x = best_hit.position.x + 0.0001 * rray.direction.x;
    rray.position.y = best_hit.position.y + 0.0001 * rray.direction.y;
    rray.position.z = best_hit.position.z + 0.0001 * rray.direction.z;
}


//refraction calculation with fresnels
void refraction(Hit &best_hit, Ray ray, Ray &refracRay, float& refCoef, float &refracCoef)
{



    float abscosi, cosi, cost, rpar, rper, refractionIndex;
    Ray I;
    Colour refracColour;
    refractionIndex = best_hit.what->material->refracI;
    I.direction.x = (ray.direction.x);
    I.direction.y = (ray.direction.y);
    I.direction.z = (ray.direction.z);
    
 
    cosi = best_hit.normal.dot(Vector(I.direction.x, I.direction.y, I.direction.z));
    
    if (cosi < 0.0f)
    {
        cosi = -cosi;
    }
    else
    {

        best_hit.normal.x = -(best_hit.normal.x);
        best_hit.normal.y = -(best_hit.normal.y);
        best_hit.normal.z = -(best_hit.normal.z);
        refractionIndex = (1.0f / best_hit.what->material->refracI);

    }

    if ((1.0f - (1.0f / pow(refractionIndex, 2)) * (1 - pow(cosi, 2))) < 0)
    {
        refracRay.direction.x = 0;
        refracRay.direction.y = 0;
        refracRay.direction.z = 0;
        refracRay.position.x = 0;
        refracRay.position.y = 0;
        refracRay.position.z = 0;
    }
    else
    {
        cost = sqrt(1.0f - (1.0f / pow(refractionIndex, 2)) * (1 - pow(cosi, 2)));
        abscosi = abs(cosi);
        if (abscosi > 1.0f)
        {
            abscosi = 1.0f;
        }
        //fresnels terms
        rpar = ((refractionIndex * abscosi) - cost) / ((refractionIndex * abscosi) + cost);
        rper = (abscosi - (refractionIndex * cost)) / (abscosi + (refractionIndex * cost));
        refCoef = (((rpar * rpar) + (rper * rper)) / 2.0f);
        refracCoef = 1.0f - refCoef;
        refracRay.direction.x = ((1 / refractionIndex) * I.direction.x) - ((cost - (1 / refractionIndex) * cosi) * best_hit.normal.x);
        refracRay.direction.y = ((1 / refractionIndex) * I.direction.y) - ((cost - (1 / refractionIndex) * cosi) * best_hit.normal.y);
        refracRay.direction.z = ((1 / refractionIndex) * I.direction.z) - ((cost - (1 / refractionIndex) * cosi) * best_hit.normal.z);
        refracRay.position.x = best_hit.position.x + (0.0001 * refracRay.direction.x);
        refracRay.position.y = best_hit.position.y + (0.0001 * refracRay.direction.y);
        refracRay.position.z = best_hit.position.z + (0.0001 * refracRay.direction.z);
    }
    return;
}

//raytracing function with detection of photons using sampling sphere
void findphotons(Ray ray, Object* objects, Light* lights, Colour& colour, int depth, KD::Tree<CORE>& kdtree, KD::Tree<CORE>& kdtree_caustic)
{
    // first step, find the closest primitive

    Hit shadow_hit;
    Hit best_hit;
    object_test(ray, objects, best_hit);

    // if we found a primitive then compute the colour we should see
    if (best_hit.flag)
    {
        best_hit.what->material->compute_base_colour(colour);
        //Photon calculations
        
        const float PI = 3.14159265359f;
        int numPhotons, numPhotons_caustic;
        int numShadowPhotons = 0;
        int numDirectPhotons = 0;
        float radiance[3] = { 0,0,0 };
        float radiance_caustic[3] = { 0,0,0 };
        float radius;
        float radius_caustic = 0.1;
        vector<Point> photons, photons_caustic;
        photons = kdtree.nearest(Point(best_hit.position.x, best_hit.position.y, best_hit.position.z), 500);
        photons_caustic = kdtree_caustic.nearest(Point(best_hit.position.x, best_hit.position.y, best_hit.position.z),500);
      
        numPhotons = photons.size();
        numPhotons_caustic = photons_caustic.size();
        
        for (int i = 0; i < numPhotons; i++)
        {
            if (numPhotons - i == 1)
            {
                Vector length;
                length.x = best_hit.position.x - photons[i].pos[0];
                length.y = best_hit.position.y - photons[i].pos[1];
                length.z = best_hit.position.z - photons[i].pos[2];
                radius = sqrt(length.len_sqr());
            }

            if (photons[i].type == 'D')
            {
                numDirectPhotons += 1;
            }
            else 
            {
                radiance[0] += photons[i].power[0];
                radiance[1] += photons[i].power[1];
                radiance[2] += photons[i].power[2];
            }

            if (photons[i].type == 'S')
            {
                numShadowPhotons += 1;
            }

        }
        
        for (int i = 0; i < numPhotons_caustic; i++)
        {
            if (numPhotons_caustic - i == 1)
            {
                Vector length;
                length.x = best_hit.position.x - photons_caustic[i].pos[0];
                length.y = best_hit.position.y - photons_caustic[i].pos[1];
                length.z = best_hit.position.z - photons_caustic[i].pos[2];
                radius_caustic = sqrt(length.len_sqr());
            }

            if (photons_caustic[i].type == 'D')
            {
                numDirectPhotons += 1;
            }
            else 
            {
                radiance_caustic[0] += photons_caustic[i].power[0];
                radiance_caustic[1] += photons_caustic[i].power[1];
                radiance_caustic[2] += photons_caustic[i].power[2];
            }

            if (photons_caustic[i].type == 'S')
            {
                numShadowPhotons += 1;
            }

        }

      
        radiance[0] = radiance[0] /(PI * (radius * radius));
        radiance[1] = radiance[1] /(PI * (radius * radius));
        radiance[2] = radiance[2] /(PI * (radius * radius));
        colour.r += radiance[0] + (radiance_caustic[0] / (PI * (radius_caustic * radius_caustic)));
        colour.g += radiance[1] + (radiance_caustic[1] / (PI * (radius_caustic * radius_caustic)));
        colour.b += radiance[2] + (radiance_caustic[2] / (PI * (radius_caustic * radius_caustic)));
       
       
        //standard raytracing from here
        
        Light* light = lights;

        while (light != (Light*)0)
        {
            Vector viewer;
            Vector ldir;

            viewer.x = -best_hit.position.x;
            viewer.y = -best_hit.position.y;
            viewer.z = -best_hit.position.z;
            viewer.normalise();

            bool lit;
            float r2;
            lit = light->get_direction_point(best_hit.position, ldir, r2);

            if (ldir.dot(best_hit.normal) > 0)
            {
                lit = false;//light is facing wrong way.
            }

            if (numShadowPhotons != 0 || numDirectPhotons == 0)
            {

                Ray shadow_ray;

                shadow_ray.direction.x = -ldir.x;
                shadow_ray.direction.y = -ldir.y;
                shadow_ray.direction.z = -ldir.z;
                shadow_ray.position.x = best_hit.position.x + (0.0001f * shadow_ray.direction.x);
                shadow_ray.position.y = best_hit.position.y + (0.0001f * shadow_ray.direction.y);
                shadow_ray.position.z = best_hit.position.z + (0.0001f * shadow_ray.direction.z);



                object_test(shadow_ray, objects, shadow_hit);

                if (shadow_hit.flag == true)
                {
                    float r = sqrt(r2);
                    if (shadow_hit.t < 1000000000.0f && r > shadow_hit.t)
                    {
                        lit = false; //there's a shadow so no lighting, if realistically close
                    }
                }
            }

            if (lit)
            {
                Colour intensity;
                Colour scaling;

                light->get_intensity(best_hit.position, scaling, r2);

                best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);

                intensity.scale(scaling);

                colour.add(intensity);
            }

            light = light->next;
        }
             
        if (best_hit.what->material->reflection > 0.0f)
        {
            Colour rcolour;

            if (depth == 0)
            {
                rcolour.r = 0.0f;
                rcolour.g = 0.0f;
                rcolour.b = 0.0f;
            }
            else
            {
                Ray rray;
                reflection(best_hit, ray, rray);
                findphotons(rray, objects, lights, rcolour, depth - 1, kdtree, kdtree_caustic);
                rcolour.r = best_hit.what->material->reflection * rcolour.r;
                rcolour.g = best_hit.what->material->reflection * rcolour.g;
                rcolour.b = best_hit.what->material->reflection * rcolour.b;
            }
            colour.add(rcolour);
        }

        if (best_hit.what->material->refraction > 0.0f)
        {
            Colour refracColour;
            if (depth == 0)
            {
                refracColour.r = 0.0f;
                refracColour.g = 0.0f;
                refracColour.b = 0.0f;
            }
            else
            {
                Ray refracRay;
                refraction(best_hit, ray, refracRay, best_hit.what->material->reflection, best_hit.what->material->refraction);
                findphotons(refracRay, objects, lights, refracColour, depth - 1, kdtree, kdtree_caustic);
                refracColour.r = best_hit.what->material->refraction * refracColour.r;
                refracColour.g = best_hit.what->material->refraction * refracColour.g;
                refracColour.b = best_hit.what->material->refraction * refracColour.b;
            }
            colour.add(refracColour);
            
        }
       
        
        
        
    }
    
    else
    {
        colour.r = 0.0f;
        colour.g = 0.0f;
        colour.b = 0.0f;
    }
}


//finds other intersections along photon trajectory and inserts shadow photon into kdtree
void shadowtrace(Ray ray, Object* objects, Light* lights, KD::Tree<CORE>& kdtree, Hit &best_hit)
{
    Hit shadow_hit;
    Ray shadow_ray;
    shadow_ray.direction = ray.direction;
    shadow_ray.position.x = best_hit.position.x + 0.0001f * ray.direction.x;
    shadow_ray.position.y = best_hit.position.y + 0.0001f * ray.direction.y;
    shadow_ray.position.z = best_hit.position.z + 0.0001f * ray.direction.z;
    object_test(shadow_ray, objects, shadow_hit);
    if (shadow_hit.flag)
    {
        if (shadow_hit.t < 1000000000.0f) 
        {
            kdtree.insert(Point(shadow_hit.position.x, shadow_hit.position.y, shadow_hit.position.z, 0.0f, 0.0f, 0.0f, shadow_hit.normal, 'S'));
            shadowtrace(ray, objects, lights, kdtree, shadow_hit);
        }
    }
}

//traces photons and inserts themn into the kd tree using russian roulette
void photontrace(Ray ray, Object* objects, Light* lights, float(&power)[3], KD::Tree<CORE> &kdtree, int &depth)
{
    Hit shadow_hit;
    Hit best_hit;
    Ray shadow_ray;
    char russian;
    object_test(ray, objects, best_hit);

    if (best_hit.flag)
    {
        
        if (depth==1)
        { 
            //First hit stored as Direct photon
            depth += 1;
            kdtree.insert(Point(best_hit.position.x, best_hit.position.y, best_hit.position.z, power[0], power[1], power[2], best_hit.normal, 'D'));
            shadowtrace(ray, objects, lights, kdtree, best_hit);
            
        }
        else 
        {
            //All subsequent hits stored as Indirect photons
            kdtree.insert(Point(best_hit.position.x, best_hit.position.y, best_hit.position.z, power[0], power[1], power[2], best_hit.normal, 'I'));
        }


        // russian roulette 
        if (best_hit.what->material->refraction > 0.0f)
        {
            best_hit.what->material->compute_russian_caustic(power, russian);
        }
        else
        {
            best_hit.what->material->compute_russian(power, russian);
        }
        
        
        // corresponding calculations if the roulette decides diffuse reflection, specular, or transmission
        if (russian == 'd')
        {
            //diffuse bounce is in random direction in hemisphere where normal is facing
            float x, y, z;
            Ray diffuseray;
            x = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
            y = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
            z = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
            while (x * x + y * y + z * z > 1.0f)
            {
                x = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
                y = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
                z = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
            }
            diffuseray.direction.x = x;
            diffuseray.direction.y = y;
            diffuseray.direction.z = z;
            if (diffuseray.direction.dot(best_hit.normal) < 0.0f)
            {
                diffuseray.direction.negate();
            }
            diffuseray.position.x = best_hit.position.x + (0.0001 * diffuseray.direction.x);
            diffuseray.position.y = best_hit.position.y + (0.0001 * diffuseray.direction.y);
            diffuseray.position.z = best_hit.position.z + (0.0001 * diffuseray.direction.z);
            photontrace(diffuseray, objects, lights, power, kdtree, depth);
   
        }
        else if (russian == 's')
        {
            //specular is exact reflection
            Ray specray;
            reflection(best_hit, ray, specray);
            photontrace(specray, objects, lights, power, kdtree, depth);

        }
        else if (russian == 't')
        { 
            
            //refraction calling the function
            Ray refracRay;
            refraction(best_hit, ray, refracRay, best_hit.what->material->reflection, best_hit.what->material->refraction);
            photontrace(refracRay, objects, lights, power, kdtree, depth);
        }
    }
}



//emitting photons into the scene for global map
void emit_photons(Object* objects, Light* lights, KD::Tree<CORE> &kdtree, int total)
{
    


    float x, y, z;
    Ray photonRay;
    Colour colour;
    float ne = 0; //number of photons
    Hit best_hit;
    Light* light = lights;
    
    

    while (light != (Light*)0)
    {
        //Photons emitted randomly in all directions
        while (ne < total)
        {
            x = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
            y = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
            z = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
            while (x * x + y * y + z * z > 1.0f)
            {
                x = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
                y = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
                z = (((float)rand() / (RAND_MAX)) * 2.0f) - 1.0f;
            }
            photonRay.direction = Vector(x, y, z);
            lights->get_position(photonRay.position);
            int depth = 1;
            float power[3] = { 1.0f*12.0f / total , 1.0f * 12.0f / total  , 1.0f * 12.0f / total };
            photontrace(photonRay, objects, lights, power, kdtree, depth);
        
            ne += 1;
        }
        light = light->next;
    }
    cout << "done global map" << endl;
}

// emitting photons into the scene for caustics
void emit_photons_caustic(Object* objects, Light* lights, KD::Tree<CORE> &kdtree_caustic, int total)
{



    float x, y, z, scale, correctx, correcty, correctz, refractionCoeff;
    Ray photonRay;
    Colour colour;
    float ne = 0; //number of photons
    Hit best_hit;
    Light* light = lights;
    //scalers to zone in on transmissive object
    scale = 2.0f;
    correctx = 1.0f;
    correcty = 1.0f;
    correctz = 1.0f;
    
    while (light != (Light*)0)
    {
        while (ne < total)
        {
            //Photons emitted randomly in all directions until it hits transmissive object
            // at which it will begin trying to send out more photons in that area
            x = (((float)rand() / (RAND_MAX)) * scale) - correctx;
            y = (((float)rand() / (RAND_MAX)) * scale) - correcty;
            z = (((float)rand() / (RAND_MAX)) * scale) - correctz;
            while (x * x + y * y + z * z > 1.0f)
            {
                x = (((float)rand() / (RAND_MAX)) * scale) - correctx;
                y = (((float)rand() / (RAND_MAX)) * scale) - correcty;
                z = (((float)rand() / (RAND_MAX)) * scale) - correctz;
            }


            photonRay.direction = Vector(x, y, z);
            lights->get_position(photonRay.position);
            
            float power[3] = { 1.0f / total, 1.0f / total , 1.0f/ total };

            object_test(photonRay, objects, best_hit);

            //to fire next photon close to transmissive object
            if (best_hit.flag && best_hit.what->material->refraction > 0.0f)
            {

                scale = 0.1f;
                correctx = -(x - 0.05f);
                correcty = -(y - 0.05f);
                correctz = -(z - 0.05f);
                int depth = 1;
                photontrace(photonRay, objects, lights, power, kdtree_caustic, depth);
                ne += 1;
            }

            

        }
        light = light->next;
    }
    cout << "done caustics" << endl;
}

//Depth of field function that will jitter our ray
void depthOfField(Ray &ray,Ray &depthRay,float focalLength, float maxAperture)
{
    Vertex focalPoint;
    Vertex centre;
    float aperture;
    float theta;
    float t_imagePlane;
    Vector normal;
    normal = Vector(0, 0, -1); 
    //calculating distance to image plane
    t_imagePlane = -(normal.dot(Vector(ray.position.x, ray.position.y, ray.position.z)) - 1) / (normal.dot(ray.direction)); 
    //centre is where the pixel is, we will perform jittering around this point
    centre.x = ray.position.x + (t_imagePlane * ray.direction.x);
    centre.y = ray.position.y + (t_imagePlane * ray.direction.y);
    centre.z = ray.position.z + (t_imagePlane * ray.direction.z);
    //calculating focal point using focal length
    focalPoint.x = ray.position.x + (focalLength * ray.direction.x);
    focalPoint.y = ray.position.y + (focalLength * ray.direction.y);
    focalPoint.z = ray.position.z + (focalLength * ray.direction.z);
    //random point on the aperture
    aperture = ((double)rand() / (RAND_MAX)) * maxAperture;
    //random angle
    theta = ((double)rand() / (RAND_MAX)) * 2 * 3.14159265359f;
    //depthRay position is now our jittered point
    depthRay.position.x = centre.x + aperture * cos(theta);
    depthRay.position.y = centre.y + aperture * sin(theta);
    depthRay.position.z = 1;
    depthRay.direction.x = focalPoint.x - depthRay.position.x;
    depthRay.direction.y = focalPoint.y - depthRay.position.y;
    depthRay.direction.z = focalPoint.z - depthRay.position.z;
    depthRay.direction.normalise();
}



int main(int argc, char *argv[])
{
  int width = 512;
  int height = 512;
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
				       0.0f, 0.0f, 1.0f, -2.7f,
				       0.0f, 1.0f, 0.0f, 7.0f,
				       0.0f, 0.0f, 0.0f, 1.0f);

  Transform* transform2 = new Transform(1.0f, 0.0f, 0.0f, 0.0f,
      0.0f, 1.0f, 0.0f, 0.0f,
      0.0f, 0.0f, 1.0f, 0.0f,
      0.0f, 0.0f, 0.0f, 1.0f);



  //  Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);

  PolyMesh* floorpm = new PolyMesh((char*)"floor.ply", transform2);
  PolyMesh* leftWall = new PolyMesh((char*)"leftWallRed.ply", transform2);
  PolyMesh* rightWall = new PolyMesh((char*)"rightWallGreen.ply", transform2);
  PolyMesh* backWall = new PolyMesh((char*)"backWallWhite.ply", transform2);
  PolyMesh* ceiling = new PolyMesh((char*)"ceiling.ply", transform2);


  Vertex v;
  v.x = -3.0f;
  v.y = -1.9f;
  v.z = 5.3f;
  
  Sphere *sphere = new Sphere(v, 1.0f);
  
  Vertex v2;
  v2.x = 3.0f;
  v2.y = -1.4f;
  v2.z = 5.3f;
  
  Sphere *sphere2 = new Sphere(v2,1.0f);

  sphere->next = pm;

  Ray ray;
  ray.position.x = 0.0001f;
  ray.position.y = -0.5f;
  ray.position.z = 0.0f;

  DirectionalandPointLight* dl = new DirectionalandPointLight(Vector(1.01f, -1.0f, 1.0f), Colour(1.0f, 1.0f, 1.0f, 0.0f), Vertex(0.01f, 6.0f, 7.0f));
  
  Phong bp1;

	bp1.ambient.r = 0.0f;
	bp1.ambient.g = 0.0f;
	bp1.ambient.b = 0.2f;
	bp1.diffuse.r = 0.0f;
	bp1.diffuse.g = 0.0f;
	bp1.diffuse.b = 0.8f;
	bp1.specular.r = 0.8f;
	bp1.specular.g = 0.8f;
	bp1.specular.b = 0.8f;
	bp1.power = 40.0f;

	pm->material = &bp1;
    pm->material->reflection = 0.7f;
    pm->material->refraction = 0.0f;
    pm->material->refracI = 0.0f;
  
  Phong floor;
    floor.ambient.r = 0.2f;
    floor.ambient.g = 0.2f;
    floor.ambient.b = 0.2f;
    floor.diffuse.r = 1.0f;
    floor.diffuse.g = 1.0f;
    floor.diffuse.b = 1.0f;
    floor.specular.r = 0.2f;
    floor.specular.g = 0.2f;
    floor.specular.b = 0.2f;
    floor.power = 20.0f;
    floorpm->material = &floor;
    floorpm->material->reflection = 0.0f;
    floorpm->material->refraction = 0.0f;

    Phong left;
    left.ambient.r = 0.2f;
    left.ambient.g = 0.0f;
    left.ambient.b = 0.0f;
    left.diffuse.r = 1.0f;
    left.diffuse.g = 0.0f;
    left.diffuse.b = 0.0f;
    left.specular.r = 0.0f;
    left.specular.g = 0.2f;
    left.specular.b = 0.0f;
    left.power = 20.0f;
    leftWall->material = &left;
    leftWall->material->reflection = 0.0f;
    leftWall->material->refraction = 0.0f;

    Phong right;
    right.ambient.r = 0.0f;
    right.ambient.g = 0.2f;
    right.ambient.b = 0.0f;
    right.diffuse.r = 0.0f;
    right.diffuse.g = 1.0f;
    right.diffuse.b = 0.0f;
    right.specular.r = 0.0f;
    right.specular.g = 0.2f;
    right.specular.b = 0.0f;
    right.power = 20.0f;
    rightWall->material = &right;
    rightWall->material->reflection = 0.0f;
    rightWall->material->refraction = 0.0f;


    backWall->material = &floor;
    backWall->material->reflection = 0.0f;
    backWall->material->refraction = 0.0f;


    ceiling->material = &floor;
    ceiling->material->reflection = 0.0f;
    ceiling->material->refraction = 0.0f;

  Phong bp3;

    bp3.ambient.r = 0.0f;
	bp3.ambient.g = 0.0f;
	bp3.ambient.b = 0.0f;
	bp3.diffuse.r = 0.0f;
	bp3.diffuse.g = 0.0f;
	bp3.diffuse.b = 0.0f;
	bp3.specular.r = 0.8f;
	bp3.specular.g = 0.8f;
	bp3.specular.b = 0.8f;
	bp3.power = 180.0f;

	sphere->material = &bp3;
    sphere->material->reflection = 1.0f;
    sphere->material->refraction = 0.0f;

 Phong bp2;

    bp2.ambient.r = 0.0f;
    bp2.ambient.g = 0.0f;
    bp2.ambient.b = 0.0f;
    bp2.diffuse.r = 0.0f;
    bp2.diffuse.g = 0.0f;
    bp2.diffuse.b = 0.0f;
    bp2.specular.r = 0.1f;
    bp2.specular.g = 0.1f;
    bp2.specular.b = 0.1f;
    bp2.power = 180.0f;
 	sphere2->material = &bp2;
    sphere2->material->reflection = 0.0f;
    sphere2->material->refraction = 1.0f;
    sphere2->material->refracI = 1.5f;

	//pm->next = floorpm;
    //sphere3->next = sphere;
    pm->next = sphere2;
    sphere2->next = floorpm;
    floorpm->next = leftWall;
    leftWall->next = rightWall;
    rightWall->next = backWall;
    backWall->next = ceiling;

    
    Point min(0, 0, 0);
    Point max(30, 30, 30);
    KD::Tree<CORE> kdtree(min, max);
    KD::Tree<CORE> kdtree_caustic(min, max);
    emit_photons(sphere, dl, kdtree, 10000);
    emit_photons_caustic(sphere, dl, kdtree_caustic, 40000);
  for (int y = 0; y < height; y += 1)
  {
    for (int x = 0; x < width; x += 1)
    {
      float fx = (float)x/(float)width;
      float fy = (float)y/(float)height;

      Vector direction;

      ray.direction.x = (fx-0.5f);
      ray.direction.y = (0.5f-fy);
      ray.direction.z = 0.5f;
      ray.direction.normalise();

      Colour colour, colourDepth;
      Ray depthRay;
      int depth = 6;
      //depth of field
      
      for (int i = 0; i < 26; i++)
      {
          depthOfField(ray, depthRay, 5.2, 0.1);
          findphotons(depthRay, sphere, dl, colour, 6, kdtree, kdtree_caustic);
          colourDepth.r += colour.r;
          colourDepth.g += colour.g;
          colourDepth.b += colour.b;
          //cout << colourDepth.r << colourDepth.g << colourDepth.b << flush;
      }
      colourDepth.r = colourDepth.r / 17;
      colourDepth.g = colourDepth.g / 17;
      colourDepth.b = colourDepth.b / 17;
      fb->plotPixel(x, y, std::min(1.0f, colourDepth.r), std::min(1.0f, colourDepth.g), std::min(1.0f, colourDepth.b));
      

      //findphotons(ray, sphere, dl, colour, 6, kdtree, kdtree_caustic);

      //fb->plotPixel(x, y, std::min(1.0f, colour.r), std::min(1.0f, colour.g), std::min(1.0f, colour.b));
      //fb->plotDepth(x,y, depth);
    }

    cerr << "*" << flush;
  }
  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");
  //  fb->writeDepthFile((char *)"depth.ppm");
  return 0;
  
}
