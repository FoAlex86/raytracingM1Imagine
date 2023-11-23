#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include "BoundingBox.h"


#include <GL/glut.h>


enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

    std::vector< BoundingBox > boundingBox;

    //paramètre pour depthOfField
    float focus_distance = 0; // distance de mise au point en mètres
    float aperture_size = 0; // taille de l'ouverture en millimètres
    float maxClarity = 0; //permet de contrôler la distance maximale qui est focus poru faire un flou de fon
    bool depthOfField = false; // depthOfField
    bool backDepthOfField = false; // back focus limited

public:


    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }

        /* //draw boundingbox
        for(size_t It = 0; It < boundingBox.size(); ++It)
        {
            boundingBox[It].draw(boundingBox[It]);
        }*/
    }




    RaySceneIntersection computeIntersection(Ray const & ray, float znear, float zfar) {
        RaySceneIntersection result;
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche
        result.intersectionExists = false;

        size_t size = meshes.size();
        size_t boundingBoxSize = boundingBox.size();

        for(size_t boundingBoxIndex = 0; boundingBoxIndex < boundingBoxSize; boundingBoxIndex++)
        {
            if(boundingBox[boundingBoxIndex].intersects(ray))
            {
                for (size_t i = 0; i < size; i++)
                {
                    RayTriangleIntersection mesh = meshes[i].intersect(ray);
                    if(mesh.intersectionExists)
                    {
                        if(mesh.t < result.t && mesh.t < zfar && mesh.t > znear)
                        {
                            result.intersectionExists = mesh.intersectionExists;
                            result.typeOfIntersectedObject = 0;
                            result.objectIndex = i;
                            result.t = mesh.t;
                            result.rayMeshIntersection = mesh;

                        }
                    }
                }
            }
        }

        size = spheres.size();
        for(size_t i = 0; i < size; i++)
        {
            RaySphereIntersection sphere = spheres[i].intersect(ray);
            if(sphere.intersectionExists)
            {
                if(sphere.t < result.t && sphere.t < zfar)
                {
                    result.intersectionExists = sphere.intersectionExists;
                    result.typeOfIntersectedObject = 1;
                    result.objectIndex = i;
                    result.t = sphere.t;
                    result.raySphereIntersection = sphere;
                }     
            }
        }

        size = squares.size();
        for(size_t i = 0; i < size; i++)
        {
                
            RaySquareIntersection square = squares[i].intersect(ray);
            if(square.intersectionExists)
            {
                if(square.t < result.t && square.t > znear && square.t < zfar)
                {
                    result.intersectionExists = square.intersectionExists;
                    result.typeOfIntersectedObject = 2;
                    result.objectIndex = i;
                    result.t = square.t;
                    result.raySquareIntersection = square; 
                }
            }
        }
        
        return result;
    }

    bool computeShadowIntersection(Ray const & ray)
    {
        RaySceneIntersection result;
        bool blocked = false;
        float epsilon = 0.001;
        float bias= 1.;
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche
        size_t size = spheres.size();

        for(size_t i = 0; i < size; i++)
        {
            RaySphereIntersection sphere = spheres[i].intersect(ray);
            if(sphere.intersectionExists && sphere.t > epsilon && sphere.t < bias)
            {
                if(spheres[i].material.type != Material_Glass)
                {
                   blocked = true;
                    break; 
                }
            }
        }

        if(!blocked)
        {
            size = squares.size();
            for(size_t i = 0; i < size; i++)
            {       
                RaySquareIntersection square = squares[i].intersect(ray);
                if(square.intersectionExists && square.t > epsilon && square.t < bias)
                {
                    if(squares[i].material.type != Material_Glass)
                    {
                        blocked = true;
                        break;
                    }
                }
            }
        }

        /*
        if(!blocked)
        {
            size_t boundingBoxSize = boundingBox.size();
            for(size_t boundingBoxIndex = 0; boundingBoxIndex < boundingBoxSize; boundingBoxIndex++)
            {
                if(boundingBox[boundingBoxIndex].intersects(ray))
                {
                    size = meshes.size();
                    for (size_t i = 0; i < size; i++)
                    {
                        RayTriangleIntersection mesh = meshes[i].intersect(ray);
                        if(mesh.intersectionExists)
                        {
                            if(mesh.t > epsilon && mesh.t < bias)
                            {
                                if(meshes[i].material.type != Material_Glass)
                                {
                                    blocked = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            } 
        }*/
        
        return blocked;
    }

    /*float isShadowed(Ray const & shadow){
        int cptShade = 0;
        int cptTotalShade = 0;
        bool blocked = false;

        //ici faire boucle for pour parcouris les multiples rayons pour ombre douces
         https://www.cs.unc.edu/~dm/UNC/COMP236/LECTURES/SoftShadows.pdf
        for(int i = 0; i < light.quad.size; i++)
        {
            for(int j = 0; j < light.quad.size; j++)
            {
                float x = (float)rand()
            }

        }
        cptTotalShade++;
        blocked = computeShadowIntersection(shadow);

        if(!blocked)
        {
            cptShade++;
        }
        //fin de la boucle
        // 1-
        return ((float)cptShade)/((float)cptTotalShade);
    }*/

    float isShadowedSoft(int nLight, int cptTotalShade, Vec3 objectIntersect)
    {
        int cptShade = 0;

        float lightRadius = lights[nLight].radius/2;
        float x = lights[nLight].pos[0] - lightRadius;
        float y = lights[nLight].pos[1] - lightRadius;
        float z = lights[nLight].pos[2] - lightRadius;

        float diam = lights[nLight].radius;

        for(int i = 0; i < cptTotalShade; i++)
        {
            float pX = (float)(rand()/(float)(RAND_MAX / diam));
            float pZ = (float)(rand()/(float)(RAND_MAX / diam));

            Vec3 L = Vec3(x + pX, y, z + pZ) - objectIntersect;
            L.normalize();

            Ray shadow = Ray(objectIntersect, L);
            if(!computeShadowIntersection(shadow))
            {
                cptShade++;
            }
        }
        
        return ((float)cptShade)/((float)cptTotalShade);
    }


    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces, float znear, float zfar ) {

        //TODO RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        RaySceneIntersection raySceneIntersection = computeIntersection(ray, znear, zfar);
        Vec3 color(0., 0., 0.); 
        if(raySceneIntersection.intersectionExists)
        {
            Vec3 objectColor, objectIntersect, objectNormal, objectAmbientMaterial, objectDifuseMaterial, objectSpecularMaterial;
            double objectShininess;
            float objectIndexMedium, objectTransparency;
            MaterialType objectType;

            switch ( raySceneIntersection.typeOfIntersectedObject )
            {  
                case 0:
                {
                    objectIntersect = raySceneIntersection.rayMeshIntersection.intersection;
                    objectNormal = raySceneIntersection.rayMeshIntersection.normal;
                    objectAmbientMaterial = meshes[raySceneIntersection.objectIndex].material.ambient_material;
                    objectDifuseMaterial = meshes[raySceneIntersection.objectIndex].material.diffuse_material;
                    objectSpecularMaterial = meshes[raySceneIntersection.objectIndex].material.specular_material;
                    objectShininess = meshes[raySceneIntersection.objectIndex].material.shininess;
                    objectType = meshes[raySceneIntersection.objectIndex].material.type;
                    if(objectType != Material_Diffuse_Blinn_Phong)
                    {
                        objectTransparency = meshes[raySceneIntersection.objectIndex].material.transparency;
                        objectIndexMedium = meshes[raySceneIntersection.objectIndex].material.index_medium;
                    }
                    break;
                }
                case 1:
                {
                    objectIntersect = raySceneIntersection.raySphereIntersection.intersection;
                    objectNormal = raySceneIntersection.raySphereIntersection.normal;
                    objectAmbientMaterial = spheres[raySceneIntersection.objectIndex].material.ambient_material;
                    objectDifuseMaterial = spheres[raySceneIntersection.objectIndex].material.diffuse_material;
                    objectSpecularMaterial = spheres[raySceneIntersection.objectIndex].material.specular_material;
                    objectShininess = spheres[raySceneIntersection.objectIndex].material.shininess;
                    objectType = spheres[raySceneIntersection.objectIndex].material.type;
                    if(objectType != Material_Diffuse_Blinn_Phong)
                    {
                        objectTransparency = spheres[raySceneIntersection.objectIndex].material.transparency;
                        objectIndexMedium = spheres[raySceneIntersection.objectIndex].material.index_medium;
                    }
                    break;
                }
                case 2:
                {
                    objectIntersect = raySceneIntersection.raySquareIntersection.intersection;
                    objectNormal = raySceneIntersection.raySquareIntersection.normal; 
                    objectAmbientMaterial = squares[raySceneIntersection.objectIndex].material.ambient_material;
                    objectDifuseMaterial = squares[raySceneIntersection.objectIndex].material.diffuse_material;
                    objectSpecularMaterial = squares[raySceneIntersection.objectIndex].material.specular_material;
                    objectShininess = squares[raySceneIntersection.objectIndex].material.shininess;
                    objectType = squares[raySceneIntersection.objectIndex].material.type;
                    if(objectType != Material_Diffuse_Blinn_Phong)
                    {
                        objectTransparency = squares[raySceneIntersection.objectIndex].material.transparency;
                        objectIndexMedium = squares[raySceneIntersection.objectIndex].material.index_medium;
                    }
                    break;
                }
                default:
                {
                    break;
                }
            }


            switch ( objectType )
            {  
                case Material_Diffuse_Blinn_Phong:
                {
                    //ia
                        for(int rgb = 0; rgb < 3; rgb++)
                        {
                            color[rgb] = objectAmbientMaterial[rgb] * lights[0].material[rgb];
                        }         

                        
                        for(size_t l = 0; l < lights.size(); l++)
                        {
                            Vec3 L = lights[l].pos - objectIntersect;
                            L.normalize();
                            float theta = std::max((Vec3::dot(L, objectNormal)), (float)0);
                            Vec3 R = 2*theta*objectNormal - L;
                            R.normalize();
                            Vec3 V = ray.origin() - objectIntersect;
                            V.normalize();
                            float alpha = std::max(Vec3::dot(R, V), (float)0);
                            
                            for(int rgb = 0; rgb < 3; rgb++)
                            {
                                float id = objectDifuseMaterial[rgb] * lights[l].material[rgb] * theta;
                                float is = objectSpecularMaterial[rgb] * lights[l].material[rgb] * pow(alpha, objectShininess);
                                color[rgb] += id + is;
                            }
                        }
                    break;
                }
                case Material_Glass:
                {
                    if(NRemainingBounces > 0)
                    {
                        float dot = Vec3::dot(objectNormal, ray.direction());
                        float x = (objectIndexMedium * objectIndexMedium) * (dot * dot);

                        Vec3 rayDir;
                        if(x >= 0.0){
                            rayDir = objectIndexMedium * ray.direction() + ((objectIndexMedium * (dot + sqrt(x))) * objectNormal);
                            rayDir.normalize();
                        }
                                
                        Ray newRay = Ray(objectIntersect, rayDir);
                        color = rayTraceRecursive(newRay, NRemainingBounces-1, 0.f, zfar);
                    }
                    break;
                }
                case Material_Mirror:
                {
                    if(NRemainingBounces > 0)
                    {
                        Vec3 Vm = ray.direction();
                        Vm.normalize();
                        float cos_teta2 = Vec3::dot(Vm, objectNormal);
                        Vec3 vecteurMirror = (2 * cos_teta2 * objectNormal) - Vm;
                        vecteurMirror.normalize();
                        Ray rayMirror(objectIntersect, -1*vecteurMirror);
                        color = rayTraceRecursive(rayMirror, NRemainingBounces - 1, 0.00001f, zfar);
                    }
                    break;
                }
                default:
                {
                    break;
                }
            }
            
            for(size_t l = 0; l < lights.size(); l++)
            {
                //SHADOW
                //Hard
                //Ray shadow = Ray(objectIntersect, L);
                //float percentShade = isShadowed(shadow);
                
                //Soft
                float percentShade = isShadowedSoft(l, 10, objectIntersect);
                color = color * percentShade;
            }
        }

        

        return color;
    }

    float rayTraceDepthOfFieldIntersect( Ray ray , int NRemainingBounces, float znear, float zfar )
    {
        RaySceneIntersection raySceneIntersection = computeIntersection(ray, znear, zfar);
        float objectIntersect;

        if(raySceneIntersection.intersectionExists)
        {
            switch ( raySceneIntersection.typeOfIntersectedObject )
            {  
                case 0:
                {
                    objectIntersect = raySceneIntersection.rayMeshIntersection.t;
                    break;
                }
                case 1:
                {
                    objectIntersect = raySceneIntersection.raySphereIntersection.t;
                    break;
                }
                case 2:
                {
                    objectIntersect = raySceneIntersection.raySquareIntersection.t;
                    break;
                }
                default:
                {
                    break;
                }
            }
        }

        return objectIntersect;
    }

    Vec3 appliesDepthOfField( Ray const & rayStart, int NRemainingBounces, float znear, float zfar)
    {
        Vec3 color = rayTraceRecursive(rayStart, NRemainingBounces, znear, zfar);
        float blur_radius = (1.0 / aperture_size) * focus_distance; // rayon de confusion en mètres

        float result = rayTraceDepthOfFieldIntersect(rayStart, NRemainingBounces, znear, zfar);
        float distance_to_focus = abs(result - focus_distance);
        float o_distance_to_focus = abs(result + focus_distance);

        if (distance_to_focus < blur_radius || (backDepthOfField && o_distance_to_focus > blur_radius + maxClarity))
        {
            Vec3 colorN = Vec3(0.,0.,0.);
            //color = colorN;
            
            for(size_t neighbor = 0; neighbor < 8; neighbor++)
            {
                Vec3 nDir = nDir.nRandom(rayStart.direction());
                Ray rayN = Ray(rayStart.origin(), nDir);
                colorN += rayTraceRecursive(rayN, NRemainingBounces, znear, zfar);
            }

            color += colorN;
            color/=9;
        }

        return color;
    }

    Vec3 rayTrace( Ray const & rayStart, float znear, float zfar) {
        //TODO appeler la fonction recursive
        Vec3 color;
        if(depthOfField)
        {
            color = appliesDepthOfField(rayStart, 5, znear, zfar);
        }
        else
        {
            color = rayTraceRecursive(rayStart, 16, znear, zfar);
        }
        
        return color;
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 0.8,0.8,0.8 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 0.,0.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3(1.,1.,.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }

    void setup_single_mesh() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            m.loadOFF("data/suzanne.off");
            m.centerAndScaleToUnit();
            m.material.type = Material_Diffuse_Blinn_Phong;
            m.material.diffuse_material = Vec3(1., 0., 0.);
            m.material.specular_material = Vec3(0.2, 0.2, 0.2);
            m.material.shininess = 20;
            m.build_arrays();

            BoundingBox boundingBoxMesh;
            std::pair<std::array<float, 3>, std::array<float, 3>> boundsMesh = boundingBoxMesh.getBounds(m);
            std::array<float, 3> min = boundsMesh.first;
            std::array<float, 3> max = boundsMesh.second;
            BoundingBox boundingBoxBound(min, max);
            boundingBoxMesh.expand(boundingBoxBound);            
            boundingBox.push_back(boundingBoxMesh);
        }
    }

    void setup_cornell_mesh()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 0.,0.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3(1.,1.,.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        {
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            m.loadOFF("data/suzanne.off");
            m.centerAndScaleToUnit();
            //m.scale(Vec3(0.5, 0.5, 0.5));
            //m.translate(Vec3(0., 0., 0.));
            m.material.type = Material_Diffuse_Blinn_Phong;
            m.material.diffuse_material = Vec3(1., 0., 0.);
            m.material.specular_material = Vec3(1., 1., 1.);
            m.material.shininess = 16;
            m.build_arrays();

            BoundingBox boundingBoxMesh;
            std::pair<std::array<float, 3>, std::array<float, 3>> boundsMesh = boundingBoxMesh.getBounds(m);
            std::array<float, 3> min = boundsMesh.first;
            std::array<float, 3> max = boundsMesh.second;
            BoundingBox boundingBoxBound(min, max);
            boundingBoxMesh.expand(boundingBoxBound);            
            boundingBox.push_back(boundingBoxMesh);
        }
    }

    void setup_cornell_box_focal()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        depthOfField = true;
        backDepthOfField = true;
        focus_distance = 3;
        aperture_size = 1;
        maxClarity = 7.5;

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 0.,0.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3(1.,1.,.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            //s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }

    void setup_default_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

    }

    void setup_cornell_finale()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,0.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1.,1.,.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        {
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            m.loadOFF("data/suzanne.off");
            m.centerAndScaleToUnit();
            //m.scale(Vec3(0.5, 0.5, 0.5));
            //m.translate(Vec3(0., 0., 0.));
            m.material.type = Material_Diffuse_Blinn_Phong;
            m.material.diffuse_material = Vec3(1., 0., 0.);
            m.material.specular_material = Vec3(1., 1., 1.);
            m.material.shininess = 16;
            m.build_arrays();

            BoundingBox boundingBoxMesh;
            std::pair<std::array<float, 3>, std::array<float, 3>> boundsMesh = boundingBoxMesh.getBounds(m);
            std::array<float, 3> min = boundsMesh.first;
            std::array<float, 3> max = boundsMesh.second;
            BoundingBox boundingBoxBound(min, max);
            boundingBoxMesh.expand(boundingBoxBound);            
            boundingBox.push_back(boundingBoxMesh);
        }
    }

};



#endif
