#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
    RaySquareIntersection() : intersectionExists(false), t(FLT_MAX) {}
};



class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    Vec3 bottomLeft() const { return this->vertices[0].position; }
    Vec3 bottomRight() const { return this->vertices[1].position; }
    Vec3 upRight() const { return this->vertices[2].position; }
    Vec3 upLeft() const { return this->vertices[3].position; }
    Vec3 normal() const { return Vec3::cross((bottomRight() - bottomLeft()) , (upLeft() - bottomLeft())); }


    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;


    }

    RaySquareIntersection intersect(const Ray &ray) const {
        RaySquareIntersection intersection;
        intersection.intersectionExists = false;
        //TODO calculer l'intersection rayon quad
        Vec3 p = bottomLeft() - ray.origin();
        Vec3 norm = normal();

        float t = Vec3::dot(p, norm)/Vec3::dot(ray.direction(), norm);

        if(t>0)
        {
            Vec3 intersection_p = ray.origin()+ t*ray.direction();
            Vec3 vertice0, vertice1, vertice2,vertice3;
            vertice0 = Vec3::cross(bottomRight() - bottomLeft(), intersection_p - bottomLeft());
            vertice1 = Vec3::cross(upRight() - bottomRight(), intersection_p - bottomRight());
            vertice2 = Vec3::cross(upLeft() - upRight(), intersection_p - upRight());
            vertice3 = Vec3::cross(bottomLeft() - upLeft(), intersection_p - upLeft());
            
            bool signe = Vec3::dot(vertice0, norm) > 0;

            if( ((Vec3::dot(vertice1, norm) > 0) == signe) && ((Vec3::dot(vertice2, norm) > 0) == signe) && ((Vec3::dot(vertice3, norm) > 0) == signe))
            {
                intersection.intersectionExists = true;
                intersection.t = t;
                //intersection.u = ;
                //intersection.v = ;
                intersection.intersection = intersection_p;
                norm.normalize();
                intersection.normal = norm;
            }
        }


        return intersection;
    }
};
#endif // SQUARE_H
