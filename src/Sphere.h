#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection{
    bool intersectionExists;
    float t;
    float theta,phi;
    float u, v;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
    RaySphereIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};

static
Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return ThetaPhiR[2] * Vec3( cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[1]) );
}
static
Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
    return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
}

static
Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
    float R = xyz.length();
    float phi = asin( xyz[2] / R );
    float theta = atan2( xyz[1] , xyz[0] );
    return Vec3( theta , phi , R );
}

static void get_sphere_uv(const Vec3& p, float& u, float& v)
{
    // p: a given point on the sphere of radius one, centered at the origin.
    // u: returned value [0,1] of angle around the Y axis from X=-1.
    // v: returned value [0,1] of angle from Y=-1 to Y=+1.
    //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
    //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
    //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

    auto theta = acos(-p[1]);
    auto phi = atan2(-p[2], p[0]) + M_PI;

    u = phi / (2*M_PI);
    v = theta / M_PI;
}

class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c , float r) : Mesh() , m_center(c) , m_radius(r) {}

    void build_arrays(){
        unsigned int nTheta = 20 , nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi );
        normalsArray.resize(3 * nTheta * nPhi );
        uvs_array.resize(2 * nTheta * nPhi );
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta ; ++thetaIt ) {
            float u = (float)(thetaIt) / (float)(nTheta-1);
            float theta = u * 2 * M_PI;
            for( unsigned int phiIt = 0 ; phiIt < nPhi ; ++phiIt ) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi-1);
                float phi = - M_PI/2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean( theta , phi );
                positions_array[ 3 * vertexIndex + 0 ] = m_center[0] + m_radius * xyz[0];
                positions_array[ 3 * vertexIndex + 1 ] = m_center[1] + m_radius * xyz[1];
                positions_array[ 3 * vertexIndex + 2 ] = m_center[2] + m_radius * xyz[2];
                normalsArray[ 3 * vertexIndex + 0 ] = xyz[0];
                normalsArray[ 3 * vertexIndex + 1 ] = xyz[1];
                normalsArray[ 3 * vertexIndex + 2 ] = xyz[2];
                uvs_array[ 2 * vertexIndex + 0 ] = u;
                uvs_array[ 2 * vertexIndex + 1 ] = v;
            }
        }
        triangles_array.clear();
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta - 1 ; ++thetaIt ) {
            for( unsigned int phiIt = 0 ; phiIt < nPhi - 1 ; ++phiIt ) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt+1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt+1) * nTheta;
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuV );
            }
        }
    }

    RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection raySphereIntersection;

        float a, b, c;
        // Delta au carré et les deux racines
        float delta; 
        float root1;
        float root2;
        
        // D'après la formule du cours de l'intersection rayon-sphère page 42, on commence par déterminer les coefficients calculables pour exprimer en fonction de t
        a = Vec3::dot(ray.direction(),ray.direction());
        b = 2 * Vec3::dot(ray.direction(),(ray.origin() - this->m_center));
        c = (ray.origin() - this->m_center).squareNorm() - pow(this->m_radius,2);

        delta = pow(b,2) - (4 * a * c);
        
        if(delta > 0) {
            root1 = ( -b + sqrt(delta) ) / ( 2 * a );
            root2 = ( -b - sqrt(delta) ) / ( 2 * a );

            if(root1 < root2) {
                raySphereIntersection.t = root1;

            } else {
                raySphereIntersection.t = root2;
            }

            if(raySphereIntersection.t < 0.00001) {
                raySphereIntersection.t = FLT_MAX;
            }
            else{
                raySphereIntersection.intersectionExists = true;
                raySphereIntersection.intersection = ray.origin() + raySphereIntersection.t * ray.direction();
                raySphereIntersection.secondintersection = ray.origin() + root2 * ray.direction();
                raySphereIntersection.normal = raySphereIntersection.intersection - this->m_center;
                raySphereIntersection.normal.normalize();
                get_sphere_uv(raySphereIntersection.normal, raySphereIntersection.u, raySphereIntersection.v);
            }
        } 

        return raySphereIntersection;
    }
    /*RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection intersection;
        //TODO calcul l'intersection rayon sphere
        Vec3 oc = ray.origin() - m_center;
        float a = ray.direction().squareLength();
        float half_b = Vec3::dot(oc, ray.direction());
        float c = oc.squareLength() - m_radius * m_radius;

        float delta = half_b * half_b - a*c;

        if( delta < 0.0)
        {
            intersection.intersectionExists = false;
        }
        else
        {
            float sqrtd = sqrt(delta);
            float root;
            float rootBis = -1;
            if(delta == 0.0)
            {
                intersection.intersectionExists = true;
                root = (-half_b - sqrtd) / 2*a;
            }
            else
            {
                float t1 = (-half_b - sqrtd) / a;
                float t2 = (-half_b + sqrtd) / a;
                if(t1 < 0 && t2 < 0)
                {
                    intersection.intersectionExists = false;
                }
                else
                {
                    intersection.intersectionExists = true;
                    if(t1 < 0)
                    {
                        root = t2;
                    }
                    else
                    {
                        if(t2 < 0)
                        {
                            root = t1;
                        }
                        else
                        {
                            root = std::min(t1, t2);
                            rootBis = std::max(t1, t2);
                        }
                    }
                }
            }
            if(root < 0.00001) {
                intersection.t = FLT_MAX;
                intersection.intersectionExists = false;
            }

            if(intersection.intersectionExists)
            {
                    intersection.t = root;
                    intersection.intersection = ray.origin() + (root*ray.direction());
                    if(rootBis != -1)
                    {
                        intersection.secondintersection = ray.origin() + (rootBis*ray.direction());
                    }
                    Vec3 norm = intersection.intersection - m_center;
                    norm.normalize();
                    intersection.normal = norm;
                    //Vec3 sphericalCoordinates = EuclideanCoordinatesToSpherical( intersection.intersection );
                    //intersection.theta = sphericalCoordinates[0];
                    //intersection.phi = sphericalCoordinates[1];
            }
            
        }
        return intersection;
    }*/
};
#endif
