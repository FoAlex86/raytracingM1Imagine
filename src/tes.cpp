#include <vector>
#include <iostream>
#include "Vec3.h"

int main () {
    Vec3 L = Vec3(-1,2,-1);
    Vec3 N = Vec3(0,1,0);
    float theta = Vec3::dot(L, N);

    Vec3 R = 2*theta*N - L;
    R.normalize();

    Vec3 V = Vec3(1,1.5,0.5);

    float alpha = Vec3::dot(R, V);
                            /*
    for(int rgb = 0; rgb < 3; rgb++)
    {
        float id = objectDifuseMaterial[rgb] * lights[l].material[rgb] * theta;
        float is = objectSpecularMaterial[rgb] * lights[l].material[rgb] * pow(alpha, objectShininess);
        color[rgb] += id + is;
    }*/

    std::cout << R[0] << " , " << R[1] << " , " << R[2] << std::endl;
    std::cout << alpha << std::endl;
    std::cout << 1*0.9*pow(alpha,5);

    return 0;
}