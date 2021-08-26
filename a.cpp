
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

typedef struct double3
{
    double x,y,z;
} double3;

typedef struct int3
{
    int x,y,z;
} int3;

typedef struct int4
{
    int x,y,z,w;
} int4;

double my_tetrahedron_vertices[12] = {
	0.0, 0.0 ,0.0,
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0
};

/// https://github.com/mapengfei-nwpu/outward_normal/blob/master/outward_normal.cpp
double3 triangle_normal(const double3* points, int opposite){
    /// for example, opposite is 1, the face will be {0,3,2}, 
    /// face[3] will never be used, points[1] will be used to decide normal direction.
    int face[4] = {0,1,2,3};
    face[opposite] = face[3];

    /// 
    const double* a = (double*)&(points[face[0]]);
    const double* b = (double*)&(points[face[1]]);
    const double* c = (double*)&(points[face[2]]);
    const double* d = (double*)&(points[opposite]);
    
    double3 result;
    double* normal = (double*)&result;

    double ab[3], ac[3], ad[3];

    ab[0] = b[0] - a[0];
    ab[1] = b[1] - a[1];
    ab[2] = b[2] - a[2];

    ac[0] = c[0] - a[0];
    ac[1] = c[1] - a[1];
    ac[2] = c[2] - a[2];

    ad[0] = d[0] - a[0];
    ad[1] = d[1] - a[1];
    ad[2] = d[2] - a[2];

    normal[0] = ab[1] * ac[2] - ab[2] * ac[1];
    normal[1] = ab[2] * ac[0] - ab[0] * ac[2];
    normal[2] = ab[0] * ac[1] - ab[1] * ac[0];

    /// decide which side is the outside
    double norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    double inner = normal[0] * ad[0] + normal[1] * ad[1] + normal[2] * ad[2];
    
    norm = inner > 0 ? norm : -norm;
    normal[0] /= norm;
    normal[1] /= norm;
    normal[2] /= norm;

    return result;
}

/// 
double tetrahedron_volume(const double3* points)
{
    // Check that we get a tetrahedr
    // Get the coordinates of the four vertices
    const double* x0 = (double*)&(points[0]);
    const double* x1 = (double*)&(points[1]);
    const double* x2 = (double*)&(points[2]);
    const double* x3 = (double*)&(points[3]);

    // Formula for volume from http://mathworld.wolfram.com
    // I see this formula in /dolfin/mesh/TetrahedronCell.cpp which is a part of fenics.
    const double v = (x0[0]*(x1[1]*x2[2] + x3[1]*x1[2] + x2[1]*x3[2]
                           - x2[1]*x1[2] - x1[1]*x3[2] - x3[1]*x2[2])
                    - x1[0]*(x0[1]*x2[2] + x3[1]*x0[2] + x2[1]*x3[2]
                             - x2[1]*x0[2] - x0[1]*x3[2] - x3[1]*x2[2])
                    + x2[0]*(x0[1]*x1[2] + x3[1]*x0[2] + x1[1]*x3[2]
                             - x1[1]*x0[2] - x0[1]*x3[2] - x3[1]*x1[2]) -
                    x3[0]*(x0[1]*x1[2] + x1[1]*x2[2] + x2[1]*x0[2]
                           - x1[1]*x0[2] - x2[1]*x1[2] - x0[1]*x2[2]));

  return std::abs(v)/6.0;
}



int main(){
    double3* points = (double3*)my_tetrahedron_vertices;
    int4 cell = {0,1,2,3};
    int3 boundary = {0,3,5000};
    /// 1. calculate the volumes
    double volume = tetrahedron_volume(points);

    /// 2. calculate the out normal
    double3 normal = triangle_normal(points,boundary.y);

    printf("volume : %f\n", volume);
    printf("%f   %f   %f\n",normal.x, normal.y, normal.z);
    return 0;
}
