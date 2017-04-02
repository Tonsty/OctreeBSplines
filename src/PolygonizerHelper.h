#ifndef POLYGONIZER_HELPER_INCLUDED
#define POLYGONIZER_HELPER_INCLUDED

#include "Function.h"
#include "Geometry.h"

class PolygonizerHelper
{
public:
	static void polygonize(Function *f,float isovalue,float cellsize,float seedx,float seedy,float seedz);
	static void save(const char* filename,const float scale,const Point3D<float> translate);
};

#endif