#include "PolygonizerHelper.h"
#include "polygonizer.h"
#include <stdio.h>

extern "C" int triangle2(int i1,int i2,int i3,VERTICES vertices);
extern "C" char* polygonize(double (*function)(double,double,double),double size,int bounds,double x,double y,double z,int (*triproc)(int,int,int,VERTICES),int mode);
extern "C" int gntris;
extern "C" VERTICES gvertices;
extern "C" TRIANGLES gtriangles;

Function *gf;
float gisovalue=0.0;

double function(double x,double y,double z) 
{
	float pos[3];
	pos[0]=(float)x;
	pos[1]=(float)y;
	pos[2]=(float)z;
	return (double)(gf->eval(pos)-gisovalue);
}

void PolygonizerHelper::polygonize(Function *f,float isovalue,float cellsize,float seedx,float seedy,float seedz) 
{
	::gf=f;
	gisovalue=isovalue;
	char *err=::polygonize(::function,cellsize,(int)(1.0/cellsize),seedx,seedy,seedz,triangle2,0);
	if(err) 
	{
		printf("%s\n",err);
		exit(1);
	}
}

void PolygonizerHelper::save(const char* filename,const float scale,const Point3D<float> translate)
{
	FILE*file=fopen(filename,"w");
	fprintf(file,"ply\n" \
		"format ascii 1.0\n" \
		"comment polygonizer generated\n" \
		"element vertex %d\n" \
		"property float x\n" \
		"property float y\n" \
		"property float z\n" \
		"property float nx\n" \
		"property float ny\n" \
		"property float nz\n" \
		"element face %d\n" \
		"property list uchar int vertex_indices\n" \
		"end_header\n",::gvertices.count,::gntris);
	for (int i=0;i<::gvertices.count;i++) {
		VERTEX v;
		v=::gvertices.ptr[i];
		Point3D<float> point={(float)v.position.x,(float)v.position.y,(float)v.position.z};
		point=point/scale-translate;
		Point3D<float> normal={(float) v.normal.x,(float) v.normal.y,(float) v.normal.z};
		fprintf(file,"%f %f %f %f %f %f\n",point[0], point[1],point[2],normal[0],normal[1],normal[2]);
	}
	for (int i=0;i<::gtriangles.count;i++) {
		TRIANGLE t;
		t=::gtriangles.ptr[i];
		fprintf(file,"3 %d %d %d\n",t.i1,t.i3,t.i2);
	}
}