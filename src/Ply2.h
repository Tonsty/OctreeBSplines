#ifndef PLY2_INCLUDED
#define PLY2_INCLUDED

class PlyVertexWithNormal
{
public:
	const static int Components=6;
	static PlyProperty Properties[];

	Point3D<float> point;
	Point3D<float> normal;

	operator Point3D<float>& ()					{return point;}
	operator const Point3D<float>& () const		{return point;}
	PlyVertexWithNormal(void)					{point[0]=point[1]=point[2]=0; normal[0]=normal[1]=normal[2]=0;}
	PlyVertexWithNormal(const Point3D<float>& p)	{point=p; normal[0]=normal[1]=normal[2]=0;}
};
PlyProperty PlyVertexWithNormal::Properties[]=
{
	{"x", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertexWithNormal,point.coords[0])), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertexWithNormal,point.coords[1])), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertexWithNormal,point.coords[2])), 0, 0, 0, 0},
	{"nx", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertexWithNormal,normal.coords[0])), 0, 0, 0, 0},
	{"ny", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertexWithNormal,normal.coords[1])), 0, 0, 0, 0},
	{"nz", PLY_FLOAT, PLY_FLOAT, int(offsetof(PlyVertexWithNormal,normal.coords[2])), 0, 0, 0, 0}
};

#endif