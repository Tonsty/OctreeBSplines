/*
Copyright (c) 2007, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/
#ifndef _WIN32
#ifndef stdext
#define stdext __gnu_cxx
#endif	
#endif	
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include "IsoOctree.h"
#include "OctreeBspline.h"
#include "CmdLineParser.h"
#include "Ply.h"
#include "Time.h"
#include "VertexData.h"
#include "Geometry.h"
#include "MAT.h"
#include "PolygonizerHelper.h"

typedef PlyVertexWithNormal InPlyVertex;
typedef PlyVertex OutPlyVertex;

template<class VertexData,class Real>
class MyNodeData
{
public:
	int mcIndex;
	Point3D<Real> center;
	VertexData v;
};

template<class NodeData,class VertexData,class Real>
int IsClippable(const IsoOctree<NodeData,Real,VertexData>& isoTree,
				OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIndex,
				stdext::hash_map<long long,std::pair<Point3D<Real>,Real> >& flatness,
				const Real& clipValue,const int& forceConforming)
{
	if(forceConforming)
	{
		for(int i=0;i<Cube::FACES;i++)
			if(IsoOctree<NodeData,Real,VertexData>::HasFaceGrandChildren(node->faceNeighbor(i),Cube::FaceReflectFaceIndex(i,i)))
				return 0;
		for(int i=0;i<Cube::EDGES;i++)
			if(IsoOctree<NodeData,Real,VertexData>::HasEdgeGrandChildren(node->edgeNeighbor(i),Cube::EdgeReflectEdgeIndex(i)))
				return 0;
	}
	long long key=OctNode<NodeData,Real>::CenterIndex(nIndex,isoTree.maxDepth);
	if(flatness.find(key)!=flatness.end())
	{
		float normalSize=Length(flatness[key].first);
		float areaSize=flatness[key].second;
		if(normalSize/areaSize>clipValue)
			return 1;
		else
			return 0;
	}
	else
		return 1;
}
template<class Vertex,class Real>
void PolygonToTriangleMesh(const std::vector<Vertex>& vertices,const std::vector< std::vector<int> >& polygons,
						   std::vector<std::vector<int> >& triangles)
{
	MinimalAreaTriangulation<Real> mat;
	triangles.clear();
	for(size_t i=0;i<polygons.size();i++)
	{
		std::vector<Point3D<Real> > loop;
		std::vector<int> vertexMap;
		std::vector<TriangleIndex> tgl;
		loop.resize(polygons[i].size());
		vertexMap.resize(polygons[i].size());
		for(size_t j=0;j<polygons[i].size();j++)
		{
			loop[j]=vertices[polygons[i][j]];
			vertexMap[j]=polygons[i][j];
		}
		mat.GetTriangulation(loop,tgl);

		size_t tSize=triangles.size();
		triangles.resize(tSize+tgl.size());
		for(size_t j=0;j<tgl.size();j++)
		{
			triangles[tSize+j].resize(3);
			for(int k=0;k<3;k++)
				triangles[tSize+j][k]=vertexMap[tgl[j].idx[k]];
		}
	}
}
template<class Vertex,class Real>
void PolygonToManifoldTriangleMesh( std::vector<Vertex>& vertices , const std::vector< std::vector<int> >& polygons ,
								    std::vector<std::vector<int> >& triangles )
{
	std::vector< int > t;
	t.resize( 3 );
	triangles.clear();
	for( int i=0 ; i<polygons.size() ; i++ )
	{
		if( polygons[i].size()==3 )
			triangles.push_back( polygons[i] );
		else if( polygons[i].size()>3 )
		{
			Point3D< Real > center;
			center *= 0;
			for( int j=0 ; j<polygons[i].size() ; j++ ) center += vertices[ polygons[i][j] ].point;
			center /= (Real) polygons[i].size();

			int idx = (int) vertices.size();
			vertices.push_back( center );
			t[2] = idx;
			for( int j=0 ; j<polygons[i].size() ; j++ )
			{
				t[0] = polygons[i][j];
				t[1] = polygons[i][(j+1)%polygons[i].size()];
				triangles.push_back( t );
			}
		}
	}
}

void ShowUsage(char* ex)
{
	printf("Usage: %s\n",ex);
	printf("\t--in  <input data>\n");
	printf("\t\tInput oriented points (.ply) used to generate the EDT.\n\n");

	printf("\t--out <ouput data>\n");
	printf("\t\tOutput mesh (.ply)\n");

	printf("\t--maxDepth <max depth value of distance tree>\n");
	printf("\t\tThis specifies the maximum depth of the\n");
	printf("\t\tgenerated distance tree.\n\n");

	printf("\t[--flatness <flatness cut-off value>]\n");
	printf("\t\tThis flag forces the octree to be clipped so the octree\n");
	printf("\t\tis not refined around planar regions. (In pracice,\n");
	printf("\t\ta value of about .99 works well.)\n\n");

	printf("\t[--curvature <curvature cut-off value>]\n");
	printf("\t\tThis flag forces the generated distance/bslines octree\n");
	printf("\t\tto be adaptive to the local point curvatures. (In pracice,\n");
	printf("\t\ta value of about .5 works well.)\n\n");

	printf("\t[--fullCaseTable]\n");
	printf("\t\tIf this flag is enabled, the full marching cubes table is\n");
	printf("\t\tused, diambiguating based on the esimated value at the\n");
	printf("\t\tcenter of faces.\n\n");

	printf("\t[--triangleMesh]\n");
	printf("\t\tIf this flag is enabled and the output is a mesh, \n");
	printf("\t\tthe mesh will be triangulated by computing the minimal\n");
	printf("\t\tarea triangulation of the polygons in the mesh\n\n");

	printf("\t[--manifold]\n");
	printf("\t\tIf this flag is enabled and the output is a mesh, the mesh\n");
	printf("\t\twill be triangulated by adding the barycenter to each \n");
	printf("\t\tpolygon with more than three vertices.\n\n");

	printf("\t[--bspline <max depth of b-splines tree>]\n");
	printf("\t\tThis flag forces the distance field to be fitted by\n");
	printf("\t\ta hierarchical b-splines. The default max depth is\n");
	printf("\t\tequal to the max depth of distance tree\n\n");

	printf("\t[--volume <grid resolution>]\n");
	printf("\t\tThis flag tell the program to output a signed distance volume\n");
	printf("\t\t(volume.vti). Generally, we set grid resolution to 128.\n\n");

	printf("\t[--smooth <the smooth weight for fitting>]\n");
	printf("\t\tThe default smooth weight is set to 0.001.\n\n");

	printf("\t[--interpolate <the interpolate weight for fitting>]\n");
	printf("\t\tThe default interpolate weight is set to 0.0.\n\n");

	printf("\t[--minDepthTree]\n");
	printf("\t\tThis flag forces the leaf node of distance tree to at least minDepthTree.\n");
	printf("\t\tNote the B-Splines tree is not forced.\n\n");

	printf("\t[--minDepthMC <minimum depth value for MC>]\n");
	printf("\t\tThis flag forces the MC leaf node to at least minDepthMC.\n\n");

	printf("\t[--splat <splat factor>]\n");
	printf("\t\tThe default splat factor is set 1.0.\n\n");

	printf("\t[--noFit]\n");
	printf("\t\tIf this flag is set, the isosurface is directly extracted\n");
	printf("\t\tfrom the adaptive signed distance field without fitting.\n\n");

	printf("\t[--octree]\n");
	printf("\t\tThis flag tell the program to output\n");
	printf("\t\tthe octree grid (octree.vtk).\n\n");

	printf("\t[--normal]\n");
	printf("\t\tThis flag tell the program to output mesh with normal\n\n");

	printf("\t[--sphereTest]\n");
	printf("\t\tThis flag tell the program to perform sphereTest\n\n");

	printf("\t[--isoValue <isovalue>]\n");
	printf("\t\tThis flag tell the program the isoValue\n\n");	

	printf("\t[--bloomenthal <bloomenthal>]\n");
	printf("\t\tThis flag tell the program the extract Bloomenthal's iso-surface\n\n");	
}

int main(int argc,char* argv[])
{
	MarchingCubes::SetCaseTable();
	MarchingCubes::SetFullCaseTable();

	typedef OctNode<MyNodeData<VertexValue<float>,float>,float> MyOctNode;

	cmdLineString In, Out;
	cmdLineReadable Conforming,FullCaseTable,TriangleMesh,Dual,Manifold,NoFit,Normal,SphereTest;
	cmdLineFloat Flatness(-1),Curvature(-1),Smooth(0.001),Interpolate(0.0),Splat(1.0),IsoValue(0.0);
	cmdLineInt MaxDepth,Bspline(-1),Volume(128),MinDepthMC(-1),Bloomenthal(128),MinDepthTree(-1);
	char* paramNames[]=
	{
		"in","out","flatness","curvature","conforming","fullCaseTable","maxDepth","triangleMesh","dual","manifold",
		"bspline","smooth","interpolate","minDepthTree","minDepthMC","volume","splat","noFit","normal","sphereTest",
		"isoValue","bloomenthal"
	};
	cmdLineReadable* params[]= 
	{
		&In,&Out,&Flatness,&Curvature,&Conforming,&FullCaseTable,&MaxDepth,&TriangleMesh,&Dual,&Manifold,
		&Bspline,&Smooth,&Interpolate,&MinDepthTree,&MinDepthMC,&Volume,&Splat,&NoFit,&Normal,&SphereTest,
		&IsoValue,&Bloomenthal
	};
	int paramNum=sizeof(paramNames)/sizeof(char*);
	cmdLineParse(argc-1,&argv[1],paramNames,paramNum,params,0);

	if((!SphereTest.set && !In.set) || !Out.set || !MaxDepth.set)
	{
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}

	if(Bspline.set && (Bspline.value<=0 || Bspline.value>MaxDepth.value)) Bspline.value=MaxDepth.value;

	Point3D<float> translate;
	float scale=1.f;
	translate[0]=translate[1]=translate[2]=0;

	typedef OctreeBspline<MyNodeData<VertexValue<float>,float>,float,VertexValue<float>> OctBspline;
	OctBspline octreeBspline;
	octreeBspline.maxBsplineDepth=Bspline.value;

	int ft;
	std::vector<InPlyVertex> inVertices;
	std::vector<std::vector<int> > polygons;
	double t;

	if(SphereTest.set)
	{
		ft=PLY_ASCII;
		octreeBspline.set4(MaxDepth.value,translate,scale);
	}
	else
	{
		t=Time();
		printf("Loading data ...\n");
		PlyReadPolygons(In.value,inVertices,polygons,ft);
		printf("Got data in: %f\n", Time()-t);

		printf("Establishing signed distance field ...\n");
		printf("maxDepth: %d\n", MaxDepth.value);
		printf("maxBsplineDepth: %d\n", Bspline.value);
		t=Time();
		octreeBspline.set3(inVertices,polygons,MaxDepth.value,Dual.set,Flatness.value,Curvature.value,Splat.value,MinDepthTree.value,translate,scale,0,!NoFit.set);
		printf("Got signed distance field in: %f\n", Time()-t);
		printf("Nodes In: %d / %d\n",octreeBspline.IsoOctree::tree.nodes(),octreeBspline.IsoOctree::tree.leaves());
		printf("Values In: %d\n",octreeBspline.IsoOctree::cornerValues.size());
		printf("Scale : %f\n",scale);
		printf("Translate : %f %f %f\n",translate[0],translate[1],translate[2]);

		std::vector<InPlyVertex> emptyVertices;
		std::vector<std::vector<int> > emptyPolygons;

		inVertices.swap(emptyVertices);
		polygons.swap(emptyPolygons);
	}

	//octreeBspline.exportOctreeGrid(scale, translate);

	if(!NoFit.set && Bspline.set && Bspline.value>0) 
	{
		t=Time();
		printf("Fitting data ...\n");
		//octreeBspline.directBsplineFitting(Smooth.value,Interpolate.value);
		octreeBspline.multigridBsplineFitting(Smooth.value,Interpolate.value);
		printf("Got fitted in: %f\n", Time()-t);
	}

	std::vector<OutPlyVertex> outVertices;
	printf("Estracting iso-surface ...\n");
	t=Time();
	if(!NoFit.set && Bspline.set && Bspline.value>0) octreeBspline.updateCornerValues();
	if(!NoFit.set && Bspline.set && Bspline.value>0 && MinDepthMC.set) octreeBspline.setMinDepthMCLeafNode(0,MinDepthMC.value,FullCaseTable.set);
	octreeBspline.getIsoSurface(IsoValue.value,outVertices,polygons,FullCaseTable.set);
	printf("Got iso-surface in: %f\n",Time()-t);

	for(size_t i=0;i<outVertices.size();i++)
		outVertices[i].point=outVertices[i].point/scale-translate;

	std::vector<std::vector<int> > triangles;
	if(Manifold.set)
	{
		double t=Time();
		PolygonToManifoldTriangleMesh<OutPlyVertex,float>(outVertices,polygons,triangles);
		printf("Converted polygons to triangles in: %f\n",Time()-t);
		PlyWritePolygons(Out.value,outVertices,triangles,ft);
		printf("Vertices: %d\n",outVertices.size());
		printf("Triangles: %d\n",triangles.size());
	}
	else if(TriangleMesh.set)
	{
		double t=Time();
		PolygonToTriangleMesh<OutPlyVertex,float>(outVertices,polygons,triangles);
		printf("Converted polygons to triangles in: %f\n",Time()-t);
		PlyWritePolygons(Out.value,outVertices,triangles,ft);
		printf("Vertices: %d\n",outVertices.size());
		printf("Triangles: %d\n",triangles.size());
	}
	else
	{
		PlyWritePolygons(Out.value,outVertices,polygons,ft);
		printf("Vertices: %d\n",outVertices.size());
		printf("Polygons: %d\n",polygons.size());
	}

	if(Bspline.set && Bspline.value>0)
	{
		if(Volume.set && Volume.value>0) octreeBspline.exportVolumeData(scale,translate,Volume.value);
		if(SphereTest.set)
		{
			printf("Extracting SphereTest iso-surface ...\n");
			t=Time();
			PolygonizerHelper::polygonize((Function*)(&octreeBspline),0.0f,1.0f/64,0.5f,0.5f,0.5f);
			printf("Got iso-surface in: %f\n", Time()-t);
			PolygonizerHelper::saveSphereTest("sphereTest.ply",scale,translate);
		}
	}

	if(Bloomenthal.set)
	{
		printf("Extracting Bloomenthal's iso-surface ...\n");
		t=Time();
		PolygonizerHelper::polygonize((Function*)(&octreeBspline),IsoValue.value,1.0f/Bloomenthal.value,0.5f,0.5f,0.5f);
		printf("Got Bloomenthal iso-surface in: %f\n", Time()-t);
		PolygonizerHelper::save("bloomenthal.ply",scale,translate);		
	}

	if(Normal.set)
	{
		std::vector<PlyVertexWithNormal> outVertexWithNormals(outVertices.size());

		for(size_t i=0;i<outVertexWithNormals.size();i++)
			outVertexWithNormals[i].point=(outVertices[i].point+translate)*scale;

		for(size_t i=0;i<outVertexWithNormals.size();i++)
		{
			octreeBspline.gradient(outVertexWithNormals[i].point.coords,outVertexWithNormals[i].normal.coords);
			outVertexWithNormals[i].normal/=Length(outVertexWithNormals[i].normal);
		}

		for(size_t i=0;i<outVertexWithNormals.size();i++)
			outVertexWithNormals[i].point=outVertexWithNormals[i].point/scale-translate;

		if(Manifold.set || TriangleMesh.set)
		{
			PlyWritePolygons("result_with_normal.ply",outVertexWithNormals,triangles,ft);
		}
		else
		{
			PlyWritePolygons("result_with_normal.ply",outVertexWithNormals,polygons,ft);
		}
	}

	return EXIT_SUCCESS;
}
/*

	if(UseNormals.set)	return Process<VertexValueAndNormal<float>,float>(argc,argv);
	else				return Process<VertexValue<float>,float>(argc,argv);
}
*/