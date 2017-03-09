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
#include "Ply2.h"
#include "PolygonizerHelper.h"

typedef PlyVertexWithNormal MyPlyVertex;

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
				stdext::hash_map<long long,std::pair<Point3D<Real>,Real>>& flatness,
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
	printf("\t\tInput mesh (.ply) used to generate the EDT.\n");

	printf("\t--out <ouput data>\n");
	printf("\t\tOutput mesh (.ply)\n");

	printf("\t--maxDepth\n");
	printf("\t\tIf the octree is generated from a mesh, this specifies the\n");
	printf("\t\tmaximum depth of the generated tree.\n");

	printf("\t[--curvature <cut-off value>]\n");
	printf("\t\tThis flag forces the octree to be clipped so the octree\n");
	printf("\t\tis only refined around high-curvature regions. (In pracice,\n");
	printf("\t\ta value of about .99 works well.)\n");

	printf("\t[--conforming]\n");
	printf("\t\tIf this flag is enabled, the octree satisfies the condition\n");
	printf("\t\tthe depth disparity between adjacent leaf nodes is never\n");
	printf("\t\tgreater than one.\n");

	printf("\t[--dual]\n");
	printf("\t\tIf this flag is enabled, the iso-surface is extracted using\n");
	printf("\t\tthe method of [Schaeffer et al. 04] by sampling the EDT at the\n");
	printf("\t\tcenters of the leaf nodes.\n");

	printf("\t[--fullCaseTable]\n");
	printf("\t\tIf this flag is enabled, the full marching cubes table is\n");
	printf("\t\tused, diambiguating based on the esimated value at the\n");
	printf("\t\tcenter of faces.\n");

	printf("\t[--triangleMesh]\n");
	printf("\t\tIf this flag is enabled and the output is a mesh, the mesh\n");
	printf("\t\twill be triangulated by computing the minimal area triangulation of\n");
	printf("\t\teach of the polygons in the mesh\n");

	printf( "\t[--manifold]\n");
	printf( "\t\tIf this flag is enabled and the output is a mesh, the mesh\n");
	printf( "\t\twill be triangulated by adding the barycenter to each polygon\n");
	printf( "\t\t\to each polygon with more than three vertices.\n");

	printf( "\t[--bspline <max depth value of the hierarchical B-Splines>]\n");
	printf( "\t\tThis flag fits the adaptive distance field by\n");
	printf( "\t\ta hierarchical implicit B-Splines function.\n");
}

int main(int argc,char* argv[])
{
	MarchingCubes::SetCaseTable();
	MarchingCubes::SetFullCaseTable();

	typedef OctNode<MyNodeData<VertexValue<float>,float>,float> MyOctNode;

	cmdLineString In, Out;
	cmdLineReadable Conforming,FullCaseTable,TriangleMesh,Dual,Manifold;
	cmdLineFloat Curvature(-1);
	cmdLineInt MaxDepth,Bspline(-1);
	char* paramNames[]=
	{
		"in","out","curvature","conforming","fullCaseTable","maxDepth","triangleMesh","dual","manifold","bspline" 
	};
	cmdLineReadable* params[]= 
	{
		&In,&Out,&Curvature,&Conforming,&FullCaseTable,&MaxDepth,&TriangleMesh,&Dual,&Manifold,&Bspline
	};
	int paramNum=sizeof(paramNames)/sizeof(char*);
	cmdLineParse(argc-1,&argv[1],paramNames,paramNum,params,0);

	if(!In.set || !Out.set || !MaxDepth.set)
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
	//VertexValue<float>::f=(Function*)(&octreeBspline);

	int ft;
	std::vector<MyPlyVertex> vertices;
	std::vector<std::vector<int> > polygons;

	double t=Time();

	printf("Loading data ...\n");
	PlyReadPolygons(In.value,vertices,polygons,ft);
	printf("Got data in: %f\n", Time()-t);

	printf("Establishing signed distance field ...\n");
	printf("maxDepth: %d\n", MaxDepth.value);
	printf("maxBsplineDepth: %d\n", Bspline.value);
	t=Time();
	octreeBspline.set2(vertices,polygons,MaxDepth.value,Dual.set,Curvature.value,translate,scale,0);
	printf("Got signed distance field in: %f\n", Time()-t);
	printf("Nodes In: %d / %d\n",octreeBspline.tree.nodes(),octreeBspline.tree.leaves());
	printf("Values In: %d\n",octreeBspline.cornerValues.size());

	if(Bspline.set && Bspline.value>0) 
	{
		t=Time();
		printf("Fitting data ...\n");
		//octreeBspline.directBsplineFitting();
		octreeBspline.multigridBsplineFitting();
		printf("Got fitted in: %f\n", Time()-t);
		//octreeBspline.exportVTKData(scale,translate);
	}

	vertices.clear();
	polygons.clear();

	printf("Estracting iso-surface ...\n");
	t=Time();
	//octreeBspline.setMCLeafNodeToMaxDepth(0,FullCaseTable.set);
	octreeBspline.getIsoSurface(0,vertices,polygons,FullCaseTable.set);
	printf("Got iso-surface in: %f\n",Time()-t);

	for(size_t i=0;i<vertices.size();i++)
		vertices[i].point=vertices[i].point/scale-translate;

	if(Manifold.set)
	{
		std::vector<std::vector<int> > triangles;
		double t=Time();
		PolygonToManifoldTriangleMesh<MyPlyVertex,float>(vertices,polygons,triangles);
		printf("Converted polygons to triangles in: %f\n",Time()-t);
		PlyWritePolygons(Out.value,vertices,triangles,ft);
		printf("Vertices: %d\n",vertices.size());
		printf("Triangles: %d\n",triangles.size());
	}
	else if(TriangleMesh.set)
	{
		std::vector<std::vector<int> > triangles;
		double t=Time();
		PolygonToTriangleMesh<MyPlyVertex,float>(vertices,polygons,triangles);
		printf("Converted polygons to triangles in: %f\n",Time()-t);
		PlyWritePolygons(Out.value,vertices,triangles,ft);
		printf("Vertices: %d\n",vertices.size());
		printf("Triangles: %d\n",triangles.size());
	}
	else
	{
		PlyWritePolygons(Out.value,vertices,polygons,ft);
		printf("Vertices: %d\n",vertices.size());
		printf("Polygons: %d\n",polygons.size());
	}

	printf("Extracting iso-surface ...\n");
	t=Time();
	PolygonizerHelper::polygonize((Function*)(&octreeBspline),0.02f,1.0f/200,0.5f,0.5f,0.5f);
	printf("Got iso-surface in: %f\n", Time()-t);
	PolygonizerHelper::save("mesh2.ply",scale,translate);

	return EXIT_SUCCESS;
}
/*

	if(UseNormals.set)	return Process<VertexValueAndNormal<float>,float>(argc,argv);
	else				return Process<VertexValue<float>,float>(argc,argv);
}
*/