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
#ifndef ISO_OCTREE_INCLUDED
#define ISO_OCTREE_INCLUDED

#include <hash_map>
#include "MarchingCubes.h"
#include "Octree.h"

class EdgeKey
{
public:
	size_t key1,key2;
	EdgeKey(void){;}
	EdgeKey(const size_t& n1,const size_t &n2);

	EdgeKey& operator = (const EdgeKey& key);
	operator size_t () const;
	operator size_t ();
	bool operator < (const EdgeKey& key) const;

	bool operator != (const EdgeKey& key) const;
};
template<class NodeData,class Real>
class NeighborKey : public OctNode<NodeData,Real>::NeighborKey
{
	void __GetCornerNeighbors(OctNode<NodeData,Real>* node,const int& depth,const int& c,OctNode<NodeData,Real>* neighbors[Cube::CORNERS]);
	OctNode<NodeData,Real>* __FaceNeighbor(OctNode<NodeData,Real>* node,const int& depth,int dir,int off);
	OctNode<NodeData,Real>* __EdgeNeighbor(OctNode<NodeData,Real>* node,const int& depth,int o,int i1,int i2);
	OctNode<NodeData,Real>* __CornerNeighbor(OctNode<NodeData,Real>* node,const int& depth,int x,int y,int z);
public:
	void GetCornerNeighbors(OctNode<NodeData,Real>* node,const int& c,OctNode<NodeData,Real>* neighbors[Cube::CORNERS]);
	OctNode<NodeData,Real>* FaceNeighbor(OctNode<NodeData,Real>* node,int dir,int off);
	OctNode<NodeData,Real>* EdgeNeighbor(OctNode<NodeData,Real>* node,int o,int i1,int i2);
	OctNode<NodeData,Real>* CornerNeighbor(OctNode<NodeData,Real>* node,int x,int y,int z);

	static void CornerIndex(const int& c,int idx[3]);
	static void EdgeIndex(const int& c,int idx[3]);
	static void FaceIndex(const int& c,int idx[3]);
};

template<class NodeData,class Real,class VertexData>
class IsoOctree
{
	class TriangleIndex
	{
	public:
		int idx[3];
	};
	class RootInfo
	{
	public:
		const OctNode<NodeData,Real>* node;
		int edgeIndex;
		long long key;
		typename OctNode<NodeData,Real>::NodeIndex nIdx;
	};
	template<class MeshReal>
	class MeshInfo
	{
	public:
		std::vector<Point3D<MeshReal> > vertexNormals;
		stdext::hash_map<EdgeKey,Point3D<MeshReal> > edgeNormals;
		std::vector<Point3D<MeshReal> > triangleNormals;
		std::vector<TriangleIndex> triangles;
		std::vector<Point3D<MeshReal> > vertices;

		template<class Vertex>
		void set(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,const Real& width,
			Point3D<Real>& translate,Real& scale,const int& noTransform);
	};

	template<class Vertex>
	void getRoots(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const Real& isoValue,stdext::hash_map<long long,int>& roots,std::vector<Vertex>& vertices);
	int getRootIndex(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const int& edgeIndex,RootInfo& ri);
	int getRootPosition(const OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const int& eIndex,const Real& isoValue,Point3D<Real>& position);
	long long getRootKey(const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const int& edgeIndex);

	int getRootPair(const RootInfo& root,const int& maxDepth,RootInfo& pair);
	void getIsoFaceEdges(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const int& faceIndex,std::vector<std::pair<RootInfo,RootInfo>>& edges,const int& flip,const int& useFull);
	void getIsoPolygons(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,stdext::hash_map<long long,int>& roots,std::vector<std::vector<int>>& polygons,const int& useFull);

	template<class C>
	void getEdgeLoops(std::vector<std::pair<C,C> >& edges,stdext::hash_map<C,int>& roots,std::vector<std::vector<int>>& polygons);

	template<class C>
	void getEdgeLoops(std::vector<std::pair<C,C> >& edges,std::vector<std::vector<C>>& loops);

	template<class MeshReal>
	void setDistanceAndNormal(const std::vector<int>& triangles,MeshInfo<MeshReal>& mInfo,const Point3D<Real>& p,	Real& v,Point3D<Real>& n);
	template<class MeshReal>
	void setDistanceAndNormal2(const std::vector<int>& vindices,MeshInfo<MeshReal>& mInfo,const Point3D<Real>& p,	Real& v,Point3D<Real>& n);
	template<class MeshReal>
	void setChildren(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
		const std::vector<int>& triangles,MeshInfo<MeshReal>& mInfo,const int& maxDepth,const int& setCenter,const Real& flatness,stdext::hash_map<long long,std::vector<int>*>* triangleMap=NULL);
	template<class MeshReal>
	void setChildren2(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
		const std::vector<int>& vindices,MeshInfo<MeshReal>& mInfo,const int& maxDepth,const int& setCenter,const Real& flatness,stdext::hash_map<long long,std::vector<int>*>* triangleMap=NULL);

	// Assumes NodeData::mcIndex
	void setMCIndex(const Real& isoValue,const int& useFull);

	NeighborKey<NodeData,Real> nKey;
public:
	// The maximum depth of the tree. This value must be at least as large as the true depth of the tree
	// as its used for assigning unique ids. (It can, however, be larger than the depth for uniqueness to still hold.)
	int maxDepth;
	// The octree itself
	OctNode<NodeData,Real> tree;
	// A hash-table of data associated to the corners of the octree nodes
	stdext::hash_map<long long,VertexData> cornerValues;

	// Sets an octree from a polygon mesh, generating an octree that is fully refined around the surface
	template<class Vertex>
	int set(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,const int& maxDepth,const int& setCenter,const Real& flatness,const int& noTransform);
	template<class Vertex>
	int set(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,const int& maxDepth,const int& setCenter,const Real& flatness,Point3D<Real>& translate,Real& scale,const int& noTransform);
	// Sets an octree from a polygon mesh, generating an octree that is fully refined around the surface
	template<class Vertex>
	int setConforming(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,const int& maxDepth,const int& setCenter,const Real& flatness,const int& noTransform);
	template<class Vertex>
	int setConforming(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,const int& maxDepth,const int& setCenter,const Real& flatness,Point3D<Real>& translate,Real& scale,const int& noTransform);

	// A clean-up method to remove un-needed entries in the cornerValues hash-table
	void resetValues(void);
	// Reads the octree from a file pointer
	int read(FILE* fp,int readData);
	// Writes the octree to a file pointer
	int write(FILE* fp,int writeData) const;

	void interpolateSharedValues(void);

	// Extracts an iso-surface from the octree
	template<class Vertex>
	void getIsoSurface(const Real& isoValue,std::vector<Vertex>& vertices,std::vector<std::vector<int>>& polygons,const int& useFull);
	template<class Vertex>
	void getIsoSoup(const Real& isoValue,std::vector<Vertex>& vertices,std::vector<std::vector<int>>& polygons,const int& useFull);
	template<class Vertex>
	void getDualIsoSurface(const Real& isoValue,std::vector<Vertex>& vertices,std::vector<std::vector<int>>& polygons,const int& useFull);
	// Generates a hash-table indexed by octree node, with each entry containing two pieces of data. The first is
	// mean-curvature vector associated to the intersection of the iso-surface with the node. The second is the area
	// of the intersection of the iso-surface with the node.
	void setNormalFlatness(const Real& isoValue,stdext::hash_map<long long,std::pair<Point3D<Real>,Real>>& curvature);

	// A method for determing if a node has grand-children along an edge
	static int HasEdgeGrandChildren(const OctNode<NodeData,Real>* node,int eIndex);
	// A method for determing if a node has grand-children along a face
	static int HasFaceGrandChildren(const OctNode<NodeData,Real>* node,int fIndex);
};

#include "IsoOctree.inl"

#endif // ISO_OCTREE_INCLUDED
