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
#define USE_MAX_DEPTH_SPEED_UP 1
#define NEW_MC_INDEX 1
#define NEW_ISO_POLYGON 1

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

/////////////
// EdgeKey //
/////////////
EdgeKey::EdgeKey(const size_t& n1,const size_t& n2)
{
	if(n1<n2)
	{
		key1=n1;
		key2=n2;
	}
	else
	{
		key1=n2;
		key2=n1;
	}
}
EdgeKey::operator size_t() const
{
	return (size_t)key1;
}
EdgeKey::operator size_t()
{
	return (size_t)key1;
}
EdgeKey& EdgeKey::operator = (const EdgeKey& key)
{
	key1=key.key1;
	key2=key.key2;
	return *this;
}
bool EdgeKey::operator < (const EdgeKey& key) const
{
	if(key1<key.key1)	return true;
	else if(key1==key.key1)
		if(key2<key.key2)
			return true;
		else
			return false;
	else
		return false;
}
bool EdgeKey::operator != (const EdgeKey& key) const
{
	if(key1!=key.key1 || key2!=key.key2)
		return true;
	else
		return false;
}
/////////////////
// NeighborKey //
/////////////////
template<class NodeData,class Real>
OctNode<NodeData,Real>* NeighborKey<NodeData,Real>::__FaceNeighbor(OctNode<NodeData,Real>* node,const int& depth,int dir,int off)
{
	if(!depth)	return NULL;
	int x,y,z;
	x=y=z=1;
	switch(dir)
	{
	case 0:
		x=off<<1;
		break;
	case 1:
		y=off<<1;
		break;
	case 2:
		z=off<<1;
		break;
	}
	for(int d=depth;d>=0;d--)
		if(neighbors[d].neighbors[x][y][z])
			return neighbors[d].neighbors[x][y][z];

	return NULL;
}
template<class NodeData,class Real>
OctNode<NodeData,Real>* NeighborKey<NodeData,Real>::__EdgeNeighbor(OctNode<NodeData,Real>* node,const int& depth,int o,int i1,int i2)
{
	if(!depth)	return NULL;
	int x,y,z,cIndex,xx,yy,zz;
	x=y=z=1;

	// Check if the edge-adjacent neighbor exists at the current depth
	switch(o)
	{
	case 0:
		y=i1<<1;
		z=i2<<1;
		break;
	case 1:
		x=i1<<1;
		z=i2<<1;
		break;
	case 2:
		x=i1<<1;
		y=i2<<1;
		break;
	}
	if(neighbors[depth].neighbors[x][y][z])
		return neighbors[depth].neighbors[x][y][z];

	cIndex=int(node-node->parent->children);
	Cube::FactorCornerIndex(cIndex,xx,yy,zz);

	// Check if the node is on the corresponding edge of the parent
	switch(o)
	{
	case 0:
		if(yy==i1 && zz==i2)	return __EdgeNeighbor(node->parent,depth-1,o,i1,i2);
		break;
	case 1:
		if(xx==i1 && zz==i2)	return __EdgeNeighbor(node->parent,depth-1,o,i1,i2);
		break;
	case 2:
		if(xx==i1 && yy==i2)	return __EdgeNeighbor(node->parent,depth-1,o,i1,i2);
		break;
	}

	// Check if the node is on a face of the parent containing the edge
	switch(o)
	{
	case 0:
		if(yy==i1)	return __FaceNeighbor(node->parent,depth-1,1,yy);
		if(zz==i2)	return __FaceNeighbor(node->parent,depth-1,2,zz);
		break;
	case 1:
		if(xx==i1)	return __FaceNeighbor(node->parent,depth-1,0,xx);
		if(zz==i2)	return __FaceNeighbor(node->parent,depth-1,2,zz);
		break;
	case 2:
		if(xx==i1)	return __FaceNeighbor(node->parent,depth-1,0,xx);
		if(yy==i2)	return __FaceNeighbor(node->parent,depth-1,1,yy);
		break;
	}
	fprintf(stderr,"We shouldn't be here: Edges\n");
	return NULL;
}
template<class NodeData,class Real>
OctNode<NodeData,Real>* NeighborKey<NodeData,Real>::__CornerNeighbor(OctNode<NodeData,Real>* node,const int& depth,int x,int y,int z)
{
	if(!depth)	return NULL;
	int cIndex,xx,yy,zz;

	// Check if the edge-adjacent neighbor exists at the current depth
	if(neighbors[depth].neighbors[x<<1][y<<1][z<<1])
		return neighbors[depth].neighbors[x<<1][y<<1][z<<1];

	cIndex=int(node-node->parent->children);
	Cube::FactorCornerIndex(cIndex,xx,yy,zz);

	// Check if the node is on the corresponding corner of the parent
	if(xx==x && yy==y && zz==z)	return __CornerNeighbor(node->parent,depth-1,x,y,z);

	// Check if the node is on an edge of the parent containing the corner
	if(xx==x && yy==y)			return __EdgeNeighbor(node->parent,depth-1,2,x,y);
	if(xx==x && zz==z)			return __EdgeNeighbor(node->parent,depth-1,1,x,z);
	if(yy==y && zz==z)			return __EdgeNeighbor(node->parent,depth-1,0,y,z);

	// Check if the node is on a face of the parent containing the edge
	if(xx==x)					return __FaceNeighbor(node->parent,depth-1,0,x);
	if(yy==y)					return __FaceNeighbor(node->parent,depth-1,1,y);
	if(zz==z)					return __FaceNeighbor(node->parent,depth-1,2,z);

	fprintf(stderr,"We shouldn't be here: Corners\n");
	return NULL;
}
template<class NodeData,class Real>
void NeighborKey<NodeData,Real>::__GetCornerNeighbors(OctNode<NodeData,Real>* node,const int& d,const int& c,OctNode<NodeData,Real>* neighbors[Cube::CORNERS])
{
	int x,y,z,xx,yy,zz,ax,ay,az;
	Cube::FactorCornerIndex(c,x,y,z);
	ax=x^1;
	ay=y^1;
	az=z^1;
	xx=x<<0;
	yy=y<<1;
	zz=z<<2;
	ax<<=0;
	ay<<=1;
	az<<=2;

	// Set the current node
	neighbors[ax|ay|az]=node;

	// Set the face adjacent neighbors
	neighbors[xx|ay|az]=__FaceNeighbor(node,d,0,x);
	neighbors[ax|yy|az]=__FaceNeighbor(node,d,1,y);
	neighbors[ax|ay|zz]=__FaceNeighbor(node,d,2,z);

	// Set the edge adjacent neighbors
	neighbors[ax|yy|zz]=__EdgeNeighbor(node,d,0,y,z);
	neighbors[xx|ay|zz]=__EdgeNeighbor(node,d,1,x,z);
	neighbors[xx|yy|az]=__EdgeNeighbor(node,d,2,x,y);

	// Set the corner adjacent neighbor
	neighbors[xx|yy|zz]=__CornerNeighbor(node,d,x,y,z);
}
template<class NodeData,class Real>
void NeighborKey<NodeData,Real>::GetCornerNeighbors(OctNode<NodeData,Real>* node,const int& c,OctNode<NodeData,Real>* neighbors[Cube::CORNERS])
{
	getNeighbors(node);
	memset(neighbors,NULL,sizeof(OctNode<NodeData,Real>*)*Cube::CORNERS);
	__GetCornerNeighbors(node,depth,c,neighbors);
}
template<class NodeData,class Real>
void NeighborKey<NodeData,Real>::CornerIndex(const int& c,int idx[3])
{
	int x,y,z;
	Cube::FactorCornerIndex(c,x,y,z);
	idx[0]=x<<1;
	idx[1]=y<<1;
	idx[2]=z<<1;
}
template<class NodeData,class Real>
void NeighborKey<NodeData,Real>::EdgeIndex(const int& e,int idx[3])
{
	int o,i1,i2;
	Cube::FactorEdgeIndex(e,o,i1,i2);
	idx[0]=idx[1]=idx[2]=1;
	switch(o)
	{
	case 0:
		idx[1]=i1<<1;
		idx[2]=i2<<1;
		break;
	case 1:
		idx[0]=i1<<1;
		idx[2]=i2<<1;
		break;
	case 2:
		idx[0]=i1<<1;
		idx[1]=i2<<1;
		break;
	}
}
template<class NodeData,class Real>
void NeighborKey<NodeData,Real>::FaceIndex(const int& f,int idx[3])
{
	int dir,off;
	Cube::FactorFaceIndex(f,dir,off);
	idx[0]=idx[1]=idx[2]=1;
	idx[dir]=off<<1;
}
///////////////
// IsoOctree //
///////////////
template<class NodeData,class Real,class VertexData>
int IsoOctree<NodeData,Real,VertexData>::HasEdgeGrandChildren(const OctNode<NodeData,Real>* node,int eIndex)
{
	if(!node || !node->children)
		return 0;
	int c1,c2;
	Cube::EdgeCorners(eIndex,c1,c2);
	if(node->children[c1].children || node->children[c2].children)
		return 1;
	else
		return 0;
}
template<class NodeData,class Real,class VertexData>
int IsoOctree<NodeData,Real,VertexData>::HasFaceGrandChildren(const OctNode<NodeData,Real>* node,int fIndex)
{
	if(!node || !node->children)
		return 0;
	int c1,c2,c3,c4;
	Cube::FaceCorners(fIndex,c1,c2,c3,c4);
	if(node->children[c1].children || node->children[c2].children || node->children[c3].children || node->children[c4].children)
		return 1;
	else
		return 0;
}

template<class NodeData,class Real,class VertexData>
template<class MeshReal>
template<class Vertex>
void IsoOctree<NodeData,Real,VertexData>::MeshInfo<MeshReal>::set(const std::vector<Vertex>& verts,const std::vector<std::vector<int> >& polys,const Real& width,
																  Point3D<Real>& translate,Real& scale,const int& noTransform)
{
	if(!noTransform)
	{
		Point3D<Real> min,max;
		for(size_t i=0;i<verts.size();i++)
		{
			for(int j=0;j<3;j++)
			{
				if(!i || Point3D<Real>(verts[i])[j]<min[j])	min[j]=Point3D<Real>(verts[i])[j];
				if(!i || Point3D<Real>(verts[i])[j]>max[j])	max[j]=Point3D<Real>(verts[i])[j];
			}
		}

		scale=max[0]-min[0];
		if( (max[1]-min[1])>scale )		scale=max[1]-min[1];
		if( (max[2]-min[2])>scale )		scale=max[2]-min[2];

		scale*=width;
		scale=Real(1.0)/scale;
		Point3D<Real> ctr;
		ctr[0]=ctr[1]=ctr[2]=Real(0.5);

		translate=ctr/scale-(max+min)/2;
	}
	else
	{
		translate[0]=translate[1]=translate[2]=0;
		scale=1;
	}
	vertices.resize(verts.size());
	for(size_t i=0;i<verts.size();i++)
		for(int j=0;j<3;j++)
			vertices[i][j]=MeshReal((translate[j]+(Point3D<Real>(verts[i]))[j])*scale);

	size_t pSize=0;
	size_t vSize=verts.size();
	for(size_t i=0;i<polys.size();i++)
	{
		if(polys.size()<3) continue;
		if(polys[i].size()==3)
		{
			Point3D<double> N=Normal(vertices[polys[i][0]],vertices[polys[i][1]],vertices[polys[i][2]]);
			if(Length(N)==0) continue;
			triangles.resize(pSize+1);
			for(int j=0;j<3;j++) triangles[pSize].idx[j]=polys[i][j];
			pSize++;
		}
		else
		{
			triangles.resize(pSize+polys[i].size());
			Point3D<MeshReal> ctr;
			ctr[0]=ctr[1]=ctr[2]=0;
			for(size_t j=0;j<polys[i].size();j++) ctr+=vertices[polys[i][j]];
			ctr/=Real(polys[i].size());
			vertices.push_back(ctr);
			for(size_t j=0;j<polys[i].size();j++)
			{
				triangles[pSize+j].idx[0]=polys[i][j];
				triangles[pSize+j].idx[1]=polys[i][(j+1)%polys[i].size()];
				triangles[pSize+j].idx[2]=int(vSize);
			}
			vSize++;
			pSize+=polys[i].size();
		}
	}
	vertexNormals.resize(vertices.size());
	triangleNormals.resize(triangles.size());

	for(size_t i=0;i<vertices.size();i++)
		vertexNormals[i][0]=vertexNormals[i][1]=vertexNormals[i][2]=0;

	for(size_t i=0;i<triangles.size();i++)
	{
		triangleNormals[i]=Normal(vertices[triangles[i].idx[0]],vertices[triangles[i].idx[1]],vertices[triangles[i].idx[2]]);
		triangleNormals[i]/=Length(triangleNormals[i]);
		for(int j=0;j<3;j++)
			vertexNormals[triangles[i].idx[j]]+=triangleNormals[i];
		for(int j=0;j<3;j++)
		{
			EdgeKey eKey(triangles[i].idx[j],triangles[i].idx[(j+1)%3]);
			if(edgeNormals.find(eKey)==edgeNormals.end())
				edgeNormals[eKey]=triangleNormals[i];
			else
				edgeNormals[eKey]+=triangleNormals[i];
		}
	}
	for(size_t i=0;i<vertices.size();i++)
		vertexNormals[i]/=Length(vertexNormals[i]);

	for(typename stdext::hash_map<EdgeKey,Point3D<MeshReal>,HashEdgeKey >::iterator iter=edgeNormals.begin();iter!=edgeNormals.end();iter++)
		iter->second/=Length(iter->second);

}
template<class NodeData,class Real,class VertexData>
template<class MeshReal>
template<class Vertex>
void IsoOctree<NodeData,Real,VertexData>::MeshInfo<MeshReal>::set2(const std::vector<Vertex>& verts,const std::vector<std::vector<int> >& polys,const Real& width,
	Point3D<Real>& translate,Real& scale,const int& noTransform)
{
	if(!noTransform)
	{
		Point3D<Real> min,max;
		for(size_t i=0;i<verts.size();i++)
		{
			for(int j=0;j<3;j++)
			{
				if(!i || Point3D<Real>(verts[i])[j]<min[j])	min[j]=Point3D<Real>(verts[i])[j];
				if(!i || Point3D<Real>(verts[i])[j]>max[j])	max[j]=Point3D<Real>(verts[i])[j];
			}
		}

		scale=max[0]-min[0];
		if( (max[1]-min[1])>scale )		scale=max[1]-min[1];
		if( (max[2]-min[2])>scale )		scale=max[2]-min[2];

		scale*=width;
		scale=Real(1.0)/scale;
		Point3D<Real> ctr;
		ctr[0]=ctr[1]=ctr[2]=Real(0.5);

		translate=ctr/scale-(max+min)/2;
	}
	else
	{
		translate[0]=translate[1]=translate[2]=0;
		scale=1;
	}
	vertices.resize(verts.size());
	for(size_t i=0;i<verts.size();i++)
		for(int j=0;j<3;j++)
			vertices[i][j]=MeshReal((translate[j]+(Point3D<Real>(verts[i]))[j])*scale);

	//vertices have precomputed normals, use point to point (plane) distance metric to compute signed distance
	vertexNormals.resize(vertices.size());
	vertexCurvatures.resize(vertices.size());
	for(size_t i=0;i<verts.size();i++) {
		for(int j=0;j<3;j++) vertexNormals[i][j]=verts[i].normal[j]; //MeshReal(*(Real*)((char*)&verts[i]+Vertex::Properties[3+j].offset));
		vertexNormals[i]/=Length(vertexNormals[i]);
		vertexCurvatures[i]=verts[i].curvature;//MeshReal(*(Real*)((char*)&verts[i]+Vertex::Properties[6].offset));
		vertexCurvatures[i]=(MeshReal)(vertexCurvatures[i]/scale);
	}
}

template<class NodeData,class Real,class VertexData>
template<class Vertex>
int IsoOctree<NodeData,Real,VertexData>::set(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,const int& maxDepth,const int& setCenter,const Real& flatness,const int& noTransform)
{
	Point3D<Real> t;
	Real s;
	return set(vertices,polygons,maxDepth,setCenter,flatness,t,s,noTransform);
}
template<class NodeData,class Real,class VertexData>
template<class Vertex>
int IsoOctree<NodeData,Real,VertexData>::set(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,
	const int& maxDepth,const int& setCenter,const Real& flatness,
	Point3D<Real>& translate,Real& scale,const int& noTransform)
{
	this->maxDepth=maxDepth;
	typename OctNode<NodeData,Real>::NodeIndex nIdx;

	MeshInfo<double> mInfo;
	std::vector<int> myTriangles;
	mInfo.set(vertices,polygons,Real(1.1),translate,scale,noTransform);
	myTriangles.resize(mInfo.triangles.size());
	for(int i=0;i<int(mInfo.triangles.size());i++) myTriangles[i]=i;

	cornerValues.clear();
	Real dist;
	Point3D<Real> n,p;
	for(int c=0;c<Cube::CORNERS;c++)
	{
		int x,y,z;
		Cube::FactorCornerIndex(c,x,y,z);
		p[0]=Real(x);
		p[1]=Real(y);
		p[2]=Real(z);

		setDistanceAndNormal(myTriangles,mInfo,p,dist,n);
		cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,c,maxDepth)]=VertexData(dist,n);
	}
	if(setCenter)
	{
		Real w;
		OctNode<NodeData,Real>::CenterAndWidth(nIdx,tree.nodeData.center,w);
		setDistanceAndNormal(myTriangles,mInfo,tree.nodeData.center,dist,n);
		tree.nodeData.v=VertexData(dist,n);
	}
	setChildren(&tree,nIdx,myTriangles,mInfo,maxDepth,setCenter,flatness);
	return 1;
}
template<class NodeData,class Real,class VertexData>
template<class Vertex>
int IsoOctree<NodeData,Real,VertexData>::setConforming(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,
													   const int& maxDepth,const int& setCenter,const Real& flatness,const int& noTransform)
{
	Point3D<Real> t;
	Real s;
	return setConforming(vertices,polygons,maxDepth,setCenter,flatness,t,s,noTransform);
}
template<class NodeData,class Real,class VertexData>
template<class Vertex>
int IsoOctree<NodeData,Real,VertexData>::setConforming(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,
													   const int& maxDepth,const int& setCenter,const Real& flatness,
													   Point3D<Real>& translate,Real& scale,const int& noTransform)
{
	this->maxDepth=maxDepth;
	typename OctNode<NodeData,Real>::NodeIndex nIdx;
	std::vector<int> myTriangles;
	stdext::hash_map<long long,std::vector<int>*> triangleMap;
	MeshInfo<double> mInfo;

	mInfo.set(vertices,polygons,Real(1.1),translate,scale,noTransform);
	myTriangles.resize(mInfo.triangles.size());
	for(int i=0;i<int(mInfo.triangles.size());i++)
		myTriangles[i]=i;

	cornerValues.clear();
	Real dist;
	Point3D<Real> n,p;
	for(int c=0;c<Cube::CORNERS;c++)
	{
		int x,y,z;
		Cube::FactorCornerIndex(c,x,y,z);
		p[0]=Real(x);
		p[1]=Real(y);
		p[2]=Real(z);

		setDistanceAndNormal(myTriangles,mInfo,p,dist,n);
		cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,c,maxDepth)]=VertexData(dist,n);
	}
	if(setCenter)
	{
		Real w;
		OctNode<NodeData,Real>::CenterAndWidth(nIdx,tree.nodeData.center,w);
		setDistanceAndNormal(myTriangles,mInfo,tree.nodeData.center,dist,n);
		tree.nodeData.v=VertexData(dist,n);
	}
	setChildren(&tree,nIdx,myTriangles,mInfo,maxDepth,setCenter,flatness,&triangleMap);

	// Refine non-conforming nodes
	int forceConforming=1;
	while(forceConforming)
	{
		forceConforming=0;
		nIdx=typename OctNode<NodeData,Real>::NodeIndex();
		for(OctNode<NodeData,Real>* node=tree.nextLeaf(NULL,nIdx) ; node ; node=tree.nextLeaf(node,nIdx) )
		{
			int setChildren=0;
			for(int i=0;i<Cube::FACES && !setChildren;i++)
				if(HasFaceGrandChildren(node->faceNeighbor(i),Cube::FaceReflectFaceIndex(i,i)))
					setChildren=1;
			for(int i=0;i<Cube::EDGES && !setChildren;i++)
				if(HasEdgeGrandChildren(node->edgeNeighbor(i),Cube::EdgeReflectEdgeIndex(i)))
					setChildren=1;
			if(setChildren)
			{
				OctNode<NodeData,Real>* temp=node;
				typename OctNode<NodeData,Real>::NodeIndex pIdx=nIdx;
				long long key;
				while(temp)
				{
					key=OctNode<NodeData,Real>::CenterIndex(pIdx,maxDepth);
					if(triangleMap.find(key)==triangleMap.end() || !triangleMap[key])
					{
						temp=temp->parent;
						--pIdx;
					}
					else
						break;
				}
				if(!temp)
				{
					fprintf(stderr,"Could not find ancestor with triangles\n");
					continue;
				}
				node->initChildren();
				forceConforming=1;

				for(int i=0;i<Cube::CORNERS;i++)
				{
					Point3D<Real> ctr;
					Real w;
					OctNode<NodeData,Real>::CenterAndWidth(nIdx.child(i),ctr,w);
					for(int c=0;c<Cube::CORNERS;c++)
					{
						int x,y,z;
						Cube::FactorCornerIndex(c,x,y,z);
						p[0]=ctr[0]-w/2+w*x;
						p[1]=ctr[1]-w/2+w*y;
						p[2]=ctr[2]-w/2+w*z;

						setDistanceAndNormal(*triangleMap[key],mInfo,p,dist,n);
						long long k=OctNode<NodeData,Real>::CornerIndex(nIdx.child(i),c,maxDepth);
						if(cornerValues.find(k)==cornerValues.end() || fabs(dist)<fabs(cornerValues[k].value()))
							cornerValues[k]=VertexData(dist,n);
					}
					if(setCenter)
					{
						node->children[i].nodeData.center=ctr;
						setDistanceAndNormal(*triangleMap[key],mInfo,ctr,dist,n);
						node->children[i].nodeData.v=VertexData(dist,n);
					}
				}
			}
		}
	}

	for(stdext::hash_map<long long,std::vector<int>*>::iterator iter=triangleMap.begin();iter!=triangleMap.end();iter++)
		if(iter->second)
			delete iter->second;
	return 1;
}

template<class NodeData,class Real,class VertexData>
template<class MeshReal>
void IsoOctree<NodeData,Real,VertexData>::setDistanceAndNormal(const std::vector<int>& triangles,
															   MeshInfo<MeshReal>& mInfo,
															   const Point3D<Real>& p,
															   Real& dist,
															   Point3D<Real>& n)
{
	Point3D<MeshReal> pp,t[3];
	pp[0]=p[0];
	pp[1]=p[1];
	pp[2]=p[2];
	size_t closest;
	for(size_t i=0;i<triangles.size();i++)
	{
		MeshReal temp;
		for(int j=0;j<3;j++)
		{
			t[0][j]=mInfo.vertices[mInfo.triangles[triangles[i]].idx[0]][j];
			t[1][j]=mInfo.vertices[mInfo.triangles[triangles[i]].idx[1]][j];
			t[2][j]=mInfo.vertices[mInfo.triangles[triangles[i]].idx[2]][j];
		}

		temp=DistanceToTriangle(pp,t);
		if(!i || temp<dist )
		{
			closest=i;
			dist=(Real)temp;
		}
	}
	Point3D<MeshReal> nn;
	int vFlag;

	closest=triangles[closest];

	for(int i=0;i<3;i++) t[i]=mInfo.vertices[mInfo.triangles[closest].idx[i]];

	Point3D<MeshReal> n2=NearestPointOnTriangle(pp,t,vFlag);
	n2=(pp-n2)/Distance(pp,n2);

	switch(vFlag)
	{
	case 7:
		nn=mInfo.triangleNormals[closest];
		break;
	case 1:
		nn=mInfo.vertexNormals[mInfo.triangles[closest].idx[0]];
		break;
	case 2:
		nn=mInfo.vertexNormals[mInfo.triangles[closest].idx[1]];
		break;
	case 4:
		nn=mInfo.vertexNormals[mInfo.triangles[closest].idx[2]];
		break;
	case 3:
		nn=mInfo.edgeNormals[EdgeKey(mInfo.triangles[closest].idx[0],mInfo.triangles[closest].idx[1])];
		break;
	case 5:
		nn=mInfo.edgeNormals[EdgeKey(mInfo.triangles[closest].idx[0],mInfo.triangles[closest].idx[2])];
		break;
	case 6:
		nn=mInfo.edgeNormals[EdgeKey(mInfo.triangles[closest].idx[1],mInfo.triangles[closest].idx[2])];
		break;
	}

	if(DotProduct(nn,n2)<0)
	{
		dist=-dist;
		n2*=-1;
	}
	n[0]=(Real)n2[0];
	n[1]=(Real)n2[1];
	n[2]=(Real)n2[2];
}

template<class NodeData,class Real,class VertexData>
template<class MeshReal>
void IsoOctree<NodeData,Real,VertexData>::setChildren(OctNode<NodeData,Real>* node,
													  const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
													  const std::vector<int>& triangles,MeshInfo<MeshReal>& mInfo,
													  const int& maxDepth,const int& setCenter,
													  const Real& flatness,
													  stdext::hash_map<long long,std::vector<int>*>* triangleMap)
{
	long long key;
	Real w,dist;
	Point3D<Real> ctr,n,p;
	if(!triangles.size())		return;
	if(nIdx.depth==maxDepth)	return;

	if(triangleMap)
	{
		std::vector<int>* myTriangles=new std::vector<int>();
		long long key=OctNode<NodeData,Real>::CenterIndex(nIdx,maxDepth);
		myTriangles->resize(triangles.size());
		for(size_t i=0;i<triangles.size();i++) (*myTriangles)[i]=triangles[i];
		(*triangleMap)[key]=myTriangles;
	}

	if(flatness>0)
	{
		MeshReal area;
		Point3D<MeshReal> mc;
		area=0;
		mc[0]=mc[1]=mc[2]=0;
		for(size_t i=0;i<triangles.size();i++)
		{
			area+=Length(mInfo.triangleNormals[triangles[i]]);
			mc+=mInfo.triangleNormals[triangles[i]];
		}
		if(Length(mc)/area>flatness) return;
	}

	if(!node->children)	node->initChildren();
	OctNode<NodeData,Real>::CenterAndWidth(nIdx,ctr,w);

	// Set the center
	setDistanceAndNormal(triangles,mInfo,ctr,dist,n);
	key=OctNode<NodeData,Real>::CenterIndex(nIdx,maxDepth);
	if( cornerValues.find(key)==cornerValues.end() || fabs(cornerValues[key].value())>fabs(dist) )
		cornerValues[key]=VertexData(dist,n);

	// Set the edge mid-points
	for(int i=0;i<Cube::EDGES;i++)
	{
		int o,i1,i2;
		Cube::FactorEdgeIndex(i,o,i1,i2);
		p=ctr;
		p[0]-=w/2;
		p[1]-=w/2;
		p[2]-=w/2;
		p[o]=ctr[o];
		switch(o)
		{
		case 0:
			p[1]+=w*i1;
			p[2]+=w*i2;
			break;
		case 1:
			p[0]+=w*i1;
			p[2]+=w*i2;
			break;
		case 2:
			p[0]+=w*i1;
			p[1]+=w*i2;
			break;
		}
		setDistanceAndNormal(triangles,mInfo,p,dist,n);
		key=OctNode<NodeData,Real>::EdgeIndex(nIdx,i,maxDepth);
		if( cornerValues.find(key)==cornerValues.end() || fabs(cornerValues[key].value())>fabs(dist) )
			cornerValues[key]=VertexData(dist,n);
	}

	// set the face mid-points
	for(int i=0;i<Cube::FACES;i++)
	{
		int dir,off;
		Cube::FactorFaceIndex(i,dir,off);
		p=ctr;
		p[dir]+=-w/2+w*off;
		setDistanceAndNormal(triangles,mInfo,p,dist,n);
		key=OctNode<NodeData,Real>::FaceIndex(nIdx,i,maxDepth);
		if( cornerValues.find(key)==cornerValues.end() || fabs(cornerValues[key].value())>fabs(dist) )
			cornerValues[key]=VertexData(dist,n);
	}

	int retCount=0;
	for(int i=0;i<Cube::CORNERS;i++)
	{
		std::vector<int> myTriangles;
		OctNode<NodeData,Real>::CenterAndWidth(nIdx.child(i),ctr,w);
		if(setCenter)
		{
			OctNode<NodeData,Real>::CenterAndWidth(nIdx.child(i),node->children[i].nodeData.center,w);
			setDistanceAndNormal(triangles,mInfo,node->children[i].nodeData.center,dist,n);
			node->children[i].nodeData.v=VertexData(dist,n);
		}

		for(size_t j=0;j<triangles.size();j++)
		{
			Point3D<MeshReal> t[3];
			for(int k=0;k<3;k++) t[k]=mInfo.vertices[mInfo.triangles[triangles[j]].idx[k]];
			Point3D<MeshReal> ctr2;
			MeshReal w2=w;
			ctr2[0]=ctr[0];
			ctr2[1]=ctr[1];
			ctr2[2]=ctr[2];
			if(TriangleInCube(ctr2,w2,t)) myTriangles.push_back(triangles[j]);
		}
		setChildren(&node->children[i],nIdx.child(i),myTriangles,mInfo,maxDepth,setCenter,flatness,triangleMap);
		if(myTriangles.size()) retCount++;
	}
	if(triangleMap && retCount==8)
	{
		long long key=OctNode<NodeData,Real>::CenterIndex(nIdx,maxDepth);
		if(triangleMap->find(key)!=triangleMap->end() && (*triangleMap)[key]!=NULL)
		{
			delete (*triangleMap)[key];
			triangleMap->erase(key);
		}
	}
}

template<class NodeData,class Real,class VertexData>
void IsoOctree<NodeData,Real,VertexData>::interpolateSharedValues(void)
{
	VertexData values[Cube::CORNERS];
	nKey.set(maxDepth);

	long long key;
	int c1,c2,c3,c4;
	int idx[3];
	for(int d=0;d<maxDepth;d++)
	{
		typename OctNode<NodeData,Real>::NodeIndex nIdx;
		for(OctNode<NodeData,Real>* temp=tree.nextNode(NULL,nIdx) ; temp ; temp=tree.nextNode(temp,nIdx) )
		{
			if(nIdx.depth==d && temp->children)
			{
				for(int c=0;c<Cube::CORNERS;c++)
					values[c]=cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,c,maxDepth)];

				nKey.getNeighbors(temp);
				for(int f=0;f<Cube::FACES;f++)
				{
					Cube::FaceCorners(f,c1,c2,c3,c4);
					key=OctNode<NodeData,Real>::FaceIndex(nIdx,f,maxDepth);
					NeighborKey<NodeData,Real>::FaceIndex(f,idx);
					if(	!nKey.neighbors[d].neighbors[idx[0]][idx[1]][idx[2]] ||
						!nKey.neighbors[d].neighbors[idx[0]][idx[1]][idx[2]]->children)
						cornerValues[key]=(values[c1]+values[c2]+values[c3]+values[c4]) / 4;
				}
				for(int e=0;e<Cube::EDGES;e++)
				{
					Cube::EdgeCorners(e,c1,c2);
					key=OctNode<NodeData,Real>::EdgeIndex(nIdx,e,maxDepth);

					int f1,f2;
					Cube::FacesAdjacentToEdge(e,f1,f2);
					NeighborKey<NodeData,Real>::FaceIndex(f1,idx);
					if(	!nKey.neighbors[d].neighbors[idx[0]][idx[1]][idx[2]] ||
						!nKey.neighbors[d].neighbors[idx[0]][idx[1]][idx[2]]->children)
						cornerValues[key]=(values[c1]+values[c2])/2;

					NeighborKey<NodeData,Real>::FaceIndex(f2,idx);
					if(	!nKey.neighbors[d].neighbors[idx[0]][idx[1]][idx[2]] ||
						!nKey.neighbors[d].neighbors[idx[0]][idx[1]][idx[2]]->children)
						cornerValues[key]=(values[c1]+values[c2])/2;
					NeighborKey<NodeData,Real>::EdgeIndex(e,idx);
					if(	!nKey.neighbors[d].neighbors[idx[0]][idx[1]][idx[2]] ||
						!nKey.neighbors[d].neighbors[idx[0]][idx[1]][idx[2]]->children)
						cornerValues[key]=(values[c1]+values[c2])/2;
				}
			}
		}
	}
}


template<class NodeData,class Real,class VertexData>
void IsoOctree<NodeData,Real,VertexData>::resetValues(void)
{
	stdext::hash_map<long long,VertexData> tempValues;
	VertexData values[Cube::CORNERS];
	nKey.set(maxDepth);
	typename OctNode<NodeData,Real>::NodeIndex nIdx;
	for(OctNode<NodeData,Real>* temp=tree.nextLeaf(NULL,nIdx) ; temp ; temp=tree.nextLeaf(temp,nIdx) )
		for(int i=0;i<Cube::CORNERS;i++)
		{
			long long key=OctNode<NodeData,Real>::CornerIndex(nIdx,i,maxDepth);
			tempValues[key]=cornerValues[key];
		}
	cornerValues.clear();
	for(typename stdext::hash_map<long long,VertexData>::iterator iter=tempValues.begin();iter!=tempValues.end();iter++)
		cornerValues[iter->first]=iter->second;
}
template<class NodeData,class Real,class VertexData>
int IsoOctree<NodeData,Real,VertexData>::write(FILE* fp,int writeData) const
{
	if(!tree.write(fp,writeData))
		return 0;
	fwrite(&maxDepth,sizeof(int),1,fp);
	int hash_size=int(cornerValues.size());
	fwrite(&hash_size,sizeof(int),1,fp);
	for(typename stdext::hash_map<long long,VertexData>::const_iterator iter=cornerValues.begin();iter!=cornerValues.end();iter++)
	{
		fwrite(&iter->first,sizeof(long long),1,fp);
		fwrite(&iter->second,sizeof(VertexData),1,fp);
	}
	return 1;
}
template<class NodeData,class Real,class VertexData>
int IsoOctree<NodeData,Real,VertexData>::read(FILE* fp,int readData)
{
	if(!tree.read(fp,readData))	return 0;
	fread(&maxDepth,sizeof(int),1,fp);
	int hash_size;
	fread(&hash_size,sizeof(int),1,fp);
	cornerValues.clear();
	for(int i=0;i<hash_size;i++)
	{
		long long key;
		VertexData value;
		fread(&key,sizeof(long long),1,fp);
		fread(&value,sizeof(VertexData),1,fp);
		cornerValues[key]=value;
	}
	return 1;
}

template<class NodeData,class Real,class VertexData>
int IsoOctree<NodeData,Real,VertexData>::getRootPosition(const OctNode<NodeData,Real>* node,
														 const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
														 const int& eIndex,const Real& isoValue,Point3D<Real>& position)
{
	int c0,c1;
	Cube::EdgeCorners(eIndex,c0,c1);

	if(!MarchingCubes::HasEdgeRoots(node->nodeData.mcIndex,eIndex))
		return 0;

	Real w;
	Point3D<Real> p1,p2;
	int x,y,z;

	OctNode<NodeData,Real>::CenterAndWidth(nIdx,p1,w);
	OctNode<NodeData,Real>::CenterAndWidth(nIdx,p2,w);

	Cube::FactorCornerIndex(c0,x,y,z);
	p1[0]+=-w/2+w*x;
	p1[1]+=-w/2+w*y;
	p1[2]+=-w/2+w*z;
	Cube::FactorCornerIndex(c1,x,y,z);
	p2[0]+=-w/2+w*x;
	p2[1]+=-w/2+w*y;
	p2[2]+=-w/2+w*z;
	position=VertexData::RootPosition(isoValue,p1,p2,
		cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,c0,maxDepth)],
		cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,c1,maxDepth)]);
	return 1;
}
template<class NodeData,class Real,class VertexData>
long long IsoOctree<NodeData,Real,VertexData>::getRootKey(const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const int& edgeIndex)
{
	int offset,eIndex[2],o,i1,i2;
	Cube::FactorEdgeIndex(edgeIndex,o,i1,i2);
	offset=BinaryNode<Real>::Index(nIdx.depth,nIdx.offset[o]);
	switch(o)
	{
	case 0:
		eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIdx.depth,nIdx.offset[1],i1);
		eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIdx.depth,nIdx.offset[2],i2);
		break;
	case 1:
		eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIdx.depth,nIdx.offset[0],i1);
		eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIdx.depth,nIdx.offset[2],i2);
		break;
	case 2:
		eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIdx.depth,nIdx.offset[0],i1);
		eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIdx.depth,nIdx.offset[1],i2);
		break;
	}
	return (long long)(o) | (long long)(eIndex[0])<<5 | (long long)(eIndex[1])<<25 | (long long)(offset)<<45;
}
template<class NodeData,class Real,class VertexData>
int IsoOctree<NodeData,Real,VertexData>::getRootIndex(OctNode<NodeData,Real>* node,
													  const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
													  const int& edgeIndex,RootInfo& ri){
	int c1,c2,f1,f2;
	const OctNode<NodeData,Real> *temp,*finest;
	int finestIndex;
	typename OctNode<NodeData,Real>::NodeIndex finestNIdx=nIdx;

	// The assumption is that the super-edge has a root along it. 
	if(!(MarchingCubes::HasEdgeRoots(node->nodeData.mcIndex,edgeIndex)))
		return 0;
#if USE_MAX_DEPTH_SPEED_UP
	if(nIdx.depth==maxDepth)
	{
		ri.node=node;
		ri.edgeIndex=edgeIndex;
		ri.nIdx=nIdx;
		ri.key=getRootKey(nIdx,edgeIndex);
		return 1;
	}
#endif // USE_MAX_DEPTH_SPEED_UP


	finest=node;
	finestIndex=edgeIndex;

	Cube::FacesAdjacentToEdge(edgeIndex,f1,f2);
	if(nIdx.depth<maxDepth)
	{
		if(!node->children)
		{
			temp=node->faceNeighbor(f1);
			if(temp && temp->children)
			{
				finest=temp;
				finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f1);
				int dir,off;
				Cube::FactorFaceIndex(f1,dir,off);
				if(off)
					finestNIdx.offset[dir]++;
				else
					finestNIdx.offset[dir]--;
			}
			else
			{
				temp=node->faceNeighbor(f2);
				if(temp && temp->children)
				{
					finest=temp;
					finestIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f2);
					int dir,off;
					Cube::FactorFaceIndex(f2,dir,off);
					if(off)
						finestNIdx.offset[dir]++;
					else
						finestNIdx.offset[dir]--;
				}
				else
				{
					temp=node->edgeNeighbor(edgeIndex);
					if(temp && temp->children)
					{
						finest=temp;
						finestIndex=Cube::EdgeReflectEdgeIndex(edgeIndex);
						int o,i1,i2;
						Cube::FactorEdgeIndex(edgeIndex,o,i1,i2);
						switch(o)
						{
						case 0:
							if(i1)	finestNIdx.offset[1]++;
							else	finestNIdx.offset[1]--;
							if(i2)	finestNIdx.offset[2]++;
							else	finestNIdx.offset[2]--;
							break;
						case 1:
							if(i1)	finestNIdx.offset[0]++;
							else	finestNIdx.offset[0]--;
							if(i2)	finestNIdx.offset[2]++;
							else	finestNIdx.offset[2]--;
							break;
						case 2:
							if(i1)	finestNIdx.offset[0]++;
							else	finestNIdx.offset[0]--;
							if(i2)	finestNIdx.offset[1]++;
							else	finestNIdx.offset[1]--;
							break;
						}
					}
				}
			}
		}
	}

	Cube::EdgeCorners(finestIndex,c1,c2);
	if(finest->children)
	{
		if		(getRootIndex(&finest->children[c1],finestNIdx.child(c1),finestIndex,ri))	{return 1;}
		else if	(getRootIndex(&finest->children[c2],finestNIdx.child(c2),finestIndex,ri))	{return 1;}
		else																				{
			fprintf(stderr,"Failed to find root index\n");
			return 0;
		}
	}
	else
	{
		ri.nIdx=finestNIdx;
		ri.node=finest;
		ri.edgeIndex=finestIndex;

		int o,i1,i2;
		Cube::FactorEdgeIndex(finestIndex,o,i1,i2);
		int offset,eIndex[2];
		offset=BinaryNode<Real>::Index(finestNIdx.depth,finestNIdx.offset[o]);
		switch(o)
		{
		case 0:
			eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,finestNIdx.depth,finestNIdx.offset[1],i1);
			eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,finestNIdx.depth,finestNIdx.offset[2],i2);
			break;
		case 1:
			eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,finestNIdx.depth,finestNIdx.offset[0],i1);
			eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,finestNIdx.depth,finestNIdx.offset[2],i2);
			break;
		case 2:
			eIndex[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,finestNIdx.depth,finestNIdx.offset[0],i1);
			eIndex[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,finestNIdx.depth,finestNIdx.offset[1],i2);
			break;
		}
		ri.key= (long long)(o) | (long long)(eIndex[0])<<5 | (long long)(eIndex[1])<<25 | (long long)(offset)<<45;
		return 1;
	}
}
template<class NodeData,class Real,class VertexData>
template<class Vertex>
void IsoOctree<NodeData,Real,VertexData>::getRoots(OctNode<NodeData,Real>* node,
												   const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
												   const Real& isoValue,
												   stdext::hash_map<long long,int>& roots,std::vector<Vertex>& vertices)
{
	Point3D<Real> position;

	int eIndex;
	RootInfo ri;

	if(!MarchingCubes::HasRoots(node->nodeData.mcIndex))
		return;

	for(eIndex=0;eIndex<Cube::EDGES;eIndex++)
	{
		if(!(MarchingCubes::HasEdgeRoots(node->nodeData.mcIndex,eIndex)))
			continue;

		if(getRootIndex(node,nIdx,eIndex,ri))
		{
			if(roots.find(ri.key)==roots.end()){
				getRootPosition(ri.node,ri.nIdx,ri.edgeIndex,isoValue,position);
				vertices.push_back(position);
				roots[ri.key]=int(vertices.size())-1;
			}
		}
		else
			fprintf(stderr,"Failed to get root index\n");
	}
}
template<class NodeData,class Real,class VertexData>
int IsoOctree<NodeData,Real,VertexData>::getRootPair(const RootInfo& ri,const int& maxDepth,RootInfo& pair)
{
	const OctNode<NodeData,Real>* node=ri.node;
	typename OctNode<NodeData,Real>::NodeIndex nIdx=ri.nIdx;
	int c1,c2,c;
	Cube::EdgeCorners(ri.edgeIndex,c1,c2);
	while(node->parent)
	{
		c=int(node-node->parent->children);
		if(c!=c1 && c!=c2)
			return 0;
		if(!MarchingCubes::HasEdgeRoots(node->parent->nodeData.mcIndex,ri.edgeIndex))
		{
			if(c==c1)
				return getRootIndex(&node->parent->children[c2],nIdx.parent().child(c2),ri.edgeIndex,pair);
			else
				return getRootIndex(&node->parent->children[c1],nIdx.parent().child(c2),ri.edgeIndex,pair);
		}
		node=node->parent;
		--nIdx;
	}
	return 0;
}
template<class NodeData,class Real,class VertexData>
void IsoOctree<NodeData,Real,VertexData>::getIsoFaceEdges(OctNode<NodeData,Real>* node,
														  const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
														  const int& faceIndex,std::vector<std::pair<RootInfo,RootInfo> >& edges,const int& flip,const int& useFull)
{
	int c1,c2,c3,c4;
	if(node->children)
	{
		Cube::FaceCorners(faceIndex,c1,c2,c3,c4);
		getIsoFaceEdges(&node->children[c1],nIdx.child(c1),faceIndex,edges,flip,useFull);
		getIsoFaceEdges(&node->children[c2],nIdx.child(c2),faceIndex,edges,flip,useFull);
		getIsoFaceEdges(&node->children[c3],nIdx.child(c3),faceIndex,edges,flip,useFull);
		getIsoFaceEdges(&node->children[c4],nIdx.child(c4),faceIndex,edges,flip,useFull);
	}
	else
	{
		int idx=node->nodeData.mcIndex;

		RootInfo ri1,ri2;
		const std::vector<std::vector<int> >& table=MarchingCubes::caseTable(idx,useFull);
		for(size_t i=0;i<table.size();i++)
		{
			size_t pSize=table[i].size();
			for(size_t j=0;j<pSize;j++)
			{
				if(faceIndex==Cube::FaceAdjacentToEdges(table[i][j],table[i][(j+1)%pSize]))
					if(getRootIndex(node,nIdx,table[i][j],ri1) && getRootIndex(node,nIdx,table[i][(j+1)%pSize],ri2))
						if(flip)
							edges.push_back(std::pair<RootInfo,RootInfo>(ri2,ri1));
						else
							edges.push_back(std::pair<RootInfo,RootInfo>(ri1,ri2));
					else
						fprintf(stderr,"Bad Edge 1: %lld %lld\n",ri1.key,ri2.key);
			}
		}
	}
}
template<class NodeData,class Real,class VertexData>
void IsoOctree<NodeData,Real,VertexData>::getIsoPolygons(OctNode<NodeData,Real>* node,
														 const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
														 stdext::hash_map<long long,int>& roots,
														 std::vector<std::vector<int> >& polygons,
														 const int& useFull)
{
	std::vector<std::pair<long long,long long> > edges;
	typename stdext::hash_map<long long,std::pair<RootInfo,int> >::iterator iter;
	stdext::hash_map<long long,std::pair<RootInfo,int> > vertexCount;
	std::vector<std::pair<RootInfo,RootInfo> > riEdges;

#if USE_MAX_DEPTH_SPEED_UP
	if(nIdx.depth==maxDepth) // Just run the standard marching cubes...
	{
		RootInfo ri;
		size_t pIndex=polygons.size();
		int idx=node->nodeData.mcIndex;
		const std::vector<std::vector<int> >& table=MarchingCubes::caseTable(idx,useFull);

		polygons.resize(pIndex+table.size());
		for(size_t i=0;i<table.size();i++)
		{
			polygons[pIndex+i].resize(table[i].size());
			for(size_t j=0;j<table[i].size();j++)
				if(getRootIndex(node,nIdx,table[i][j],ri))
					polygons[pIndex+i][j]=roots[ri.key];
				else
					fprintf(stderr,"Bad Edge 1: %lld\n",ri.key);
		}
		return;
	}
	else
	{
		int x[3];
		nKey.getNeighbors(node);
		x[0]=x[1]=x[2]=1;
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<2;j++)
			{
				x[i]=j<<1;
				if(!nKey.neighbors[nIdx.depth].neighbors[x[0]][x[1]][x[2]] || !nKey.neighbors[nIdx.depth].neighbors[x[0]][x[1]][x[2]]->children)
					getIsoFaceEdges(node,nIdx,Cube::FaceIndex(i,j),riEdges,0,useFull);
				else
				{
					typename OctNode<NodeData,Real>::NodeIndex idx=nIdx;
					if(j)	idx.offset[i]++;
					else	idx.offset[i]--;
					getIsoFaceEdges(nKey.neighbors[idx.depth].neighbors[x[0]][x[1]][x[2]],idx,Cube::FaceIndex(i,j^1),riEdges,1,useFull);
				}
			}
			x[i]=1;
		}
	}
	for(size_t i=0;i<riEdges.size();i++)
	{
		edges.push_back(std::pair<long long,long long>(riEdges[i].first.key,riEdges[i].second.key));
		iter=vertexCount.find(riEdges[i].first.key);
		if(iter==vertexCount.end())
		{
			vertexCount[riEdges[i].first.key].first=riEdges[i].first;
			vertexCount[riEdges[i].first.key].second=0;
		}
		iter=vertexCount.find(riEdges[i].second.key);
		if(iter==vertexCount.end())
		{
			vertexCount[riEdges[i].second.key].first=riEdges[i].second;
			vertexCount[riEdges[i].second.key].second=0;
		}
		vertexCount[riEdges[i].first.key ].second++;
		vertexCount[riEdges[i].second.key].second--;
	}
#else // !USE_MAX_DEPTH_SPEED_UP
	for(fIndex=0;fIndex<Cube::FACES;fIndex++)
	{
		ref=Cube::FaceReflectFaceIndex(fIndex,fIndex);
		temp=node->faceNeighbor(fIndex);

		riEdges.clear();
		if(temp && temp->children)
			getIsoFaceEdges(temp,ref,riEdges,1);
		else
			getIsoFaceEdges(node,fIndex,riEdges,0);
		for(int i=0;i<riEdges.size();i++)
		{
			edges.push_back(std::pair<long long,long long>(riEdges[i].first.key,riEdges[i].second.key));
			iter=vertexCount.find(riEdges[i].first.key);
			if(iter==vertexCount.end())
			{
				vertexCount[riEdges[i].first.key].first=riEdges[i].first;
				vertexCount[riEdges[i].first.key].second=0;
			}
			iter=vertexCount.find(riEdges[i].second.key);
			if(iter==vertexCount.end())
			{
				vertexCount[riEdges[i].second.key].first=riEdges[i].second;
				vertexCount[riEdges[i].second.key].second=0;
			}
			vertexCount[riEdges[i].first.key ].second++;
			vertexCount[riEdges[i].second.key].second--;
		}
	}
#endif // USE_MAX_DEPTH_SPEED_UP
	for(int i=0;i<int(edges.size());i++)
	{
		iter=vertexCount.find(edges[i].first);
		if(iter==vertexCount.end())
			printf("Could not find vertex: %lld\n",edges[i].first);
		else if(vertexCount[edges[i].first].second)
		{
			RootInfo ri;
			if(!getRootPair(vertexCount[edges[i].first].first,maxDepth,ri))
				fprintf(stderr,"Failed to get root pair 1: %lld %d\n",edges[i].first,vertexCount[edges[i].first].second);
			iter=vertexCount.find(ri.key);
			if(iter==vertexCount.end())
				printf("Vertex pair not in list\n");
			else
			{
				edges.push_back(std::pair<long long,long long>(ri.key,edges[i].first));
				vertexCount[ri.key].second++;
				vertexCount[edges[i].first].second--;
			}
		}

		iter=vertexCount.find(edges[i].second);
		if(iter==vertexCount.end())
			printf("Could not find vertex: %lld\n",edges[i].second);
		else if(vertexCount[edges[i].second].second)
		{
			RootInfo ri;
			if(!getRootPair(vertexCount[edges[i].second].first,maxDepth,ri))
				fprintf(stderr,"Failed to get root pair 2: %lld %d\n",edges[i].second,vertexCount[edges[i].second].second);
			iter=vertexCount.find(ri.key);
			if(iter==vertexCount.end())
				printf("Vertex pair not in list\n");
			else{
				edges.push_back(std::pair<long long,long long>(edges[i].second,ri.key));
				vertexCount[edges[i].second].second++;
				vertexCount[ri.key].second--;
			}
		}
	}
	getEdgeLoops(edges,roots,polygons);
}
template<class NodeData,class Real,class VertexData>
template<class C>
void IsoOctree<NodeData,Real,VertexData>::getEdgeLoops(std::vector<std::pair<C,C> >& edges,
													   stdext::hash_map<C,int>& roots,
													   std::vector<std::vector<int> >& polygons)
{
	size_t polygonSize=polygons.size();
	C frontIdx,backIdx;
	std::pair<C,C> e,temp;

	while(edges.size())
	{
		std::vector<std::pair<C,C> > front,back;
		e=edges[0];
		polygons.resize(polygonSize+1);
		edges[0]=edges[edges.size()-1];
		edges.pop_back();
		frontIdx=e.second;
		backIdx=e.first;
		for(int j=int(edges.size())-1;j>=0;j--){
			if(edges[j].first==frontIdx || edges[j].second==frontIdx){
				if(edges[j].first==frontIdx)	{temp=edges[j];}
				else							{temp.first=edges[j].second;temp.second=edges[j].first;}
				frontIdx=temp.second;
				front.push_back(temp);
				edges[j]=edges[edges.size()-1];
				edges.pop_back();
				j=int(edges.size());
			}
			else if(edges[j].first==backIdx || edges[j].second==backIdx){
				if(edges[j].second==backIdx)	{temp=edges[j];}
				else							{temp.first=edges[j].second;temp.second=edges[j].first;}
				backIdx=temp.first;
				back.push_back(temp);
				edges[j]=edges[edges.size()-1];
				edges.pop_back();
				j=int(edges.size());
			}
		}
		polygons[polygonSize].resize(back.size()+front.size()+1);
		int idx=0;
		for(int j=int(back.size())-1;j>=0;j--)	polygons[polygonSize][idx++]=roots[back[j].first];
		polygons[polygonSize][idx++]=roots[e.first];
		for(int j=0;j<int(front.size());j++)	polygons[polygonSize][idx++]=roots[front[j].first];
		polygonSize++;
	}
}
template<class NodeData,class Real,class VertexData>
template<class C>
void IsoOctree<NodeData,Real,VertexData>::getEdgeLoops(std::vector<std::pair<C,C> >& edges,
													   std::vector<std::vector<C> >& polygons)
{
	int polygonSize=polygons.size();
	C frontIdx,backIdx;
	std::pair<C,C> e,temp;

	while(edges.size())
	{
		std::vector<std::pair<C,C> > front,back;
		e=edges[0];
		polygons.resize(polygonSize+1);
		edges[0]=edges[edges.size()-1];
		edges.pop_back();
		frontIdx=e.second;
		backIdx=e.first;
		for(int j=int(edges.size())-1;j>=0;j--){
			if(edges[j].first==frontIdx || edges[j].second==frontIdx){
				if(edges[j].first==frontIdx)	{temp=edges[j];}
				else							{temp.first=edges[j].second;temp.second=edges[j].first;}
				frontIdx=temp.second;
				front.push_back(temp);
				edges[j]=edges[edges.size()-1];
				edges.pop_back();
				j=int(edges.size());
			}
			else if(edges[j].first==backIdx || edges[j].second==backIdx){
				if(edges[j].second==backIdx)	{temp=edges[j];}
				else							{temp.first=edges[j].second;temp.second=edges[j].first;}
				backIdx=temp.first;
				back.push_back(temp);
				edges[j]=edges[edges.size()-1];
				edges.pop_back();
				j=int(edges.size());
			}
		}
		polygons[polygonSize].resize(back.size()+front.size()+1);
		int idx=0;
		for(int j=int(back.size())-1;j>=0;j--)
			polygons[polygonSize][idx++]=back[j].first;
		polygons[polygonSize][idx++]=e.first;
		for(int j=0;j<int(front.size());j++)
			polygons[polygonSize][idx++]=front[j].first;
		polygonSize++;
	}
}
template<class NodeData,class Real,class VertexData>
void IsoOctree<NodeData,Real,VertexData>::setMCIndex(const Real& isoValue,const int& useFull)
{
	OctNode<NodeData,Real>* temp;
	Real cValues[Cube::CORNERS];

	// Clear the indices
	for(temp=tree.nextNode(NULL) ; temp ; temp=tree.nextNode(temp) )
		temp->nodeData.mcIndex=0;

	// Get the values at the leaf nodes and propogate up to the parents
	typename OctNode<NodeData,Real>::NodeIndex nIdx;
	for(temp=tree.nextLeaf(NULL,nIdx) ; temp ; temp=tree.nextLeaf(temp,nIdx) )
	{
		for(int i=0;i<Cube::CORNERS;i++)
		{
			if(cornerValues.find(OctNode<NodeData,Real>::CornerIndex(nIdx,i,maxDepth))==cornerValues.end())
				fprintf(stderr,"Could not find value in corner value table!\n");
			cValues[i]=cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,i,maxDepth)].value();
		}
		if(useFull)
			temp->nodeData.mcIndex=MarchingCubes::GetFullIndex(cValues,isoValue);
		else
			temp->nodeData.mcIndex=MarchingCubes::GetIndex(cValues,isoValue);
		if(temp->parent)
		{
			int cIndex=int(temp-temp->parent->children);
			int bitFlag = temp->nodeData.mcIndex & (1<<cIndex);
			if(bitFlag)
			{
				OctNode<NodeData,Real> *parent,*child;
				child=temp;
				parent=temp->parent;
				while(parent && (child-parent->children)==cIndex)
				{
					parent->nodeData.mcIndex |= bitFlag;
					child=parent;
					parent=parent->parent;
				}
			}
		}
	}
}
template<class NodeData,class Real,class VertexData>
template<class Vertex>
void IsoOctree<NodeData,Real,VertexData>::getDualIsoSurface(const Real& isoValue,
															std::vector<Vertex>& vertices,
															std::vector<std::vector<int> >& polygons,
															const int& useFull)
{
	OctNode<NodeData,Real>* temp;
	OctNode<NodeData,Real>* neighbors[Cube::CORNERS];
	stdext::hash_map<long long,char> cornerSet;
	Real cValues[Cube::CORNERS];
	VertexData cVertexData[Cube::CORNERS];
	Point3D<Real> p;
	int c1,c2;
	EdgeKey eKey;
	std::vector<std::vector<int> > polys;
	std::vector<std::pair<int,int> > eList;
	stdext::hash_map<EdgeKey,int,HashEdgeKey> vTable;

	nKey.set(maxDepth);

	for(int d=maxDepth;d>=0;d--)
	{
		typename OctNode<NodeData,Real>::NodeIndex nIdx;
		for(temp=tree.nextNode(NULL,nIdx) ; temp ; temp=tree.nextNode(temp,nIdx) )
		{
			if(nIdx.depth==d && !temp->children)
				for(int c=0;c<Cube::CORNERS;c++)
				{
					long long key=OctNode<NodeData,Real>::CornerIndex(nIdx,c,maxDepth);
					if(cornerSet.find(key)==cornerSet.end())
					{
						cornerSet[key]=1;
						nKey.GetCornerNeighbors(temp,c,neighbors);
						int count=0;
						for(int i=0;i<Cube::CORNERS;i++)
							if(neighbors[i])
							{
								cVertexData[i]=neighbors[i]->nodeData.v;
								cValues[i]=cVertexData[i].value();
								count++;
							}
						if(count==8)
						{
							int idx;
							if(useFull)
								idx=MarchingCubes::GetFullIndex(cValues,isoValue);
							else
								idx=MarchingCubes::GetIndex(cValues,isoValue);
							// Add the necessary vertices
							const std::vector<std::vector<int> >& table=MarchingCubes::caseTable(idx,useFull);
							for(int i=0;i<int(table.size());i++)
							{
								for(int j=0;j<int(table[i].size());j++)
								{
									Cube::EdgeCorners(table[i][j],c1,c2);
									eKey=EdgeKey(size_t(neighbors[c1]),size_t(neighbors[c2]));
									if(vTable.find(eKey)==vTable.end())
									{
										vTable[eKey]=int(vertices.size());
										p=VertexData::RootPosition(isoValue,neighbors[c1]->nodeData.center,neighbors[c2]->nodeData.center,
											cVertexData[c1],cVertexData[c2]);
										vertices.push_back(p);
									}
								}
							}
							// Add the loops
							for(size_t i=0;i<table.size();i++)
							{
								std::vector<int> polygon;
								for(size_t j=0;j<table[i].size();j++)
								{
									Cube::EdgeCorners(table[i][j],c1,c2);
									eKey=EdgeKey(size_t(neighbors[c1]),size_t(neighbors[c2]));
									polygon.push_back(vTable[eKey]);
								}
								polygons.push_back(polygon);
							}
						}
					}
			}
		}
	}
}

template<class NodeData,class Real,class VertexData>
template<class Vertex>
void IsoOctree<NodeData,Real,VertexData>::getIsoSurface(const Real& isoValue,
														std::vector<Vertex>& vertices,
														std::vector<std::vector<int> >& polygons,
														const int& useFull)
{
	OctNode<NodeData,Real>* temp;
	stdext::hash_map<long long,int> roots;
	nKey.set(maxDepth);

	// Set the marching cubes values
	setMCIndex(isoValue,useFull);

	// Set the iso-vertex positions
	typename OctNode<NodeData,Real>::NodeIndex nIdx;
	for(temp=tree.nextLeaf(NULL,nIdx) ; temp ; temp=tree.nextLeaf(temp,nIdx) )
		getRoots(temp,nIdx,isoValue,roots,vertices);

	nIdx=typename OctNode<NodeData,Real>::NodeIndex();
	for(temp=tree.nextLeaf(NULL,nIdx) ; temp ; temp=tree.nextLeaf(temp,nIdx) )
		getIsoPolygons(temp,nIdx,roots,polygons,useFull);
}
template<class NodeData,class Real,class VertexData>
template<class Vertex>
void IsoOctree<NodeData,Real,VertexData>::getIsoSoup(const Real& isoValue,
													 std::vector<Vertex>& vertices,
													 std::vector<std::vector<int> >& polygons,
													 const int& useFull)
{
	OctNode<NodeData,Real>* temp;

	nKey.set(maxDepth);

	// Set the marching cubes values
	setMCIndex(isoValue,useFull);

	typename OctNode<NodeData,Real>::NodeIndex nIdx;
	for(temp=tree.nextLeaf(NULL,nIdx) ; temp ; temp=tree.nextLeaf(temp,nIdx) )
	{
		int pIndex=int(polygons.size());
		int idx=temp->nodeData.mcIndex;
		const std::vector<std::vector<int> >& table=MarchingCubes::caseTable(idx,useFull);
		polygons.resize(pIndex+table.size());
		for(size_t i=0;i<table.size();i++)
		{
			polygons[pIndex+i].resize(table[i].size());
			for(size_t j=0;j<table[i].size();j++)
			{
				polygons[pIndex+i][j]=int(vertices.size());
				Point3D<Real> position;
				getRootPosition(temp,nIdx,table[i][j],isoValue,position);
				vertices.push_back(position);
			}
		}
	}
}

template<class NodeData,class Real,class VertexData>
void IsoOctree<NodeData,Real,VertexData>::setNormalFlatness(const Real& isoValue,stdext::hash_map<long long,std::pair<Point3D<Real>,Real> >& flatness)
{
	OctNode<NodeData,Real>* temp;
	stdext::hash_map<long long,int> roots;
	std::vector<Point3D<Real> > vertices;

	nKey.set(maxDepth);

	// Set the marching cubes values
	setMCIndex(isoValue,0);

	// Set the iso-vertex positions
	typename OctNode<NodeData,Real>::NodeIndex nIdx;
	for(temp=tree.nextLeaf(NULL,nIdx) ; temp ; temp=tree.nextLeaf(temp,nIdx) )
		getRoots(temp,nIdx,isoValue,roots,vertices);

	// Get the polygon normals
	nIdx=typename OctNode<NodeData,Real>::NodeIndex();
	for(temp=tree.nextLeaf(NULL,nIdx) ; temp ; temp=tree.nextLeaf(temp,nIdx) )
	{
		Point3D<Real> normal;
		Real area=0;

		normal[0]=normal[1]=normal[2]=0;

		std::vector<std::vector<int> > polygons;
		getIsoPolygons(temp,nIdx,roots,polygons,0);
		if(!polygons.size())
			continue;
		for(size_t i=0;i<polygons.size();i++)
		{
			Point3D<Real> n;
			n[0]=n[1]=n[2]=0;
			for(size_t j=0;j<polygons[i].size();j++)
			{
				Point3D<Real> temp;
				CrossProduct(vertices[polygons[i][j]],vertices[polygons[i][(j+1)%polygons[i].size()]],temp);
				n+=temp;
			}
			area+=Length(n);
			normal+=n;
		}
		long long key;

		key=OctNode<NodeData,Real>::CenterIndex(nIdx,maxDepth);
		flatness[key].first=normal;
		flatness[key].second=area;

		OctNode<NodeData,Real>* parent=temp->parent;
		typename OctNode<NodeData,Real>::NodeIndex pIdx=nIdx.parent();
		while(parent)
		{
			key=OctNode<NodeData,Real>::CenterIndex(pIdx,maxDepth);

			if(flatness.find(key)==flatness.end())
			{
				flatness[key].first=normal;
				flatness[key].second=area;
			}
			else
			{
				flatness[key].first+=normal;
				flatness[key].second+=area;
			}
			--pIdx;
			parent=parent->parent;
		}
	}
}
