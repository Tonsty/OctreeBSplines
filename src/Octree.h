/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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
#ifndef OCT_NODE_INCLUDED
#define OCT_NODE_INCLUDED

#include "BinaryNode.h"
#include "MarchingCubes.h"
#include "Geometry.h"

//#include <assert.h>
//#define MAX_OCTNODE_NUM 800000

template<class NodeData,class Real=float>
class OctNode
{
private:
	const OctNode* __faceNeighbor(const int& dir,const int& off) const;
	const OctNode* __edgeNeighbor(const int& o,const int i[2],const int idx[2]) const;
	OctNode* __faceNeighbor(const int& dir,const int& off,const int& forceChildren);
	OctNode* __edgeNeighbor(const int& o,const int i[2],const int idx[2],const int& forceChildren);
public:
	class NodeIndex
	{
	public:
		NodeIndex(void);
		int depth,offset[3];
		NodeIndex child(const int& cIndex) const;
		NodeIndex parent(void) const;
		NodeIndex& operator += (const int& cIndex);
		NodeIndex& operator -- (void);
	};

	OctNode* parent;
	OctNode* children;
	NodeData nodeData;

	OctNode(void);
	~OctNode(void);
	int initChildren(void);
	void deleteChildren(void);

	static inline void CenterAndWidth(const NodeIndex& nIndex,Point3D<Real>& center,Real& width);

	int leaves(void) const;
	//int maxDepthLeaves(const int& maxDepth) const;
	int nodes(void) const;
	int maxDepth(void) const;

	const OctNode* root(void) const;

	const OctNode* nextLeaf(const OctNode* currentLeaf) const;
	OctNode* nextLeaf(OctNode* currentLeaf);
	const OctNode* nextNode(const OctNode* currentNode) const;
	OctNode* nextNode(OctNode* currentNode);
	const OctNode* nextBranch(const OctNode* current) const;
	OctNode* nextBranch(OctNode* current);

	const OctNode* nextLeaf(const OctNode* currentLeaf,NodeIndex &nIndex) const;
	OctNode* nextLeaf(OctNode* currentLeaf,NodeIndex &nIndex);
	const OctNode* nextNode(const OctNode* currentNode,NodeIndex &nIndex) const;
	OctNode* nextNode(OctNode* currentNode,NodeIndex &nIndex);
	const OctNode* nextBranch(const OctNode* current,NodeIndex &nIndex) const;
	OctNode* nextBranch(OctNode* current,NodeIndex &nIndex);

	void setFullDepth(const int& maxDepth);
	void printLeaves(void) const;
	void printRange(void) const;

	OctNode* faceNeighbor(const int& faceIndex,const int& forceChildren=0);
	const OctNode* faceNeighbor(const int& faceIndex) const;
	OctNode* edgeNeighbor(const int& edgeIndex,const int& forceChildren=0);
	const OctNode* edgeNeighbor(const int& edgeIndex) const;
	OctNode* cornerNeighbor(const int& cornerIndex,const int& forceChildren=0);
	const OctNode* cornerNeighbor(const int& cornerIndex) const;


	int write(const char* fileName,int writeData) const;
	int write(FILE* fp,int writeData) const;
	int read(const char* fileName,int readData);
	int read(FILE* fp,int readData);

	static long long CenterIndex(const NodeIndex &nIndex,const int& maxDepth);
	static long long CenterIndex(const NodeIndex &nIndex,const int& maxDepth,int index[3]);
	static long long CornerIndex(const NodeIndex &nIndex,const int& cIndex,const int& maxDepth);
	static long long CornerIndex(const NodeIndex &nIndex,const int& cIndex,const int& maxDepth,int index[3]);
	static long long EdgeIndex(const NodeIndex &nIndex,const int& eIndex,const int& maxDepth);
	static long long EdgeIndex(const NodeIndex &nIndex,const int& eIndex,const int& maxDepth,int idx[3]);
	static long long FaceIndex(const NodeIndex &nIndex,const int& fIndex,const int& maxDepth);
	static long long FaceIndex(const NodeIndex &nIndex,const int& fIndex,const int& maxDepth,int idx[3]);

	class Neighbors{
	public:
		OctNode* neighbors[3][3][3];
		NodeIndex nIndex;
		Neighbors(void);
		void clear(void);
	};
	class NeighborKey
	{
		int _depth;
		Neighbors& _setNeighbors(OctNode* node,const int& d);
		Neighbors& _getNeighbors(OctNode* node,const int& d);
	public:
		Neighbors* neighbors;
		int depth;

		NeighborKey(void);
		~NeighborKey(void);

		void set(const int& depth);
		Neighbors& setNeighbors(OctNode* node);
		Neighbors& getNeighbors(OctNode* node);
	};

	//void* operator new(size_t size)
	//{
	//	assert(freeList.head!=0 && freeList.cnt<MAX_OCTNODE_NUM*sizeof(OctNode));
	//	size_t temp=freeList.cnt;
	//	freeList.cnt+=size;
	//	return freeList.head+temp;
	//}

	//void* operator new[](size_t size)
	//{
	//	assert(freeList.head!=0 && freeList.cnt<MAX_OCTNODE_NUM*sizeof(OctNode));
	//	size_t temp=freeList.cnt;
	//	freeList.cnt+=size;
	//	return freeList.head+temp;
	//}

	//void operator delete[](void*)
	//{
	//	//printf("delete eight children\n");
	//}

	//void operator delete(void*)
	//{
	//	//printf("delete a node\n");
	//}

	//class FreeOctNodeList
	//{
	//public:
	//	FreeOctNodeList():head(0),cnt(0) 
	//	{
	//		head=new char[sizeof(OctNode)*MAX_OCTNODE_NUM];
	//	}
	//	~FreeOctNodeList()
	//	{
	//		delete head;
	//	}
	//	char* head;
	//	size_t cnt;
	//};
	//static FreeOctNodeList freeList;
};

#include "Octree.inl"

#endif // OCT_NODE
