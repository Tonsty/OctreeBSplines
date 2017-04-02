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
#include <stdlib.h>
#include <math.h>
#include <algorithm>

////////////////////////
// OctNode::NodeIndex //
////////////////////////
template<class NodeData,class Real>
OctNode<NodeData,Real>::NodeIndex::NodeIndex(void)
{
	depth=offset[0]=offset[1]=offset[2]=0;
}
template<class NodeData,class Real>
typename OctNode<NodeData,Real>::NodeIndex OctNode<NodeData,Real>::NodeIndex::child(const int& cIndex) const
{
	int x,y,z;
	NodeIndex idx;
	Cube::FactorCornerIndex(cIndex,x,y,z);
	idx.depth=depth+1;
	idx.offset[0]=(offset[0]<<1)|x;
	idx.offset[1]=(offset[1]<<1)|y;
	idx.offset[2]=(offset[2]<<1)|z;
	return idx;
}
template<class NodeData,class Real>
typename OctNode<NodeData,Real>::NodeIndex OctNode<NodeData,Real>::NodeIndex::parent(void) const
{
	NodeIndex idx;
	idx.depth=depth-1;
	idx.offset[0]=offset[0]>>1;
	idx.offset[1]=offset[1]>>1;
	idx.offset[2]=offset[2]>>1;
	return idx;
}
template<class NodeData,class Real>
typename OctNode<NodeData,Real>::NodeIndex& OctNode<NodeData,Real>::NodeIndex::operator += (const int& cIndex)
{
	int x,y,z;
	Cube::FactorCornerIndex(cIndex,x,y,z);
	depth++;
	offset[0]=(offset[0]<<1)|x;
	offset[1]=(offset[1]<<1)|y;
	offset[2]=(offset[2]<<1)|z;
	return *this;
}
template<class NodeData,class Real>
typename OctNode<NodeData,Real>::NodeIndex& OctNode<NodeData,Real>::NodeIndex::operator -- (void)
{
	depth--;
	offset[0]>>=1;
	offset[1]>>=1;
	offset[2]>>=1;
	return *this;
}

/////////////
// OctNode //
/////////////

template <class NodeData,class Real>
OctNode<NodeData,Real>::OctNode(void)
{
	parent=children=NULL;
}

template <class NodeData,class Real>
OctNode<NodeData,Real>::~OctNode(void){
	if(children){delete[] children;}
	parent=children=NULL;
}
template <class NodeData,class Real>
void OctNode<NodeData,Real>::setFullDepth(const int& maxDepth){
	if(maxDepth){
		if(!children){initChildren();}
		for(int i=0;i<8;i++){children[i].setFullDepth(maxDepth-1);}
	}
}
template <class NodeData,class Real>
void OctNode<NodeData,Real>::deleteChildren(void)
{
	if(children)
		delete[] children;
	children=NULL;
}
template <class NodeData,class Real>
int OctNode<NodeData,Real>::initChildren(void){
	int i,j,k;

	if(children){delete[] children;}
	children=NULL;
	children=new OctNode[Cube::CORNERS];
	if(!children){
		fprintf(stderr,"Failed to initialize children in OctNode::initChildren\n");
		exit(0);
		return 0;
	}
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			for(k=0;k<2;k++){
				int idx=Cube::CornerIndex(i,j,k);
				children[idx].parent=this;
				children[idx].children=NULL;
			}
		}
	}
	return 1;
}
template <class NodeData,class Real>
inline void OctNode<NodeData,Real>::CenterAndWidth(const NodeIndex &nIndex,Point3D<Real>& center,Real& width)
{
	width=Real(1.0/(1<<nIndex.depth));
	for(int dim=0;dim<3;dim++)
		center[dim]=Real(0.5+nIndex.offset[dim])*width;
}

template <class NodeData,class Real>
int OctNode<NodeData,Real>::maxDepth(void) const{
	if(!children){return 0;}
	else{
		int c,d;
		for(int i=0;i<Cube::CORNERS;i++){
			d=children[i].maxDepth();
			if(!i || d>c){c=d;}
		}
		return c+1;
	}
}
template <class NodeData,class Real>
int OctNode<NodeData,Real>::nodes(void) const{
	if(!children){return 1;}
	else{
		int c=0;
		for(int i=0;i<Cube::CORNERS;i++){c+=children[i].nodes();}
		return c+1;
	}
}
template <class NodeData,class Real>
int OctNode<NodeData,Real>::leaves(void) const{
	if(!children){return 1;}
	else{
		int c=0;
		for(int i=0;i<Cube::CORNERS;i++){c+=children[i].leaves();}
		return c;
	}
}
// template<class NodeData,class Real>
// int OctNode<NodeData,Real>::maxDepthLeaves(const int& maxDepth) const{
// 	if(depth()>maxDepth){return 0;}
// 	if(!children){return 1;}
// 	else{
// 		int c=0;
// 		for(int i=0;i<Cube::CORNERS;i++){c+=children[i].maxDepthLeaves(maxDepth);}
// 		return c;
// 	}
// }
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::root(void) const{
	const OctNode* temp=this;
	while(temp->parent){temp=temp->parent;}
	return temp;
}


template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextBranch(const OctNode* current) const{
	if(!current->parent || current==this){return NULL;}
	if(current-current->parent->children==Cube::CORNERS-1){return nextBranch(current->parent);}
	else{return current+1;}
}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextBranch(OctNode* current){
	if(!current->parent || current==this){return NULL;}
	if(current-current->parent->children==Cube::CORNERS-1){return nextBranch(current->parent);}
	else{return current+1;}
}
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextLeaf(const OctNode* current) const{
	if(!current){
		const OctNode<NodeData,Real>* temp=this;
		while(temp->children){temp=&temp->children[0];}
		return temp;
	}
	if(current->children){return current->nextLeaf(NULL);}
	const OctNode* temp=nextBranch(current);
	if(!temp){return NULL;}
	else{return temp->nextLeaf(NULL);}
}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextLeaf(OctNode* current){
	if(!current){
		OctNode<NodeData,Real>* temp=this;
		while(temp->children){temp=&temp->children[0];}
		return temp;
	}
	if(current->children){return current->nextLeaf(NULL);}
	OctNode* temp=nextBranch(current);
	if(!temp){return NULL;}
	else{return temp->nextLeaf(NULL);}
}

template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextNode(const OctNode* current) const{
	if(!current){return this;}
	else if(current->children){return &current->children[0];}
	else{return nextBranch(current);}
}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextNode(OctNode* current){
	if(!current){return this;}
	else if(current->children){return &current->children[0];}
	else{return nextBranch(current);}
}
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextBranch(const OctNode* current,NodeIndex& nIndex) const
{
	if(!current->parent || current==this)
		return NULL;
	int c=current-current->parent->children;
	nIndex--;
	if(c==Cube::CORNERS-1)
		return nextBranch(current->parent,nIndex);
	else
	{
		nIndex+=c;
		return current+1;
	}
}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextBranch(OctNode* current,NodeIndex& nIndex)
{
	if(!current->parent || current==this)
		return NULL;
	int c=int(current-current->parent->children);
	--nIndex;
	if(c==Cube::CORNERS-1)
		return nextBranch(current->parent,nIndex);
	else
	{
		nIndex+=c+1;
		return current+1;
	}
}
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextLeaf(const OctNode* current,NodeIndex& nIndex) const
{
	if(!current)
	{
		const OctNode<NodeData,Real>* temp=this;
		while(temp->children)
		{
			nIndex+=0;
			temp=&temp->children[0];
		}
		return temp;
	}
	if(current->children)
		return current->nextLeaf(NULL,nIndex);
	const OctNode* temp=nextBranch(current,nIndex);
	if(!temp)
		return NULL;
	else
		return temp->nextLeaf(NULL,nIndex);
}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextLeaf(OctNode* current,NodeIndex& nIndex)
{
	if(!current){
		OctNode<NodeData,Real>* temp=this;
		while(temp->children)
		{
			nIndex+=0;
			temp=&temp->children[0];
		}
		return temp;
	}
	if(current->children)
		return current->nextLeaf(NULL,nIndex);
	OctNode* temp=nextBranch(current,nIndex);
	if(!temp)
		return NULL;
	else
		return temp->nextLeaf(NULL,nIndex);
}
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextNode(const OctNode* current,NodeIndex& nIndex) const
{
	if(!current)
		return this;
	else if(current->children)
	{
		nIndex+=0;
		return &current->children[0];
	}
	else
		return nextBranch(current,nIndex);
}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::nextNode(OctNode* current,NodeIndex& nIndex)
{
	if(!current)
		return this;
	else if(current->children)
	{
		nIndex+=0;
		return &current->children[0];
	}
	else
		return nextBranch(current,nIndex);
}

template <class NodeData,class Real>
void OctNode<NodeData,Real>::printRange(void) const{
	Point3D<Real> center;
	Real width;
	centerAndWidth(center,width);
	for(int dim=0;dim<3;dim++){
		printf("%[%f,%f]",center[dim]-width/2,center[dim]+width/2);
		if(dim<3-1){printf("x");}
		else printf("\n");
	}
}

template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::faceNeighbor(const int& faceIndex,const int& forceChildren){return __faceNeighbor(faceIndex>>1,faceIndex&1,forceChildren);}
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::faceNeighbor(const int& faceIndex) const {return __faceNeighbor(faceIndex>>1,faceIndex&1);}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::__faceNeighbor(const int& dir,const int& off,const int& forceChildren){
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	pIndex^=(1<<dir);
	if((pIndex & (1<<dir))==(off<<dir)){return &parent->children[pIndex];}
//	if(!(((pIndex>>dir)^off)&1)){return &parent->children[pIndex];}
	else{
		OctNode* temp=parent->__faceNeighbor(dir,off,forceChildren);
		if(!temp){return NULL;}
		if(!temp->children){
			if(forceChildren){temp->initChildren();}
			else{return temp;}
		}
		return &temp->children[pIndex];
	}
}
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::__faceNeighbor(const int& dir,const int& off) const {
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	pIndex^=(1<<dir);
	if((pIndex & (1<<dir))==(off<<dir)){return &parent->children[pIndex];}
//	if(!(((pIndex>>dir)^off)&1)){return &parent->children[pIndex];}
	else{
		const OctNode* temp=parent->__faceNeighbor(dir,off);
		if(!temp || !temp->children){return temp;}
		else{return &temp->children[pIndex];}
	}
}

template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::edgeNeighbor(const int& edgeIndex,const int& forceChildren){
	int idx[2],o,i[2];
	Cube::FactorEdgeIndex(edgeIndex,o,i[0],i[1]);
	switch(o){
		case 0:	idx[0]=1;	idx[1]=2;	break;
		case 1:	idx[0]=0;	idx[1]=2;	break;
		case 2:	idx[0]=0;	idx[1]=1;	break;
	};
	return __edgeNeighbor(o,i,idx,forceChildren);
}
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::edgeNeighbor(const int& edgeIndex) const {
	int idx[2],o,i[2];
	Cube::FactorEdgeIndex(edgeIndex,o,i[0],i[1]);
	switch(o){
		case 0:	idx[0]=1;	idx[1]=2;	break;
		case 1:	idx[0]=0;	idx[1]=2;	break;
		case 2:	idx[0]=0;	idx[1]=1;	break;
	};
	return __edgeNeighbor(o,i,idx);
}
template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::__edgeNeighbor(const int& o,const int i[2],const int idx[2]) const{
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	int aIndex,x[3];

	Cube::FactorCornerIndex(pIndex,x[0],x[1],x[2]);
	aIndex=(~((i[0] ^ x[idx[0]]) | ((i[1] ^ x[idx[1]])<<1))) & 3;
	pIndex^=(7 ^ (1<<o));
	if(aIndex==1)	{	// I can get the neighbor from the parent's face adjacent neighbor
		const OctNode* temp=parent->__faceNeighbor(idx[0],i[0]);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==2)	{	// I can get the neighbor from the parent's face adjacent neighbor
		const OctNode* temp=parent->__faceNeighbor(idx[1],i[1]);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==0)	{	// I can get the neighbor from the parent
		return &parent->children[pIndex];
	}
	else if(aIndex==3)	{	// I can get the neighbor from the parent's edge adjacent neighbor
		const OctNode* temp=parent->__edgeNeighbor(o,i,idx);
		if(!temp || !temp->children){return temp;}
		else{return &temp->children[pIndex];}
	}
	else{return NULL;}
}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::__edgeNeighbor(const int& o,const int i[2],const int idx[2],const int& forceChildren){
	if(!parent){return NULL;}
	int pIndex=int(this-parent->children);
	int aIndex,x[3];

	Cube::FactorCornerIndex(pIndex,x[0],x[1],x[2]);
	aIndex=(~((i[0] ^ x[idx[0]]) | ((i[1] ^ x[idx[1]])<<1))) & 3;
	pIndex^=(7 ^ (1<<o));
	if(aIndex==1)	{	// I can get the neighbor from the parent's face adjacent neighbor
		OctNode* temp=parent->__faceNeighbor(idx[0],i[0],0);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==2)	{	// I can get the neighbor from the parent's face adjacent neighbor
		OctNode* temp=parent->__faceNeighbor(idx[1],i[1],0);
		if(!temp || !temp->children){return NULL;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==0)	{	// I can get the neighbor from the parent
		return &parent->children[pIndex];
	}
	else if(aIndex==3)	{	// I can get the neighbor from the parent's edge adjacent neighbor
		OctNode* temp=parent->__edgeNeighbor(o,i,idx,forceChildren);
		if(!temp){return NULL;}
		if(!temp->children){
			if(forceChildren){temp->initChildren();}
			else{return temp;}
		}
		return &temp->children[pIndex];
	}
	else{return NULL;}
}

template <class NodeData,class Real>
const OctNode<NodeData,Real>* OctNode<NodeData,Real>::cornerNeighbor(const int& cornerIndex) const {
	int pIndex,aIndex=0;
	if(!parent){return NULL;}

	pIndex=int(this-parent->children);
	aIndex=(cornerIndex ^ pIndex);	// The disagreement bits
	pIndex=(~pIndex)&7;				// The antipodal point
	if(aIndex==7){					// Agree on no bits
		return &parent->children[pIndex];
	}
	else if(aIndex==0){				// Agree on all bits
		const OctNode* temp=((const OctNode*)parent)->cornerNeighbor(cornerIndex);
		if(!temp || !temp->children){return temp;}
		else{return &temp->children[pIndex];}
	}
	else if(aIndex==6){				// Agree on face 0
		const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(0,cornerIndex & 1);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==5){				// Agree on face 1
		const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(1,(cornerIndex & 2)>>1);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==3){				// Agree on face 2
		const OctNode* temp=((const OctNode*)parent)->__faceNeighbor(2,(cornerIndex & 4)>>2);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==4){				// Agree on edge 2
		const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(8 | (cornerIndex & 1) | (cornerIndex & 2) );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==2){				// Agree on edge 1
		const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(4 | (cornerIndex & 1) | ((cornerIndex & 4)>>1) );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==1){				// Agree on edge 0
		const OctNode* temp=((const OctNode*)parent)->edgeNeighbor(((cornerIndex & 2) | (cornerIndex & 4))>>1 );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else{return NULL;}
}
template <class NodeData,class Real>
OctNode<NodeData,Real>* OctNode<NodeData,Real>::cornerNeighbor(const int& cornerIndex,const int& forceChildren){
	int pIndex,aIndex=0;
	if(!parent){return NULL;}

	pIndex=int(this-parent->children);
	aIndex=(cornerIndex ^ pIndex);	// The disagreement bits
	pIndex=(~pIndex)&7;				// The antipodal point
	if(aIndex==7){					// Agree on no bits
		return &parent->children[pIndex];
	}
	else if(aIndex==0){				// Agree on all bits
		OctNode* temp=((OctNode*)parent)->cornerNeighbor(cornerIndex,forceChildren);
		if(!temp){return NULL;}
		if(!temp->children){
			if(forceChildren){temp->initChildren();}
			else{return temp;}
		}
		return &temp->children[pIndex];
	}
	else if(aIndex==6){				// Agree on face 0
		OctNode* temp=((OctNode*)parent)->__faceNeighbor(0,cornerIndex & 1,0);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==5){				// Agree on face 1
		OctNode* temp=((OctNode*)parent)->__faceNeighbor(1,(cornerIndex & 2)>>1,0);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==3){				// Agree on face 2
		OctNode* temp=((OctNode*)parent)->__faceNeighbor(2,(cornerIndex & 4)>>2,0);
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==4){				// Agree on edge 2
		OctNode* temp=((OctNode*)parent)->edgeNeighbor(8 | (cornerIndex & 1) | (cornerIndex & 2) );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==2){				// Agree on edge 1
		OctNode* temp=((OctNode*)parent)->edgeNeighbor(4 | (cornerIndex & 1) | ((cornerIndex & 4)>>1) );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else if(aIndex==1){				// Agree on edge 0
		OctNode* temp=((OctNode*)parent)->edgeNeighbor(((cornerIndex & 2) | (cornerIndex & 4))>>1 );
		if(!temp || !temp->children){return NULL;}
		else{return & temp->children[pIndex];}
	}
	else{return NULL;}
}
////////////////////////
// OctNodeNeighborKey //
////////////////////////
template<class NodeData,class Real>
OctNode<NodeData,Real>::Neighbors::Neighbors(void){clear();}
template<class NodeData,class Real>
void OctNode<NodeData,Real>::Neighbors::clear(void){
	for(int i=0;i<3;i++){for(int j=0;j<3;j++){for(int k=0;k<3;k++){neighbors[i][j][k]=NULL;}}}
}
template<class NodeData,class Real>
OctNode<NodeData,Real>::NeighborKey::NeighborKey(void){neighbors=NULL;}
template<class NodeData,class Real>
OctNode<NodeData,Real>::NeighborKey::~NeighborKey(void){
	if(neighbors){delete[] neighbors;}
	neighbors=NULL;
}

template<class NodeData,class Real>
void OctNode<NodeData,Real>::NeighborKey::set(const int& d){
	if(neighbors){delete[] neighbors;}
	neighbors=NULL;
	if(d<0){return;}
	_depth=d;
	neighbors=new Neighbors[d+1];
}
template<class NodeData,class Real>
typename OctNode<NodeData,Real>::Neighbors& OctNode<NodeData,Real>::NeighborKey::setNeighbors(OctNode<NodeData,Real>* node){
	OctNode<NodeData,Real>* temp=node;
	depth=0;
	while(temp->parent)
	{
		depth++;
		temp=temp->parent;
	}
	if(node!=neighbors[depth].neighbors[1][1][1])
		for(int i=depth;i<=_depth;i++)
			neighbors[i].clear();
	
	return _setNeighbors(node,depth);
}
template<class NodeData,class Real>
typename OctNode<NodeData,Real>::Neighbors& OctNode<NodeData,Real>::NeighborKey::_setNeighbors(OctNode<NodeData,Real>* node,const int& d)
{
	if(node!=neighbors[d].neighbors[1][1][1]){
		neighbors[d].clear();

		if(!node->parent)
		{
			neighbors[d].nIndex=OctNode<NodeData,Real>::NodeIndex();
			neighbors[d].neighbors[1][1][1]=node;
		}
		else
		{
			Neighbors& temp=_setNeighbors(node->parent,d-1);
			int i,j,k,x1,y1,z1,x2,y2,z2;
			int idx=int(node-node->parent->children);
			neighbors[d].nIndex=neighbors[d-1].nIndex.child(idx);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
			for(i=0;i<2;i++){
				for(j=0;j<2;j++){
					for(k=0;k<2;k++){
						neighbors[d].neighbors[x2+i][y2+j][z2+k]=&node->parent->children[Cube::CornerIndex(i,j,k)];
					}
				}
			}

			// Set the neighbors from across the faces
			i=x1<<1;
			if(temp.neighbors[i][1][1]){
				if(!temp.neighbors[i][1][1]->children){temp.neighbors[i][1][1]->initChildren();}
				for(j=0;j<2;j++){for(k=0;k<2;k++){neighbors[d].neighbors[i][y2+j][z2+k]=&temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];}}
			}
			j=y1<<1;
			if(temp.neighbors[1][j][1]){
				if(!temp.neighbors[1][j][1]->children){temp.neighbors[1][j][1]->initChildren();}
				for(i=0;i<2;i++){for(k=0;k<2;k++){neighbors[d].neighbors[x2+i][j][z2+k]=&temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];}}
			}
			k=z1<<1;
			if(temp.neighbors[1][1][k]){
				if(!temp.neighbors[1][1][k]->children){temp.neighbors[1][1][k]->initChildren();}
				for(i=0;i<2;i++){for(j=0;j<2;j++){neighbors[d].neighbors[x2+i][y2+j][k]=&temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];}}
			}

			// Set the neighbors from across the edges
			i=x1<<1;	j=y1<<1;
			if(temp.neighbors[i][j][1]){
				if(!temp.neighbors[i][j][1]->children){temp.neighbors[i][j][1]->initChildren();}
				for(k=0;k<2;k++){neighbors[d].neighbors[i][j][z2+k]=&temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];}
			}
			i=x1<<1;	k=z1<<1;
			if(temp.neighbors[i][1][k]){
				if(!temp.neighbors[i][1][k]->children){temp.neighbors[i][1][k]->initChildren();}
				for(j=0;j<2;j++){neighbors[d].neighbors[i][y2+j][k]=&temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];}
			}
			j=y1<<1;	k=z1<<1;
			if(temp.neighbors[1][j][k]){
				if(!temp.neighbors[1][j][k]->children){temp.neighbors[1][j][k]->initChildren();}
				for(i=0;i<2;i++){neighbors[d].neighbors[x2+i][j][k]=&temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];}
			}

			// Set the neighbor from across the corner
			i=x1<<1;	j=y1<<1;	k=z1<<1;
			if(temp.neighbors[i][j][k]){
				if(!temp.neighbors[i][j][k]->children){temp.neighbors[i][j][k]->initChildren();}
				neighbors[d].neighbors[i][j][k]=&temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
			}
		}
	}
	return neighbors[d];
}
template<class NodeData,class Real>
typename OctNode<NodeData,Real>::Neighbors& OctNode<NodeData,Real>::NeighborKey::getNeighbors(OctNode<NodeData,Real>* node){
	depth=0;
	OctNode<NodeData,Real>* temp=node;
	while(temp->parent)
	{
		depth++;
		temp=temp->parent;
	}
	if(node!=neighbors[depth].neighbors[1][1][1])
		for(int i=depth;i<=_depth;i++)
			neighbors[i].clear();
	
	return _getNeighbors(node,depth);
}
template<class NodeData,class Real>
typename OctNode<NodeData,Real>::Neighbors& OctNode<NodeData,Real>::NeighborKey::_getNeighbors(OctNode<NodeData,Real>* node,const int& d)
{
	if(node!=neighbors[d].neighbors[1][1][1]){
		neighbors[d].clear();

		if(!node->parent)
		{
			neighbors[d].neighbors[1][1][1]=node;
			neighbors[d].nIndex=OctNode<NodeData,Real>::NodeIndex();
		}
		else
		{
			Neighbors& temp=_getNeighbors(node->parent,d-1);
			int i,j,k,x1,y1,z1,x2,y2,z2;
			int idx=int(node-node->parent->children);
			neighbors[d].nIndex=neighbors[d-1].nIndex.child(idx);
			Cube::FactorCornerIndex(  idx   ,x1,y1,z1);
			Cube::FactorCornerIndex((~idx)&7,x2,y2,z2);
			for(i=0;i<2;i++){
				for(j=0;j<2;j++){
					for(k=0;k<2;k++){
						neighbors[d].neighbors[x2+i][y2+j][z2+k]=&node->parent->children[Cube::CornerIndex(i,j,k)];
					}
				}
			}

			// Set the neighbors from across the faces
			i=x1<<1;
			if(temp.neighbors[i][1][1] && temp.neighbors[i][1][1]->children){
				for(j=0;j<2;j++){for(k=0;k<2;k++){neighbors[d].neighbors[i][y2+j][z2+k]=&temp.neighbors[i][1][1]->children[Cube::CornerIndex(x2,j,k)];}}
			}
			j=y1<<1;
			if(temp.neighbors[1][j][1] && temp.neighbors[1][j][1]->children){
				for(i=0;i<2;i++){for(k=0;k<2;k++){neighbors[d].neighbors[x2+i][j][z2+k]=&temp.neighbors[1][j][1]->children[Cube::CornerIndex(i,y2,k)];}}
			}
			k=z1<<1;
			if(temp.neighbors[1][1][k] && temp.neighbors[1][1][k]->children){
				for(i=0;i<2;i++){for(j=0;j<2;j++){neighbors[d].neighbors[x2+i][y2+j][k]=&temp.neighbors[1][1][k]->children[Cube::CornerIndex(i,j,z2)];}}
			}

			// Set the neighbors from across the edges
			i=x1<<1;	j=y1<<1;
			if(temp.neighbors[i][j][1] && temp.neighbors[i][j][1]->children){
				for(k=0;k<2;k++){neighbors[d].neighbors[i][j][z2+k]=&temp.neighbors[i][j][1]->children[Cube::CornerIndex(x2,y2,k)];}
			}
			i=x1<<1;	k=z1<<1;
			if(temp.neighbors[i][1][k] && temp.neighbors[i][1][k]->children){
				for(j=0;j<2;j++){neighbors[d].neighbors[i][y2+j][k]=&temp.neighbors[i][1][k]->children[Cube::CornerIndex(x2,j,z2)];}
			}
			j=y1<<1;	k=z1<<1;
			if(temp.neighbors[1][j][k] && temp.neighbors[1][j][k]->children){
				for(i=0;i<2;i++){neighbors[d].neighbors[x2+i][j][k]=&temp.neighbors[1][j][k]->children[Cube::CornerIndex(i,y2,z2)];}
			}

			// Set the neighbor from across the corner
			i=x1<<1;	j=y1<<1;	k=z1<<1;
			if(temp.neighbors[i][j][k] && temp.neighbors[i][j][k]->children){
				neighbors[d].neighbors[i][j][k]=&temp.neighbors[i][j][k]->children[Cube::CornerIndex(x2,y2,z2)];
			}
		}
	}
	return neighbors[d];
}

template <class NodeData,class Real>
int OctNode<NodeData,Real>::write(const char* fileName,int writeData) const{
	FILE* fp=fopen(fileName,"wb");
	if(!fp){return 0;}
	int ret=write(fp,writeData);
	fclose(fp);
	return ret;
}
template <class NodeData,class Real>
int OctNode<NodeData,Real>::write(FILE* fp,int writeData) const
{
	if(writeData)
		fwrite(this,sizeof(OctNode<NodeData,Real>),1,fp);
	else
	{
		fwrite(&parent,sizeof(OctNode<NodeData,Real>*),1,fp);
		fwrite(&children,sizeof(OctNode<NodeData,Real>*),1,fp);
	}
	if(children){for(int i=0;i<Cube::CORNERS;i++){children[i].write(fp,writeData);}}
	return 1;
}
template <class NodeData,class Real>
int OctNode<NodeData,Real>::read(const char* fileName,int readData){
	FILE* fp=fopen(fileName,"rb");
	if(!fp){return 0;}
	int ret=read(fp,readData);
	fclose(fp);
	return ret;
}
template <class NodeData,class Real>
int OctNode<NodeData,Real>::read(FILE* fp,int readData){
	if(readData)
		fread(this,sizeof(OctNode<NodeData,Real>),1,fp);
	else
	{
		fread(&parent,sizeof(OctNode<NodeData,Real>*),1,fp);
		fread(&children,sizeof(OctNode<NodeData,Real>*),1,fp);
	}
	parent=NULL;
	if(children)
	{
		children=NULL;
		initChildren();
		for(int i=0;i<Cube::CORNERS;i++){
			children[i].read(fp,readData);
			children[i].parent=this;
		}
	}
	return 1;
}
////////////////
// VertexData //
////////////////
template<class NodeData,class Real>
long long OctNode<NodeData,Real>::CornerIndex(const NodeIndex& nIndex,const int& cIndex,const int& maxDepth)
{
	int idx[3];
	return CornerIndex(nIndex,cIndex,maxDepth,idx);
}
template<class NodeData,class Real>
long long OctNode<NodeData,Real>::CornerIndex(const NodeIndex& nIndex,const int& cIndex,const int& maxDepth,int idx[3]){
	int x[3];
	Cube::FactorCornerIndex(cIndex,x[0],x[1],x[2]);
	for(int i=0;i<3;i++)
		idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[i],x[i]);
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
template<class NodeData,class Real>
long long OctNode<NodeData,Real>::CenterIndex(const NodeIndex& nIndex,const int& maxDepth)
{
	int idx[3];
	return CenterIndex(nIndex,maxDepth,idx);
}

template<class NodeData,class Real>
long long OctNode<NodeData,Real>::CenterIndex(const NodeIndex& nIndex,const int& maxDepth,int idx[3])
{
	for(int i=0;i<3;i++)
		idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth+1,nIndex.offset[i]<<1,1);
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
template<class NodeData,class Real>
long long OctNode<NodeData,Real>::EdgeIndex(const NodeIndex& nIndex,const int& eIndex,const int& maxDepth)
{
	int idx[3];
	return EdgeIndex(nIndex,eIndex,maxDepth,idx);
}
template<class NodeData,class Real>
long long OctNode<NodeData,Real>::EdgeIndex(const NodeIndex& nIndex,const int& eIndex,const int& maxDepth,int idx[3])
{
	int o,i1,i2;
	for(int i=0;i<3;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth+1,nIndex.offset[i]<<1,1);}
	Cube::FactorEdgeIndex(eIndex,o,i1,i2);
	switch(o){
		case 0:
			idx[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[1],i1);
			idx[2]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[2],i2);
			break;
		case 1:
			idx[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[0],i1);
			idx[2]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[2],i2);
			break;
		case 2:
			idx[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[0],i1);
			idx[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[1],i2);
			break;
	};
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
template<class NodeData,class Real>
long long OctNode<NodeData,Real>::FaceIndex(const NodeIndex& nIndex,const int& fIndex,const int& maxDepth)
{
	int idx[3];
	return FaceIndex(nIndex,fIndex,maxDepth,idx);
}
template<class NodeData,class Real>
long long OctNode<NodeData,Real>::FaceIndex(const NodeIndex& nIndex,const int& fIndex,const int& maxDepth,int idx[3])
{
	int dir,offset;
	Cube::FactorFaceIndex(fIndex,dir,offset);
	for(int i=0;i<3;i++)
		idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth+1,nIndex.offset[i]<<1,1);
	idx[dir]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[dir],offset);
	return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
