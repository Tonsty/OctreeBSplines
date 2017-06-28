#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <omp.h>
#include "Kdtree.h"
#include "vtkHelper.h"

template<class NodeData,class Real,class VertexData>
int OctreeBspline<NodeData,Real,VertexData>::set4(const int& maxDepth,Point3D<Real>& translate,Real& scale)
{
	this->maxDepth=maxDepth;

	OctNode<NodeData,Real>* temp;
	OctNode<NodeData,Real>::NodeIndex nIndex;

	compressedKeys.resize(maxBsplineDepth+1);
	for(int depth=maxDepth-1;depth<=maxDepth-1;depth++)
	{
		tree.setFullDepth(depth);
		temp=tree.nextLeaf(NULL,nIndex);
		while(temp)
		{
			int bsplineOffset[3];
			int bsplineDepth=nIndex.depth+1;
			bsplineOffset[0]=nIndex.offset[0]<<1;
			bsplineOffset[1]=nIndex.offset[1]<<1;
			bsplineOffset[2]=nIndex.offset[2]<<1;
			long long compressedKey=(long long)(bsplineOffset[0]) | (long long)(bsplineOffset[1])<<15 | (long long)(bsplineOffset[2])<<30;
			compressedKeys[bsplineDepth].push_back(compressedKey);
			temp=tree.nextLeaf(temp,nIndex);
		}
	}

	cornerValues.clear();

	tree.setFullDepth(maxDepth);
	temp=tree.nextLeaf(NULL,nIndex);
	while(temp)
	{
		for(int c=0;c<Cube::CORNERS;c++)
			cornerValues[OctNode<NodeData,Real>::CornerIndex(nIndex,c,maxDepth)]=VertexData(); //allocate hash-table entry
		temp=tree.nextLeaf(temp,nIndex);
	}

	translate[0]=translate[1]=translate[2]=0.5;
	scale=1.0;

	float unitLen=(float)1.0/(1<<(maxDepth+1));
	Point3D<Real> center;
	center[0]=(Real)0.5;
	center[1]=(Real)0.5;
	center[2]=(Real)0.5;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++)
	{
		long long cornerKey=it->first;
		Point3D<Real> pos;
		getPosFromCornerKey(cornerKey,unitLen,pos.coords);
		
		Real dist=Distance(pos,center)-(Real)0.5;
		Point3D<Real> n=(pos-center)/Distance(pos,center);
		Real w=(Real)1.0;
		it->second=VertexData(dist,n,w);
	}

	return 1;
}

template<class NodeData,class Real,class VertexData>
template<class Vertex>
int OctreeBspline<NodeData,Real,VertexData>::set3(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,
	const int& maxDepth,const int& setCenter,const Real& flatness,const Real& curvature,const Real& splat,const int& maxDepthTree,
	Point3D<Real>& translate,Real& scale,const int& noTransform,const int&noFit)
{
	this->maxDepth=maxDepth;
	OctNode<NodeData,Real>::NodeIndex nIdx;

	MeshInfo<float> &mInfo=mInfoGlobal;
	mInfo.set2(vertices,polygons,Real(1.1),translate,scale,noTransform);

	std::vector<int> myVertices;
	for(int i=0;i<mInfo.vertices.size();i++) myVertices.push_back(i);

	cornerValues.clear();
	for(int c=0;c<Cube::CORNERS;c++)
		cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,c,maxDepth)]=VertexData(); //allocate hash-table entry

	printf("Allocating cornerValues ...\n");
	compressedKeys.resize(maxBsplineDepth+1);
	setChildren3(&tree,nIdx,myVertices,mInfo,maxDepth,setCenter,flatness,curvature,splat,maxDepthTree,NULL,!noFit);
	printf("Finished allocating\n");

	KDTree<float> kdtree;
	kdtree.setInputPoints((float*)mInfo.vertices.data(),mInfo.vertices.size());
	float unitLen=(float)1.0/(1<<(maxDepth+1));
	std::vector<Point3D<float> > queries;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++)
	{
		long long cornerKey=it->first;
		Point3D<float> point;
		getPosFromCornerKey(cornerKey,unitLen,point.coords);
		queries.push_back(point);
	}

	std::vector<std::vector<int> > indices; 
	std::vector<std::vector<float> > dists;
	kdtree.KnnSearch((float*)queries.data(),queries.size(),indices,dists,1);

	Real dist,w;
	Point3D<Real> n;
	int i=0;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++,i++)
	{
		Point3D<Real> point=queries[i];
		Point3D<Real> point2;
		for(int j=0;j<3;j++) point2[j]=(Real)mInfo.vertices[indices[i][0]][j];
		Point3D<Real> normal2;
		for(int j=0;j<3;j++) normal2[j]=(Real)mInfo.vertexNormals[indices[i][0]][j];
		setDistanceAndNormal3(point,point2,normal2,dist,n,w);
		it->second=VertexData(dist,n,w);
	}

	return 1;
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::setDistanceAndNormal3(const Point3D<Real>& p,const Point3D<Real>& p2,const Point3D<Real>& n2,Real& dist,Point3D<Real>& n,Real& w)
{
	dist=Distance(p,p2);
	n=(p-p2)/dist;

	w=DotProduct(n,n2);

	//only points locating inside narrow band are computed for point to plane distance
	Real narrow_band=(Real)0.5;
	if (dist<narrow_band) 
	{
		//point to plane distance
		dist=DotProduct(p-p2,n2); 
		//point to plane gradient
		n=n2; 
	}
	else 
	{
		if(DotProduct(p-p2,n2)<0)
		{
			dist=-dist;
			n*=-1;
		}
	}
}

template<class NodeData,class Real,class VertexData>
template<class MeshReal>
int OctreeBspline<NodeData,Real,VertexData>::setChildren3(OctNode<NodeData,Real>* node,
	const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
	std::vector<int>& vindices,MeshInfo<MeshReal>& mInfo,
	const int& maxDepth,const int& setCenter,const Real& flatness,
	const Real&curvature,const Real& splat,const int& maxDepthTree,
	stdext::hash_map<long long,std::vector<int>*>* triangleMap,
	int bFlag)
{
	long long key;
	Real w;
	Point3D<Real> ctr;
	if(vindices.size()<=1)		return -1;
	if(nIdx.depth==maxDepth)	return -1;

	OctNode<NodeData,Real>::CenterAndWidth(nIdx,ctr,w);

	if(flatness>0)
	{
		MeshReal area;
		Point3D<MeshReal> mc;
		area=0;
		mc[0]=mc[1]=mc[2]=0;
		for(size_t i=0;i<vindices.size();i++)
		{
			area+=Length(mInfo.vertexNormals[vindices[i]]);
			mc+=mInfo.vertexNormals[vindices[i]];
		}
		if(Length(mc)/area>flatness) return -1;
	}

	if(bFlag)
	{
		if(curvature>0)
		{
			std::vector<int> vindicesTemp;
			for(size_t i=0;i<vindices.size();i++)
			{
				if(curvature/mInfo.vertexCurvatures[vindices[i]]<w)
					vindicesTemp.push_back(vindices[i]);
			}
			if(!vindicesTemp.size()) 
			{
				if(maxDepthTree) bFlag=0;
				else return -1;
			}
			//vindices.swap(vindicesTemp);
		}
	}

	if(!node->children)	node->initChildren();

	// Set the center
	key=OctNode<NodeData,Real>::CenterIndex(nIdx,maxDepth);
	cornerValues[key]=VertexData();

	// Set the edge mid-points
	for(int i=0;i<Cube::EDGES;i++)
	{
		key=OctNode<NodeData,Real>::EdgeIndex(nIdx,i,maxDepth);
		cornerValues[key]=VertexData();
	}

	// set the face mid-points
	for(int i=0;i<Cube::FACES;i++)
	{
		key=OctNode<NodeData,Real>::FaceIndex(nIdx,i,maxDepth);
		cornerValues[key]=VertexData();
	}

	int retCount=0;
	for(int i=0;i<Cube::CORNERS;i++)
	{
		OctNode<NodeData,Real>::CenterAndWidth(nIdx.child(i),ctr,w);
		std::vector<int> myVertices;
		for(int j=0;j<vindices.size();j++)
		{
			Point3D<MeshReal> t=mInfo.vertices[vindices[j]];
			Point3D<MeshReal> ctr2;
			MeshReal w2=w;
			ctr2[0]=ctr[0];
			ctr2[1]=ctr[1];
			ctr2[2]=ctr[2];
			if(PointInCube(ctr2,w2*(Real)splat,t)) myVertices.push_back(vindices[j]);
		}
		if(myVertices.size())
		{
			if(setChildren3(&node->children[i],nIdx.child(i),myVertices,mInfo,maxDepth,setCenter,flatness,curvature,splat,maxDepthTree,triangleMap,bFlag)>0) retCount++;
		}
	}

	//std::vector<int> myVerticess[Cube::CORNERS];
	//for(size_t i=0;i<vindices.size();i++)
	//{
	//	Point3D<MeshReal> t=mInfo.vertices[vindices[i]];
	//	int offset[3];
	//	offset[0]=t[0]<ctr[0]?0:1;
	//	offset[1]=t[1]<ctr[1]?0:1;
	//	offset[2]=t[2]<ctr[2]?0:1;
	//	int cIndex= (offset[0]<<0) | (offset[1]<<1) | (offset[2]<<2);
	//	myVerticess[cIndex].push_back(vindices[i]);
	//}
	//int retCount=0;
	//for(int i=0;i<Cube::CORNERS;i++)
	//	if(myVerticess[i].size() && setChildren3(&node->children[i],nIdx.child(i),myVerticess[i],mInfo,maxDepth,setCenter,flatness,triangleMap)>0) retCount++;

	//if(maxBsplineDepth>0 && ((nIdx.depth==maxBsplineDepth-1) || (nIdx.depth<maxBsplineDepth-1 && retCount<8)))
	if(bFlag && maxBsplineDepth>0 && (nIdx.depth<=maxBsplineDepth-1))
	{
		int bsplineOffset[3];
		int bsplineDepth=nIdx.depth+1;
		bsplineOffset[0]=nIdx.offset[0]<<1;
		bsplineOffset[1]=nIdx.offset[1]<<1;
		bsplineOffset[2]=nIdx.offset[2]<<1;

		long long compressedKey=(long long)(bsplineOffset[0]) | (long long)(bsplineOffset[1])<<15 | (long long)(bsplineOffset[2])<<30;
		compressedKeys[bsplineDepth].push_back(compressedKey);
	}
	return 1;
}

template<class NodeData,class Real,class VertexData>
template<class Vertex>
int OctreeBspline<NodeData,Real,VertexData>::set2(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,
	const int& maxDepth,const int& setCenter,const Real& flatness,
	Point3D<Real>& translate,Real& scale,const int& noTransform)
{
	this->maxDepth=maxDepth;
	OctNode<NodeData,Real>::NodeIndex nIdx;

	MeshInfo<float> &mInfo=mInfoGlobal;
	std::vector<int> myVertices;

	mInfo.set2(vertices,polygons,Real(1.1),translate,scale,noTransform);
	myVertices.resize(mInfo.vertices.size());
	for(int i=0;i<int(mInfo.vertices.size());i++) myVertices[i]=i;

	KDTree<float> kdtree;
	kdtree.setInputPoints(mInfo.vertices);
	std::vector<Point3D<Real> > queries;
	queries.push_back(RandomBallPoint<Real>());
	std::vector<std::vector<int> > indices; 
	std::vector<std::vector<Real> > dists;
	kdtree.KnnSearch(queries,indices,dists,1);

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

		setDistanceAndNormal2(myVertices,mInfo,p,dist,n);
		cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,c,maxDepth)]=VertexData(dist,n);
	}
	if(setCenter)
	{
		Real w;
		OctNode<NodeData,Real>::CenterAndWidth(nIdx,tree.nodeData.center,w);
		setDistanceAndNormal2(myVertices,mInfo,tree.nodeData.center,dist,n);
		tree.nodeData.v=VertexData(dist,n);
	}
	compressedKeys.resize(maxBsplineDepth+1);
	setChildren2(&tree,nIdx,myVertices,mInfo,maxDepth,setCenter,flatness);
	return 1;
}

template<class NodeData,class Real,class VertexData>
template<class MeshReal>
void OctreeBspline<NodeData,Real,VertexData>::setDistanceAndNormal2(const std::vector<int>& vindices,
	MeshInfo<MeshReal>& mInfo,
	const Point3D<Real>& p,
	Real& dist,
	Point3D<Real>& n)
{
	Point3D<MeshReal> pp,t;
	pp[0]=p[0];
	pp[1]=p[1];
	pp[2]=p[2];
	size_t closest;
	MeshReal minSqDist;
	for(size_t i=0;i<vindices.size();i++) 
	{
		MeshReal temp;
		for(int j=0;j<3;j++) t[j]=mInfo.vertices[vindices[i]][j];
		temp=SquareDistance(pp,t); 
		if(!i || temp<minSqDist ) 
		{
			closest=i;
			minSqDist=temp;
		}
	}
	Point3D<MeshReal> nn;
	closest=vindices[closest];

	for(int i=0;i<3;i++) 
	{
		t[i]=mInfo.vertices[closest][i];
		nn[i]=mInfo.vertexNormals[closest][i];
	}

	Real w;
	setDistanceAndNormal3(pp,t,nn,dist,n,w);
}

template<class NodeData,class Real,class VertexData>
template<class MeshReal>
int OctreeBspline<NodeData,Real,VertexData>::setChildren2(OctNode<NodeData,Real>* node,
	const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
	const std::vector<int>& vindices,MeshInfo<MeshReal>& mInfo,
	const int& maxDepth,const int& setCenter,
	const Real& flatness,
	stdext::hash_map<long long,std::vector<int>*>* triangleMap)
{
	long long key;
	Real w,dist;
	Point3D<Real> ctr,n,p;
	if(!vindices.size())		return -1;
	if(nIdx.depth==maxDepth)	return -1;

	if(flatness>0)
	{
		MeshReal area;
		Point3D<MeshReal> mc;
		area=0;
		mc[0]=mc[1]=mc[2]=0;
		for(size_t i=0;i<vindices.size();i++)
		{
			area+=Length(mInfo.vertexNormals[vindices[i]]);
			mc+=mInfo.vertexNormals[vindices[i]];
		}
		if(Length(mc)/area>flatness) return -1;
	}

	if(!node->children)	node->initChildren();
	OctNode<NodeData,Real>::CenterAndWidth(nIdx,ctr,w);

	// Set the center
	setDistanceAndNormal2(vindices,mInfo,ctr,dist,n);
	key=OctNode<NodeData,Real>::CenterIndex(nIdx,maxDepth);
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
		key=OctNode<NodeData,Real>::EdgeIndex(nIdx,i,maxDepth);
		auto it=cornerValues.find(key);
		if(it==cornerValues.end())
		{
			setDistanceAndNormal2(vindices,mInfo,p,dist,n);
			cornerValues[key]=VertexData(dist,n);
		}
		//setDistanceAndNormal2(vindices,mInfo,p,dist,n);
		//key=OctNode<NodeData,Real>::EdgeIndex(nIdx,i,maxDepth);
		//auto it=cornerValues.find(key);
		//if(it==cornerValues.end() || fabs(it->second.value())>fabs(dist))
		//	cornerValues[key]=VertexData(dist,n);
	}

	// set the face mid-points
	for(int i=0;i<Cube::FACES;i++)
	{
		int dir,off;
		Cube::FactorFaceIndex(i,dir,off);
		p=ctr;
		p[dir]+=-w/2+w*off;
		key=OctNode<NodeData,Real>::FaceIndex(nIdx,i,maxDepth);
		auto it=cornerValues.find(key);
		if(it==cornerValues.end())
		{
			setDistanceAndNormal2(vindices,mInfo,p,dist,n);
			cornerValues[key]=VertexData(dist,n);
		}
		//setDistanceAndNormal2(vindices,mInfo,p,dist,n);
		//key=OctNode<NodeData,Real>::FaceIndex(nIdx,i,maxDepth);
		//auto it=cornerValues.find(key);
		//if(it==cornerValues.end() || fabs(it->second.value())>fabs(dist))
		//	cornerValues[key]=VertexData(dist,n);
	}

	int retCount=0;
	for(int i=0;i<Cube::CORNERS;i++)
	{
		OctNode<NodeData,Real>::CenterAndWidth(nIdx.child(i),ctr,w);
		if(setCenter)
		{
			OctNode<NodeData,Real>::CenterAndWidth(nIdx.child(i),node->children[i].nodeData.center,w);
			setDistanceAndNormal2(vindices,mInfo,node->children[i].nodeData.center,dist,n);
			node->children[i].nodeData.v=VertexData(dist,n);
		}

		std::vector<int> myVerticesExp;
		std::vector<int> myVertices;
		for(size_t j=0;j<vindices.size();j++)
		{
			Point3D<MeshReal> t=mInfo.vertices[vindices[j]];
			Point3D<MeshReal> ctr2;
			MeshReal w2=w;
			ctr2[0]=ctr[0];
			ctr2[1]=ctr[1];
			ctr2[2]=ctr[2];
			if(PointInCube(ctr2,w2*4,t))
			{
				myVerticesExp.push_back(vindices[j]);
				if(PointInCube(ctr2,w2,t)) myVertices.push_back(vindices[j]);
			}
		}
		if(myVertices.size() && setChildren2(&node->children[i],nIdx.child(i),myVerticesExp,mInfo,maxDepth,setCenter,flatness,triangleMap)>0) retCount++;
	}
	//if(maxBsplineDepth>0 && ((nIdx.depth==maxBsplineDepth-1) || (nIdx.depth<maxBsplineDepth-1 && retCount<8)))
	if(maxBsplineDepth>0 && (nIdx.depth<=maxBsplineDepth-1))
	{
		int bsplineOffset[3];
		int bsplineDepth=nIdx.depth+1;
		bsplineOffset[0]=nIdx.offset[0]<<1;
		bsplineOffset[1]=nIdx.offset[1]<<1;
		bsplineOffset[2]=nIdx.offset[2]<<1;

		long long compressedKey=(long long)(bsplineOffset[0]) | (long long)(bsplineOffset[1])<<15 | (long long)(bsplineOffset[2])<<30;
		compressedKeys[bsplineDepth].push_back(compressedKey);
	}
	return 1;
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getBposIntU(const int bsplineDepth,const float pos[3],int bposInt[3],float u[3])
{
	float delta=(float)1.0/(1<<bsplineDepth);
	float bposReal[3]={pos[0]/delta+(float)1.5,pos[1]/delta+(float)1.5,pos[2]/delta+(float)1.5};
	bposInt[0]=(int)bposReal[0]; 
	bposInt[1]=(int)bposReal[1]; 
	bposInt[2]=(int)bposReal[2];
	u[0]=bposReal[0]-bposInt[0];
	u[1]=bposReal[1]-bposInt[1];
	u[2]=bposReal[2]-bposInt[2];
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getPosFromCornerKey(const long long cornerKey,const float unitLen,float pos[3])
{
	int idx[3];
	idx[0]=cornerKey&0x7fff;
	idx[1]=(cornerKey>>15)&0x7fff;
	idx[2]=(cornerKey>>30)&0x7fff;
	pos[0]=unitLen*idx[0];
	pos[1]=unitLen*idx[1];
	pos[2]=unitLen*idx[2];
}

template<class NodeData,class Real,class VertexData>
float OctreeBspline<NodeData,Real,VertexData>::eval(const float pos[3])
{
	Eigen::Matrix<float,3,3> B;
	B<<1,-2,1,1,2,-2,0,0,1;
	B/=2;

	std::vector<int> coeffIndices;
	std::vector<float> coeffWeights;
	float value;
	getCoeffIndicesWeightsValueFromPos(B,1,maxBsplineDepth,pos,coeffIndices,coeffWeights,value);
	return value;
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getCoeffIndicesWeightsValueFromPos(
	const Eigen::Matrix<float,3,3>& B,const int minBsplineDepth,const int maxBsplineDepth,const float pos[3], 
	std::vector<int>& coeffIndices,std::vector<float>& coeffWeights,float& value)
{
	value=0;
	coeffIndices.clear();
	coeffWeights.clear();
	if(maxBsplineDepth<minBsplineDepth || minBsplineDepth<1) return;
	for(int bsplineDepth=minBsplineDepth;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
	{
		if(!coeffValues[bsplineDepth].size()) continue;

		int bposInt[3];
		float u[3];
		getBposIntU(bsplineDepth,pos,bposInt,u);

		Eigen::Matrix<float,3,3> u_mat;
		u_mat<<1,1,1,u[0],u[1],u[2],u[0]*u[0],u[1]*u[1],u[2]*u[2];
		Eigen::Matrix<float,3,3> Bu=B*u_mat;

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				for(int k=0;k<3;k++)
				{
					int bsplineOffset[3];
					bsplineOffset[0]=bposInt[0]-2+i;
					bsplineOffset[1]=bposInt[1]-2+j;
					bsplineOffset[2]=bposInt[2]-2+k;
					if(bsplineOffset[0]>=0 && bsplineOffset[1]>=0 && bsplineOffset[2]>=0)
					{
						long long key=(long long)(bsplineOffset[0]) | (long long)(bsplineOffset[1])<<15 | (long long)(bsplineOffset[2])<<30;
						auto it2=coeffValues[bsplineDepth].find(key);
						if(it2!=coeffValues[bsplineDepth].end())
						{
							int coeffIndex=it2->second.first;
							float coeffValue=it2->second.second;
							float coeffWeight=Bu(i,0)*Bu(j,1)*Bu(k,2);
							coeffIndices.push_back(coeffIndex);
							coeffWeights.push_back(coeffWeight);
							value+=coeffValue*coeffWeight;
						}
					}
				}
	}
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::setCoeffValuesFromCompressedKeys(
	const int bsplineDepth,const std::vector<std::vector<long long> >& compressedKeys,
	std::vector<stdext::hash_map<long long,std::pair<int,Real> > >& coeffValues)
{
	for(auto it=compressedKeys[bsplineDepth].begin();it!=compressedKeys[bsplineDepth].end();it++)
	{
		long long compressedKey=*it;
		int bsplineOffset[3];
		bsplineOffset[0]=compressedKey&0x7fff;
		bsplineOffset[1]=(compressedKey>>15)&0x7fff;
		bsplineOffset[2]=(compressedKey>>30)&0x7fff;
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				for(int k=0;k<3;k++)
				{
					long long coeffKey=(long long)(bsplineOffset[0]+i) | (long long)(bsplineOffset[1]+j)<<15 | (long long)(bsplineOffset[2]+k)<<30;
					std::pair<int,Real>& coeffValue=coeffValues[bsplineDepth][coeffKey]; //allocate hash-table entry
					coeffValue.first=(int)0;
					coeffValue.second=(Real)0.0;
				}
	}
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::directBsplineFitting(const Real &smooth,const Real &interpolate)
{
	int maxDepth=this->maxDepth;
	int maxBsplineDepth=this->maxBsplineDepth;

	coeffValues.resize(maxBsplineDepth+1);
	int r=0;
	for(int bsplineDepth=1;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
	{
		setCoeffValuesFromCompressedKeys(bsplineDepth,compressedKeys,coeffValues);
		r+=(int)compressedKeys[bsplineDepth].size();
		std::cout<<"Num. of compressedKeys at bsplineDepth "<<bsplineDepth<<" is "<<compressedKeys[bsplineDepth].size()<<std::endl;
	}
	std::cout<<"Num. of compressedKeys is "<<r<<std::endl;

	int p=0;
	for(int bsplineDepth=1;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
	{
		for(auto it=coeffValues[bsplineDepth].begin();it!=coeffValues[bsplineDepth].end();it++,p++)it->second.first=p;
		std::cout<<"Num. of coeffValues at bsplineDepth "<<bsplineDepth<<" is "<<coeffValues[bsplineDepth].size()<<std::endl;
	}
	std::cout<<"Num. of coeffValues is "<<p<<std::endl;

	Vector b(cornerValues.size());
	std::vector<Triplet> A_triplets;

	Eigen::Matrix<float,3,3> B;
	B<<1,-2,1,1,2,-2,0,0,1;
	B/=2;
	float unitLen=(float)1.0/(1<<(maxDepth+1));
	int q=0;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++,q++)
	{
		b(q)=it->second.v*it->second.w;
		long long cornerKey=it->first;
		float pos[3];
		getPosFromCornerKey(cornerKey,unitLen,pos);
		//std::cout<<q<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
		std::vector<int> coeffIndices;
		std::vector<float> coeffWeights;
		float value;
		getCoeffIndicesWeightsValueFromPos(B,1,maxBsplineDepth,pos,coeffIndices,coeffWeights,value);
		for(int i=0;i<coeffIndices.size();i++) A_triplets.push_back(Triplet(q,coeffIndices[i],coeffWeights[i]*it->second.w));
	}
	std::cout<<"Num. of cornerValues is "<<q<<std::endl;

	SparseMatrix A(q, p);
	A.setFromTriplets(A_triplets.begin(),A_triplets.end());

	//Eigen::ConjugateGradient<SparseMatrix> cg;
	//cg.compute(A.transpose()*A);
	//x=cg.solve(b);
	Eigen::LeastSquaresConjugateGradient<SparseMatrix> lscg;
	lscg.compute(A);
	Vector x=lscg.solve(b);
	std::cout<<"finised solving"<<std::endl;

	for(int bsplineDepth=1;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
		for(auto it=coeffValues[bsplineDepth].begin();it!=coeffValues[bsplineDepth].end();it++)
		{
			int &index=it->second.first;
			Real &value=it->second.second;
			value=x(index);
		}

	Vector y=A*x;
	q=0;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++,q++) it->second.v=y(q);
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getCornerKeysFromCompressedKeys(const int maxDepth,
	const int bsplineDepth,const std::vector<std::vector<long long> >& compressedKeys,
	std::vector<long long>& cornerKeys)
{
	stdext::hash_set<long long> cornerKeys_hash;
	for(auto it=compressedKeys[bsplineDepth].begin();it!=compressedKeys[bsplineDepth].end();it++)
	{
		long long compressedKey=*it;
		int bsplineOffset[3];
		bsplineOffset[0]=compressedKey&0x7fff;
		bsplineOffset[1]=(compressedKey>>15)&0x7fff;
		bsplineOffset[2]=(compressedKey>>30)&0x7fff;

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				for(int k=0;k<3;k++)
				{
					int idx[3];
					idx[0]=(bsplineOffset[0]+i)<<(maxDepth+1-bsplineDepth);
					idx[1]=(bsplineOffset[1]+j)<<(maxDepth+1-bsplineDepth);
					idx[2]=(bsplineOffset[2]+k)<<(maxDepth+1-bsplineDepth);
					long long cornerKey=(long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
					cornerKeys_hash.insert(cornerKey);
				}
	}
	for(auto it=cornerKeys_hash.begin();it!=cornerKeys_hash.end();it++) cornerKeys.push_back(*it);
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::multigridBsplineFitting(const Real& smooth,const Real& interpolate)
{
	int maxDepth=this->maxDepth;
	int maxBsplineDepth=this->maxBsplineDepth;

	Eigen::Matrix<float,3,3> B;
	B<<1,-2,1,1,2,-2,0,0,1;
	B/=2;
	Eigen::Matrix<float,3,3> dB;
	dB<<-2,2,0,2,-4,0,0,2,0;
	dB/=2;
	Eigen::Matrix<float,3,3> ddB;
	ddB<<2,0,0,-4,0,0,2,0,0;
	ddB/=2;

	coeffValues.resize(maxBsplineDepth+1);
	int r=0,p=0;
	for(int bsplineDepth=1;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
	{
		double t=Time();
		setCoeffValuesFromCompressedKeys(bsplineDepth,compressedKeys,coeffValues);
		r+=(int)compressedKeys[bsplineDepth].size();
		p+=(int)coeffValues[bsplineDepth].size();
		std::cout<<"Num. of compressedKeys at bsplineDepth "<<bsplineDepth<<" is "<<compressedKeys[bsplineDepth].size()<<std::endl;
		std::cout<<"Num. of coeffValues at bsplineDepth "<<bsplineDepth<<" is "<<coeffValues[bsplineDepth].size()<<std::endl;
		printf("Set coeffValues in: %f\n", Time()-t);

		if(coeffValues[bsplineDepth].size()>0 && compressedKeys[bsplineDepth].size()>0)
		{
			int pp=0;
			for(auto it=coeffValues[bsplineDepth].begin();it!=coeffValues[bsplineDepth].end();it++,pp++)it->second.first=pp;
			int K=coeffValues[bsplineDepth].size();

			t=Time();
			SparseMatrix fitMatrix(K,K);
			Vector fitVector(K);
			int N;
			getFitMatrixVector(B,bsplineDepth,bsplineDepth,fitMatrix,fitVector,N);
			printf("Set fit matrix and vector in: %f\n", Time()-t);
			
			t=Time();
			SparseMatrix interpolateMatrix(K,K);
			Vector interpolateVector(K);
			int M;
			getInterpolateMatrixVector(B,bsplineDepth,bsplineDepth,interpolateMatrix,interpolateVector,M);
			float w1=0;//(float)interpolate*N/M;
			printf("Set interpolate matrix and vector in: %f\n", Time()-t);

			t=Time();
			SparseMatrix smoothMatrix(K,K);
			getSmoothMatrix(B,dB,ddB,bsplineDepth,bsplineDepth,smoothMatrix);
			float w2=(float)smooth*N*(1<<bsplineDepth);
			printf("Set smooth matrix and vector in: %f\n", Time()-t);
			
			t=Time();
			SparseMatrix globalMatrix=fitMatrix+w1*interpolateMatrix+w2*smoothMatrix;
			Vector globalVector=fitVector+w1*interpolateVector;
			printf("Set global matrix and vector in: %f\n", Time()-t);
			
			//printf("Eigen use %d threads\n", Eigen::nbThreads( ));
			t=Time();
			Eigen::ConjugateGradient<SparseMatrix> solver;
			solver.compute(globalMatrix);
			Vector x=solver.solve(globalVector);
			printf("Solved in: %f\n", Time()-t);

			//Eigen::JacobiSVD<Matrix> svd(fitMatrix+w1*interpolateMatrix+w2*smoothMatrix);
			//float cond=svd.singularValues()(0)/svd.singularValues()(svd.singularValues().size()-1);
			//printf("condition number = %f\n",cond);

			for(auto it=coeffValues[bsplineDepth].begin();it!=coeffValues[bsplineDepth].end();it++)
			{
				int &index=it->second.first;
				Real &value=it->second.second;
				value=x(index);
			}
		}
	}
	std::cout<<"Num. of compressedKeys is "<<r<<std::endl;
	std::cout<<"Num. of coeffValues is "<<p<<std::endl;
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::updateCornerValues()
{
	float unitLen=(float)1.0/(1<<(maxDepth+1));
	
	std::vector<long long> cornerKeys;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++)
	{
		long long cornerKey=it->first;
		cornerKeys.push_back(cornerKey);
	}
	int nprocs=omp_get_num_procs();
	omp_set_num_threads(nprocs);
	//printf("using %d threads\n", nprocs);

	std::vector<float> cornerVs(cornerKeys.size());
#pragma omp parallel for schedule (dynamic,1000)
	for(int i=0;i<cornerKeys.size();i++)
	{
		long long cornerKey=cornerKeys[i];
		float pos[3];
		getPosFromCornerKey(cornerKey,unitLen,pos);
		cornerVs[i]=eval(pos);
	}
	int i=0;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++,i++)
	{
		it->second.v=cornerVs[i];
	}
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::exportVolumeData(const float scale, const Point3D<float> translate, const int d)
{
#ifndef _DEBUG
	const int dim[3] = {d,d,d};
#else
	const int dim[3] = {8,8,8};
#endif
	Vector volume(dim[0]*dim[1]*dim[2]);
	for(int zi = 0; zi < dim[2]; zi++) 
		for(int yi = 0; yi < dim[1]; yi++) 
			for (int xi = 0; xi < dim[0]; xi++) 
			{
				float pos[3]={(float)xi/dim[0],(float)yi/dim[1],(float)zi/dim[2]};
				int index=xi+yi*dim[0]+zi*dim[0]*dim[1];
				volume(index)=eval(pos);
			}
	Eigen::Matrix<double,4,4> transform=Eigen::Matrix<double,4,4>::Identity();
	transform(0,0)=transform(1,1)=transform(2,2)=(float)1.0/scale;
	transform(0,3)=-translate[0];
	transform(1,3)=-translate[1];
	transform(2,3)=-translate[2];
	sameVTIVolumeAndVTKMesh(volume,dim,"volume.vti",transform,"mesh.ply");
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::exportOctreeGrid(const float scale, const Point3D<float> translate)
{
	std::vector<Point3D<float> > vertices;
	std::vector<std::pair<int,int> > edges;
	stdext::hash_set<long long> edgesSet;
	stdext::hash_map<long long,int> vertexIndexMap;

	float unitLen=(float)1.0/(1<<(maxDepth+1));

	OctNode<NodeData,Real>* temp;
	OctNode<NodeData,Real>::NodeIndex nIndex;
	temp=tree.nextLeaf(NULL,nIndex);
	while(temp)
	{
		for(int e=0;e<Cube::EDGES;e++)
		{
			int idx[3];
			int o,i1,i2;
			for(int i=0;i<3;i++){idx[i]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth+1,nIndex.offset[i]<<1,1);}
			Cube::FactorEdgeIndex(e,o,i1,i2);
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
			long long edgeKey=(long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
			if(edgesSet.find(edgeKey)!=edgesSet.end()) continue;
			edgesSet.insert(edgeKey);

			switch(o){
			case 0:
				idx[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[0],0);
				break;
			case 1:
				idx[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[1],0);
				break;
			case 2:
				idx[2]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[2],0);
				break;
			};
			long long vertexKey0=(long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
			auto it0=vertexIndexMap.find(vertexKey0);
			int vertexIndex0;
			if(it0!=vertexIndexMap.end()) vertexIndex0=it0->second;
			else
			{
				Point3D<float> vertex;
				getPosFromCornerKey(vertexKey0,unitLen,vertex.coords);
				vertexIndex0=vertices.size();
				vertexIndexMap[vertexKey0]=vertexIndex0;
				vertices.push_back(vertex);
			}

			switch(o){
			case 0:
				idx[0]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[0],1);
				break;
			case 1:
				idx[1]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[1],1);
				break;
			case 2:
				idx[2]=BinaryNode<Real>::CornerIndex(maxDepth+1,nIndex.depth,nIndex.offset[2],1);
				break;
			};
			long long vertexKey1=(long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
			auto it1=vertexIndexMap.find(vertexKey1);
			int vertexIndex1;
			if(it1!=vertexIndexMap.end()) vertexIndex1=it1->second;
			else
			{
				Point3D<float> vertex;
				getPosFromCornerKey(vertexKey1,unitLen,vertex.coords);
				vertexIndex1=vertices.size();
				vertexIndexMap[vertexKey1]=vertexIndex1;
				vertices.push_back(vertex);
			}

			edges.push_back(std::make_pair(vertexIndex0,vertexIndex1));
		}

		temp=tree.nextLeaf(temp,nIndex);
	}

	for(size_t i=0;i<vertices.size();i++)
		vertices[i]=vertices[i]/scale-translate;

	sameVTKOctree(vertices,edges,"octree.vtk");
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::setMinDepthMCLeafNode(const Real& isoValue,const int& minDepth,const int& useFull)
{
	OctNode<NodeData,Real>* temp;
	OctNode<NodeData,Real>::NodeIndex nIdx;
	temp=tree.nextLeaf(NULL,nIdx);
	while(temp)
	{
		OctNode<NodeData,Real>* temp2=temp;
		OctNode<NodeData,Real>::NodeIndex nIdx2=nIdx;
		temp=tree.nextLeaf(temp,nIdx);
		if(nIdx2.depth<minDepth) setMinDepthMCLeafNode(temp2,nIdx2,isoValue,minDepth,useFull);
	}
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::setMinDepthMCLeafNode(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const Real& isoValue,const int& minDepth,const int& useFull)
{
	if(nIdx.depth==minDepth) return;

	Real cValues[Cube::CORNERS];
	for(int i=0;i<Cube::CORNERS;i++)
	{
		if(cornerValues.find(OctNode<NodeData,Real>::CornerIndex(nIdx,i,maxDepth))==cornerValues.end())
			fprintf(stderr,"Could not find value in corner value table!\n");
		cValues[i]=cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,i,maxDepth)].value();
	}
	if(useFull)
		node->nodeData.mcIndex=MarchingCubes::GetFullIndex(cValues,isoValue);
	else
		node->nodeData.mcIndex=MarchingCubes::GetIndex(cValues,isoValue);

	if(MarchingCubes::HasRoots(node->nodeData.mcIndex))
	{
		long long key;
		Real w,dist;
		Point3D<Real> ctr,n,p;

		if(!node->children)	node->initChildren();
		OctNode<NodeData,Real>::CenterAndWidth(nIdx,ctr,w);

		// Set the center
		dist=this->eval(ctr.coords);
		key=OctNode<NodeData,Real>::CenterIndex(nIdx,maxDepth);
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
			key=OctNode<NodeData,Real>::EdgeIndex(nIdx,i,maxDepth);
			if(cornerValues.find(key)==cornerValues.end())
			{
				dist=this->eval(p.coords);
				cornerValues[key]=VertexData(dist,n);
			}
		}

		// set the face mid-points
		for(int i=0;i<Cube::FACES;i++)
		{
			int dir,off;
			Cube::FactorFaceIndex(i,dir,off);
			p=ctr;
			p[dir]+=-w/2+w*off;
			key=OctNode<NodeData,Real>::FaceIndex(nIdx,i,maxDepth);
			if(cornerValues.find(key)==cornerValues.end())
			{
				dist=this->eval(p.coords);
				cornerValues[key]=VertexData(dist,n);
			}
		}

		for(int i=0;i<Cube::CORNERS;i++)
			setMinDepthMCLeafNode(&node->children[i],nIdx.child(i),isoValue,minDepth,useFull);
	}
}

typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<float,Eigen::Dynamic,1> Vector;

Vector MULT(const Vector& U,const Vector& V) 
{
	Vector W=Vector::Zero(U.size()+V.size()-1);
	for(int i=0;i<U.size();i++) 
		for (int j=0;j<V.size();j++)
			W(i+j) += U(i)*V(j);
	return W;
}

Matrix INTofMULT(const Matrix& A, const Matrix& B) 
{
	Matrix C=Matrix::Zero(A.rows(),B.rows());
	for (int i=0;i<A.rows();i++) 
		for (int j=0;j<B.rows();j++) 
		{
			Vector W=MULT(A.row(i),B.row(j));
			//std::cerr<<"W=\n"<<W<<std::endl; 
			float integral_sum=0;
			for(int k=0;k<W.size();k++)
				integral_sum += W(k)/(k+1);
			C(i,j) = integral_sum;
		}
	return C;
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getSmoothMatrix(const Eigen::Matrix<float,3,3>& B,const Eigen::Matrix<float,3,3>& dB,const Eigen::Matrix<float,3,3>& ddB,
	const int minBsplineDepth,const int maxBsplineDepth, SparseMatrix& smoothMatrix)
{
	assert(minBsplineDepth==maxBsplineDepth);

	Matrix IB2,IdB2,IddB2;
	IB2=INTofMULT(B,B);
	IdB2=INTofMULT(dB,dB);
	IddB2=INTofMULT(ddB,ddB);

	Eigen::Matrix<float,5,1> hashtableIB2,hashtableIdB2,hashtableIddB2;
	hashtableIB2=Eigen::Matrix<float,5,1>::Zero();
	hashtableIdB2=Eigen::Matrix<float,5,1>::Zero();
	hashtableIddB2=Eigen::Matrix<float,5,1>::Zero();
	for(int i=2;i>=-2;i--)
		for(int j=0;j<3;j++)
			if(i+j>=0 && i+j<3)
			{
				hashtableIB2(2-i)+=IB2(j,i+j);
				hashtableIdB2(2-i)+=IdB2(j,i+j);
				hashtableIddB2(2-i)+=IddB2(j,i+j);
			}

	float hashtableWeights[5][5][5];
	for(int i=0;i<5;i++)
		for(int j=0;j<5;j++)
			for(int k=0;k<5;k++)
			{
				hashtableWeights[i][j][k]=hashtableIddB2(i)*hashtableIB2(j)*hashtableIB2(k)
					+hashtableIB2(i)*hashtableIddB2(j)*hashtableIB2(k)
					+hashtableIB2(i)*hashtableIB2(j)*hashtableIddB2(k)
					+hashtableIdB2(i)*hashtableIdB2(j)*hashtableIB2(k)*2
					+hashtableIdB2(i)*hashtableIB2(j)*hashtableIdB2(k)*2
					+hashtableIB2(i)*hashtableIdB2(j)*hashtableIdB2(k)*2;
			}

	std::vector<Triplet> smoothMatrix_triplets;
	int bsplineDepth=minBsplineDepth;

	std::vector<long long> coeffKeys;
	std::vector<int> coeffIndices;
	for(auto it=coeffValues[bsplineDepth].begin();it!=coeffValues[bsplineDepth].end();it++) 
	{
		long long coeffKey=it->first;
		coeffKeys.push_back(coeffKey);
		coeffIndices.push_back(it->second.first);
	}

	int nprocs=omp_get_num_procs();
	omp_set_num_threads(nprocs);
	//printf("Get smooth matrix using %d threads\n", nprocs);

	std::vector<std::vector<Triplet> > smoothMatrix_triplets_buffers(nprocs);
#pragma omp parallel for schedule (dynamic,1000)
	for(int p=0;p<coeffKeys.size();p++)
	{
		long long coeffKey=coeffKeys[p];
		int bsplineOffset[3];
		bsplineOffset[0]=coeffKey&0x7fff;
		bsplineOffset[1]=(coeffKey>>15)&0x7fff;
		bsplineOffset[2]=(coeffKey>>30)&0x7fff;

		int thread_id=omp_get_thread_num();

		for(int i=0;i<5;i++)
			for(int j=0;j<5;j++)
				for(int k=0;k<5;k++)
				{
					int bsplineOffset2[3];
					bsplineOffset2[0]=bsplineOffset[0]-2+i;
					bsplineOffset2[1]=bsplineOffset[1]-2+j;
					bsplineOffset2[2]=bsplineOffset[2]-2+k;
					if(bsplineOffset2[0]>=0 && bsplineOffset2[1]>=0 && bsplineOffset2[2]>=0)
					{
						long long coeffKey2=(long long)(bsplineOffset2[0]) | (long long)(bsplineOffset2[1])<<15 | (long long)(bsplineOffset2[2])<<30;
						auto it2=coeffValues[bsplineDepth].find(coeffKey2);
						if(it2!=coeffValues[bsplineDepth].end())
						{
							int row=coeffIndices[p];
							int col=it2->second.first;
							if(col>=row)
							{
								float weight=hashtableWeights[i][j][k];
								smoothMatrix_triplets_buffers[thread_id].push_back(Triplet(row,col,weight));
								if(col>row) 
									smoothMatrix_triplets_buffers[thread_id].push_back(Triplet(col,row,weight));
							}
						}
					}
				}
	}
	coeffKeys.swap(std::vector<long long>());
	coeffIndices.swap(std::vector<int>());
	for(int i=0;i<nprocs;i++) 
	{
		smoothMatrix_triplets.insert(smoothMatrix_triplets.end(), smoothMatrix_triplets_buffers[i].begin(), smoothMatrix_triplets_buffers[i].end());
		smoothMatrix_triplets_buffers[i].swap(std::vector<Triplet>());
	}

	smoothMatrix.resize(coeffValues[bsplineDepth].size(),coeffValues[bsplineDepth].size());
	smoothMatrix.setFromTriplets(smoothMatrix_triplets.begin(),smoothMatrix_triplets.end());

	float unitLen=(float)1.0/(1<<bsplineDepth);
	smoothMatrix*=(unitLen*unitLen*unitLen);

	//std::cout<<"smoothMatrix=\n"<<smoothMatrix.toDense()<<std::endl;
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getFitMatrixVector(
	const Eigen::Matrix<float,3,3>& B,
	const int minBsplineDepth,
	const int maxBsplineDepth, 
	SparseMatrix& fitMatrix, 
	Vector& fitVector, 
	int& N)
{
	assert(minBsplineDepth==maxBsplineDepth);

	int bsplineDepth=minBsplineDepth;
	float unitLen=(float)1.0/(1<<(maxDepth+1));
	int K=coeffValues[bsplineDepth].size();

	std::vector<long long> cornerKeys;
	getCornerKeysFromCompressedKeys(maxDepth,bsplineDepth,compressedKeys,cornerKeys);
	N=cornerKeys.size();

	int nprocs=omp_get_num_procs();
	omp_set_num_threads(nprocs);
	//printf("Get fit matrix using %d threads\n", nprocs);

	std::vector<Triplet> fitMatrix_temp_triplets;
	Vector fitVector_temp=Vector::Zero(N);
	std::vector<std::vector<Triplet> > fitMatrix_temp_triplets_buffers(nprocs);
#pragma omp parallel for schedule (dynamic,1000)
	for(int i=0;i<N;i++)
	{
		int thread_id=omp_get_thread_num();

		long long cornerKey=cornerKeys[i];
		float pos[3];
		getPosFromCornerKey(cornerKey,unitLen,pos);

		std::vector<int> coeffIndices;
		std::vector<float> coeffWeights;
		float value;
		getCoeffIndicesWeightsValueFromPos(B,1,bsplineDepth-1,pos,coeffIndices,coeffWeights,value);
		VertexData cornerValue=cornerValues[cornerKey];
		float residual=cornerValue.v-value;
		fitVector_temp(i)=residual*cornerValue.w;

		getCoeffIndicesWeightsValueFromPos(B,bsplineDepth,bsplineDepth,pos,coeffIndices,coeffWeights,value);
		for(int j=0;j<coeffIndices.size();j++) 
		{
			fitMatrix_temp_triplets_buffers[thread_id].push_back(Triplet(i,coeffIndices[j],coeffWeights[j]*cornerValue.w));
		}
	}
	SparseMatrix fitMatrix_temp(N,K);
	for(int i=0;i<nprocs;i++) 
	{
		fitMatrix_temp_triplets.insert(fitMatrix_temp_triplets.end(), fitMatrix_temp_triplets_buffers[i].begin(), fitMatrix_temp_triplets_buffers[i].end());
		fitMatrix_temp_triplets_buffers[i].swap(std::vector<Triplet>());
	}
	fitMatrix_temp.setFromTriplets(fitMatrix_temp_triplets.begin(),fitMatrix_temp_triplets.end());
	fitMatrix=fitMatrix_temp.transpose()*fitMatrix_temp;
	fitVector=fitMatrix_temp.transpose()*fitVector_temp;

	//printf("Creating A_mapTriplets ...\n");
	//stdext::hash_map<long long, float> A_mapTriplets;
	//Vector b=Vector::Zero(coeffValues[bsplineDepth].size());
	//for(auto it=cornerKeys.begin();it!=cornerKeys.end();it++)
	//{
	//	long long cornerKey=*it;
	//	float pos[3];
	//	getPosFromCornerKey(cornerKey,unitLen,pos);

	//	std::vector<int> coeffIndices;
	//	std::vector<float> coeffWeights;
	//	float value;
	//	getCoeffIndicesWeightsValueFromPos(B,1,bsplineDepth-1,pos,coeffIndices,coeffWeights,value);
	//	VertexData cornerValue=cornerValues[cornerKey];
	//	float residual=cornerValue.v-value;

	//	getCoeffIndicesWeightsValueFromPos(B,bsplineDepth,bsplineDepth,pos,coeffIndices,coeffWeights,value);
	//	float w2=cornerValue.w*cornerValue.w;
	//	for(int i=0;i<coeffIndices.size();i++) 
	//	{
	//		for(int j=0;j<coeffIndices.size();j++)
	//		{
	//			float wiwj=coeffWeights[i]*coeffWeights[j]*w2;
	//			long long key=(long long)coeffIndices[i]<<0 | (long long)coeffIndices[j]<<31;
	//			auto it2=A_mapTriplets.find(key);
	//			if(it2==A_mapTriplets.end())
	//				A_mapTriplets[key]=wiwj;
	//			else
	//				it2->second+=wiwj;
	//		}
	//		b(coeffIndices[i])+=coeffWeights[i]*residual*w2;
	//	}
	//}
	//printf("Finished creating A_mapTriplets\n");

	//printf("Creating A_triplets ...\n");
	//std::vector<Triplet> A_triplets;
	//for(auto it=A_mapTriplets.begin();it!=A_mapTriplets.end();it++)
	//{
	//	long long key=it->first;
	//	int i=(key>>0)&0x7fffffff;
	//	int j=(key>>31)&0x7fffffff;
	//	A_triplets.push_back(Triplet(i,j,it->second));
	//}
	//A_mapTriplets.swap(stdext::hash_map<long long, float>());
	//printf("Finished creating A_triplets\n");

	//printf("Creating A ...\n");
	//SparseMatrix A(coeffValues[bsplineDepth].size(), coeffValues[bsplineDepth].size());
	//A.setFromTriplets(A_triplets.begin(),A_triplets.end());
	//A_triplets.swap(std::vector<Triplet>());
	//printf("Finished Creating A\n");

	//Eigen::ConjugateGradient<SparseMatrix> solver;
	//solver.compute(A+0.0001*bsplineDepth*bsplineDepth*cornerKeys.size()*smoothMatrix);
	//Vector x=solver.solve(b);
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::sortPointsByCurvatures(std::vector<int>& beginIndices)
{
	//for(int i=0;i<mInfoGlobal.vertices.size();i++)
	//{
	//	if(mInfoGlobal.vertexCurvatures[i]>(float)(1<<maxBsplineDepth) && mInfoGlobal.vertexCurvatures[i]<(float)(1<<(maxBsplineDepth+1)))
	//	{
	//	}
	//}
}


template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getInterpolateMatrixVector(
	const Eigen::Matrix<float,3,3>& B,
	const int minBsplineDepth,
	const int maxBsplineDepth, 
	SparseMatrix& interpolateMatrix, 
	Vector& interpolateVector,
	int& M)
{
	assert(minBsplineDepth==maxBsplineDepth);

	int bsplineDepth=minBsplineDepth;
	float unitLen=(float)1.0/(1<<(maxDepth+1));
	int K=coeffValues[bsplineDepth].size();

	//M=0;
	//for(int i=0;i<mInfoGlobal.vertices.size();i++)
	//{
	//	if(mInfoGlobal.vertexCurvatures[i]>(float)(1<<maxBsplineDepth) && mInfoGlobal.vertexCurvatures[i]<(float)(1<<(maxBsplineDepth+1)))
	//	{
	//		float pos[3];
	//		pos[0]=mInfoGlobal.vertices[i].coords[0];
	//		pos[1]=mInfoGlobal.vertices[i].coords[1];
	//		pos[2]=mInfoGlobal.vertices[i].coords[2];

	//		std::vector<int> coeffIndices;
	//		std::vector<float> coeffWeights;
	//		float value;
	//		getCoeffIndicesWeightsValueFromPos(B,1,bsplineDepth-1,pos,coeffIndices,coeffWeights,value);
	//		float residual=(float)0-value;
	//		d_vec.push_back(residual*(float)1.0f);

	//		getCoeffIndicesWeightsValueFromPos(B,bsplineDepth,bsplineDepth,pos,coeffIndices,coeffWeights,value);
	//		for(int i=0;i<coeffIndices.size();i++) C_triplets.push_back(Triplet(q3,coeffIndices[i],coeffWeights[i]*(float)1.0f));
	//		q3++;
	//	}
	//}

	//SparseMatrix CtC(p2,p2);
	//Vector Ctd=Vector::Zero(p2);
	//std::vector<Triplet> C_triplets;
	//std::vector<float> d_vec;
	//int q3=0;

	//printf("q3=%d\n",q3);
	//Vector d(q3);
	//for(int i=0;i<q3;i++) d(i)=d_vec[i];
	//SparseMatrix C(q3,p2);
	//C.setFromTriplets(C_triplets.begin(),C_triplets.end());
	//CtC=C.transpose()*C;
	//Ctd=C.transpose()*d;
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::gradient(const float pos[3], float gradient[3])
{
	//float delta=(float)1.0/(1<<maxDepth)/100;
	//float posX[3]={pos[0]+delta,pos[1],pos[2]};
	//float posY[3]={pos[0],pos[1]+delta,pos[2]};
	//float posZ[3]={pos[0],pos[1],pos[2]+delta};

	//float value=eval(pos);
	//float valueX=eval(posX);
	//float valueY=eval(posY);
	//float valueZ=eval(posZ);

	//gradient[0]=(valueX-value)/delta;
	//gradient[1]=(valueY-value)/delta;
	//gradient[2]=(valueZ-value)/delta;

	Eigen::Matrix<float,3,3> B;
	B<<1,-2,1,
		1,2,-2,
		0,0,1;
	B/=2;
	Eigen::Matrix<float,3,3> dB;
	dB<<-2,2,0,
		2,-4,0,
		0,2,0;
	dB/=2;

	std::vector<int> coeffIndices;
	std::vector<float> coeffWeights;

	getCoeffIndicesWeightsValueFromPos(dB,B,B,1,maxBsplineDepth,pos,coeffIndices,coeffWeights,gradient[0]);
	getCoeffIndicesWeightsValueFromPos(B,dB,B,1,maxBsplineDepth,pos,coeffIndices,coeffWeights,gradient[1]);
	getCoeffIndicesWeightsValueFromPos(B,B,dB,1,maxBsplineDepth,pos,coeffIndices,coeffWeights,gradient[2]);
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getCoeffIndicesWeightsValueFromPos(
	const Eigen::Matrix<float,3,3>& Bx,const Eigen::Matrix<float,3,3>& By,const Eigen::Matrix<float,3,3>& Bz, 
	const int minBsplineDepth,const int maxBsplineDepth,const float pos[3], 
	std::vector<int>& coeffIndices,std::vector<float>& coeffWeights,float& value)
{
	value=0;
	coeffIndices.clear();
	coeffWeights.clear();
	if(maxBsplineDepth<minBsplineDepth || minBsplineDepth<1) return;
	for(int bsplineDepth=minBsplineDepth;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
	{
		if(!coeffValues[bsplineDepth].size()) continue;

		int bposInt[3];
		float u[3];
		getBposIntU(bsplineDepth,pos,bposInt,u);

		Eigen::Matrix<float,3,3> u_mat;
		u_mat<< 1,1,1,u[0],u[1],u[2],u[0]*u[0],u[1]*u[1],u[2]*u[2];
		Eigen::Matrix<float,3,3> Bu;
		Bu.col(0)=Bx*u_mat.col(0);
		Bu.col(1)=By*u_mat.col(1);
		Bu.col(2)=Bz*u_mat.col(2);

		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				for(int k=0;k<3;k++)
				{
					int bsplineOffset[3];
					bsplineOffset[0]=bposInt[0]-2+i;
					bsplineOffset[1]=bposInt[1]-2+j;
					bsplineOffset[2]=bposInt[2]-2+k;
					if(bsplineOffset[0]>=0 && bsplineOffset[1]>=0 && bsplineOffset[2]>=0)
					{
						long long key=(long long)(bsplineOffset[0]) | (long long)(bsplineOffset[1])<<15 | (long long)(bsplineOffset[2])<<30;
						auto it2=coeffValues[bsplineDepth].find(key);
						if(it2!=coeffValues[bsplineDepth].end())
						{
							int coeffIndex=it2->second.first;
							float coeffValue=it2->second.second;
							float coeffWeight=Bu(i,0)*Bu(j,1)*Bu(k,2);
							coeffWeight*= (1<<bsplineDepth); //Note this scale is very important!
							coeffIndices.push_back(coeffIndex);
							coeffWeights.push_back(coeffWeight);
							value+=coeffValue*coeffWeight;
						}
					}
				}
	}
}