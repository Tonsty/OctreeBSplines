#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Kdtree.h"
#include "vtkHelper.h"

template<class NodeData,class Real,class VertexData>
template<class Vertex>
int OctreeBspline<NodeData,Real,VertexData>::set3(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,
	const int& maxDepth,const int& setCenter,const Real& flatness,
	Point3D<Real>& translate,Real& scale,const int& noTransform)
{
	this->maxDepth=maxDepth;
	OctNode<NodeData,Real>::NodeIndex nIdx;

	MeshInfo<Real> &mInfo=mInfoGlobal;
	mInfo.set2(vertices,polygons,Real(1.1),translate,scale,noTransform);
	std::vector<int> myVertices;
	for(int i=0;i<mInfo.vertices.size();i++) myVertices.push_back(i);

	cornerValues.clear();
	for(int c=0;c<Cube::CORNERS;c++)
		cornerValues[OctNode<NodeData,Real>::CornerIndex(nIdx,c,maxDepth)]=VertexData(); //allocate hash-table entry

	printf("Allocating cornerValues ...\n");
	compressedKeys.resize(maxBsplineDepth+1);
	setChildren3(&tree,nIdx,myVertices,mInfo,maxDepth,setCenter,flatness);
	printf("Finished allocating\n");

	KDTree<Real> kdtree;
	kdtree.setInputPoints(mInfo.vertices);
	Real unitLen=(Real)1.0/(1<<(maxDepth+1));
	std::vector<Point3D<Real> > queries;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++)
	{
		long long cornerKey=it->first;
		Point3D<Real> point;
		getPosFromCornerKey(cornerKey,unitLen,point.coords);
		queries.push_back(point);
	}

	std::vector<std::vector<int> > indices; 
	std::vector<std::vector<Real> > dists;
	kdtree.KnnSearch(queries,indices,dists,1);

	Real dist;
	Point3D<Real> n;
	int i=0;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++,i++)
	{
		Point3D<Real> point=queries[i];
		Point3D<Real> point2=mInfo.vertices[indices[i][0]];
		Point3D<Real> normal2=mInfo.vertexNormals[indices[i][0]];
		setDistanceAndNormal3(point,point2,normal2,dist,n);
		it->second=VertexData(dist,n);
	}

	//scale=1.f;
	//translate[0]=translate[1]=translate[2]=-0.5;
	//Point3D<Real> center;
	//center[0]=(Real)0.5;
	//center[1]=(Real)0.5;
	//center[2]=(Real)0.5;
	//for(auto it=cornerValues.begin();it!=cornerValues.end();it++)
	//{
	//	long long cornerKey=it->first;
	//	int idx[3];
	//	idx[0]=cornerKey&0x7fff;
	//	idx[1]=(cornerKey>>15)&0x7fff;
	//	idx[2]=(cornerKey>>30)&0x7fff;
	//	idx[0]>>=(maxDepth+1-maxDepth);
	//	idx[1]>>=(maxDepth+1-maxDepth);
	//	idx[2]>>=(maxDepth+1-maxDepth);
	//	Point3D<Real> pos;
	//	pos[0]=(Real)idx[0]/(1<<(maxDepth));
	//	pos[1]=(Real)idx[1]/(1<<(maxDepth));
	//	pos[2]=(Real)idx[2]/(1<<(maxDepth));
	//	it->second.v=Distance(pos,center)-(Real)0.5;
	//}

	return 1;
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::setDistanceAndNormal3(const Point3D<Real>& p,const Point3D<Real>& p2,const Point3D<Real>& n2,Real& dist,Point3D<Real>& n)
{
	dist=Distance(p,p2);
	n=(p-p2)/dist;

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
	const std::vector<int>& vindices,MeshInfo<MeshReal>& mInfo,
	const int& maxDepth,const int& setCenter,
	const Real& flatness,
	stdext::hash_map<long long,std::vector<int>*>* triangleMap)
{
	long long key;
	Real w;
	Point3D<Real> ctr;
	if(!vindices.size())		return -1;
	if(nIdx.depth==maxDepth)	return -1;

	if(flatness>0 && nIdx.depth>0)
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
		for(size_t j=0;j<vindices.size();j++)
		{
			Point3D<MeshReal> t=mInfo.vertices[vindices[j]];
			Point3D<MeshReal> ctr2;
			MeshReal w2=w;
			ctr2[0]=ctr[0];
			ctr2[1]=ctr[1];
			ctr2[2]=ctr[2];
			if(PointInCube(ctr2,w2,t)) myVertices.push_back(vindices[j]);
		}
		if(myVertices.size() && setChildren2(&node->children[i],nIdx.child(i),myVertices,mInfo,maxDepth,setCenter,flatness,triangleMap)>0) retCount++;
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

	setDistanceAndNormal3(pp,t,nn,dist,n);
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

typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<float,Eigen::Dynamic,1> Vector;
typedef Eigen::SparseMatrix<float> SparseMatrix;
typedef Eigen::Triplet<float> Triplet;

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::directBsplineFitting()
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
		b(q)=it->second.v;
		long long cornerKey=it->first;
		float pos[3];
		getPosFromCornerKey(cornerKey,unitLen,pos);
		//std::cout<<q<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
		std::vector<int> coeffIndices;
		std::vector<float> coeffWeights;
		float value;
		getCoeffIndicesWeightsValueFromPos(B,1,maxBsplineDepth,pos,coeffIndices,coeffWeights,value);
		for(int i=0;i<coeffIndices.size();i++) A_triplets.push_back(Triplet(q,coeffIndices[i],coeffWeights[i]));
	}
	std::cout<<"Num. of cornerValues is "<<q<<std::endl;

	SparseMatrix A(q, p);
	A.setFromTriplets(A_triplets.begin(),A_triplets.end());
//#ifdef _DEBUG
//	std::cout<<"A=\n"<<A.toDense()<<std::endl;
//	std::cout<<"b=\n"<<b<<std::endl; 
//#endif

	Vector x;
	//Eigen::ConjugateGradient<SparseMatrix> cg;
	//cg.compute(A.transpose()*A);
	//x=cg.solve(b);
	Eigen::LeastSquaresConjugateGradient<SparseMatrix> lscg;
	lscg.compute(A);
	x=lscg.solve(b);
	std::cout<<"finised solving"<<std::endl;

	for(int bsplineDepth=1;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
		for(auto it=coeffValues[bsplineDepth].begin();it!=coeffValues[bsplineDepth].end();it++)
		{
			int &index=it->second.first;
			Real &value=it->second.second;
			value=x(index);
		}

	Vector y=A*x;
//#ifdef _DEBUG
//	std::cout<<"y=A*x=\n"<<y<<std::endl; 
//	std::cout<<"y-b=\n"<<y-b<<std::endl; 
//#endif
	q=0;
	for(auto it=cornerValues.begin();it!=cornerValues.end();it++,q++) it->second.v=y(q);
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::getCornerKeysFromCompressedKeys(const int maxDepth,
	const int bsplineDepth,const std::vector<std::vector<long long> >& compressedKeys,
	stdext::hash_set<long long>& cornerKeys)
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
					int idx[3];
					idx[0]=(bsplineOffset[0]+i)<<(maxDepth+1-bsplineDepth);
					idx[1]=(bsplineOffset[1]+j)<<(maxDepth+1-bsplineDepth);
					idx[2]=(bsplineOffset[2]+k)<<(maxDepth+1-bsplineDepth);
					long long cornerKey=(long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
					cornerKeys.insert(cornerKey);
				}
	}
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::multigridBsplineFitting()
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

	float unitLen=(float)1.0/(1<<(maxDepth+1));

	coeffValues.resize(maxBsplineDepth+1);
	int r=0,p=0,q=0;
	for(int bsplineDepth=1;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
	{
		setCoeffValuesFromCompressedKeys(bsplineDepth,compressedKeys,coeffValues);
		if(bsplineDepth==1) continue;
		r+=(int)compressedKeys[bsplineDepth].size();
		p+=(int)coeffValues[bsplineDepth].size();
		std::cout<<"Num. of compressedKeys at bsplineDepth "<<bsplineDepth<<" is "<<compressedKeys[bsplineDepth].size()<<std::endl;
		std::cout<<"Num. of coeffValues at bsplineDepth "<<bsplineDepth<<" is "<<coeffValues[bsplineDepth].size()<<std::endl;

		if(coeffValues[bsplineDepth].size()>0 && compressedKeys[bsplineDepth].size()>0)
		{
			int p2=0;
			for(auto it=coeffValues[bsplineDepth].begin();it!=coeffValues[bsplineDepth].end();it++,p2++)it->second.first=p2;

			stdext::hash_set<long long> cornerKeys;
			getCornerKeysFromCompressedKeys(maxDepth,bsplineDepth,compressedKeys,cornerKeys);
			q+=(int)cornerKeys.size();

			std::vector<Triplet> A_triplets;
			Vector b(cornerKeys.size());
			int q2=0;
			for(auto it=cornerKeys.begin();it!=cornerKeys.end();it++,q2++)
			{
				long long cornerKey=*it;
				float pos[3];
				getPosFromCornerKey(cornerKey,unitLen,pos);

				std::vector<int> coeffIndices;
				std::vector<float> coeffWeights;
				float value;

				getCoeffIndicesWeightsValueFromPos(B,1,bsplineDepth-1,pos,coeffIndices,coeffWeights,value);
				b(q2)=cornerValues[cornerKey].v-value;

				getCoeffIndicesWeightsValueFromPos(B,bsplineDepth,bsplineDepth,pos,coeffIndices,coeffWeights,value);
				for(int i=0;i<coeffIndices.size();i++) A_triplets.push_back(Triplet(q2,coeffIndices[i],coeffWeights[i]));
			}

			SparseMatrix A(q2, p2);
			A.setFromTriplets(A_triplets.begin(),A_triplets.end());
			Vector x;
			//solveAbx(A,b,x);
			//std::cout<<"finised pre-solving at bsplineDepth "<<bsplineDepth<<std::endl;

			SparseMatrix smoothMatrix;
			getSmoothMatrix(B,dB,ddB,bsplineDepth,bsplineDepth,smoothMatrix);
			Eigen::ConjugateGradient<SparseMatrix> solver;
			solver.compute(A.transpose()*A+0.0*bsplineDepth*bsplineDepth*q2*smoothMatrix);
			x=solver.solve(A.transpose()*b);

			//std::cout<<"finised solving at bsplineDepth "<<bsplineDepth<<std::endl;

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

	for(auto it=cornerValues.begin();it!=cornerValues.end();it++,q++)
	{
		long long cornerKey=it->first;
		float pos[3];
		getPosFromCornerKey(cornerKey,unitLen,pos);
		//std::cout<<q<<" "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
		it->second.v=eval(pos);
	}
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::exportVTKData(const float scale, const Point3D<float> translate)
{
	Eigen::Matrix<float,3,3> B;
	B<<1,-2,1,1,2,-2,0,0,1;
	B/=2;

#ifndef _DEBUG
	const int dim[3] = {128,128,128};
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
	sameVTKVolumeAndMesh(volume,dim,"volume.vti",transform,"mesh.ply");
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::setMCLeafNodeToMaxDepth(const Real& isoValue,const int& useFull)
{
	OctNode<NodeData,Real>* temp;
	Real cValues[Cube::CORNERS];

	OctNode<NodeData,Real>::NodeIndex nIdx;
	temp=tree.nextLeaf(NULL,nIdx);
	while(temp)
	{
		OctNode<NodeData,Real>* temp2=temp;
		OctNode<NodeData,Real>::NodeIndex nIdx2=nIdx;
		temp=tree.nextLeaf(temp,nIdx);
		if(nIdx2.depth<maxDepth) setMCLeafNodeToMaxDepth(temp2,nIdx2,isoValue,useFull);
	}
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::setMCLeafNodeToMaxDepth(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const Real& isoValue,const int& useFull)
{
	if(nIdx.depth==maxDepth) return;

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
			setMCLeafNodeToMaxDepth(&node->children[i],nIdx.child(i),isoValue,useFull);
	}
}

Vector MULT(const Vector &U,const Vector &V) 
{
	Vector W=Vector::Zero(U.size()+V.size()-1);
	for(int i=0;i<U.size();i++) 
		for (int j=0;j<V.size();j++)
			W(i+j) += U(i)*V(j);
	return W;
}

Matrix INTofMULT(const Matrix &A, const Matrix &B) 
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
	const int minBsplineDepth,const int maxBsplineDepth, Eigen::SparseMatrix<float> &smoothMatrix)
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
	for(auto it=coeffValues[bsplineDepth].begin();it!=coeffValues[bsplineDepth].end();it++)
	{
		long long coeffKey=it->first;
		int bsplineOffset[3];
		bsplineOffset[0]=coeffKey&0x7fff;
		bsplineOffset[1]=(coeffKey>>15)&0x7fff;
		bsplineOffset[2]=(coeffKey>>30)&0x7fff;

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
							int row=it->second.first;
							int col=it2->second.first;
							if(col>=row)
							{
								float weight=hashtableWeights[i][j][k];
								smoothMatrix_triplets.push_back(Triplet(row,col,weight));
								if(col>row)
									smoothMatrix_triplets.push_back(Triplet(col,row,weight));
							}
						}
					}
				}
	}
	smoothMatrix.resize(coeffValues[bsplineDepth].size(),coeffValues[bsplineDepth].size());
	smoothMatrix.setFromTriplets(smoothMatrix_triplets.begin(),smoothMatrix_triplets.end());

	float unitLen=(float)1.0/(1<<bsplineDepth);
	smoothMatrix*=(unitLen*unitLen*unitLen);

	//std::cout<<"smoothMatrix=\n"<<smoothMatrix.toDense()<<std::endl;
}