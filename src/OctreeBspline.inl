#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "vtkHelper.h"

template<class NodeData,class Real,class VertexData>
template<class Vertex>
int OctreeBspline<NodeData,Real,VertexData>::set2(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,
	const int& maxDepth,const int& setCenter,const Real& flatness,
	Point3D<Real>& translate,Real& scale,const int& noTransform)
{
	this->maxDepth=maxDepth;
	OctNode<NodeData,Real>::NodeIndex nIdx;

	MeshInfo<double> &mInfo=mInfoGlobal;
	std::vector<int> myVertices;

	mInfo.set2(vertices,polygons,Real(1.05),translate,scale,noTransform);
	myVertices.resize(mInfo.vertices.size());
	for(int i=0;i<int(mInfo.vertices.size());i++) myVertices[i]=i;

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
	for(size_t i=0;i<vindices.size();i++) 
	{
		MeshReal temp;
		for(int j=0;j<3;j++) t[j]=mInfo.vertices[vindices[i]][j];

		//point to point distance
		temp=Distance(pp,t); 
		if(!i || temp<dist )
		{
			closest=i;
			dist=(Real)temp;
		}
	}
	Point3D<MeshReal> nn;
	int vFlag;

	closest=vindices[closest];

	for(int i=0;i<3;i++) 
	{
		t[i]=mInfo.vertices[closest][i];
		nn[i]=mInfo.vertexNormals[closest][i];
	}

	//only points locating inside narrow band consider point to plane distance
	Real narrow_band=(Real)0.1;
	Point3D<MeshReal> n2;
	if (dist<narrow_band) {
		n2=NearestPointOnPlane(pp,t,nn,vFlag);
		//point to plane distance
		Real dist2=(Real)Distance(pp,n2); 
		dist = dist2;
		//point to plane orientation
		n2=(pp-n2)/dist2; 
	}
	else 
	{
		//point to point orientation
		n2 = (pp-t)/Distance(pp,t);
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
		if(cornerValues.find(key)==cornerValues.end())
		{
			setDistanceAndNormal2(vindices,mInfo,p,dist,n);
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
			setDistanceAndNormal2(vindices,mInfo,p,dist,n);
			cornerValues[key]=VertexData(dist,n);
		}
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
		for(size_t j=0;j<vindices.size();j++)
		{
			Point3D<MeshReal> t = mInfo.vertices[vindices[j]];
			Point3D<MeshReal> ctr2;
			MeshReal w2=w;
			ctr2[0]=ctr[0];
			ctr2[1]=ctr[1];
			ctr2[2]=ctr[2];
			if(PointInCube(ctr2,w2,t)) myVerticesExp.push_back(vindices[j]);
		}
		//std::vector<int> myVertices;
		//for (size_t j=0;j<myVerticesExp.size();j++)
		//{
		//	Point3D<MeshReal> t = mInfo.vertices[myVerticesExp[j]];
		//	Point3D<MeshReal> ctr2;
		//	MeshReal w2=w;
		//	ctr2[0]=ctr[0];
		//	ctr2[1]=ctr[1];
		//	ctr2[2]=ctr[2];
		//	if(PointInCube(ctr2,w2,t)) myVertices.push_back(myVerticesExp[j]);
		//}
		
		if(myVerticesExp.size() && setChildren2(&node->children[i],nIdx.child(i),myVerticesExp,mInfo,maxDepth,setCenter,flatness,triangleMap)>0) retCount++;
	}
	if(maxBsplineDepth>0 && ((nIdx.depth==maxBsplineDepth-1) || (nIdx.depth<maxBsplineDepth-1 && retCount<8)))
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

		for(int i=0;i<=2;i++)
			for(int j=0;j<=2;j++)
				for(int k=0;k<=2;k++)
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
void OctreeBspline<NodeData,Real,VertexData>::solveAbx(
	const Eigen::SparseMatrix<float>& A, 
	const Eigen::Matrix<float,Eigen::Dynamic,1>& b, 
	Eigen::Matrix<float,Eigen::Dynamic,1>& x)
{
	typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matrix;
	typedef Eigen::SparseMatrix<float> SparseMatrix;

	if(A.cols()<500) 
	{
		Eigen::ColPivHouseholderQR<Matrix> qr;
		//qr.compute((A.transpose()*A).toDense());
		//x=qr.solve(A.transpose()*b);
		qr.compute(A.toDense());
		x=qr.solve(b);
	}
	else
	{
		Eigen::LeastSquaresConjugateGradient<SparseMatrix> cg;
		//cg.compute(A.transpose()*A);
		//x=cg.solve(A.transpose()*b);
		cg.compute(A);
		x=cg.solve(b);
		std::cout << "cg.iterations = " << cg.iterations() << std::endl;
		std::cout << "cg.maxIterations = " << cg.maxIterations() << std::endl;
		std::cout << "cg.tolerance = " << cg.tolerance() << std::endl;
	}
	//std::cout<<"x=\n"<<x<<std::endl; 
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::directBsplineFitting()
{
	typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matrix;
	typedef Eigen::Matrix<float,Eigen::Dynamic,1> Vector;
	typedef Eigen::SparseMatrix<float> SparseMatrix;
	typedef Eigen::Triplet<float> Triplet;

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
	solveAbx(A,b,x);
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
	typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matrix;
	typedef Eigen::Matrix<float,Eigen::Dynamic,1> Vector;
	typedef Eigen::SparseMatrix<float> SparseMatrix;
	typedef Eigen::Triplet<float> Triplet;

	int maxDepth=this->maxDepth;
	int maxBsplineDepth=this->maxBsplineDepth;

	Eigen::Matrix<float,3,3> B;
	B<<1,-2,1,1,2,-2,0,0,1;
	B/=2;
	float unitLen=(float)1.0/(1<<(maxDepth+1));

	coeffValues.resize(maxBsplineDepth+1);
	int r=0,p=0,q=0;
	for(int bsplineDepth=1;bsplineDepth<=maxBsplineDepth;bsplineDepth++)
	{
		setCoeffValuesFromCompressedKeys(bsplineDepth,compressedKeys,coeffValues);
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
			solveAbx(A,b,x);
			std::cout<<"finised solving at bsplineDepth "<<bsplineDepth<<std::endl;

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
void OctreeBspline<NodeData,Real,VertexData>::exportVTKData(const float scale, const Point3D<float> translate)
{
	typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matrix;
	typedef Eigen::Matrix<float,Eigen::Dynamic,1> Vector;

	Eigen::Matrix<float,3,3> B;
	B<<1,-2,1,1,2,-2,0,0,1;
	B/=2;

#ifndef _DEBUG
	const int dim[3] = {256,256,256};
#else
	const int dim[3] = {8,8,8};
#endif
	Vector volume(dim[0]*dim[1]*dim[2]);
	for(int zi = 0; zi < dim[2]; zi++) 
		for(int yi = 0; yi < dim[1]; yi++) 
			for (int xi = 0; xi < dim[0]; xi++) 
			{
				float pos[3]={(float)xi/dim[0]/2+0.25,(float)yi/dim[1]/2,(float)zi/dim[2]/2+0.25};
				std::vector<int> coeffIndices;
				std::vector<float> coeffWeights;
				float value;
				getCoeffIndicesWeightsValueFromPos(B,1,maxBsplineDepth,pos,coeffIndices,coeffWeights,value);
				int index=xi+yi*dim[0]+zi*dim[0]*dim[1];
				volume(index)=value;
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

	// Get the values at the leaf nodes and propogate up to the parents
	OctNode<NodeData,Real>::NodeIndex nIdx;
	temp=tree.nextLeaf(NULL,nIdx);
	while(temp)
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

		OctNode<NodeData,Real>* temp2=temp;
		OctNode<NodeData,Real>::NodeIndex nIdx2=nIdx;
		temp=tree.nextLeaf(temp,nIdx);
		if(MarchingCubes::HasRoots(temp2->nodeData.mcIndex) && nIdx2.depth<maxDepth)
			setChildren(temp2,nIdx2,isoValue,useFull);
	}
}

template<class NodeData,class Real,class VertexData>
void OctreeBspline<NodeData,Real,VertexData>::setChildren(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const Real& isoValue,const int& useFull)
{
	if(nIdx.depth==maxDepth)
	{
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
	}
	else
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
				dist=this->eval(ctr.coords);
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
				dist=this->eval(ctr.coords);
				cornerValues[key]=VertexData(dist,n);
			}
		}

		for(int i=0;i<Cube::CORNERS;i++)
			setChildren(&node->children[i],nIdx.child(i),isoValue,useFull);
	}
}