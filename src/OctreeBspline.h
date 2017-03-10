#ifndef OCTREE_BSPLINE_INCLUDED
#define OCTREE_BSPLINE_INCLUDED

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <hash_set>
#include "IsoOctree.h"
#include "Function.h"

//Hierarchical Implicit Surface Reconstruction based on Adaptive Distance Field

//This paper presents an hierarchical implicit surface reconstruction method. 
//Our method begins by establishing an adaptive signed distance field from the input data using an octree. 
//the sampling density adapts to the local curvature of data.
//In this way, varying-scale geometric details are automatically preserved. 
//Then, the signed distance field is fitted by a hierarhcial tensor-product B-splines, 
//of which the basic spline 

template<class NodeData,class Real,class VertexData>
class OctreeBspline: public IsoOctree<NodeData,Real,VertexData>, public Function
{
	template<class MeshReal>
	void setDistanceAndNormal2(const std::vector<int>& vindices,MeshInfo<MeshReal>& mInfo,const Point3D<Real>& p,Real& v,Point3D<Real>& n);
	template<class MeshReal>
	int setChildren2(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
		const std::vector<int>& vindices,MeshInfo<MeshReal>& mInfo,const int& maxDepth,const int& setCenter,const Real& flatness,
		stdext::hash_map<long long,std::vector<int>*>* triangleMap=NULL);

	void setDistanceAndNormal3(const Point3D<Real>& p,const Point3D<Real>& p2,const Point3D<Real>& n2,Real& dist,Point3D<Real>& n);
	template<class MeshReal>
	int setChildren3(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,
		const std::vector<int>& vindices,MeshInfo<MeshReal>& mInfo,const int& maxDepth,const int& setCenter,const Real& flatness,
		stdext::hash_map<long long,std::vector<int>*>* triangleMap=NULL);

	MeshInfo<float> mInfoGlobal;

	void setCoeffValuesFromCompressedKeys(const int bsplineDepth,const std::vector<std::vector<long long> >& compressedKeys,
		std::vector<stdext::hash_map<long long,std::pair<int,Real> > >& coeffValues);
	void getCornerKeysFromCompressedKeys(const int maxDepth,const int bsplineDepth,const std::vector<std::vector<long long> >& compressedKeys,
		stdext::hash_set<long long>& cornerKeys);

	inline void getPosFromCornerKey(const long long cornerKey,const float unitLen,float pos[3]);
	inline void getBposIntU(const int bsplineDepth,const float pos[3],int bposInt[3],float u[3]);
	void getCoeffIndicesWeightsValueFromPos(const Eigen::Matrix<float,3,3>& B,const int minBsplineDepth,const int maxBsplineDepth,const float pos[3], 
		std::vector<int>& coeffIndices,std::vector<float>& coeffWeights,float& value);

	void getSmoothMatrix(const Eigen::Matrix<float,3,3>& B,const Eigen::Matrix<float,3,3>& dB,const Eigen::Matrix<float,3,3>& ddB,
		const int minBsplineDepth,const int maxBsplineDepth, Eigen::SparseMatrix<float> &smoothMatrix);
	
	// Set children according to isoValue until reaching maxDepth
	void setMCLeafNodeToMaxDepth(OctNode<NodeData,Real>* node,const typename OctNode<NodeData,Real>::NodeIndex& nIdx,const Real& isoValue,const int& useFull);

public:
	// The maximum depth of the hierarchical bsplines. The value must be at least 1. (>=1.)
	int maxBsplineDepth;
	// A vector of hash-table of the multi-level bsplines (coeffKey: coeffIndex,coeffValue)
	std::vector<stdext::hash_map<long long,std::pair<int,Real> > > coeffValues;
	// A vector of vector of the compressed keys of the multi-level bsplines
	std::vector<std::vector<long long> > compressedKeys;

	template<class Vertex>
	int set2(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,const int& maxDepth,const int& setCenter,const Real& flatness,
		Point3D<Real>& translate,Real& scale,const int& noTransform);

	template<class Vertex>
	int set3(const std::vector<Vertex>& vertices,const std::vector<std::vector<int> >& polygons,const int& maxDepth,const int& setCenter,const Real& flatness,
		Point3D<Real>& translate,Real& scale,const int& noTransform);

	void directBsplineFitting();
	void multigridBsplineFitting();
	
	void exportVTKData(const float scale, const Point3D<float> translate); 
	virtual float eval(const float pos[3]);

	void setMCLeafNodeToMaxDepth(const Real& isoValue,const int& useFull);
};

#include "OctreeBSpline.inl"
#endif