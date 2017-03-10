#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include <flann/flann.hpp>

template<typename Real>
class KDTree: public flann::Index<flann::L2<Real> >
{
	using flann::Index<flann::L2<Real> >::knnSearch;
public:
	KDTree();
	void setInputPoints(const std::vector<Point3D<Real> >& points);
	int KnnSearch(const std::vector<Point3D<Real> >& queries, 
		std::vector<std::vector<int> >& indices, 
		std::vector<std::vector<Real> >& dists, 
		const size_t knn);

};

#include "Kdtree.inl"

#endif