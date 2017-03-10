template<typename Real>
KDTree<Real>::KDTree(): Index(flann::KDTreeSingleIndexParams(10)) {}

template<typename Real>
void KDTree<Real>::setInputPoints(const std::vector<Point3D<Real> >& points)
{
	size_t npts=points.size();
	assert(npts>0);
	const flann::Matrix<Real> pts((Real*)points.data(), npts, 3);
	buildIndex(pts);
}

template<typename Real>
int KDTree<Real>::KnnSearch(const std::vector<Point3D<Real> >& queries, 
	std::vector<std::vector<int> >& indices, 
	std::vector<std::vector<Real> >& dists, 
	const size_t knn)
{
	size_t nqrs=queries.size();
	const flann::Matrix<Real> qrs((Real*)queries.data(),nqrs,3);
	flann::SearchParams searchParams(32,0,true);
	searchParams.cores=0;
	return knnSearch(qrs,indices,dists,knn,searchParams);
}

