template<typename Real>
KDTree<Real>::KDTree(): Index(flann::KDTreeSingleIndexParams(10)) {}

template<typename Real>
void KDTree<Real>::setInputPoints(Real* points, size_t npts)
{
	assert(npts>0);
	const flann::Matrix<Real> pts(points,npts,3);
	buildIndex(pts);
}

template<typename Real>
int KDTree<Real>::KnnSearch(Real* queries, size_t nqrs,
	std::vector<std::vector<int> >& indices, 
	std::vector<std::vector<Real> >& dists, 
	const size_t knn) const
{
	const flann::Matrix<Real> qrs(queries,nqrs,3);
	flann::SearchParams searchParams(32,0,true);
	searchParams.cores=0;
	return knnSearch(qrs,indices,dists,knn,searchParams);
}

