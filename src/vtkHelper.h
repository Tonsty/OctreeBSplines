#ifndef VTK_HELPER_INCLUDED
#define VTK_HELPER_INCLUDED

#include <Eigen/Dense>
#include <string>
#include "Geometry.h"

void sameVTIVolumeAndVTKMesh(const Eigen::Matrix<float,Eigen::Dynamic,1> &volume,const int dim[3],const std::string volumeFileName,
	const Eigen::Matrix<double,4,4> &transform,const std::string meshFileName);

void sameVTKOctree(const std::vector<Point3D<float> > &vertices, const std::vector<std::pair<int,int> > &edges, const std::string octreeFileName);

#endif