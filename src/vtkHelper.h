#ifndef VTK_HELPER_INCLUDED
#define VTK_HELPER_INCLUDED

#include <Eigen/Dense>
#include <string>

void sameVTKVolumeAndMesh(const Eigen::Matrix<float,Eigen::Dynamic,1> &volume,const int dim[3],const std::string volumeFileName,
	const Eigen::Matrix<double,4,4> &transform,const std::string meshFileName);

#endif