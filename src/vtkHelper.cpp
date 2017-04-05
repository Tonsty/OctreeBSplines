#include <Eigen/Dense>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkMarchingCubes.h>
#include <vtkPLYWriter.h>
#include <vtkReverseSense.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include "Geometry.h"

void sameVTIVolumeAndVTKMesh(const Eigen::Matrix<float,Eigen::Dynamic,1> &volume,const int dim[3],const std::string volumeFileName,
	const Eigen::Matrix<double,4,4> &transform,const std::string meshFileName)
{
	vtkSmartPointer<vtkImageData> volumeVTK=vtkSmartPointer<vtkImageData>::New();
	volumeVTK->SetDimensions(dim[0],dim[1],dim[2]);
	volumeVTK->SetSpacing(1.0/dim[0],1.0/dim[1],1.0/dim[2]);
	volumeVTK->SetOrigin(0,0,0);
	volumeVTK->SetNumberOfScalarComponents(1);
	volumeVTK->SetScalarTypeToFloat();
	volumeVTK->AllocateScalars();
	float *ptr=(float*)volumeVTK->GetScalarPointer();
	memcpy(ptr,volume.data(),dim[0]*dim[1]*dim[2]*sizeof(float));

	vtkSmartPointer<vtkXMLImageDataWriter> imageWriter=vtkSmartPointer<vtkXMLImageDataWriter>::New();
	imageWriter->SetFileName(volumeFileName.c_str());
	imageWriter->SetInputConnection(volumeVTK->GetProducerPort());
	imageWriter->Write();

	vtkSmartPointer<vtkMarchingCubes> mc=vtkSmartPointer<vtkMarchingCubes>::New();
	mc->SetInput(volumeVTK);
	mc->ComputeNormalsOn();
	mc->SetValue(0,0.0);

	vtkSmartPointer<vtkReverseSense> reverse=vtkSmartPointer<vtkReverseSense>::New();
	reverse->SetInputConnection(mc->GetOutputPort());
	reverse->ReverseCellsOn();

	vtkSmartPointer<vtkTransform> transformVTK=vtkSmartPointer<vtkTransform>::New();
	Eigen::Matrix<double,4,4> transform_t=transform.transpose();
	transformVTK->SetMatrix(transform_t.data());

	vtkSmartPointer<vtkTransformPolyDataFilter> filter=vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	filter->SetTransform(transformVTK);
	filter->SetInputConnection(reverse->GetOutputPort());

	vtkSmartPointer<vtkPLYWriter> plyWriter=vtkSmartPointer<vtkPLYWriter>::New();
	plyWriter->SetFileName(meshFileName.c_str());
	plyWriter->SetInputConnection(filter->GetOutputPort());
	plyWriter->Write();
}

void sameVTKOctree(const std::vector<Point3D<float> > &vertices, const std::vector<std::pair<int,int> > &edges, const std::string octreeFileName)
{
	// Open file
	FILE*vtkFile=fopen(octreeFileName.c_str(),"w");

	int npts=vertices.size();

	// Write the header information
	fprintf(vtkFile,"# vtk DataFile Version 3.0\n");
	fprintf(vtkFile,"vtk output\n");
	fprintf(vtkFile,"ASCII\n");
	fprintf(vtkFile,"DATASET POLYDATA\n");
	fprintf(vtkFile,"POINTS %d float\n",npts);

	// Iterate through the points
	for(int i=0;i<npts;++i) 
	{
		fprintf(vtkFile,"%f %f %f\n",vertices[i][0],vertices[i][1],vertices[i][2]);
	}

	int negs=edges.size();

	// Write lines
	fprintf(vtkFile,"\nLINES %d %d\n", negs, 3*negs);
	for(int i=0;i<negs;++i) fprintf(vtkFile,"2 %d %d\n", edges[i].first, edges[i].second);

	// Close file
	fclose(vtkFile);
}