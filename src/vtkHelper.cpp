#include <Eigen/Dense>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkMarchingCubes.h>
#include <vtkPLYWriter.h>
#include <vtkReverseSense.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

void sameVTKVolumeAndMesh(const Eigen::Matrix<float,Eigen::Dynamic,1> &volume,const int dim[3],const std::string volumeFileName,
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