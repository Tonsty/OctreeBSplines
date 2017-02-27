/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/
#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED
#include <math.h>
#include <vector>

template<class Real>
Real Random(void);

template<class Real>
struct Point3D{
	Real coords[3];

	Real& operator[] (const int& idx);
	const Real& operator[] (const int& idx) const;
	Point3D operator + (const Point3D& p) const;
	Point3D operator - (const Point3D& p) const;
	Point3D operator * (const Real& s) const;
	Point3D operator / (const Real& s) const;
	Point3D& operator += (const Point3D& p);
	Point3D& operator -= (const Point3D& p);
	Point3D& operator *= (const Real& s);
	Point3D& operator /= (const Real& s);
};

template<class Real>
Point3D<Real> RandomBallPoint(void);

template<class Real>
Point3D<Real> RandomSpherePoint(void);

template<class Real>
Real Length(const Point3D<Real>& p);

template<class Real>
Real SquareLength(const Point3D<Real>& p);

template<class Real>
Real DotProduct(const Point3D<Real>& p,const Point3D<Real>& q);

template<class Real>
Real Distance(const Point3D<Real>& p1,const Point3D<Real>& p2);

template<class Real>
Real SquareDistance(const Point3D<Real>& p1,const Point3D<Real>& p2);

template <class Real>
void CrossProduct(const Point3D<Real>& p1,const Point3D<Real>& p2,Point3D<Real>& p);

template <class Real>
Point3D<Real> Normal(const Point3D<Real>& p1,const Point3D<Real>& p2,const Point3D<Real>& p3);

template <class Real>
Real DistanceToEdge(const Point3D<Real>& p,const Point3D<Real> e[2]);
template <class Real>
Real DistanceToTriangle(const Point3D<Real>& p,const Point3D<Real> t[3]);

template <class Real>
Point3D<Real> NearestPointOnEdge(const Point3D<Real>& p,const Point3D<Real> e[2],int& vFlag);
template <class Real>
Point3D<Real> NearestPointOnTriangle(const Point3D<Real>& p,const Point3D<Real> t[3],int& vFlag);
template <class Real>
Point3D<Real> NearestPointOnPlane(const Point3D<Real>& p,const Point3D<Real> t,const Point3D<Real> n,int& vFlag);

template <class Real>
int OutCode(const Point3D<Real>& ctr,const Real& w,const Point3D<Real>& p);

template <class Real>
int PointInCube(const Point3D<Real>& ctr,const Real& w,const Point3D<Real>& p);
template <class Real>
int EdgeInCube(const Point3D<Real>& ctr,const Real& w,const Point3D<Real> e[2]);
template <class Real>
int TriangleInCube(const Point3D<Real>& ctr,const Real& w,const Point3D<Real> t[3]);

class TriangleIndex{
public:
	size_t idx[3];
};

#include "Geometry.inl"

#endif // GEOMETRY_INCLUDED
