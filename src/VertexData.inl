/*
Copyright (c) 2007, Michael Kazhdan
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
/////////////////
// VertexValue //
/////////////////
template<class Real>
VertexValue<Real>::VertexValue(void)
{
	v=0;
}
template<class Real>
VertexValue<Real>::VertexValue(const Real& v,const Point3D<Real>& n)
{
	this->v=v;
}
template<class Real>
Real VertexValue<Real>::value(void) const
{
	return v;
}
template<class Real>
VertexValue<Real> VertexValue<Real>::operator + (const VertexValue<Real>& v) const
{
	VertexValue<Real> vv;
	vv.v=this->v+v.v;
	return vv;
}
template<class Real>
VertexValue<Real> VertexValue<Real>::operator - (const VertexValue<Real>& v) const
{
	VertexValue<Real> vv;
	vv.v=this->v-v.v;
	return vv;
}
template<class Real>
VertexValue<Real> VertexValue<Real>::operator * (const Real& v) const
{
	VertexValue<Real> vv;
	vv.v=this->v*v;
	return vv;
}
template<class Real>
VertexValue<Real> VertexValue<Real>::operator / (const Real& v) const
{
	VertexValue<Real> vv;
	vv.v=this->v/v;
	return vv;
}
template<class Real>
Point3D<Real> VertexValue<Real>::RootPosition(const Real& isoValue,const Point3D<Real>& p1,const Point3D<Real>& p2,VertexValue v1,VertexValue v2)
{
	Real t=(v1.v-isoValue)/(v1.v-v2.v);
	return p1*(Real(1.0)-t)+p2*t;
}

