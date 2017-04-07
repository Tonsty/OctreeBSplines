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
/////////////
// Point3D //
/////////////
template<class Real>
Real& Point3D<Real>::operator[] (const int& idx)
{
	return coords[idx];
}
template<class Real>
const Real& Point3D<Real>::operator[] (const int& idx) const
{
	return coords[idx];
}
template<class Real>
Point3D<Real> Point3D<Real>::operator + (const Point3D<Real>& p) const
{
	Point3D<Real> q;
	q.coords[0]=coords[0]+p.coords[0];
	q.coords[1]=coords[1]+p.coords[1];
	q.coords[2]=coords[2]+p.coords[2];
	return q;
}
template<class Real>
Point3D<Real> Point3D<Real>::operator - (const Point3D<Real>& p) const
{
	Point3D<Real> q;
	q.coords[0]=coords[0]-p.coords[0];
	q.coords[1]=coords[1]-p.coords[1];
	q.coords[2]=coords[2]-p.coords[2];
	return q;
}
template<class Real>
Point3D<Real> Point3D<Real>::operator * (const Real& s) const
{
	Point3D<Real> q;
	q.coords[0]=coords[0]*s;
	q.coords[1]=coords[1]*s;
	q.coords[2]=coords[2]*s;
	return q;
}
template<class Real>
Point3D<Real> Point3D<Real>::operator / (const Real& s) const
{
	Point3D<Real> q;
	q.coords[0]=coords[0]/s;
	q.coords[1]=coords[1]/s;
	q.coords[2]=coords[2]/s;
	return q;
}
template<class Real>
Point3D<Real>& Point3D<Real>::operator += (const Point3D<Real>& p)
{
	coords[0]+=p.coords[0];
	coords[1]+=p.coords[1];
	coords[2]+=p.coords[2];
	return *this;
}
template<class Real>
Point3D<Real>& Point3D<Real>::operator -= (const Point3D<Real>& p)
{
	coords[0]-=p.coords[0];
	coords[1]-=p.coords[1];
	coords[2]-=p.coords[2];
	return *this;
}
template<class Real>
Point3D<Real>& Point3D<Real>::operator *= (const Real& s)
{
	coords[0]*=s;
	coords[1]*=s;
	coords[2]*=s;
	return *this;
}
template<class Real>
Point3D<Real>& Point3D<Real>::operator /= (const Real& s)
{
	coords[0]/=s;
	coords[1]/=s;
	coords[2]/=s;
	return *this;
}
template<class Real>
Real Random(void){return Real(rand())/RAND_MAX;}

template<class Real>
Point3D<Real> RandomBallPoint(void){
	Point3D<Real> p;
	while(1){
		p.coords[0]=Real(1.0-2.0*Random<Real>());
		p.coords[1]=Real(1.0-2.0*Random<Real>());
		p.coords[2]=Real(1.0-2.0*Random<Real>());
		Real l=SquareLength(p);
		if(l<=1){return p;}
	}
}
template<class Real>
Point3D<Real> RandomSpherePoint(void){
	Point3D<Real> p=RandomBallPoint<Real>();
	return p/Length(p);
}

template<class Real>
Real SquareLength(const Point3D<Real>& p)
{
	return DotProduct(p,p);
}

template<class Real>
Real DotProduct(const Point3D<Real>& p,const Point3D<Real>& q)
{
	return p.coords[0]*q.coords[0]+p.coords[1]*q.coords[1]+p.coords[2]*q.coords[2];
}

template<class Real>
Real Length(const Point3D<Real>& p)
{
	return Real(sqrt(SquareLength(p)));
}

template<class Real>
Real SquareDistance(const Point3D<Real>& p1,const Point3D<Real>& p2)
{
	return SquareLength(p1-p2);
}

template<class Real>
Real Distance(const Point3D<Real>& p1,const Point3D<Real>& p2){return Real(sqrt(SquareDistance(p1,p2)));}

template <class Real>
void CrossProduct(const Point3D<Real>& p1,const Point3D<Real>& p2,Point3D<Real>& p){
	p.coords[0]= p1.coords[1]*p2.coords[2]-p1.coords[2]*p2.coords[1];
	p.coords[1]=-p1.coords[0]*p2.coords[2]+p1.coords[2]*p2.coords[0];
	p.coords[2]= p1.coords[0]*p2.coords[1]-p1.coords[1]*p2.coords[0];
}
template <class Real>
Point3D<Real> Normal(const Point3D<Real>& p1,const Point3D<Real>& p2,const Point3D<Real>& p3)
{
	Point3D<Real> q1,q2,n;
	q1=p2-p1;
	q2=p3-p1;
	CrossProduct(q1,q2,n);
	return n;
}
template <class Real>
Real Area(const Point3D<Real>& p1,const Point3D<Real>& p2,const Point3D<Real>& p3)
{
	Point3D<Real> p12=p2-p1,p13=p3-p1;
	Point3D<Real> tmp;
	CrossProduct(p12,p13,tmp);
	return Length(tmp);
}

template <class Real>
Real DistanceToEdge(const Point3D<Real>& p,const Point3D<Real> e[2])
{
	Point3D<Real> q,v;
	Real dot;
	q=p-e[0];
	v=e[1]-e[0];
	dot=DotProduct(q,v);
	if(dot<=0)
		return Length(q);
	else if (dot>SquareLength(v))
		return Distance(p,e[1]);
	else
	{
		Real t=dot/SquareLength(v);
		v=e[0]*(Real(1.0)-t)+e[1]*t;
		return Distance(p,v);
	}
}
template <class Real>
Real DistanceToTriangle(const Point3D<Real>& p,const Point3D<Real> t[3])
{
	Point3D<Real> e[2];
	Point3D<Real> q,v,n,nn;
	
	n=Normal(t[0],t[1],t[2]);
	for(int i=0;i<3;i++)
	{
		v=t[(i+1)%3]-t[i];
		q=p-t[i];
		CrossProduct(n,v,nn);
		if(DotProduct(q,nn)<=0)
		{
			e[0]=t[i];
			e[1]=t[(i+1)%3];
			return DistanceToEdge(p,e);
		}
	}
	return fabs(DotProduct(q,n))/Length(n);
}
template <class Real>
Point3D<Real> NearestPointOnEdge(const Point3D<Real>& p,const Point3D<Real> e[2],int& vFlag)
{
	Point3D<Real> q,v;
	Real dot;

	q=p-e[0];
	v=e[1]-e[0];

	dot=DotProduct(q,v);
	if(dot<=0)
	{
		vFlag=1;
		return e[0];
	}
	else if (dot>=SquareLength(v))
	{
		vFlag=2;
		return e[1];
	}
	else
	{
		Real t=dot/Real(SquareLength(v));
		v=e[0]*(Real(1.0)-t)+e[1]*t;
		vFlag=3;
		return v;
	}
}
template <class Real>
Point3D<Real> NearestPointOnTriangle(const Point3D<Real>& p,const Point3D<Real> t[3],int& vFlag)
{
	Point3D<Real> e[2];
	Point3D<Real> q,v,n,nn,nearest;
	vFlag=0;
	
	n=Normal(t[0],t[1],t[2]);
	for(int i=0;i<3;i++)
	{
		v=t[(i+1)%3]-t[i];
		q=p-t[i];

		CrossProduct(n,v,nn);
		if(DotProduct(q,nn)<=0)
		{
			int tempFlag;
			e[0]=t[i];
			e[1]=t[(i+1)%3];
			nearest=NearestPointOnEdge(p,e,tempFlag);
			if(tempFlag&1) vFlag|=1<<i;
			if(tempFlag&2) vFlag|=1<<((i+1)%3);
			return nearest;
		}
	}

	n/=Length(n);
	v=p-n*DotProduct(q,n);
	vFlag=7;
	return v;
}
template <class Real>
Point3D<Real> NearestPointOnPlane(const Point3D<Real>& p,const Point3D<Real> t,const Point3D<Real> n,int& vFlag)
{
	Point3D<Real> q,nearest;
	vFlag=0;

	q=p-t;
	nearest=p-n*DotProduct(q,n);
	return nearest;
}

template <class Real>
int OutCode(const Point3D<Real>& ctr,const Real& w,const Point3D<Real>& p)
{
	int oc=0;
	if(p.coords[0]<ctr.coords[0]-w/2)
		oc|=1;
	if(p.coords[0]>ctr.coords[0]+w/2)
		oc|=2;
	if(p.coords[1]<ctr.coords[1]-w/2)
		oc|=4;
	if(p.coords[1]>ctr.coords[1]+w/2)
		oc|=8;
	if(p.coords[2]<ctr.coords[2]-w/2)
		oc|=16;
	if(p.coords[2]>ctr.coords[2]+w/2)
		oc|=32;
	return oc;
}
template <class Real>
int PointInCube(const Point3D<Real>& ctr,const Real& w,const Point3D<Real>& p)
{
	return !OutCode(ctr,w,p);
}
template <class Real>
int EdgeInCube(const Point3D<Real>& ctr,const Real& w,const Point3D<Real> e[2])
{
	int oc[2],dir,off;
	Real t,x;
	oc[0]=OutCode(ctr,w,e[0]);
	oc[1]=OutCode(ctr,w,e[1]);
	if(!oc[0] || !oc[1])	return 1;
	if(oc[0] & oc[1])		return 0;
#if 1
	for(dir=0;dir<3;dir++)
		if((oc[0]>>(dir<<1))&3)
		{
			off=( (oc[0]>>(dir<<1))&2) >> 1;
			t=( e[0][dir]-(ctr[dir] - w/2 + w*off) ) / (e[0][dir]-e[1][dir]);
			int inside=0;
			for(int i=1;i<3;i++)
			{
				int j=(dir+i)%3;
				x=e[0][j]*(Real(1.0)-t)+e[1][j]*t;
				if(x>=(ctr[j]-w/2) && x<=(ctr[j]+w/2))
					inside++;
			}
			if(inside==2)
				return 1;
		}
	return 0;
#else
	for(dir=0;dir<3;dir++)
		if((oc[0]>>(dir<<1))&3)
			break;
	off=( (oc[0]>>(dir<<1))&2) >> 1;
	t=( e[0][dir]-(ctr[dir] -w/2 + w*off) ) / (e[0][dir]-e[1][dir]);
	for(int i=1;i<3;i++)
	{
		int j=(dir+i)%3;
		x=e[0][j]*(Real(1.0)-t)+e[1][j]*t;
		if(x<(ctr[j]-w/2) || x>(ctr[j]+w/2))
			return 0;
	}
	return 1;
#endif
}
template <class Real>
int TriangleInCube(const Point3D<Real>& ctr,const Real& w,const Point3D<Real> t[3])
{
	Point3D<Real> e[2],n,nn[3];
	int oc[3];
	oc[0]=OutCode(ctr,w,t[0]);
	oc[1]=OutCode(ctr,w,t[1]);
	oc[2]=OutCode(ctr,w,t[2]);
	if(!oc[0] || !oc[1] || !oc[2])
		return 1;
	if(oc[0] & oc[1] & oc[2])
		return 0;

	for(int i=0;i<3;i++)
	{
		e[0]=t[i];
		e[1]=t[(i+1)%3];
		if(EdgeInCube(ctr,w,e))
			return 1;
	}

	n=Normal(t[0],t[1],t[2]);
	for(int i=0;i<3;i++)
		CrossProduct(n,t[(i+1)%3]-t[i],nn[i]);

	for(int i=0;i<Cube::EDGES;i++)
	{
		int c1,c2,x[3];
		Real tt,dot[2];
		Point3D<Real> p;
		Cube::EdgeCorners(i,c1,c2);
		Cube::FactorCornerIndex(c1,x[0],x[1],x[2]);
		for(int j=0;j<3;j++)	e[0][j]=ctr[j]-w/2+w*x[j];

		Cube::FactorCornerIndex(c2,x[0],x[1],x[2]);
		for(int j=0;j<3;j++)	e[1][j]=ctr[j]-w/2+w*x[j];

		dot[0]=DotProduct(n,e[0]-t[0]);
		dot[1]=DotProduct(n,e[1]-t[0]);
		if(dot[0]*dot[1] >=0 )	continue;
		tt=dot[0]/(dot[0]-dot[1]);
		p=e[0]*(Real(1.0)-tt)+e[1]*tt;
		if(DotProduct(p-t[0],nn[0])>0 && DotProduct(p-t[1],nn[1])>0 && DotProduct(p-t[2],nn[2])>0 )
			return 1;
	}
	return 0;
}


