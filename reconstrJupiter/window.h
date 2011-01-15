// window : this program reads in a 3-D datafile and applies a windowing

/*   This file is part of a software package written by 
     Rainer Heintzmann
     Institute of Applied Optics and Information Processing
     Albert Ueberle Strasse 3-5
     69120 Heidelberg
     Tel.: ++49 (0) 6221  549264
     No garantee, whatsoever is given about functionallaty and  savety.
     No warranty is taken for any kind of damage it may cause.
     No support for it is provided !

     THIS IS NOT FREE SOFTWARE. All rights are reserved, and the software is only
     given to a limited number of persons for evaluation purpose !
     Please do not modify and/or redistribute this file or any other files, libraries
     object files or executables of this software package !
*/


// used as : 

// Window WinObj(SizeX/2,SizeY/2,SizeZ/2,0.0,1.0);
// image.apply(WinObj.window);
#include "stlfixes.h"
#include <complex>


template<class ArrayType>
class Window {

private:
  float MidX,MidY,MidZ;
  float DistX,DistY,DistZ;
  float OutX,OutY,OutZ; // = 1.0 ; outer border as relative image width
  float InX,InY,InZ; // = 0.0  ; inner border relative

public:
  typedef typename ArrayType::ArrayBaseType ArrayBType;

  float gaussstretch;
  bool rectborder;
  bool ring;

  void InitWindow(float IX,float OX,float IY, float OY, float IZ, float OZ) 
    {
      InX=IX;InY=IY;InZ=IZ;
      OutX=OX;OutY=OY;OutZ=OZ;
    }

  Window(float IX,float OX,float IY, float OY, float IZ, float OZ, bool rectb, bool rng) : MidX(0),MidY(0),MidZ(0),gaussstretch(1.0),rectborder(rectb),ring(rng)
    {
      InitWindow(IX,OX,IY,OY,IZ,OZ);
    }
  
  /* Window(float I, float O) : MidX(0),MidY(0),MidZ(0),gaussstretch(1.0),rectborder(0),ring(0)
    {
      InitWindow(I,O,I,O,I,O);
      };*/

  inline double CLIP(double val)
    {
      if (! ring)
	return (val > 0.0) ? ((val < 1.0) ? val : 1.0)  : 0.0;
      else
	return (val > 0.0) ? ((val < 1.0) ? ((val < 0.5) ? val*2.0 : (1.0-val)*2.0) : 0.0)  : 0.0;
    }

/// returns distance to border relative to outer and inner borders.
/// returning 1.0 inside and 0.0 outside
inline float BorderDist(int x, int y ,int z)
{
  float dx=x-MidX, dy=y-MidY, dz=z-MidZ,r,minr,maxr;   // Distances to the midpoint of coordinate system
  if (MidX==0.5) dx -=0.5;  // in case size in that dimension is == 1
  if (MidY==0.5) dy -=0.5;
  if (MidZ==0.5) dz -=0.5;

  if (dx > DistX) dx = 2.0*DistX-dx;   // simulates cyclic FFT type borders
  if (dy > DistY) dy = 2.0*DistY-dy;   // simulates cyclic FFT type borders
  if (dz > DistZ) dz = 2.0*DistZ-dz;   // simulates cyclic FFT type borders

  if (rectborder)
    {
      dx = (OutX- fabs(dx)/DistX) / (OutX-InX);  // rel distance to borders
      dy = (OutY- fabs(dy)/DistY) / (OutY-InY);
      dz = (OutZ- fabs(dz)/DistZ) / (OutZ-InZ);
      if (dx < dy)
	if (dx < dz)
	  return CLIP(dx);
	else
	  return CLIP(dz);
      else
	if (dy < dz)
	  return CLIP(dy);
	else
	  return CLIP(dz);
    }
  else
    {
      float Sx=(InX+OutX)/2.0;
      float Sy=(InY+OutY)/2.0;
      float Sz=(InZ+OutZ)/2.0;
      dx =fabs(dx)/DistX/Sx;  // coordinates normed to image
      dy =fabs(dy)/DistY/Sy;
      dz =fabs(dz)/DistZ/Sz;
      r = dx*dx+dy*dy+dz*dz;
      if (r <= 0)
	if (ring)
	  return 0.0;
	else
	  return (OutX>=InX) ? 1.0 : 0.0;

      r = sqrt(r); // normed radius inner ellypse == 1
      // minr = sqrt(dx*dx*InX*InX+dy*dy*InY*InY+dz*dz*InZ*InZ);

      // if (r == 0)
	// return 1.0;

      float x,y,z;
      x = dx*OutX/Sx/r;
      y = dy*OutY/Sy/r;
      z = dz*OutZ/Sz/r;
      maxr = sqrt(x*x+y*y+z*z); // radius of Vector projected to Outer Ellypse
      x = dx*InX/Sx/r;
      y = dy*InY/Sy/r;
      z = dz*InZ/Sz/r;
      minr = sqrt(x*x+y*y+z*z); // radius of Vector projected to Inner Ellypse

      // r = (maxr- r) / (maxr-minr);
      r = (maxr- r) / (maxr-minr);

      return CLIP(r); // distance to normed circle . Inner border is 0 , outer border is 1.0
    }
}

bool IsInUpperHalf(float x,float y, float z, float DimX, float DimY, float DimZ) {
    // The scalar product between the vector characterising the diagonal plane (through 0,0,0) and the position is calculated and compared to zero
    // The planar vector is obtained by the cross-product (vectorproduct) of two characterising directions (DX,DY,0) x (DX,0,DZ)
    float prod= x*(DimY*DimZ) - y *(DimX*DimZ) - z * (DimY *DimX);
    if (prod >= 0) return true;
    else return false;
  }

 bool FlipX(int & x, int DimX)
    {
      float dx=x-MidX;
      if (MidX==0.5) dx -=0.5;  // in case size in that dimension is == 1
      if (dx > 0.0) return false; // flip only one side!
      dx = (OutX- fabs(dx)/DistX) / (OutX-InX);  // rel distance to borders
      if (CLIP(dx) > 0 && CLIP(dx) < 1.0)
	{
	  x=  DimX-1-x;
	  return true;
	}
      else return false;
    }
  
 bool FlipY(int & y,int DimY)
    {
      float dy=y-MidY;
      if (MidY==0.5) dy -=0.5;  // in case size in that dimension is == 1
      if (dy > 0.0) return false; // flip only one side!
      dy = (OutY- fabs(dy)/DistY) / (OutY-InY);  // rel distance to borders
      if (CLIP(dy) > 0 && CLIP(dy) < 1.0)
	{
	  y= DimY-1-y;
	  return true;
	}
      else return false;
    }
  
 bool FlipZ(int & z, int DimZ)
    {
      float dz=z-MidZ;
      if (MidZ==0.5) dz -=0.5;  // in case size in that dimension is == 1
      if (dz > 0.0) return false; // flip only one side!
      dz = (OutZ- fabs(dz)/DistZ) / (OutZ-InZ);  // rel distance to borders
      if (CLIP(dz) > 0 && CLIP(dz) < 1.0)
	{
	  z= DimZ-1-z;
	  return true;
	}
      else return false;
    }
  
 void ApplyMirror(ArrayType * arr)
    {
      ArrayBType weight,val,aval;
      int DimX=arr->GetSize(0),DimY=arr->GetSize(1),DimZ=arr->GetSize(2);
      int fx,fy,fz;
      int x,y,z;
      for(z=0;z<DimZ;z++)
	for(y=0;y<DimY;y++)
	  for(x=0;x<DimX;x++)
	    {
	      fx=x,fy=y,fz=z;
	      if (FlipX(fx,DimX))
		{
		  weight = CalcWindow(x,DimY/2,DimZ/2);
		  val = arr->Value(x,y,z);
		  aval = arr->Value(fx,y,z);
		  arr-> SetValue(x,y,z,(ArrayBType(1.0)-weight)*aval + weight*val); 
		  arr-> SetValue(fx,fy,fz,(ArrayBType(1.0)-weight)*val +weight*aval);
		}
	      fx=x,fy=y,fz=z;
	      if (FlipY(fy,DimY))
		{
		  weight = CalcWindow(DimX/2,y,DimZ/2);
		  val = arr->Value(x,y,z);
		  aval = arr->Value(x,fy,z);
		  arr-> SetValue(x,y,z,(ArrayBType(1.0)-weight)*aval + weight*val); 
		  arr-> SetValue(fx,fy,fz,(ArrayBType(1.0)-weight)*val +weight*aval);
		}
	      fx=x,fy=y,fz=z;
	      if (FlipZ(fz,DimZ))
		{
		  weight = CalcWindow(DimX/2,DimY/2,z);
		  val = arr->Value(x,y,z);
		  aval = arr->Value(x,y,fz);
		  arr-> SetValue(x,y,z,(ArrayBType(1.0)-weight)*aval + weight*val); 
		  arr-> SetValue(fx,fy,fz,(ArrayBType(1.0)-weight)*val +weight*aval);
		}
	    }
    }
  
 virtual ArrayBType CalcWindow(int x, int y, int z)=0;  // will be overwritten in derived classes

   /*ArrayBType CalcWindow(int x, int y, int z)  // will be overwritten in derived classes
    {
      return 1.0;
      }*/

void ApplyWindow(ArrayType * arr, bool Center=true, bool mirror=false) 
    {
      int DimX=arr->GetSize(0),DimY=arr->GetSize(1),DimZ=arr->GetSize(2);
      if (Center)
	{
	  MidX=DimX/2.0,MidY=DimY/2.0,MidZ=DimZ/2.0;
	  DistX=DimX/2.0,DistY=DimY/2.0,DistZ=DimZ/2.0;
	}
      else
	{
	  MidX=0,MidY=0,MidZ=0;
	  DistX=DimX,DistY=DimY/2.0,DistZ=DimZ/2.0;  // Distance to Border
	}
      ArrayBType weight;
      if (mirror) 
	{
	  ApplyMirror(arr);
	  return;
	}
      else
	{
	  for(int z=0;z<DimZ;z++)
	    for(int y=0;y<DimY;y++)
	      for(int x=0;x<DimX;x++)
		{
		  weight = CalcWindow(x,y,z);
		  arr-> SetValue(x,y,z,weight*arr->Value(x,y,z)); 
		}
	}
    }
      
};

template <class ArrayType>
class NearGaussWindow : public Window<ArrayType> {
 public:
  NearGaussWindow(float IX,float OX,float IY, float OY, float IZ, float OZ, bool rectb, bool ring) : Window<ArrayType> (IX,OX,IY,OY,IZ,OZ,rectb,ring)
    {
    }
  typedef typename ArrayType::ArrayBaseType ArrayBType;
  virtual ArrayBType CalcWindow(int x, int y, int z)
    {
      float bd = (1.0-Window<ArrayType>::BorderDist(x,y,z));
      if (Window<ArrayType>::BorderDist(x,y,z) > 0.0)
	return (ArrayBType) float(exp(-(1.0/(1.0-bd*bd)-1.0)));
      else
	return 0.0;
    }
};

template <class ArrayType>
class GaussWindow : public Window<ArrayType> {
 public:
  GaussWindow(float IX,float OX,float IY, float OY, float IZ, float OZ, bool rectb, bool ring) : Window<ArrayType> (IX,OX,IY,OY,IZ,OZ,rectb,ring)
    { }
  typedef typename ArrayType::ArrayBaseType ArrayBType;
  virtual ArrayBType CalcWindow(int x, int y, int z)
    {
      float bd = (1.0-Window<ArrayType>::BorderDist(x,y,z)) / Window<ArrayType>::gaussstretch;
      return (ArrayBType) float(exp(-bd*bd));
    }
};

template <class ArrayType>
class EdgeWindow : public Window<ArrayType> {
 public:
  EdgeWindow(float IX,float OX,float IY, float OY, float IZ, float OZ, bool rectb, bool ring) : Window<ArrayType> (IX,OX,IY,OY,IZ,OZ,rectb,ring)
    { }
  typedef typename ArrayType::ArrayBaseType ArrayBType;
  virtual ArrayBType CalcWindow(int x, int y, int z)
    {
      if (Window<ArrayType>::BorderDist(x,y,z) > 0.0)
	return 1.0;
      else
	return 0;
    }
};

template <class ArrayType>
class LinWindow : public Window<ArrayType> {
 public:
  LinWindow(float IX,float OX,float IY, float OY, float IZ, float OZ, bool rectb, bool ring) : Window<ArrayType> (IX,OX,IY,OY,IZ,OZ,rectb,ring)
    { }
  typedef typename ArrayType::ArrayBaseType ArrayBType;
  virtual ArrayBType CalcWindow(int x, int y, int z)
    {
      return (ArrayBType) Window<ArrayType>::BorderDist(x,y,z);
    }
};

template <class ArrayType>
class SinWindow : public Window<ArrayType> {
 public:
  SinWindow(float IX,float OX,float IY, float OY, float IZ, float OZ, bool rectb, bool ring) : Window<ArrayType> (IX,OX,IY,OY,IZ,OZ,rectb,ring)
    { }
  typedef typename ArrayType::ArrayBaseType ArrayBType;
  virtual ArrayBType CalcWindow(int x, int y, int z)
    {
      return (ArrayBType) (1-cos(M_PI*Window<ArrayType>::BorderDist(x,y,z)))/2.0f;
    }
};

template <class ArrayType>
class TunnelWindow : public Window<ArrayType> {
 public:
  TunnelWindow(float IX,float OX,float IY, float OY, float IZ, float OZ, bool rectb, bool ring) : Window<ArrayType> (IX,OX,IY,OY,IZ,OZ,rectb,ring)
    { }

  typedef typename ArrayType::ArrayBaseType ArrayBType;
  virtual ArrayBType CalcWindow(int x, int y, int z)
    {
      return 1.0;
    }

void DoTunneling(ArrayType * arr, bool Center=true, bool mirror=false) {
  int x,y,z;
  int DimX=arr->GetSize(0),DimY=arr->GetSize(1),DimZ=arr->GetSize(2);
  if (Center)
    {
      Window<ArrayType>::MidX=DimX/2.0,Window<ArrayType>::MidY=DimY/2.0,Window<ArrayType>::MidZ=DimZ/2.0;
      Window<ArrayType>::DistX=DimX/2.0,Window<ArrayType>::DistY=DimY/2.0,Window<ArrayType>::DistZ=DimZ/2.0;
    }
  else
    {
      Window<ArrayType>::MidX=0,Window<ArrayType>::MidY=0,Window<ArrayType>::MidZ=0;
      Window<ArrayType>::DistX=DimX,Window<ArrayType>::DistY=DimY/2.0,Window<ArrayType>::DistZ=DimZ/2.0;  // Distance to Border
    }

  cout << "sqrt(-1) = " <<sqrt(complex<float>(-1,0)) << "\n";
  float k,eps=Window<ArrayType>::gaussstretch,d=Window<ArrayType>::OutX*Window<ArrayType>::MidX;

  ArrayBType dk,a;
  ArrayBType i(0,1);

  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
	{
	  k=M_PI*(x-Window<ArrayType>::MidX)/Window<ArrayType>::MidX;
	  dk=sqrt(ArrayBType(k*k-eps*eps));
	  if ((k != 0) && (dk != ArrayBType(0)))
	    a=exp(-i*k*d) / (-i*(k*k+dk*dk)/2.0f /k/dk * sin (dk*d)+cos(dk*d));
	  else
	    a=0.0;
	  if (y==0)
	    cout << "x="<<x<<"  k="<<k<<"  dk = " <<dk << "  a=" << a <<"\n";
	  if (mirror)
	    arr-> SetValue(x,y,z,((1.0-a)*arr->Value(DimX-x-1,DimY-y-1,DimZ-z-1) + a*arr->Value(x,y,z))/2.0); 
	  else
	    arr-> SetValue(x,y,z,a*arr->Value(x,y,z)); 
	}
}
};
