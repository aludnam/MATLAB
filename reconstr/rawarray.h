#ifndef rawarray_h  		// -*- C++ -*- 
#define rawarray_h

// The strange adressing is not yet completely implemented !

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>

// STL includes
#include <iterator>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>     // for stl numerics like accumulate, ...

// STL fixes for missing things :
#include "stlfixes.h"
#include "vec.h"
#include "accut.h"   // defines types for accumulation operations
#include "helpcompl.h"  // helper funktions, like CastToReal

#include <sys/types.h>  // needed for getpid();
#include <unistd.h>     // needed for getpid();

using namespace std;

// STL extensions :
template <class InputIterator, class OutputValt>
OutputValt * ccopy(InputIterator first, InputIterator last, OutputValt * result) {
    while (first != last) *result++ = (OutputValt) *first++;
    return result;
}

#include "basearrays.h"      // basic array templates :  static and dynamic
#include "asi.h"
#include "lsm.h"
#include "khoros.h"
#define SQR(x) ((x)*(x))

const double MAXARRAYTYPENR=(1e6);  // Will be return in case of div by zero in ArgDivSelf !

const int     NoSel=0;
const int CenterSel=1;
const int   ClipSel=2;
const int InterpSel=4;

// The 3-D array is taken to be the most general case for now.

// Lower Dimension-Arrays are derived from it, so other templates work on them too, if TArray is required.

/// to make arbitrary numbers (complex or real) real

/// General 3d-array type for variable type ("TheArrayType" == static or dynamic). 

template<class ArrayBType> 
class TArray3d {  // This is a baseclass !
public:
  typedef ArrayBType ArrayBaseType;  // get the basis type of myself and make it accessable to the outside
  typedef typename Accu<ArrayBaseType> ::type ArrayAccuType;  // this type is used for all kinds of accumulations

protected:
DynamicArray<ArrayBaseType> TheArray;        // This class can (only) assume a flat (as dynamic) array representation.

   /// Recursive Function for Linear interpolation in dim dimensions
ArrayBaseType LInterp(CalcFloatType * D,int * Dl,int * Dh,const int dim)
{
  ArrayBaseType sum=0;
  CalcFloatType weight;
  int mDl,mDh;  

  if (dim==0) return TheArray.GetValue(Dl[0],Dl[1],Dl[2]);

  mDh=Dh[dim-1];
  mDl=Dl[dim-1];

  weight=(mDh-D[dim-1]);

  Dh[dim-1]=mDl;
  sum += (ArrayBaseType) weight*LInterp(D,Dl,Dh,dim-1);
  Dl[dim-1]=mDh;
  Dh[dim-1]=mDh;
  sum += (ArrayBaseType) (1.0-weight)*LInterp(D,Dl,Dh,dim-1);
  Dl[dim-1]=mDl;              // Back to initial state

  return (ArrayBaseType) sum;

}

  /// does a trilinear interpolation, clipping depends on the value of the "Selector".
ArrayBaseType DoInterp(const unsigned char Selector,CalcFloatType x,CalcFloatType y,CalcFloatType z) 
{
	CalcFloatType D[3];
	int Dl[3],Dh[3];

	Dh[0]=(Dl[0]=(int) (D[0]=x))+1;
	Dh[1]=(Dl[1]=(int) (D[1]=y))+1;
	Dh[2]=(Dl[2]=(int) (D[2]=z))+1;

    if (Selector & ClipSel)
      {
	D[0] = TheArray.DoClipX(D[0]);       // Maybe one Clipping (Dl) can be eliminated
        Dl[0]= TheArray.DoClipX(Dl[0]);
	Dh[0]= TheArray.DoClipX(Dh[0]);

	D[1] = TheArray.DoClipY(D[1]);
	Dl[1]= TheArray.DoClipY(Dl[1]);
	Dh[1]= TheArray.DoClipY(Dh[1]);

	D[2] = TheArray.DoClipZ(D[2]);
	Dl[2]= TheArray.DoClipZ(Dl[2]);
	Dh[2]= TheArray.DoClipZ(Dh[2]);
      }     // else Values are right

    return LInterp(D,Dl,Dh,3);
}
public:

/// empty constructor
TArray3d() : TheArray()
{
  // Clear(); // array is created with zeros
}

/// empty constructor
TArray3d(int SizeX,int SizeY=1,int SizeZ=1) : TheArray(SizeX,SizeY,SizeZ)
{
  // Clear(); // array is created with zeros
}

int GetDimensions(void)
{
  return TheArray.Dimension();
}

template <class AnyType>
bool SizesEqual(TArray3d<AnyType> * Other)
{
  if (GetDimensions() != Other->GetDimensions())
   {
     cerr << "Number of Dimensions wrong!";
     return false;
   }

  for (int i=0;i < GetDimensions();i++)
    if (GetSize(i) != Other->GetSize(i))
    {
      cerr << GetSize(i) << " is not equal to " << Other->GetSize(i) << "\n";
      return false;
    }
    else
      cerr << GetSize(i) << " is correct \n";

  return true;
}

/// rezises array to new sizes (will be cleared)
void Resize(int SizeX,int SizeY=1,int SizeZ=1)
{
  TheArray.Resize(SizeX,SizeY,SizeZ);
}

/// rezises array to same as *AP (this array will be cleared)
void Resize(TArray3d * Ap)
{
  TheArray.Resize(& Ap->TheArray);
}

TArray3d(ArrayBaseType * ExistingArray,int SizeX,int SizeY,int SizeZ) : TheArray(ExistingArray,SizeX,SizeY,SizeZ)
{}

/// Erases array. (Set to 0)
void Clear()
{
  Set((ArrayBaseType) 0.0);
}

/// copy from raw array, e.g. for extracting a slice
void Copy(const ArrayBaseType * pcpy)
{
  // transform(TheArray.begin(),TheArray.end(),pcpy,TheArray.begin(),project2nd<ArrayBaseType,ArrayBaseType>());
  copy_n(pcpy,TheArray.size(),TheArray.begin());
}

/// copy into raw array, e.g. for extracting a slice
void CopyTo(ArrayBaseType * pcpy)
{
  copy(TheArray.begin(),TheArray.end(),pcpy);
}

/// FLips the data about a specified axis
void Flip(bool flipx, bool flipy=false, bool flipz=false)
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z,x2,y2,z2;
  ArrayBaseType val,val2;
  int SizeX=DimX,SizeY=DimY,SizeZ=DimZ;

  if (flipx)
    SizeX=DimX/2;
  else if (flipy)   // important: Only one flipped size can be reduced
    SizeY=DimY/2;
  else if (flipz)
    SizeZ=DimZ/2;

  for (z=0; z < SizeZ; z++)
      for (y=0; y < SizeY; y++)
	for (x=0; x < SizeX; x++)
	  {
	    x2=x;y2=y;z2=z;
	    if (flipx)
	      x2 = DimX-x-1;
	    if (flipy)
	      y2 = DimY-y-1;
	    if (flipz)
	      z2 = DimZ-z-1;
	    val2=Value(x2,y2,z2);
	    val = Value(x,y,z);
	    SetValue(x,y,z,val2);
	    SetValue(x2,y2,z2,val);
	  }
}


/// copy from array * cpy
void Copy(TArray3d * cpy)
  {
    if (cpy->TheArray.end() - cpy->TheArray.begin() < TheArray.end() - TheArray.begin())
      copy(cpy->TheArray.begin(),cpy->TheArray.end(),
	   TheArray.begin());
    else
      copy(cpy->TheArray.begin(),cpy->TheArray.begin() + (TheArray.end() - TheArray.begin()),
	   TheArray.begin());
  }

void ResampleFrom(TArray3d * from, int BinX=1,int BinY=1, int BinZ=1)
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z,rx,ry,rz,BinsX,BinsY,BinsZ;
  ArrayBaseType sum;

  if (BinX < 1) BinsX=1,BinX=-BinX;
  else BinsX=BinX;
  if (BinY < 1) BinsY=1,BinY=-BinY;
  else BinsY=BinY;
  if (BinZ < 1) BinsZ=1,BinZ=-BinZ;
  else BinsZ=BinZ;

  for (z=0; z < DimZ; z++)
    for (y=0; y < DimY; y++)
      for (x=0; x < DimX; x++)
	{
	  sum=0;
	for (rz=0; rz < BinsZ; rz++)
	  for (ry=0; ry < BinsY; ry++)
	    for (rx=0; rx < BinsX; rx++)
	      {
		sum += from->Value(x*BinX+rx,y*BinY+ry,z*BinZ+rz);
	      }
	SetValue(x,y,z,sum);
	}

}

  /// overloaded [ operator returning a reference
ArrayBaseType & operator() (const int x,const int y=0,const int z=0)
{
  return TheArray(x,y,z);
}

/// returns pointer to position (x,y,z) in array.
ArrayBaseType * Pointer(const int x=0,const int y=0,const int z=0)
{
  return TheArray.GetPValue(x,y,z);
}

/// returns value at (x,y,z) in array.
ArrayBaseType Value(const int x,const int y=0,const int z=0)
{
  return TheArray.GetValue(x,y,z);
}

/** returns value near floats (x,y,z) in array. The results depend on the used "Selector" which sais, if a clipping will be done
  (ClipSel), interpolation will be used (InterpSel), and/or the array will be used in centered coordinates (CenterSel) (center is the
  coordinates (SizeX div 2, SizeY div 2, SizeZ div 2). Most of the switches involved (many if's) are eliminated by the optimizer ! */

ArrayBaseType Value(const unsigned char Selector,const CalcFloatType x,const CalcFloatType y=0,const CalcFloatType z=0)
  {
    if (Selector & CenterSel)            // Hopefully the compiler will do a good optimizing job
      if (Selector & InterpSel)
        return DoInterp(Selector,TheArray.DoCenterX(x),TheArray.DoCenterY(y),TheArray.DoCenterZ(z));
      else 
	if (Selector & ClipSel)
	  return TheArray.GetValue(TheArray.DoClipX(RINT(TheArray.DoCenterX(x))),TheArray.DoClipY(RINT(TheArray.DoCenterY(y))),TheArray.DoClipZ(RINT(TheArray.DoCenterZ(z))));
        else
	  return TheArray.GetValue(RINT(TheArray.DoCenterX(x)),RINT(TheArray.DoCenterY(y)),RINT(TheArray.DoCenterZ(z)));
    else
      if (Selector & InterpSel)
        return DoInterp(Selector,x,y,z);
      else 
	if (Selector & ClipSel)
	  return TheArray.GetValue(TheArray.DoClipX(RINT(x)),TheArray.DoClipY(RINT(y)),TheArray.DoClipZ(RINT(z)));
        else
	  return TheArray.GetValue(RINT(x),RINT(y),RINT(z));
  }

/// returns value in the array at vectorial position (*pPos1 - *Ppos2). For "Selector" see other value function.
ArrayBaseType Value(const unsigned char Selector,Vector * pPos1, Vector * pPos2) 
  {
    return Value(Selector,pPos1->comp(0)-pPos2->comp(0),pPos1->comp(1)-pPos2->comp(1),pPos1->comp(2)-pPos2->comp(2));
  }

/// returns value in the array at vectorial position fac*(*pPos1 - *Ppos2). For "Selector" see other value function.
ArrayBaseType Value(const unsigned char Selector,Vector * pPos1, Vector * pPos2, CalcFloatType fac) 
  {
    return Value(Selector,fac*(pPos1->comp(0)-pPos2->comp(0)),fac*(pPos1->comp(1)-pPos2->comp(1)),fac*(pPos1->comp(2)-pPos2->comp(2)));
  }

/// returns value in the array at vectorial position (*pPos1). For "Selector" see other value function.
ArrayBaseType Value(const unsigned char Selector,Vector * pPos1) 
  {
    return Value(Selector,pPos1->comp(0),pPos1->comp(1),pPos1->comp(2));
  }

/// returns size of array in the dimension dim (0==x, 1==y, 2==z) in units of Numbers inside (ArrayBaseType's)
int GetSize(int dim)
{
  return TheArray.DimVec[dim];
}

/// calculates the complete sum of the array
ArrayAccuType  Integral(void)
{
  // return accumulate(TheArray.begin(),TheArray.end(), ArrayAccuType(0));
  ArrayAccuType Sum = 0;
  ArrayBaseType * myp1= TheArray.EndPos();
  ArrayBaseType maxVal = * myp1;
  while (myp1>=TheArray.StartPos())
      {
	  Sum += * myp1;
          myp1--;
      }
  return Sum;
}

ArrayAccuType Mean(void)
{
  return Integral() / ArrayAccuType(GetSize(0)) / ArrayAccuType(GetSize(1)) / ArrayAccuType(GetSize(2));  
}

/// does only count first n points
ArrayAccuType Mean(int num)
{
  return accumulate(TheArray.begin(),TheArray.begin()+num, ArrayAccuType(0))/num;
}

/// does only count pixels, that are in the range [clip0 .. clip1]
ArrayAccuType Mean(ArrayBaseType clip0, ArrayBaseType clip1)
{
  ArrayBaseType * myp1= TheArray.EndPos();
  ArrayAccuType integral=0;
  int num=0;

  while (myp1>=TheArray.StartPos())
    if ((CastToReal(*myp1) >= CastToReal(clip0)) && (CastToReal(*myp1) <= CastToReal(clip1)))
      {
	integral += * (myp1--);
	num++;
      }
    else
      myp1--;
 
  return integral/num;

}


/// returns maximal value in array
ArrayBaseType Maximum(void)
{
  ArrayBaseType * myp1= TheArray.EndPos();
  ArrayBaseType maxVal = * myp1;
  while (myp1>=TheArray.StartPos())
      {
                if (CastToReal(*myp1) >  CastToReal(maxVal))
                        maxVal= *myp1;
                myp1--;
      }
  return maxVal;
    // Used to work, but now does not compile:
   // return * max_element(TheArray.begin(),TheArray.end());
}

/// returns minimal value in array
ArrayBaseType Minimum(void)
{
  ArrayBaseType * myp1= TheArray.EndPos();
  ArrayBaseType minVal = * myp1;
  while (myp1>=TheArray.StartPos())
      {
                if (CastToReal(*myp1) < CastToReal(minVal))
                        minVal= *myp1;
                myp1--;
      }
  return minVal;
// Used to work, but now does not compile:
//  return * min_element(TheArray.begin(),TheArray.end());
}

/// normalizes Y-direction, resulting in a first value of one
void FirstYNormalize(ArrayAccuType val)  
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;
  ArrayBaseType norm;
  for (z=0; z < DimZ; z++)
  for (x=0; x < DimX; x++)
    {
      norm=1.0/Value(x,0,z)*val;
      for (y=0; y < DimY; y++)
	SetValue(x,y,z,Value(x,y,z)*norm);
    }
}

/// normalizes the sum to be "Integ" (if possible) along only one direction
void DirNormalize(ArrayAccuType Integ,int dir)  
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;
  ArrayBaseType sum;

  if (dir == 0)
  for (z=0; z < DimZ; z++)
  for (y=0; y < DimY; y++)
    {
      sum=0.0;
      for (x=0; x < DimX; x++)
	sum += Value(x,y,z);
      for (x=0; x < DimX; x++)
	SetValue(x,y,z,Value(x,y,z)/sum);
    }
  else if (dir == 1)
  for (z=0; z < DimZ; z++)
  for (x=0; x < DimX; x++)
    {
      sum=0.0;
      for (y=0; y < DimY; y++)
	sum += Value(x,y,z);
      for (y=0; y < DimY; y++)
	SetValue(x,y,z,Value(x,y,z)/sum);
    }
  else if (dir == 2)
  for (x=0; x < DimX; x++)
  for (y=0; y < DimY; y++)
    {
      sum=0.0;
      for (z=0; z < DimZ; z++)
	sum += Value(x,y,z);
      for (z=0; z < DimZ; z++)
	SetValue(x,y,z,Value(x,y,z)/sum);
    }
}

/// normalizes the sum to be "Integ" (if possible)
ArrayAccuType Normalize(ArrayAccuType Integ)    // normalize integral to Integ
{
  ArrayAccuType integral= Integral();

  ArrayBaseType * myp1= TheArray.EndPos();

 if (integral == (ArrayAccuType) 0.0) cerr << "ERROR Cannot normalize Array, Integral is Zero !\n",exit(-1);

  myp1= TheArray.EndPos();
  while (myp1>=TheArray.StartPos())
    {
      (* (myp1)) = (ArrayBaseType) ((* myp1) * Integ/integral);
      myp1--;
    }
  return Integ/integral;
}

/// normalizes the maximum to be "Max" (if possible)
void NormalizeToMax(CalcFloatType Max)  
{
  ArrayBaseType max= Maximum();

  ArrayBaseType * myp1= TheArray.EndPos();

 if (max == 0.0) cerr << "ERROR Cannot normalize Array, maximum is Zero !\n",exit(-1);

  myp1= TheArray.EndPos();
  while (myp1>=TheArray.StartPos())
    {
      (* (myp1)) = (ArrayBaseType) ((* myp1) * Max/max);
      myp1--;
    }
}

int CntBelowToMax(ArrayBaseType Val, TArray3d * ArgArr, ArrayBaseType NewVal)   // Counts voxel below threshold and sets other array at these positions to NewVal
{
  ArrayBaseType * myp1= TheArray.EndPos();
  ArrayBaseType * myp2= ArgArr->TheArray.EndPos();
  int res=0;
  while (myp1>=TheArray.StartPos())
      {
                if (CastToReal(*myp1) < Val)
			{
                        res++;
			(* myp2) = NewVal;
			}
                myp1--;
                myp2--;
      }
  // cerr << "CntBelow " << Val << ", " << res << ", " << other <<"\n";
  return res;
}

/// sets all values below "Val" to "Val"
void ClipAt(ArrayBaseType Val)
{
  replace_if(TheArray.begin(),TheArray.end(),bind2nd(less<ArrayBaseType>(),Val),Val);
}
/// sets all values below "Val" to "Val2"
void ClipAt(ArrayBaseType Val,ArrayBaseType Val2)
{
  replace_if(TheArray.begin(),TheArray.end(),bind2nd(less<ArrayBaseType>(),Val),Val2);
}

/// sets all values below "Val" to "Val2"
void ClipBetween(ArrayBaseType Val,ArrayBaseType ValM, ArrayBaseType Val2)
{
  replace_if(TheArray.begin(),TheArray.end(),bind2nd(less<ArrayBaseType>(),Val),Val2);
}

/// sets all values below "Val" to "Val2" and above or equal, to Val3
void ThreshAt(ArrayBaseType Val,ArrayBaseType Val2, ArrayBaseType Val3)
{
  replace_if_else(TheArray.begin(),TheArray.end(),bind2nd(less<ArrayBaseType>(),Val),Val2,Val3);
}

/// sets all values below "Val" to "Val2" and above or equal, to Val3
void ThreshBetween(ArrayBaseType Val, ArrayBaseType ValM, ArrayBaseType Val2, ArrayBaseType Val3)
{
  replace_if_else(TheArray.begin(),TheArray.end(),
	between<ArrayBaseType>(Val,ValM),
	Val3,Val2);
}

/// sets all values in array to "Val"
void Set(ArrayBaseType Val)
{
  fill(TheArray.begin(),TheArray.end(),Val);
}


typedef ArrayBaseType (* FillFktType) (Vector * Pos);        // FillFunktion Type
/// fills the array with values calculated by the "FillFkt" invoked with a vector pointing to the actual array-position.
void Fill(FillFktType FillFkt)
{
  int i=0,j=0;
  Vector Tmp(TheArray.Dimension()); // allocate a vector of length Dimension
  ArrayBaseType * myp1= TheArray.StartPos();

  Tmp.SetZero();

  while(myp1<=TheArray.EndPos())
    {
    (* myp1) = FillFkt(& Tmp);

    Tmp.ChangeComp(0,Tmp.comp(0)+1.0); // increment X-coordinate by 1
    i=0;
    while ((i < TheArray.Dimension() -1) && (Tmp.comp(i) >= GetSize(i)))
      {
	Tmp.ChangeComp(i,0);
	Tmp.ChangeComp(i+1,Tmp.comp(i+1)+1.0); // increment next dimension
	i++; // Now next dimension
      }
    j++;
    myp1++;
    }
}

/** returns true (1) if voxel at (x,y,z) is an outer voxel above the "Threshhold". Outer voxel means a neighbour is below
"Threshhold" */ 
int CheckOuterThreshVoxel(CalcFloatType Threshhold, int x, int y=0, int z=0)
{
  int dx,dy,dz;

  if (CastToReal(Value(ClipSel,x,y,z)) < Threshhold)
    return 0;
  else
  {
    for (dz=-1;dz<=1;dz++)         // if one of the 27 neigbours is outside, then voxel is valid
      for (dy=-1;dy<=1;dy++)
	for (dx=-1;dx<=1;dx++)
	  if (CastToReal(Value(ClipSel,x-dx,y-dy,z-dz)) < Threshhold)
	    return 1;
  }
  return 0;
}

/** returns true (1) if voxel at (x,y,z) is an outer voxel above the "Threshhold". Outer voxel means a neighbour is below
"Threshhold" */ 
int Check6OuterThreshVoxel(CalcFloatType Threshhold, int x, int y=0, int z=0)
{
  int dx,dy,dz;

  if (CastToReal(Value(ClipSel,x,y,z)) < Threshhold)
    return 0;
  else
  {
    dx=-1,dy=0,dz=0;
    if (CastToReal(Value(ClipSel,x-dx,y-dy,z-dz)) < Threshhold)
      return 1;
    dx=1,dy=0,dz=0;
    if (CastToReal(Value(ClipSel,x-dx,y-dy,z-dz)) < Threshhold)
      return 1;
    dx=0,dy=-1,dz=0;
    if (CastToReal(Value(ClipSel,x-dx,y-dy,z-dz)) < Threshhold)
      return 1;
    dx=0,dy=1,dz=0;
    if (CastToReal(Value(ClipSel,x-dx,y-dy,z-dz)) < Threshhold)
      return 1;
    dx=0,dy=0,dz=-1;
    if (CastToReal(Value(ClipSel,x-dx,y-dy,z-dz)) < Threshhold)
      return 1;
    dx=0,dy=0,dz=1;
    if (CastToReal(Value(ClipSel,x-dx,y-dy,z-dz)) < Threshhold)
      return 1;
  }
  return 0;
}


/// calculates the COM of a voxel range x+/- range  ,...
void RangedCOM(CalcFloatType & mx, CalcFloatType & my, CalcFloatType & mz, int rangex, int rangey, int rangez)
{
  int x,y,z;
  CalcFloatType sx=0,sy=0,sz=0,min,Sum=0,val,max;
  // int rangex=IntoRange(range,GetSize(0));
  // int rangey=IntoRange(range,GetSize(1));
  // int rangez=IntoRange(range,GetSize(2));
	
  min = CastToReal(Value(int(mx),int(my),int(mz)));
  max = min;
  for (z=int(mz)-rangez;z<=int(mz)+rangez;z++)              // find minimum
    for (y=int(my)-rangey;y<=int(my)+rangey;y++)
      for (x=int(mx)-rangex;x<=int(mx)+rangex;x++)
	{
	  val = CastToReal(Value(IntoRange(x,GetSize(0)),IntoRange(y,GetSize(1)),IntoRange(z,GetSize(2))));
	  if (val < min) min=val;
	}
  for (z=int(mz)-rangez;z<=int(mz)+rangez;z++)              // apply range and find COM
    for (y=int(my)-rangey;y<=int(my)+rangey;y++)
      for (x=int(mx)-rangex;x<=int(mx)+rangex;x++)
	{
	  val = CastToReal(Value(IntoRange(x,GetSize(0)),IntoRange(y,GetSize(1)),IntoRange(z,GetSize(2))));
	  if (val==val) // If not NaN
	    {
	      sx += (val-min)*x;
	      sy += (val-min)*y;
	      sz += (val-min)*z;
	      Sum+= (val-min);
	    }
	  // cout << "x = " << x << " y = " << y << " z = " << z << " -> val = " << val << " min = " << min << "\n";
	}
  if (Sum != 0)   // meaning min == max ->  just return the original spot
    {
      mx = IntoRange(sx/Sum,(CalcFloatType) GetSize(0));  // otherwise it might fall out of range !
      my = IntoRange(sy/Sum,(CalcFloatType) GetSize(1));
      mz = IntoRange(sz/Sum,(CalcFloatType) GetSize(2));
    }
  return;
}

ArrayBaseType MeanRegionIntensity(int mx, int my, int mz, int rangex, int rangey, int rangez)
{
        ArrayBaseType val=0;
    for (int z=mz-rangez;z<=mz+rangez;z++)              // find minimum
      for (int y=my-rangey;y<=my+rangey;y++)
        for (int x=mx-rangex;x<=mx+rangex;x++)
	  val += Value(IntoRange(x,GetSize(0)),IntoRange(y,GetSize(1)),IntoRange(z,GetSize(2)));
    return val / (double(rangex*2+1)*(rangey*2+1)*(rangez*2+1));
}

/** determines the center of Mass of a spot with 
  values above "ClipVal" 
  Function returns mean value of Spot. */
void ClippedCOM(ArrayBaseType ClipVal, Vector * Result)
{
  int i=0,j=0;
  CalcFloatType Voxels=0.0;
  CalcFloatType WSum=0.0,Val;

  Vector Tmp(TheArray.Dimension()); // allocate a vector of length Dimension pointing to array-elements
  Vector Tmp2(TheArray.Dimension());
  Vector COM(TheArray.Dimension());

  ArrayBaseType * myp1= TheArray.StartPos();

  while(myp1<=TheArray.EndPos())
    {
      Val = CastToReal(* myp1) - CastToReal(ClipVal);   // ehemals faelschlicherweise abs(..)

      if (Val > 0.0)
	{
	  Tmp2.copy(&Tmp);
	  Tmp2.Mul(Val);
	  COM.Add(&Tmp2);
	  WSum += Val;
	  Voxels++;
	}

      // Now advance to next coordinate
      Tmp.ChangeComp(0,Tmp.comp(0)+1.0);
      i=0;
      while ((i < TheArray.Dimension()-1) && (Tmp.comp(i) >= GetSize(i)))
	{
	  Tmp.ChangeComp(i,0);
	  Tmp.ChangeComp(i+1,Tmp.comp(i+1)+1.0);
	  i++; // Now next dimension
	}
      j++;
      
      myp1++;
    }
  COM.Div(WSum);
  Result->copy(& COM);
  return ; // WSum/Voxels;
} 

/// replaces all occurences of "ValEquals" with value "ReplaceWith" in array.
void ReplaceVal(ArrayBaseType ValEquals,ArrayBaseType ReplaceWith)
{
  replace_if(TheArray.begin(),TheArray.end(),bind2nd(equal_to<ArrayBaseType>(),ValEquals),ReplaceWith);
}

// min_element ... accumulate ...

/// add another array (ArgArr) to this one
void Add(TArray3d * ArgArr)
{
  transform(TheArray.begin(),TheArray.end(),ArgArr->TheArray.begin(),TheArray.begin(),
	    plus<ArrayBaseType>()); 
}


/// subtract another array (ArgArr) from this one
void Sub(TArray3d * ArgArr)
{
  transform(TheArray.begin(),TheArray.end(),ArgArr->TheArray.begin(),TheArray.begin(),
	    minus<ArrayBaseType>()); 
}

/// subtract this array from "ArgArr" and overwrite this array.
void SubFrom(TArray3d * ArgArr)
{
  transform(ArgArr->TheArray.begin(),ArgArr->TheArray.end(),TheArray.begin(),TheArray.begin(),
	    minus<ArrayBaseType>()); 
}

/// multiply this array voxel by voxel with "ArgArr"
void Mul(TArray3d * ArgArr)
{
  transform(TheArray.begin(),TheArray.end(),ArgArr->TheArray.begin(),TheArray.begin(),
	    multiplies<ArrayBaseType>()); // deprecated ! (now : multiplies) before : times
}

/// multiply this array by a number "Arg"
void Mul(ArrayBaseType Arg)
{
  transform(TheArray.begin(),TheArray.end(),TheArray.begin(),
	    bind2nd(multiplies<ArrayBaseType>(),Arg)); // deprecated ! (now : multiplies) before : times
}

/// raise the array the a power of exponent
void Pow(ArrayBaseType exponent)   // Negative values will be set to zero!
{
   transform(TheArray.begin(),TheArray.end(),TheArray.begin(),
      bind2nd(raisetopower<ArrayBaseType>(),exponent)); // see definition in stlfixes, will clip at zero!
}

void PlaneMul(TArray3d * ArgArr,int Z)  // Multiplies whole Z-stack with the Z's plane of ArgArr
{
  int x,y,z;
  for (z=0;z<GetSize(2);z++) 
    for (y=0;y<GetSize(1);y++)
      for (x=0;x<GetSize(0);x++)
	{
	  SetValue(x,y,z,Value(x,y,z)*ArgArr->Value(x,y,Z));
	}

}


/// Add  a number "Arg"
void Add(ArrayBaseType Arg)
{
  transform(TheArray.begin(),TheArray.end(),TheArray.begin(),
	    bind2nd(plus<ArrayBaseType>(),Arg));
}

/// Sub  a number "Arg"
void Sub(ArrayBaseType Arg)
{
  transform(TheArray.begin(),TheArray.end(),TheArray.begin(),
	    bind2nd(minus<ArrayBaseType>(),Arg));
}

void SubFrom(ArrayBaseType Arg)
{
  transform(TheArray.begin(),TheArray.end(),TheArray.begin(),
	    bind1st(minus<ArrayBaseType>(),Arg));
}


ArrayBaseType TheMedian(int mx,int my, int mz,int SX,int SY, int SZ,ArrayBaseType* rawarray)
{
  int x,y,z,i=0;
  int AllSize=(2*SX+1)*(2*SY+1)*(2*SZ+1);

    for (z=int(mz)-SX;z<=int(mz)+SX;z++)              // copy into tmpArray
    for (y=int(my)-SY;y<=int(my)+SY;y++)
      for (x=int(mx)-SZ;x<=int(mx)+SZ;x++)
	{
	  rawarray[i] = Value(IntoRange(x,GetSize(0)),
				      IntoRange(y,GetSize(1)),
				      IntoRange(z,GetSize(2)));
	  i++;
	}
    sort(rawarray,rawarray+AllSize);
    return (rawarray[(AllSize-1)/2]+rawarray[AllSize/2])/2.0;  // In case of even size : take average
}


bool DIterationStep(int d0,int d1, int d2, bool fwd, TArray3d & Penalty,TArray3d & DistX,TArray3d & DistY,TArray3d & DistZ)
{
  int Dim0=GetSize(d0),Dim1=GetSize(d1),Dim2=GetSize(d2),i[3];
  int start0=0,end0=Dim0,step0=1;
  bool updated = false;
  float BestDist,OldPenalty,Best[3],NewDist,NewDV[3];

  if (! fwd) {step0 = -1; start0 = Dim0-1;end0 = -1;}    
  
  for (i[d2]=0; i[d2] < Dim2; i[d2]++)
  for (i[d1]=0; i[d1] < Dim1; i[d1]++)
  {
    i[d0]=start0;
    Best[0] = DistX(i[0],i[1],i[2]);
    Best[1] = DistY(i[0],i[1],i[2]);
    Best[2] = DistZ(i[0],i[1],i[2]);
    OldPenalty = Penalty(i[0],i[1],i[2]);
    Best[d0] += step0;  // make a move
    i[d0]+=step0;
  for (; i[d0] != end0; i[d0]+=step0)
    {
    BestDist = OldPenalty + Best[0]*Best[0]+Best[1]*Best[1]+Best[2]*Best[2];
      // cout << "("<<i[0]<<","<<i[1]<<","<<i[2]<<") "<< flush;
      NewDV[0] = DistX(i[0],i[1],i[2]);
      NewDV[1] = DistY(i[0],i[1],i[2]);
      NewDV[2] = DistZ(i[0],i[1],i[2]);
      NewDist = Penalty(i[0],i[1],i[2]) + NewDV[0]*NewDV[0]+NewDV[1]*NewDV[1]+NewDV[2]*NewDV[2];
      if (BestDist < NewDist)  // a shorter path was found
        {
          Penalty(i[0],i[1],i[2]) = OldPenalty;  // write a better estimate to the data-fields
          DistX(i[0],i[1],i[2]) = Best[0];
          DistY(i[0],i[1],i[2]) = Best[1];
          DistZ(i[0],i[1],i[2]) = Best[2];
          updated = true;
        }
        else
        {
                                                // get a better estimate from the data fields
          Best[0]=NewDV[0];
          Best[1]=NewDV[1];
          Best[2]=NewDV[2];
          OldPenalty = Penalty(i[0],i[1],i[2]);
        }
    Best[d0] += step0;  // make a move
    }
  }
  return updated;
}

void DTransFinalize(TArray3d & Penalty,TArray3d & DistX,TArray3d & DistY,TArray3d & DistZ)
{
  int Dim0=GetSize(0),Dim1=GetSize(1),Dim2=GetSize(2),x,y,z;
  for (z=0; z < Dim2; z++)
  for (y=0; y < Dim1; y++)
  for (x=0; x < Dim0; x++)
    Penalty(x,y,z) = sqrt(SQR(Penalty(x,y,z))+SQR(DistX(x,y,z))+SQR(DistY(x,y,z))+SQR(DistZ(x,y,z)));
}

/// 3D Distance transform algorithm
void DTrans(TArray3d & Penalty,TArray3d & DistX,TArray3d & DistY,TArray3d & DistZ)
{
  bool updated=false;
  int i=1;

  do {
    updated = false;
    updated |= DIterationStep(0,1,2,true,Penalty,DistX,DistY,DistZ);
    updated |= DIterationStep(0,1,2,false,Penalty,DistX,DistY,DistZ);
    updated |= DIterationStep(1,0,2,true,Penalty,DistX,DistY,DistZ);
    updated |= DIterationStep(1,0,2,false,Penalty,DistX,DistY,DistZ);
    updated |= DIterationStep(2,1,0,true,Penalty,DistX,DistY,DistZ);
    updated |= DIterationStep(2,1,0,false,Penalty,DistX,DistY,DistZ);
    cout << "Did update nr. " << i << "\n";
    i++;
  } while (updated); 
}
             
             
void Median(int SX,int SY, int SZ,TArray3d * Out)
{
  int x,y,z;
  int AllSize=(2*SX+1)*(2*SY+1)*(2*SZ+1);
  ArrayBaseType* rawarray = new ArrayBaseType[AllSize];
  if (! rawarray) cerr << "Fatal Error : Could not allocate median filter heap !\n",exit(-1);
  for (z=0;z<GetSize(2);z++)              // find minimum
    for (y=0;y<GetSize(1);y++)
      for (x=0;x<GetSize(0);x++)
	{
	  Out->SetValue(x,y,z,TheMedian(x,y,z,SX,SY,SZ,rawarray));
	}
  delete(rawarray);

}

void Sort()
{
  sort(TheArray.begin(),TheArray.end());
} 

void Reverse()
{
  reverse(TheArray.begin(),TheArray.end());
} 

template <class ArgArray>
void extract(int dim, int nr, ArgArray * Arg)  // writes only into the 0th layer
{
  int x,y,z,sx=0,ex=GetSize(0),sy=0,ey=GetSize(1),sz=0,ez=GetSize(2);

  if (dim == 0)
    sx=nr,ex=nr+1;

  if (dim == 1)
    sy=nr,ey=nr+1;

  if (dim == 2)
    sz=nr,ez=nr+1;

  for (z=sz;z<ez;z++)
    for (y=sy;y<ey;y++)
      for (x=sx;x<ex;x++)
	{
	  if (dim == 0)
	    SetValue(0,y,z,(ArrayBaseType) Arg->Value(x,y,z));
	  else if (dim == 1)
	    SetValue(x,0,z,(ArrayBaseType) Arg->Value(x,y,z));
	  else if (dim == 2)
	    SetValue(x,y,0,(ArrayBaseType) Arg->Value(x,y,z));
	}
}

template <class ArgArray>
void insert(int dim, int nr, ArgArray * Arg)  // writes 0th plane somewhere into stack
{
  int x,y,z,sx=0,ex=GetSize(0),sy=0,ey=GetSize(1),sz=0,ez=GetSize(2);

  if (dim == 0)
    sx=nr,ex=nr+1;

  if (dim == 1)
    sy=nr,ey=nr+1;

  if (dim == 2)
    sz=nr,ez=nr+1;

  for (z=sz;z<ez;z++)
    for (y=sy;y<ey;y++)
      for (x=sx;x<ex;x++)
	{
	  if (dim == 0)
	    SetValue(x,y,z,(ArrayBaseType) Arg->Value(0,y,z));
	  else if (dim == 1)
	    SetValue(x,y,z,(ArrayBaseType) Arg->Value(x,0,z));
	  else if (dim == 2)
	    SetValue(x,y,z,(ArrayBaseType) Arg->Value(x,y,0));
	}
}


typedef ArrayBaseType (* PCompareFkt) (ArrayBaseType,ArrayBaseType);

/** apply a function (having two arguments (ArraybaseType)  and returning ArrayBaseType) to 
   every voxel */
void Apply(PCompareFkt  pFkt, TArray3d * ArgArr)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();

  while(myp1>=TheArray.StartPos())
    {
      (* myp1) = pFkt(* myp2,* myp1);           // The order is very important !!
      myp1--;
      myp2--;
    }
}
typedef ArrayBaseType (* PCorAppFkt) (int,int,int,ArrayBaseType);
void Apply(PCorAppFkt  pFkt)
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;
  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
          SetValue(x,y,z,(*pFkt)(x,y,z,Value(x,y,z)));
}

typedef ArrayBaseType (* PAppFkt) (ArrayBaseType);
/** apply a function (having one argument) to every voxel */
void Apply(PAppFkt  pFkt)
{
  ArrayBaseType * myp1 = TheArray.EndPos();

  while(myp1>=TheArray.StartPos())
    {
      (* myp1) = pFkt(* myp1); 
      myp1--;
    }
}

/// add another array squared (ArgArr) to this one
void SqrAdd(TArray3d * ArgArr)   // x += arg^2
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();

   while(myp1>=TheArray.StartPos())
    {
      (* myp1) += (* myp2) * (* myp2);   
      myp1--;
      myp2--;
    }
}

/// Calculates the Variance and Mean from Sum, SumSqr and Number of counts
void CalcVarMean(TArray3d * Sum, int Counts)   //  (SumSqr-sum*sum/Num)/Num
{
  ArrayBaseType * myp1 = TheArray.EndPos(),  // Contains sum of Squares
    * myp2 = Sum->TheArray.EndPos();

    double rCounts = double(Counts);

   while(myp1>=TheArray.StartPos())
    {
      (* myp1) = ((* myp1) - (* myp2)*(* myp2) / rCounts) / rCounts;   
      (* myp2) /= rCounts;   
      myp1--;
      myp2--;
    }
}

void MulMinusSelfAdd(TArray3d * ArgArr, CalcFloatType overrelax, double numElem)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();

   while(myp1>=TheArray.StartPos())
    {
      (* myp1) += overrelax*((* myp2) * (* myp1) - numElem * (* myp1));   // This function is for deconvolution : ReconImg * CorrectionVal - ReconImg
      if ((* myp1) < 0.000001) (* myp1)= 0.000001;    // Constraint to prevent impossible overrelaxiations
      myp1--;
      myp2--;
    }
}

void MulSelfAdd(TArray3d * ArgArr, CalcFloatType overrelax)  // for applying the correction data in ML-deconv (when 1.0 was subtracted)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();

   while(myp1>=TheArray.StartPos())
    {
      (* myp1) *= (1.0 + overrelax*(* myp2));   // This function is for deconvolution : ReconImg * CorrectionVal - ReconImg
      if ((* myp1) < 0.000001) (* myp1)= 0.000001;    // Constraint to prevent impossible overrelaxiations
      myp1--;
      myp2--;
    }
}


/// add another array (ArgArr) to this one with weight arg.
void AddMul(TArray3d * ArgArr, ArrayBaseType Arg)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();

  while(myp1>=TheArray.StartPos())
    {
      (* myp1) += (* myp2) * Arg;   // This will NOT work properly in int-type-arrays !
      myp1--;
      myp2--;
    }
}
		 
/// replace array voxel by voxel by "ArgArr" divided by this array values. 0/0 -> 0, a/0 -> MAXARRAYTYPENR 
void ArgDivSelf(TArray3d * ArgArr)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();

  while(myp1>=TheArray.StartPos())
    {
#ifdef NegIsInvalid
      if (CastToReal(* myp2) >= 0.0)
#endif
      	if (CastToReal(* myp1) > 0.0)
		(* myp1) = (* myp2) / (* myp1);   // This will NOT work properly in int-type-arrays !
      	else
        	if (CastToReal(* myp2) == 0)
	  	(* myp1) = 0;  // Maybe no good guess ?
        	else
	  	(* myp1) = (ArrayBaseType) MAXARRAYTYPENR;
#ifdef NegIsInvalid
      else
	(* myp1) = (ArrayBaseType) 1.0;
#endif
      myp1--;
      myp2--;
    }
}
		 
		 
/// replace array voxel by voxel by "ArgArr" divided by this array values and subtracts one. 0/0 -> 0, a/0 -> MAXARRAYTYPENR 
void ArgDivSelfM1(TArray3d * ArgArr)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();

  while(myp1>=TheArray.StartPos())
    {
#ifdef NegIsInvalid
      if (CastToReal(* myp2) >= 0.0)
#endif
      if (CastToReal(* myp1) > 0.0)
	(* myp1) = (* myp2) / (* myp1) - 1.0;   // This will NOT work properly in int-type-arrays !
      else
        if (CastToReal(* myp2) == 0)
	  (* myp1) = 0 - 1.0;  // Maybe no good guess ?
        else
	  (* myp1) = (ArrayBaseType) MAXARRAYTYPENR;
#ifdef NegIsInvalid
      else
	(* myp1) = (ArrayBaseType) 0.0;
#endif
      myp1--;
      myp2--;
    }
}

/// replace array voxel by voxel by "ArgArr" divided by this array values. 0/0 -> 0, a/0 -> MAXARRAYTYPENR 
// computes also the log likelihood sum(Mi ln Ei - Ei) omitting  - ln (Mi !)
double ArgDivSelfLogLikelihood(TArray3d * ArgArr)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();
  double loglikelihood=0.0;

  while(myp1>=TheArray.StartPos())
    {
#ifdef NegIsInvalid
      if (CastToReal(* myp2) > 0.0)  // negative values in the data do not count in the likelihood computations
#endif
      if ((* myp1) > 0.0)  // negativities can also arise
	{
	  loglikelihood += (* myp2) * log((* myp1)) - (* myp1);  // computes also the log likelihood sum(Mi ln Ei - Ei) omitting  - ln (Mi !)
	  (* myp1) = (* myp2) / (* myp1);   // This will NOT work properly in int-type-arrays !
	}
      else
        if ((* myp2) == 0)
	  (* myp1) = 0;  // Maybe no good guess ?
        else
	  (* myp1) = (ArrayBaseType) MAXARRAYTYPENR;
      myp1--;
      myp2--;
    }
  return loglikelihood;
}


/// replace array voxel by voxel by "ArgArr" divided by this array values. 0/0 -> 0, a/0 -> MAXARRAYTYPENR 
// computes also the log likelihood sum(Mi ln Ei - Ei) omitting  - ln (Mi !)
double ArgDivSelfM1LogLikelihood(TArray3d * ArgArr)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();
  double loglikelihood=0.0;

  while(myp1>=TheArray.StartPos())
    {
#ifdef NegIsInvalid
      if (CastToReal(* myp2) > 0.0)  // negative values in the data do not count in the likelihood computations
#endif
      if ((* myp1) > 0.0)  // negativities can also arise
	{
	  loglikelihood += (* myp2) * log((* myp1)) - (* myp1);  // computes also the log likelihood sum(Mi ln Ei - Ei) omitting  - ln (Mi !)
	  (* myp1) = (* myp2) / (* myp1) - 1.0;   // This will NOT work properly in int-type-arrays !
	}
      else
        if ((* myp2) == 0)
	  (* myp1) = 0 -1.0;  // Maybe no good guess ?
        else
	  (* myp1) = (ArrayBaseType) MAXARRAYTYPENR;
      myp1--;
      myp2--;
    }
  return loglikelihood;
}

/// replace array voxel by voxel by number "Arg" divided by this array values. 0/0 -> 0, Arg/0 -> MAXARRAYTYPENR 
void ArgDivSelf(ArrayBaseType Arg)
{
  ArrayBaseType * myp1 = TheArray.EndPos();
 
  while(myp1>=TheArray.StartPos())
    {
      if ((* myp1) != 0)
	(* myp1) = Arg / (* myp1);   // This will NOT work properly in int-type-arrays !
      else
        if ((Arg) == 0)
	  (* myp1) = 0;  // Maybe no good guess ?
        else
	  (* myp1) = (ArrayBaseType) MAXARRAYTYPENR;
      myp1--;
    }
}

void Sqrt(void)
{
  ArrayBaseType * myp1 = TheArray.EndPos();
  while(myp1>=TheArray.StartPos())
    {
      (* myp1) = (ArrayBaseType) sqrt((* myp1));
      myp1--;
    }
}

void Sqr(void)
{
  ArrayBaseType * myp1 = TheArray.EndPos();
  while(myp1>=TheArray.StartPos())
    {
      (* myp1) = (ArrayBaseType) ((* myp1)*(* myp1));
      myp1--;
    }
}

/// replace voxels in this by (ArgArr - this) / Sqrt(ArgArr + this)
void SubsDivSqrtAdd(TArray3d * ArgArr)
{
  ArrayBaseType * myp1 = TheArray.EndPos(),
                * myp2 = ArgArr->TheArray.EndPos();

  while(myp1>=TheArray.StartPos())
    {
      (* myp1) = (ArrayBaseType) (((* myp2) - (* myp1))/sqrt((* myp2) + (* myp1)));
      myp1--;
      myp2--;
    }
}

  /// Calculates self as the numerical derivative of "ArgArr" in direction "direction"
void DeriveFrom(TArray3d * ArgArr,const int direction)
{
  int x,y,z;
  ArrayBaseType v;

  if (this == ArgArr) cerr << "ERROR : Derivefrom mustnt be called on the same array ! \n";

  for (z=0;z<GetSize(2);z++)
    for (y=0;y<GetSize(1);y++)
      for (x=0;x<GetSize(0);x++)
	{
	v=(ArrayBaseType) ((ArgArr->Value(ClipSel,x+(direction==0),y+(direction==1),z+(direction==2)) -
	   ArgArr->Value(ClipSel,x-(direction==0),y-(direction==1),z-(direction==2)))/2.0);
	TheArray.SetValue(x,y,z,v);
	}
}

int IntoRange(int pos,int Max)
{
  if (pos < 0)
      return (Max-((-pos) % Max)) % Max;
  else
    return (pos % Max);
}

CalcFloatType FMod(CalcFloatType Val,CalcFloatType Div)
{
  return Val-int(Val/Div)*Div;
}

CalcFloatType IntoRange(CalcFloatType pos,CalcFloatType Max)
{
  if (pos < 0)
      return FMod(Max - FMod((-pos) , Max) , Max);
  else
    return FMod(pos , Max);
}

/// range : apply a cyclic range and use a COM based determination of the Position (cyclic).
CalcFloatType MaxIntensityPosition (CalcFloatType * px, CalcFloatType * py, CalcFloatType * pz,int range=0, bool forceplanes=false)
{
  int x,y,z,XM=GetSize(0),YM=GetSize(1),ZM=GetSize(2);
  CalcFloatType mx=0,my=0,mz=0;
  CalcFloatType max=abs(Value(0,0,0)); // CastToReal(Value(0,0,0));
  int rangez=range;

  if (forceplanes) 
    {
      ZM = 1;  // restrict search to z = 0
      rangez=0;
    }

  for (z=0;z<ZM;z++)
    for (y=0;y<YM;y++)
      for (x=0;x<XM;x++)
	{
	  if (abs(Value(x,y,z)) > max) 
	    {
	      max=CastToReal(Value(x,y,z));
	      mx=x;
	      my=y;
	      mz=z;
	    }
	}
  cout << "Size : " << XM << " " << YM << " " << ZM << "\n";
  cout << "Max at : " << mx << " " << my << " " << mz << " is " << max << "\n";
  
  if (range > 0)  // calculate mean
    {
      RangedCOM(mx,my,mz,range,range,rangez);
    }


  (*px) = mx;
  (*py) = my;
  (*pz) = mz;
  return max;
}

/// range : apply a cyclic range and use a Maximum based determination of the Position (cyclic).
CalcFloatType MaxIntensityNeighbour (int & px, int & py, int & pz,int range)
{
  int x,y,z,
    mx=IntoRange(int(px)-range,GetSize(0)),
    my=IntoRange(int(py)-range,GetSize(1)),
    mz=IntoRange(int(pz)-range,GetSize(2));
  CalcFloatType max=CastToReal(Value(mx,my,mz)),val;
  

  for (z=int(pz)-range;z<=int(pz)+range;z++)              // find minimum
    for (y=int(py)-range;y<=int(py)+range;y++)
      for (x=int(px)-range;x<=int(px)+range;x++)
	{
	  val = CastToReal(Value(IntoRange(x,GetSize(0)),IntoRange(y,GetSize(1)),IntoRange(z,GetSize(2))));
	  if (val > max) 
	    {
	      max=val;
	      mx=IntoRange(x,GetSize(0));
	      my=IntoRange(y,GetSize(1));
	      mz=IntoRange(z,GetSize(2));
	    }
	}
  px=mx,py=my,pz=mz;
  
  return max;
}


/// The range is from - Size/2 to Size/2
CalcFloatType WrappedMaxIntensityPosition (CalcFloatType * px, CalcFloatType * py, CalcFloatType * pz,int range=0, bool forceplanes=false)
{
  CalcFloatType max=MaxIntensityPosition(px,py,pz,range, forceplanes);
  if ((*px) > GetSize(0)/2)  // to get positive and negative values
    (*px) = (*px)-GetSize(0);

  if ((*py) > GetSize(1)/2)
    (*py) = (*py)-GetSize(1);

  if ((*pz) > GetSize(2)/2)
    (*pz) = (*pz)-GetSize(2);
  return max;
}

typedef ArrayBaseType (* WeightFunctor) (ArrayBaseType,CalcFloatType);

static ArrayBaseType LinearWeightor(ArrayBaseType val, CalcFloatType reldist) 
{
  if (reldist > 1.0) 
    return 0.0;
  else
    return val*(ArrayBaseType) (1.0-reldist);
}

static ArrayBaseType GaussWeightor(ArrayBaseType val, CalcFloatType reldist) 
{
    return val * (ArrayBaseType) exp(-(reldist*reldist));
}

/// The range is from - Size/2 to Size/2
void ApplyDistWeight(double MaxRelDist, double MaxRelDistZ, WeightFunctor MyWeightFunctor = &LinearWeightor)
{
  double px,py,pz,weight;
  int SX=GetSize(0);
  int SY=GetSize(1);
  int SZ=GetSize(2);
  for (int z=0;z<SZ;z++)              // find minimum
    for (int y=0;y<SY;y++)
      for (int x=0;x<SX;x++)
	{
	  px=x;py=y;pz=z;
	  if (px > SX/2)  // to get positive and negative values
	    px = px-SX;

	  if (py > SY/2)
	    py = py-SY;
	  
	  if (pz > SZ/2)
	    pz = pz-SZ;

	  if (MaxRelDist > 0.0)
	    {
	      px /= (MaxRelDist*SX / 2.0);
	      py /= (MaxRelDist*SY / 2.0);
	    }
	  else
	    {
	      px=0;py=0;
	    }
	  
	  if (MaxRelDistZ > 0.0)
	    pz /= (MaxRelDistZ*SZ / 2.0);
	  else
	    pz=0;
	  
	  weight = sqrt(px*px+py*py+pz*pz);
	  TheArray.SetValue(x,y,z,MyWeightFunctor(TheArray.GetValue(x,y,z),weight));
	}
}


/// calculates "MaxVector" pointing to position in array with maximal intensity (sets "* MaxInten")
void MaxIntensityVec(Vector * MaxVector,CalcFloatType * MaxInten,int range=0)
{
  CalcFloatType mx=0,my=0,mz=0;
  
  (* MaxInten) = MaxIntensityPosition(&mx,&my,&mz,range);
  MaxVector->copy((CalcFloatType)mx,(CalcFloatType)my,(CalcFloatType)mz);
}

void EdgeWindow(int distance) // sets outer "distance" pixels to 0
{
  int x,y,z;
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2);
  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
	{
	  if (((x < distance) || (x >= GetSize(0)-distance))
	    || ((y < distance) || (y >= GetSize(1)-distance))
	    || ((z < distance) || (z >= GetSize(2)-distance)))
          SetValue(x,y,z,0);
	}
}


void CircWindow(CalcFloatType radius) // sets outer "distance" pixels to 0
{
  int x,y,z;
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2);
  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
	if ((SQR((x - (DimX-1)/2.0)/(DimX/2.0)) + SQR((y-(DimY-1)/2.0)/(DimY/2.0)) + SQR((z-(DimZ-1)/2.0)/(DimZ/2.0))) >= SQR(radius) )
          SetValue(x,y,z,0);
}

/// sets the border-pixel (coordinate of dim == 0) to zero
void RemoveEdge(int dim)
{
  int x,y,z;
  int DimX=GetSize(0), DimY=GetSize(1), DimZ=GetSize(2);
  
  if (dim == 0) DimX=1;
  if (dim == 1) DimY=1;
  if (dim == 2) DimZ=1;

  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	{
	  SetValue(x,y,z,0);
	}
}

double CreateHisto(TArray3d * Histo, ArrayBaseType min, ArrayBaseType max, TArray3d * ToHisto = 0, double binsize = 0.0)
{
  int DimX=GetSize(0), DimY=GetSize(1), DimZ=GetSize(2);
  int HDimX=Histo->GetSize(0),hpi,x,y,z;
  double hp;

  if (binsize == 0.0)
    binsize = (max-min) / (HDimX -1);
  else
  {
    min = 0.0;
    cout << "WARNING: BinSize was given, histogram will start at zero!\n";
  }

  Histo->Clear();

  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	{
	  hp = (Value(x,y,z)- min) / binsize;
	  hpi = int ( hp+0.5);
	  if (hpi < 0) hpi=0;
	  if (hpi >= HDimX) hpi=HDimX-1;

	  if (! ToHisto)
	    Histo->SetValue(hpi,0,0,Histo->Value(hpi,0,0)+1);
	  else
	    Histo->SetValue(hpi,0,0,Histo->Value(hpi,0,0)+ToHisto->Value(x,y,z));
	}
  return binsize;
}

/// Creates a histogram using the sum of square array and the sum array (self)
double CreateVarHisto(TArray3d * VarHisto, TArray3d * Histo, TArray3d * SumHisto, ArrayBaseType min, ArrayBaseType max, 
		      TArray3d * SSQArray, TArray3d * SumArray, int NumsSelf, int Nums, double binsize=0.0)
{
  int DimX=GetSize(0), DimY=GetSize(1), DimZ=GetSize(2);
  int HDimX=Histo->GetSize(0),hpi,x,y,z;
  double hp,num,sum,ssq;

  if (binsize == 0.0)
    binsize = (max-min) / (HDimX -1);
  else
  {
    min = 0.0;
    cout << "WARNING: BinSize was given, histogram will start at zero!\n";
  }

  Histo->Clear();
  VarHisto->Clear();
  SumHisto->Clear();

  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	{
	  hp = ((Value(x,y,z) / NumsSelf) - min) / binsize;
	  hpi = int ( hp+0.5);
	  if (hpi < 0) hpi=0;
	  if (hpi >= HDimX) hpi=HDimX-1;

	  Histo->SetValue(hpi,0,0,Histo->Value(hpi,0,0)+1);   // Just count how many SSQs have been summed in this bin
	  VarHisto->SetValue(hpi,0,0,VarHisto->Value(hpi,0,0)+SSQArray->Value(x,y,z));  // Collect SSQs here
	  SumHisto->SetValue(hpi,0,0,SumHisto->Value(hpi,0,0)+SumArray->Value(x,y,z));  // Collect sums here
	}

  for (x=0;x < HDimX;x++)
    {
      num = Histo->Value(x,0,0) * Nums;
      ssq = VarHisto->Value(x,0,0);
      sum = SumHisto->Value(x,0,0);
      if (num)
	{
	  VarHisto->SetValue(x,0,0, (ssq - sum*sum/num)/num ) ; // Formula for the variance of the whole dataset
	}
    }
  return binsize;
}

/// determines the mean of the individual rings
double  NormedRingSum(TArray3d * Histo)  
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2);
  int hsize=Histo->GetSize(0);
  vector<int> cnt(hsize);
  int x,y,z,ind;
  CalcFloatType r2,x2,y2,z2,scale;

  Histo->Set(0);
  cnt.clear();

  scale=CalcFloatType((DimX > 1)+(DimY > 1)+(DimZ > 1));

  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
	{
	  x2=(x-DimX/2)/CalcFloatType(DimX/2.0);
	  y2=(y-DimY/2)/CalcFloatType(DimY/2.0);
	  z2=(z-DimZ/2)/CalcFloatType(DimZ/2.0);
	  r2 = x2*x2+y2*y2+z2*z2;
	  ind =int(sqrt(r2/scale)*(hsize-1)+0.5);
	  Histo->SetValue(ind,0,0,Histo->Value(ind,0,0)+Value(x,y,z));
	  cnt[ind]++;
	}

  for (ind=0;ind< hsize;ind++)
    Histo->SetValue(ind,0,0,Histo->Value(ind,0,0)/cnt[ind]);

  return sqrt(1.0/scale)/CalcFloatType(DimX/2.0)*(hsize-1);  // returns the scaling for X only
}


// **************** FFT related stuff follows here **********************


/// The following Macro does the StrangeAddressing in three dimension 

ArrayBaseType * STRADDR3DPointer(int X,int Y,int Z)
{
  int DimensionX=GetSize(0);
  int DimensionY=GetSize(1);
  int DimensionZ=GetSize(2);

  if (Z >= DimensionZ/2)
    if (Y >= DimensionY/2)
      if  (X >= DimensionX/2)
        return TheArray.GetPValue(X-DimensionX/2,Y-DimensionY/2,Z-DimensionZ/2);
      else
        return TheArray.GetPValue(X+DimensionX/2,Y-DimensionY/2,Z-DimensionZ/2);
    else
      if (X >= DimensionX/2)
	return TheArray.GetPValue(X-DimensionX/2,Y+DimensionY/2,Z-DimensionZ/2);
      else
	return TheArray.GetPValue(X+DimensionX/2,Y+DimensionY/2,Z-DimensionZ/2);
  else
    if (Y >= DimensionY/2)
      if  (X >= DimensionX/2)
        return TheArray.GetPValue(X-DimensionX/2,Y-DimensionY/2,Z+DimensionZ/2);
      else
        return TheArray.GetPValue(X+DimensionX/2,Y-DimensionY/2,Z+DimensionZ/2);
    else
      if (X >= DimensionX/2)
	return TheArray.GetPValue(X-DimensionX/2,Y+DimensionY/2,Z+DimensionZ/2);
      else
	return TheArray.GetPValue(X+DimensionX/2,Y+DimensionY/2,Z+DimensionZ/2);
}

ArrayBaseType STRADDR3DValue(int X,int Y,int Z)
{
  return * STRADDR3DPointer(X,Y,Z);
}

/// maps val into range [0.. Size-1]
int intoRange(int val,int Size)
{
  return (val >= 0) ? val%Size : Size-((-val)%Size);
}

/// rotates datastack cyclic into positive direction in all coordinates by (dx,dy,dz) voxels
void rotate(int dx, int dy, int dz)
{
  int DimX=GetSize(0), DimY=GetSize(1), DimZ=GetSize(2),
    nx=0,ny=0,nz=0,i,xb=0,yb=0,zb=0,px=0,py=0,pz=0;
  ArrayBaseType swp=Value(0,0,0),nswp;

  for (i=0;i<DimX*DimY*DimZ;i++) 
    {

      nx = intoRange((nx+dx),DimX);
      ny = intoRange((ny+dy),DimY);
      nz = intoRange((nz+dz),DimZ);
      nswp=Value(nx,ny,nz);
      SetValue(nx,ny,nz,swp);
      // TEST : SetValue(nx,ny,nz,py*10000+px);
      swp=nswp;
      if ((ny==yb)&&(nz==zb))  // counts number of pixels that were changed in this row .
	px ++;
      if ((nx==0)&&(nz==zb))
	py ++;
      if ((nx==0)&&(ny==0))
	pz ++;

      if ((nx == xb) && (ny == yb) && (nz == zb)) // reached startingpoint
	{
	  if (px < DimX) // in this case the right neighbour must be free (only equal steps along line are possible)
	    {
	      xb ++;
	    }
	  else if (py < DimY) // x-line is full. So are all other lines reached by algorithm
	    {
	      xb=0;
	      yb ++;
	      px=0;
	    }
	  else // if (pz < DimZ) // xy-plane is full. All other planes reached are also full
	    {
	      xb=0;
	      yb=0;
	      zb ++;
	      px=0; py=0;
	    }
	  xb=intoRange(xb,DimX);
	  yb=intoRange(yb,DimY);
	  zb=intoRange(zb,DimZ);
	  nx=xb;
	  ny=yb;
	  nz=zb;
	  swp=Value(xb,yb,zb);
	} 
    }
}

void Conjugate(void)
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;

  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
	SetValue(x,y,z,Conj(Value(x,y,z)));
}

/// sets a value at (x,y,z) to v
void SetValue(int x,int y,int z,ArrayBaseType v) {TheArray.SetValue(x,y,z,v);}


/** reads an array from a file in raw-format. Stream file has to be set to right position before. */
void  Read(ifstream * file)
{
  
  file->read((char *) TheArray.StartPos(),TheArray.Size());

  if (file->eof())
    {
      cerr << "Reached EOF unexpectedly !\n";
      exit(-1);
    }
}

protected:  // the IO functions below shall not be used by every user

/** reads an array from a file in raw-format. FILE * file has to be set to right position before. */
void  Read(FILE * file)
{
  
  if (fread(TheArray.StartPos(),TheArray.Size(),1,file) != (unsigned) 1)
    {
      cerr << "Reached EOF unexpectedly !\n";
      exit(-1);
    }
}

/// loads a file ("filename") from disk in raw format.
void Load(const char * filename,int offset=0,char * DataIsType="Byte",int ElementNr=0)
{
  if ( * filename == 0)  // empty string
    { cerr << "ERROR: No filename given in KLoad! \n";
      exit (-1);
     }

  ifstream file(filename);

  if (! file) {
    cerr << "couldn't open datafile " << filename << "\n";
    exit (-1);
    }

  char dummy;

  if ( * filename == 0)  // empty string
     return;             // dont do anything

  cerr << "TArray LoadfromDisk : " << filename << ", AllDim :" << TheArray.Size() << "\n";

  if (offset > 0)
    file.seekg(offset,ios::beg);

  ArrayBType Tdummy;
  const char * ReqType=TypeString(Tdummy);

  ReadAndConvert(file, ReqType, DataIsType, ElementNr);

  //if (ElementNr== ValueDimE-1)
  //  {
      file.read(&dummy,sizeof(char));
      if (! file.eof()) 
	cerr << "Warning ! Datafile " << filename << " longer than expected !\n";
  //  }

  file.close();
}

ArrayBaseType BaseTypeUC(unsigned char a) { return (ArrayBaseType) a;}
ArrayBaseType BaseTypeC(char a) { return (ArrayBaseType) a;}
ArrayBaseType BaseTypeS(short a) { return (ArrayBaseType) a;}
ArrayBaseType BaseTypeF(float a) { return (ArrayBaseType) a;}
ArrayBaseType BaseTypeD(double a) { return (ArrayBaseType) a;}

/** Reads the data and converts to the array base type if necessary */
void ReadAndConvert(ifstream & file, const char * ReqType, char * TypeString, int ElementNr=0,int normalize=0)
{
  int ValueDimX=GetSize(0),ValueDimY=GetSize(1),ValueDimZ=GetSize(2);

  if (strcmp(ReqType,TypeString) != 0) 
    {
      cerr << "Kdf has wrong type ! Requested : " << ReqType << ", but has type : " << TypeString << "\n";
      // exit(-1);
      cerr << "Trying a type-conversion\n";
      if (strcmp("Unsigned Byte",TypeString) == 0) 
	{
	  typedef TArray3d<unsigned char>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  
	  // file.rdbuf()->seekoff(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur,ios::in);
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
	  tmpArray->Read(&file);
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
      else if (strcmp("Byte",TypeString) == 0) 
	{
	  typedef TArray3d<char>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
	  tmpArray->Read(&file);
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
      else if (strcmp("Short",TypeString) == 0) 
	{
	  typedef TArray3d<short>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
	  tmpArray->Read(&file);
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
      else if (strcmp("Long",TypeString) == 0) 
	{
	  typedef TArray3d<long>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
	  tmpArray->Read(&file);
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
      else if (strcmp("Unsigned Short",TypeString) == 0) 
	{
	  typedef TArray3d<unsigned short>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
	  tmpArray->Read(&file);
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
      else if (strcmp("Unsigned Long",TypeString) == 0) 
	{
	  typedef TArray3d<unsigned long>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
	  tmpArray->Read(&file);
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
      else if (strcmp("Integer",TypeString) == 0) 
	{
	  typedef TArray3d<int>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
	  tmpArray->Read(&file);
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
      else if (strcmp("Float",TypeString) == 0) 
	{
	  typedef TArray3d<float>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
                tmpArray->Read(&file);
	  if (normalize)
	    tmpArray->NormalizeToMax(255.0);
	    
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
      else if (strcmp("Double",TypeString) == 0) 
	{
	  typedef TArray3d<double>  TA;
	  TA * tmpArray = new TA (GetSize(0),GetSize(1),GetSize(2));
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(TA::ArrayBaseType),ios::cur);
	  tmpArray->Read(&file);
	  if (normalize)
	    tmpArray->NormalizeToMax(255.0);
	  ccopy(tmpArray->Pointer(0,0,0),tmpArray->Pointer(ValueDimX-1,ValueDimY-1,ValueDimZ-1)+1,TheArray.begin());
	  delete tmpArray;
	}
       else cerr << "unable to convert type " << TypeString << " to " << ReqType << " ! Bailing out \n",exit(-1);
    }
  else
    {
	  file.seekg(GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(ArrayBaseType),ios::cur);
      // cout << "seeking position:" << ElementNr << " at " <<  GetSize(0)*GetSize(1)*GetSize(2)*ElementNr*sizeof(ArrayBaseType) << "\n";
      Read(& file);
    }
}

/// loads a file ("filename") from disk in Khoros data format (KDF). Return Number of Elements in file.
int KLoad(const char * filename,const char * ReqType,int GetElem=-1,int normalize=0)        
{
  char dummy;
  int ValueDimX,ValueDimY,ValueDimZ,ValueDimE;
  char TypeString[100];
  int ElementNr=GetElem;

  if (ElementNr < 0)
    ElementNr=0;

  if ( * filename == 0)  // empty string
    { cerr << "ERROR: No filename given in KLoad! \n";
      exit (-1);
     }

  ifstream file(filename);

  if (! file) {
    cerr << "couldn't open datafile " << filename << "\n";
    exit (-1);
    }

  ReadKhorosHeader(& file,TypeString,ValueDimX,ValueDimY,ValueDimZ,ValueDimE);

  cerr << "TArray LoadfromDisk (Kdf) : " << filename << " Type :" << TypeString << "\n";


  if (ValueDimX != GetSize(0)
      || ValueDimY != GetSize(1)
      || ValueDimZ != GetSize(2))
    {
      cerr << "KLoad has to resize array to " << ValueDimX <<" x " << ValueDimY << " x " << ValueDimZ <<"\n";
      Resize(ValueDimX,ValueDimY,ValueDimZ);
    }

  if (ElementNr>=ValueDimE) 
    cerr << "Error : KLoad requested element " << ElementNr << " out of " << ValueDimE << " elements in file " << filename << "\n",exit(-1);

  ReadAndConvert(file,ReqType,TypeString,ElementNr,normalize);

  if (ElementNr== ValueDimE-1)
    {
      file.read(&dummy,sizeof(char));
      if (! file.eof()) 
	cerr << "Warning ! Datafile " << filename << " longer than expected !\n";
    }

  if (GetElem == -1) // warning only if elementnr was not explicitely given
    cerr << "Warning ! Read only element " << ElementNr << " out of " << ValueDimE << " elements in file " << filename << "\n";

  /*  if (dummy != '<')
    cerr << "Warning ! Datafile " << filename << " longer than expected !\n";

  file.read(&dummy,sizeof(char));
  if (dummy != '>')
    cerr << "Warning ! Datafile " << filename << " longer than expected !\n"; */

  file.close();
  // cout << "File loaded\n" << flush;
  return ValueDimE;
}

/** writes array to a file in raw-format. FILE * file has to be set to right position before. */
void  Write(FILE * file)
{
  if (fwrite(TheArray.StartPos(),TheArray.Size(),1,file) != 1) 
    cerr << "Error writing file, no space on disk \n";
}

/// saves array to disk under name ("filename") in raw format.
void  Save(const char * filename)       // Saves Imaginary part as well, if it exists (real1,imag1,real2,imag2,...)
{
  cerr << "TArray SaveToDisk : " << filename << "\n";

  ofstream to(filename);
  if (! to) 
    {
      cerr << "Couldn't open file " << filename << " for writing !!\n";
      exit(-1);
    }

  Write(& to);

  to.close();
}

void KHeader(ofstream * to,int Elem=1,int SizeT=1)  // write khoros header with appropriate type
{
  if (! to) 
    {
      cerr << "File not open  for writing !!\n";
      exit(-1);
    }
  int SizeX=GetSize(0);
  int SizeY=GetSize(1);
  int SizeZ=GetSize(2);
   ArrayBaseType MyBaseType=0;

  if (SizeZ % SizeT !=0)
    {cerr << "Error Writing Header: Internal Z-size " << SizeZ << " is not compatible with time-size " << SizeT << "\n";
     exit(-1);
    }

   // cerr << "Writing Header a: " << SizeZ << ", " << SizeT << "\n";
   /*if (SizeT > 1 )
     {
     if (SizeZ > SizeT)  // If SizeZ==1 write it as apure Z-Stack
        SizeZ /= SizeT;  // correct for incorrect assumption during load
     else
        SizeT = 1;  // write as a pure Z-Stack
     }*/
   // cerr << "Writing Header b: " << SizeZ << ", " << SizeT << "\n";
     
   WriteKhorosHeader(to,"Generated by Reconstruction Set 2007",TypeString(MyBaseType),SizeX,SizeY,SizeZ/SizeT,Elem,SizeT);
}

ofstream * KWOpen(const char * filename,int Elem=1) // Opens an outputfile in KDF for writing
{
  ofstream * to=WOpen(filename);
  KHeader(to,Elem);
  return to;
}

/// saves array to disk under name ("filename") in KDF format.
void  KSave(const char * filename)       // Saves data in Khoros Data format 1,3, 1994 Version 2
{
  ofstream * to = KWOpen(filename);

  Write(to);

  to->close();
  delete to;
} 

ofstream * WOpen(const char * filename) // Opens an outputfile in KDF for writing
{
  cerr << "TArray OpenOutputFile(KDF): " << filename << "\n";
  ofstream * file= new ofstream(filename);
  if (! file) 
    {
      cerr << "Couldn't open file " << filename << " for writing !!\n";
      exit(-1);
    }
  return file;
}

public:   // these are legal to use :

ofstream* DWOpen(int kflag,const char * filename,int Elem=1)
{
if (kflag)
    return KWOpen(filename,Elem);
else
    return WOpen(filename);
}
 

/** writes array to a file in raw-format. Stream "to" has to be set to right position before. */
void  Write(ofstream * to)
{
  if (! to)
    cerr << "Fatal Error : Output stream not existent for writing !\n",exit(-1);

  to->write((char *) TheArray.StartPos(),TheArray.Size());
  if (to->bad()) cerr << "Error writing file, no space on disk \n";
}

void DHeader(int kflag,ofstream & to,int Elements=1,int SizeT=1) // write header depending on the k-flag
{
  if (! to)
    {
      cerr << "File not open  for writing !!\n";
      exit(-1);
    }

  if (kflag)
    {
      KHeader(&to,Elements,SizeT);
    }
}

void DSave(int kflag,const char * filename) // Save data in format dependent on kflag
{
if (kflag)
  {
    KSave(filename);
  }
else
  {
    Save(filename);
  }
}


/// loads a file ("filename") from disk in LSM data format.
int LSMLoad(const char * filename, int &SizeX,int &SizeY, int &SizeZ, int &SizeT,
            int ElementNr=0, int & StartZ=0, int & StopZ=-1,int stepZ=1,
            TArray3d<double> * TimeStamps=0, bool to12bits= false, TArray3d * ROIs=0, TArray3d * DatROIs=0, bool NoDisableROIs=false,
	    int & StartTime=0,int & StopTime=-1,int stepTime=1)  // ROI image will be generated if supplied
{
  int SizeE=1;
  if ( * filename == 0)  // empty string
    { cerr << "No filename given! \n";
      exit (-1);
     }

  ifstream file(filename);

  if (! file) {
    cerr << "couldn't open datafile " << filename << "\n";
    exit (-1);
    } 

  LSMParser mylsm(& file, NoDisableROIs);
  // LSMParser mylsm(filename);

  mylsm.ParseHeader();  // also parses first directory
  SizeX=mylsm.SizeX;
  SizeY=mylsm.SizeY;
  SizeT=mylsm.SizeT;
  SizeE=mylsm.SizeE;

  if (StartZ <0) StartZ=0;
  if (StartTime <0) StopTime=0;

  if (StopZ < 0) 
    StopZ=mylsm.SizeZ - 1;
  if (StopTime < 0) 
    StopTime = mylsm.SizeT-1;
  SizeZ = (StopZ - StartZ)/ stepZ  + 1;
  SizeT = (StopTime - StartTime)/ stepTime + 1;

  if (stepZ <1) stepZ=1;
  if (stepTime <1) stepTime=1;

  if (SizeX != GetSize(0)
      || SizeY != GetSize(1)
      || SizeZ*SizeT != GetSize(2))
    {
      cerr << "LSMLoad has to resize array to " << SizeX <<" x " << SizeY << " x " << SizeZ*SizeT <<"\n";
      Resize(SizeX,SizeY,SizeZ*SizeT);
    }

  if (ElementNr>=SizeE)
  {
    cerr << "Error : LSMLoad requested element " << ElementNr << " out of " << SizeE << " elements in file " << filename << "\n";
    exit(-1);
    }

  cout << "\nLSMLoad loading element " << ElementNr << " out of " << SizeE << " elements in file " << filename << "\n";

  char ppid[10];
  string tmpname("/tmp/buffer");
  sprintf(ppid,"%d",getpid());
  tmpname += ppid;
  tmpname += ".tif";

  for (int j=0;j<StartTime;j++)
   for (int jj=0;jj<mylsm.SizeZ;jj++)
     mylsm.ReadDir();  // parse next "step" directories to reach next valid slice for import

  for (int t=0;t<SizeT;t++) //Size Z * Time points
  {
  for (int j=0;j<StartZ;j++)
     mylsm.ReadDir();  // ignore next "StartFrame" data-directories to reach first valid slice for import
  for (int z=0;z<SizeZ;z++) //Size Z * Time points
  {
    //cout << "Reading Slice Nr. " << z << " (of " << SizeZ << "), timepoint " << t << " (of " << SizeT << ") \n" << flush;
    mylsm.SliceToBuffer(ElementNr);
    mylsm.WriteBufferToTif(tmpname.c_str(),ElementNr);
    mylsm.ReadTiff(tmpname.c_str(), (unsigned char *) Pointer(0,0,z+t*SizeZ),sizeof(ArrayBaseType));
    for (int j=0;j<stepZ;j++)
      mylsm.ReadDir();  // parse next "step" directories to reach next valid slice for import
  }
    for (int j=0;j<(stepTime-1);j++)
     for (int jj=0;jj<mylsm.SizeZ;jj++)
       mylsm.ReadDir();  // parse next "step" directories to reach next valid slice for import
  }

  if (TimeStamps)
    {
      //cout << "resizing\n" << SizeZ/mylsm.SizeZ << flush;
      TimeStamps->Resize(SizeT,1,1);
      //cout << "resized\n" << flush;
      mylsm.GetTimeStamps(TimeStamps->Pointer(0,0,0),TimeStamps->GetSize(0), StartTime, stepTime, SizeE);
    }

  // cout << "printing notes\n" << flush;
  mylsm.PrintNotes();
  // cout << "printed notes\n" << flush;
  if (to12bits && mylsm.DatType == 1)  // 8-bit data has to be rescaled to 12 bit
    Mul(ArrayBaseType(16));

  // cout << "generating ROIs\n" << flush;
  if (ROIs != 0)
  {
    int BRois=mylsm.GetBleachROICnt();
    cout << "generating ROIs\n" << "x,y,z: " << SizeX << ", " << SizeY << ", " << BRois << "\n" << flush;
    if (BRois<1) {
	BRois=1;
	cerr << "WARNING: No BleachROI found\n";
	}
    ROIs->Resize(SizeX,SizeY,BRois);
    for (int n=0;n<mylsm.GetBleachROICnt();n++)
     for (int y=0;y<SizeY;y++)
      for (int x=0;x <SizeX;x++)
      {
        if (mylsm.IsInBleachROI(x,y,n))
          ROIs->SetValue(x,y,n,1);
        else
          ROIs->SetValue(x,y,n,0);
      }
  }
  // cout << "generating data ROIs\n" << flush;
  if (DatROIs != 0)
  {
    int Rois=mylsm.GetDataROICnt();
    cout << "generating ROIs\n" << "x,y,z: " << SizeX << ", " << SizeY << ", " << Rois << "\n" << flush;
    if (Rois<1) {
	Rois=1;
	cerr << "WARNING: No DataROI found\n";
	}
    DatROIs->Resize(SizeX,SizeY,Rois);
    for (int n=0;n<mylsm.GetDataROICnt();n++)
     for (int y=0;y<SizeY;y++)
      for (int x=0;x <SizeX;x++)
      {
        if (mylsm.IsInDataROI(x,y,n))
          DatROIs->SetValue(x,y,n,1);
        else
          DatROIs->SetValue(x,y,n,0);
      }
  }
  // cout << "generated data ROIs\n" << flush;
  system(("rm "+tmpname).c_str());
  // cout << "removed tmp file:" << ("rm "+tmpname).c_str() << "\n" << flush;
  file.close();
  return SizeE;
}

/// loads a file ("filename") from disk in ASI data format.
void ASILoad(int filetype , const char * filename, int posx=-1, int posy=-1, int posz=-1, int headersize=-1)        
{
  char dummy;
  int ValueDimX,ValueDimY,ValueDimZ;  // Z is the color direction

  if ( * filename == 0)  // empty string
     return;           // dont do anything

  cerr << "TArray LoadfromDisk (ASI) : " << filename << "\n";

  ifstream file(filename);

  if (! file) {
    cerr << "couldn't open datafile " << filename << "\n";
    exit (-1);
    }

  if (filetype <= 1)
    ReadRAWHeader(& file,ValueDimX,ValueDimY,ValueDimZ);
  else if (filetype == 2)
    ReadSDIHeader(& file,ValueDimX,ValueDimY,ValueDimZ,0,posx,posy,posz,headersize);
  else if (filetype == 3)
    ReadASIHeader(& file,ValueDimX,ValueDimY,ValueDimZ);
  else 
    ReadSDIHeader(& file,ValueDimX,ValueDimY,ValueDimZ,1,posx,posy,posz,headersize);

  if (ValueDimX != GetSize(0)
      || ValueDimY != GetSize(1)
      || ValueDimZ != GetSize(2))
    {
      cerr << "KLoad has to resize array to " << ValueDimX <<" x " << ValueDimY << " x " << ValueDimZ <<"\n";
      Resize(ValueDimX,ValueDimY,ValueDimZ);
    }

  Read(& file);

  file.read(&dummy,sizeof(char));
  if (! file.eof()) 
    cerr << "Warning ! Datafile " << filename << " longer than expected !\n";

  file.close();
  return;
}


/** loads a file ("filename") from disk in data format (KDF or RAW) depending on kflag. Return Number of Elements in file.
    This may change the Sizes (SizeX, SizeY,SizeZ), or rezise the array according to them. **/ 
int DLoad(int kflag,const char * filename,const char * ReqType,int * SizeX,int * SizeY=0, int * SizeZ=0,int ElementNr=-1,int offset=0)
{
int rtv=1;
 if (kflag)
   {
        rtv=KLoad(filename,ReqType,ElementNr);  // ReqType is the datatype of this array
        (* SizeX) = GetSize(0);  

        if (SizeY != 0)
          (* SizeY) = GetSize(1);
        if (SizeZ != 0)
          (* SizeZ) = GetSize(2);
   }
 else
   {
        Resize(* SizeX,* SizeY,* SizeZ);
        Load(filename,offset,(char *) ReqType,ElementNr);  // The meaning is changed here! ReqType is the type of the data to read
   }
 return rtv;
} 

/** A DLoad version without the requirement for a typestring. The required type is the datatype to convert to */
int DLoad(int kflag,const char * filename,int * SizeX,int * SizeY=0, int * SizeZ=0,int ElementNr=-1,int offset=0)
{
  ArrayBType dummy;
  return DLoad(kflag,filename,TypeString(dummy),SizeX,SizeY,SizeZ,ElementNr,offset);
}

void LoadSlice(const char * filename,int slicenr,int offset=0)  // Only in Raw Format so far
{ 
  FILE *file;

  char dummy;

  if ( * filename == 0)  // empty string
     return;             // dont do anything

  cerr << "TArray LoadfromDisk : " << filename << ", AllDim :" << TheArray.Size() << "\n";

  if ((file = fopen(filename,"rb")) == NULL) {
    cerr << "couldn't open datafile " << filename << "\n";
    exit (-1);
    }

  if (offset > 0)
    fseek(file,offset,SEEK_SET);


  int SliceSize = GetSize(0)*GetSize(1);

  if (fread(TheArray.StartPos()+slicenr * SliceSize,SliceSize*sizeof(ArrayBaseType),1,file) != (unsigned) 1)
    {
      cerr << "Reached EOF unexpectedly !\n";
      exit(-1);
    }

  if (fread(&dummy,sizeof(char),1,file) == 1)
    cerr << "Warning ! Datafile " << filename << " longer than expected !\n";

  fclose(file);
}


};



// ++++++++++++++++++++++++++++++++++++++++++++++ Derived Classes +++++++++++++++++++++++++++++++++++++++++++++++

/// class with a possibility for simple iteration. (First, Next, Change, Add)
template<class TheArrayType> class TArrayIter3d : public TArray3d<TheArrayType> {
public:
  typedef typename TArray3d<TheArrayType>::ArrayBaseType ArrayBaseType;  // get the basis type of myself

protected:
  int LinIndexX,LinIndexY,LinIndexZ;

  ArrayBaseType * ActPtr;

public:


TArrayIter3d() : TArray3d<TheArrayType>()
{}

TArrayIter3d(int SizeX,int SizeY,int SizeZ) : TArray3d<TheArrayType>(SizeX,SizeY,SizeZ)
{}

TArrayIter3d(ArrayBaseType * ExistingArray,int SizeX,int SizeY,int SizeZ) : 
  TArray3d<TheArrayType>(ExistingArray,SizeX,SizeY,SizeZ)
{}

int GetSize(int dim)
{return TArray3d<TheArrayType>::GetSize(dim);}


  /// gives First Array value, inits counters, returns 0 Vec
ArrayBaseType FirstValue()
{
  LinIndexX = GetSize(0);
  LinIndexY = GetSize(1);
  LinIndexZ = GetSize(2);
  ActPtr =  TArray3d<TheArrayType>::TheArray.StartPos();
  return * ActPtr;
}

  /// Gives the next Value and the Relative Position to the last one. Returns 1 if no next value exists.
int NextValue(ArrayBaseType * ValToChange)
{
  ActPtr++;  // Step to next Voxel

  if (--LinIndexX == 0)
      if (--LinIndexY == 0)
	if (--LinIndexZ == 0)
	  return -1;   // ready
        else
	  {
	    LinIndexY=GetSize(1);
	    LinIndexX=GetSize(0);
	    (* ValToChange) = (* ActPtr);
	    return 2;
	  }
      else
	{
	  LinIndexX=GetSize(0);
          (* ValToChange) = (* ActPtr);
          return 1;
        }
	  
   else
     {
       (* ValToChange) = (* ActPtr);
       return 0;
     }
}

/// Changes actual value to "NewVal"
void ChangeValue(ArrayBaseType NewVal) {  (* ActPtr)=NewVal; }
/// Changes "NewVal" to actual value
void AddValue(ArrayBaseType NewVal) {  (* ActPtr) += NewVal; }

};

/// code dublication !! Complex to real conversion. Arguments must be pointers !
template <class FirstArr, class SecondArr>
void Copy(FirstArr * arr1,SecondArr * arr2, int SrcXAdjust=0,int DestXAdjust=0) // copy part of array centered
{
  int Dim1X=arr1->GetSize(0) + SrcXAdjust,  // if the logical size is smaller 
      Dim1Y=arr1->GetSize(1),
      Dim1Z=arr1->GetSize(2);

  int Dim2X=arr2->GetSize(0) + DestXAdjust,
      Dim2Y=arr2->GetSize(1),
      Dim2Z=arr2->GetSize(2);

  int dx1=0,dy1=0,dz1=0,dx2=0,dy2=0,dz2=0,x,y,z,SizeX=Dim2X,SizeY=Dim2Y,SizeZ=Dim2Z;

  if (Dim1X < Dim2X)
    SizeX=Dim1X,dx1=Dim2X/2-Dim1X/2;
  else
    dx2=Dim1X/2-Dim2X/2;

  if (Dim1Y < Dim2Y)
    SizeY=Dim1Y,dy1=Dim2Y/2-Dim1Y/2;
  else
    dy2=Dim1Y/2-Dim2Y/2;

  if (Dim1Z < Dim2Z)
    SizeZ=Dim1Z,dz1=Dim2Z/2-Dim1Z/2;
  else
    dz2=Dim1Z/2-Dim2Z/2;

  for (z=0;z<SizeZ;z++)
    for (y=0;y<SizeY;y++)
      for (x=0;x<SizeX;x++)
	arr2->SetValue(x+dx1,y+dy1,z+dz1,arr1->Value(x+dx2,y+dy2,z+dz2));
}


/// copies by cyclic mapping starting from edge
template <class FirstArr, class SecondArr, class AFunctor>
void CyclicCopy(FirstArr & arr1,SecondArr & arr2,
		int shiftX=0,
		int shiftY=0,
		int shiftZ=0,
		AFunctor Fun = project1st<typename FirstArr::ArrayBaseType, typename SecondArr::ArrayBaseType>()
		) // copy with cyclic wrap-around , if needed, apply a function
{
  int Dim1X=arr1.GetSize(0),
      Dim1Y=arr1.GetSize(1),
      Dim1Z=arr1.GetSize(2);

  int Dim2X=arr2.GetSize(0),
      Dim2Y=arr2.GetSize(1),
      Dim2Z=arr2.GetSize(2);

  int x,y,z;

  for (z=0;z<Dim2Z;z++)
    for (y=0;y<Dim2Y;y++)
      for (x=0;x<Dim2X;x++)
              arr2.SetValue(x,y,z,
		       Fun(arr1(((x+shiftX)%Dim1X),((y+shiftY)%Dim1Y),((z+shiftZ)%Dim1Z)),
			   arr2(x,y,z))
		       );
}

/// copies by left alignment, the rest will be left as it is
template <class PFirstArr, class PSecondArr>
void LCopy(PFirstArr arr1,PSecondArr arr2) // copy with cyclic wrap-around 
{
  int DimX= (arr1->GetSize(0) < arr2->GetSize(0)) ? arr1->GetSize(0) : arr2->GetSize(0),   // always take the smaller one
      DimY= (arr1->GetSize(1) < arr2->GetSize(1)) ? arr1->GetSize(1) : arr2->GetSize(1),
      DimZ= (arr1->GetSize(2) < arr2->GetSize(2)) ? arr1->GetSize(2) : arr2->GetSize(2);

  int x,y,z;

  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	arr2->SetValue(x,y,z,arr1->Value(x,y,z));
}

/// copies by cyclic mapping starting from edge
template <class PFirstArr, class PSecondArr>
void BlendCopy(PFirstArr arr1,PSecondArr arr2) // copy with blending the last slice
{
  int Dim1X=arr1->GetSize(0),
      Dim1Y=arr1->GetSize(1),
      Dim1Z=arr1->GetSize(2);

  int Dim2X=arr2->GetSize(0),
      Dim2Y=arr2->GetSize(1),
      Dim2Z=arr2->GetSize(2);

  int x,y,z,xb,yb,zb;
  CalcFloatType xweight,yweight,zweight;

  for (z=0;z<Dim2Z;z++)
    for (y=0;y<Dim2Y;y++)
      for (x=0;x<Dim2X;x++)
	{
	  if (x >= Dim1X) xb=Dim1X-1,xweight=(x-(Dim1X-1))/CalcFloatType(1+Dim2X-Dim1X);
	  else xb=x,xweight=0.0;
	  if (y >= Dim1Y) yb=Dim1Y-1,yweight=(y-(Dim1Y-1))/CalcFloatType(1+Dim2Y-Dim1Y);
	  else yb=y,yweight=0.0;
	  if (z >= Dim1Z) zb=Dim1Z-1,zweight=(z-(Dim1Z-1))/CalcFloatType(1+Dim2Z-Dim1Z);
	  else zb=z,zweight=0.0;
 
	  arr2->SetValue(x,y,z,
			 arr1->Value(xb,yb,zb)*(1.0-xweight)*(1.0-yweight)*(1.0-zweight)+
			 arr1->Value(0,yb,zb)*xweight*(1.0-yweight)*(1.0-zweight)+
			 arr1->Value(xb,0,zb)*yweight*(1.0-xweight)*(1.0-zweight)+
			 arr1->Value(xb,yb,0)*zweight*(1.0-xweight)*(1.0-yweight)
			 );
	}
}



/// Adds Real part of First Array into Second Array
template <class FirstArr, class SecondArr>
void AddReal(FirstArr * arr1,SecondArr *  arr2)
{
  int Dim1X=arr1->GetSize(0),
      Dim1Y=arr1->GetSize(1),
      Dim1Z=arr1->GetSize(2);

  int Dim2X=arr2->GetSize(0),
      Dim2Y=arr2->GetSize(1),
      Dim2Z=arr2->GetSize(2);

  int x,y,z;

  if ((Dim1X != Dim2X) ||
      (Dim1Y != Dim2Y) ||
      (Dim1Z != Dim2Z))
    cerr << "Error Add(Arr1,Arr2) wrong sizes !\n" << flush, exit(-1);

  for (z=0;z<Dim1Z;z++)
    for (y=0;y<Dim1Y;y++)
      for (x=0;x<Dim1X;x++)
	arr2->SetValue(x,y,z,arr2->Value(x,y,z)+CastToReal(arr1->Value(x,y,z)));
}

/// code dublication ?! Complex to real conversion. Array is taken to be centered
template <class FirstArr, class SecondArr>
void CopyCtoR(FirstArr * arr1,SecondArr * arr2, int SrcXAdjust=0,int DestXAdjust=0)
{
  int Dim1X=arr1->GetSize(0) + SrcXAdjust,
      Dim1Y=arr1->GetSize(1),
      Dim1Z=arr1->GetSize(2);

  int Dim2X=arr2->GetSize(0) + DestXAdjust,
      Dim2Y=arr2->GetSize(1),
      Dim2Z=arr2->GetSize(2);

  int dx1=0,dy1=0,dz1=0,dx2=0,dy2=0,dz2=0,x,y,z,SizeX=Dim2X,SizeY=Dim2Y,SizeZ=Dim2Z;

  if (Dim1X < Dim2X)
    SizeX=Dim1X,dx1=(Dim2X-Dim1X)/2;
  else
    dx2=(Dim1X-Dim2X)/2;

  if (Dim1Y < Dim2Y)
    SizeY=Dim1Y,dy1=(Dim2Y-Dim1Y)/2;
  else
    dy2=(Dim1Y-Dim2Y)/2;

  if (Dim1Z < Dim2Z)
    SizeZ=Dim1Z,dz1=(Dim2Z-Dim1Z)/2;
  else
    dz2=(Dim1Z-Dim2Z)/2;

  for (z=0;z<SizeZ;z++)
    for (y=0;y<SizeY;y++)
      for (x=0;x<SizeX;x++)
	arr2->SetValue(x+dx1,y+dy1,z+dz1,CastToReal(arr1->Value(x+dx2,y+dy2,z+dz2)));   // abs is not correct because of negative real values !

  return;
}

/// copies from Matlab style complex into a complex valued array, Array is taken to be centered
template <class BaseT, class SecondArr>
void CopyMCtoC(BaseT * pR,BaseT * pI,SecondArr * arr2,int * dims1=0, int SrcXAdjust=0,int DestXAdjust=0)
{
  int Dim1X,Dim1Y,Dim1Z;
  if (dims1==0)
      {Dim1X=arr2->GetSize(0) + SrcXAdjust;
      Dim1Y=arr2->GetSize(1);
      Dim1Z=arr2->GetSize(2);}
  else {
      Dim1X=dims1[0] + SrcXAdjust;
      Dim1Y=dims1[1];
      Dim1Z=dims1[2]; }

  int Dim2X=arr2->GetSize(0) + DestXAdjust,
      Dim2Y=arr2->GetSize(1),
      Dim2Z=arr2->GetSize(2);

  int dx1=0,dy1=0,dz1=0,dx2=0,dy2=0,dz2=0,x,y,z,SizeX=Dim2X,SizeY=Dim2Y,SizeZ=Dim2Z;

  if (Dim1X < Dim2X)
    SizeX=Dim1X,dx1=(Dim2X-Dim1X)/2;
  else
    dx2=(Dim1X-Dim2X)/2;

  if (Dim1Y < Dim2Y)
    SizeY=Dim1Y,dy1=(Dim2Y-Dim1Y)/2;
  else
    dy2=(Dim1Y-Dim2Y)/2;

  if (Dim1Z < Dim2Z)
    SizeZ=Dim1Z,dz1=(Dim2Z-Dim1Z)/2;
  else
    dz2=(Dim1Z-Dim2Z)/2;

  for (z=0;z<SizeZ;z++)
    for (y=0;y<SizeY;y++)
      for (x=0;x<SizeX;x++)
	arr2->SetValue(x+dx1,y+dy1,z+dz1,complex<BaseT>(*(pR++),*(pI++)));

  return;
}

/// copies into Matlab style complex from a complex valued array
template <class BaseT, class SecondArr>
void CopyCtoMC(SecondArr * arr2, BaseT * pR,BaseT * pI)
{
  int DimX=arr2->GetSize(0), //  + DestXAdjust
      DimY=arr2->GetSize(1),
      DimZ=arr2->GetSize(2);

  int x,y,z,SizeX=DimX,SizeY=DimY,SizeZ=DimZ;

  for (z=0;z<SizeZ;z++)
    for (y=0;y<SizeY;y++)
      for (x=0;x<SizeX;x++)
	{
	*(pR++)=real(arr2->Value(x,y,z));
	*(pI++)=imag(arr2->Value(x,y,z));
	}

  return;
}

// #define DynArray(Type)  TArray3d <DynamicArray< Type > >

// #define StatArray(Type,X,Y,Z) TArray3d <StaticArray3d< Type,X,Y,Z>,Clipperator >

#endif

