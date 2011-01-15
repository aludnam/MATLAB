#ifndef basearrays_h	// -*- C++ -*- 
#define basearrays_h

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

#include <assert.h>

template<class T>
T CLIPBTW(T Low,T High,T Val) 
{return (((Val) < (Low)) ? (Low) : ((Val) > (High)) ? (High) : (Val));}  // Proves to be faster

template<class T>
T CLIPZERO(T Val) {return ((Val > 0)*Val);}
// #define CLIPZUN(High,Val) (High - (Val < High) * (High-CLIPZERO(Val)))

template<class T> 
T CLIPZUN(T High,T Val) {return (((Val) < (0)) ? (0) : ((Val) > (High)) ? (High) : (Val));}  // Proves to be faster

template<class T>
int RINT(T X) {return ((int) ((X)+0.5));}       // ONLY FOR POSITIVE NUMBERS !!!

static const int MAXARRAYDIM=3;

/// Basetype of every array, static or dynamic. 
class ArrayAnchestor {        

protected:
int Dimensions;

void checkbounds(int x, int y, int z)
    {
      assert (x>=0);
      assert (y>=0);
      assert (z>=0);
      assert (x<DimVec[0]);
      assert (y<DimVec[1]);
      assert (z<DimVec[2]);
    }

public:
int DimVec[MAXARRAYDIM+1];

ArrayAnchestor()
  {
    Dimensions=MAXARRAYDIM;
  }

int Dimension(void)
{
  return Dimensions;
}


double DoCenterX(double Val)   // CalcFloatType ??
{
  return (Val+DimVec[0]/2.0);
}

double DoCenterY(double Val)
{
  return (Val+DimVec[1]/2.0);
}

double DoCenterZ(double Val)
{
  return (Val+DimVec[2]/2.0);
}


double DoClipX(double Val) {return CLIPZUN(DimVec[0]-1.0,Val);}
double DoClipY(double Val) {return CLIPZUN(DimVec[1]-1.0,Val);}
double DoClipZ(double Val) {return CLIPZUN(DimVec[2]-1.0,Val);}

int DoClipX(int Val) {return CLIPZUN(DimVec[0]-1,Val);}
int DoClipY(int Val) {return CLIPZUN(DimVec[1]-1,Val);}
int DoClipZ(int Val) {return CLIPZUN(DimVec[2]-1,Val);}

};


// Now Templates for PrimitiveArrayType


/** Type for a 3-dimensional static array. The sizes have to be known in advance and must be
 specified as template arguments */

template<class ArrayBaseType, const int SIZEX, const int SIZEY, const int SIZEZ> class StaticArray3d : public ArrayAnchestor {
public:
  typedef ArrayBaseType BaseType;  // export this typename to users. They can refer by classname::BaseType
  void checkpointer(ArrayBaseType * p)
    {
      assert(p > StartPos());
      assert(p <= EndPos());
    }

  // static  const int SizeX=SIZEX;   // makes Template-Parameters accessable
  // static  const int SizeY=SIZEY;
  // static  const int SizeZ=SIZEZ;   // some compilers dont like this !

protected:
ArrayBaseType TheArray[SIZEZ][SIZEY][SIZEX];

  ///  clears the array completely (by set to 0)
void Clr(void)
  {
    int x,y,z;
    for (z=0;z<SIZEZ;z++)
    for (y=0;y<SIZEY;y++)
    for (x=0;x<SIZEX;x++)
      {
	TheArray[z][y][x]=0;
      }
  }

public:

  /// empty (normal) constructor
StaticArray3d()
  {
    DimVec[0]=SIZEX;
    DimVec[1]=SIZEY;
    DimVec[2]=SIZEZ;
    DimVec[3]=0;
    Clr();
  }

  /// dummy constructor, will yield error
StaticArray3d(int SX,int SY,int SZ)
  {
    cout << "ERROR : Static Array of size " << SX << "x" << SY << "x" << SZ << " allocated, giving Dynamic Sizes for it.\n";
    exit(-1);
  }

  /// dummy constructor, will yield error
StaticArray3d(ArrayBaseType * ExistingArray,int SX,int SY, int SZ)
  {
    cout << "ERROR : Tried to init ( " << SX << "x" << SY << "x" << SZ << " ) Static Array  from existing array at "<<ExistingArray <<".\n";
    exit(-1);
  }

  /// dummy "resize"-fuction, will yield error
void Resize(int SX,int SY,int SZ)
  {
    SX=0,SY=0,SZ=0; // to prevent "unused" warning
    cerr << "ERROR : Tried to resize static array \n";
    exit(-1);
  }

  /// dummy "resize"-fuction, will yield error
void Resize(ArrayAnchestor * Ap)
  {
    Ap=0; // to prevent "unused" warning
    cerr << "ERROR : Tried to resize static array \n";
    exit(-1);
  }


  /// returns size of array in chars
size_t Size(void)
  {
    return sizeof(TheArray);
  }

/// overloaded [ operator returning a reference
ArrayBaseType & operator() (int x,int y=0,int z=0)
{ checkbounds(x,y,z); return TheArray[z][y][x];}

  /// Sets a value in the array at (x,y,z) to "v"
void SetValue(int x,int y,int z,ArrayBaseType v) {checkbounds(x,y,z); TheArray[z][y][x] = v;}
  /// returns value of the array at (x,y,z) 
ArrayBaseType GetValue(int x,int y,int z) {checkbounds(x,y,z); return TheArray[z][y][x];}
  /// returns poiter to value of the array at (x,y,z) 
ArrayBaseType * GetPValue(int x,int y,int z) {checkbounds(x,y,z); return & TheArray[z][y][x];}

  /// returns poiter to first element of array 
ArrayBaseType * StartPos(void) {return & TheArray[0][0][0];}
  /// returns poiter to last element of array 
ArrayBaseType * EndPos(void) {return & TheArray[SIZEZ-1][SIZEY-1][SIZEX-1];}

ArrayBaseType * begin(void) {return StartPos();}
ArrayBaseType * end(void) {return EndPos()+1;}
size_t size(void) {return end()-begin();}

};

/// as static 3darray but with some special functions
template<class ArrayBaseType,int SizeX,int SizeY> class StaticArray2d : public StaticArray3d<ArrayBaseType,SizeX,SizeY,1> { 
public:
StaticArray2d() : StaticArray3d<ArrayBaseType,SizeX,SizeY,1>() {}
StaticArray2d(int SiX,int SiY,int SiZ) : StaticArray3d<ArrayBaseType,SizeX,SizeY,1>(SiX,SiY,SiZ) {};
StaticArray2d(ArrayBaseType * ExistingArray,int SiX,int SiY, int SiZ) : StaticArray3d<ArrayBaseType,SizeX,SizeY,1>(ExistingArray,SiX,SiY,SiZ) {};
};

/// as static 2darray but with some special functions
template<class ArrayBaseType,int SizeX> class StaticArray1d : public StaticArray2d<ArrayBaseType,SizeX,1> { 
public:
StaticArray1d() : StaticArray2d<ArrayBaseType,SizeX,1>() {}
StaticArray1d(int SiX,int SiY,int SiZ) : StaticArray2d<ArrayBaseType,SizeX,1>(SiX,SiY,SiZ) {};
StaticArray1d(ArrayBaseType * ExistingArray,int SiX,int SiY, int SiZ) : StaticArray2d<ArrayBaseType,SizeX,1>(ExistingArray,SiX,SiY,SiZ) {};
};


/** Type for a 3-dimensional dynamic array. The sizes have to be known when constructing it and must be
 specified as constructor arguments. This class is only very slightly slower than the static class */

template<class ArrayBaseType> class DynamicArray : public ArrayAnchestor {

public:
  typedef ArrayBaseType BaseType;  // export this typename to users. They can refer by classname::BaseType

  void checkpointer(ArrayBaseType * p)
    {
      assert(p >= StartPos());
      assert(p <= EndPos());
    }

protected:

 /// here the array will be stored
  ArrayBaseType * TheArray;

int AllSize;
ArrayBaseType * LastElem;
int DidAllocate;

 /// will set the array contens to 0
void Clr(void)
  {
    int i;
    for (i=0;i<AllSize;i++)
      {
	TheArray[i]=0;
      }
  }

public:

 /// an empty array (size==0) is created, which can be resized later on
DynamicArray() : DidAllocate(0)
  {
    DimVec[0]=0;
    DimVec[1]=0;
    DimVec[2]=0;
    DimVec[3]=0;
    TheArray=(ArrayBaseType *) 0;     // Resize has to be called later !
 //   cout << "ERROR : Dynamic Array allocated without giving Sizes for it.\n";
 //   exit(-1);
  }

 /// resizes the array. The contens will be set to zero.
void Resize(int SizeX,int SizeY,int SizeZ)     // contens is lost !
  {
    if ((! DidAllocate) || (DimVec[0] != SizeX) ||
	(DimVec[1] != SizeY) ||
	(DimVec[2] != SizeZ))
      {
	if (DidAllocate != 0 && TheArray != 0) delete[] TheArray;

    AllSize=SizeX*SizeY*SizeZ;
	if (! (TheArray=new ArrayBaseType[AllSize]))
	    cerr << "Fatal Error : New failed in resize !\n",exit(-1);
	DidAllocate=1;
	DimVec[0]=SizeX;
	DimVec[1]=SizeY;
	DimVec[2]=SizeZ;
	DimVec[3]=0;
	LastElem= & TheArray[AllSize-1]; 
      }
    Clr();   // in any case
  }

 /// resizes the array. The contens will be set to zero. The other array is used as template for sizes.
void Resize(ArrayAnchestor * AnotherArray)     // contens is lost !
  {
    Resize(AnotherArray->DimVec[0],AnotherArray->DimVec[1],AnotherArray->DimVec[2]);
  }

  /// constructs array of the named sizes
DynamicArray(int SizeX,int SizeY,int SizeZ) : DidAllocate(0)
  {
    TheArray=(ArrayBaseType *) 0;     // Resize has to be called later !
    Resize(SizeX,SizeY,SizeZ);
  }

  /// constructs array of the named sizes, but doesn't allocate any memory and incorporates given pointer instead
DynamicArray(ArrayBaseType * ExistingArray,int SizeX,int SizeY,int SizeZ)
  {
    AllSize=SizeX*SizeY*SizeZ;
    TheArray=ExistingArray;
    DidAllocate=0;
    DimVec[0]=SizeX;
    DimVec[1]=SizeY;
    DimVec[2]=SizeZ;
    DimVec[3]=0;
    LastElem= & TheArray[AllSize-1];
  }

  /// destructor
~DynamicArray()
  {
    if (DidAllocate) 
	{
		delete[] TheArray;
		DidAllocate = 0;
	}
  }

  /// returns size of array in chars
size_t Size(void)
  {
    return sizeof(ArrayBaseType)*AllSize;
  }

 /// overloaded [ operator returning a reference
ArrayBaseType & operator() (const int x,const int y=0,const int z=0)
{ checkbounds(x,y,z); return TheArray[z*DimVec[0]*DimVec[1]+y*DimVec[0]+x];}

  /// Sets a value in the array at (x,y,z) to "v"
void SetValue(const int x,const int y,const int z,ArrayBaseType v) {checkbounds(x,y,z); TheArray[z*DimVec[0]*DimVec[1]+y*DimVec[0]+x] = v;}
  /// returns value of the array at (x,y,z) 
ArrayBaseType   GetValue(const int x,const int y=0,const int z=0) {checkbounds(x,y,z); return TheArray[z*DimVec[0]*DimVec[1]+y*DimVec[0]+x];}
  /// returns poiter to value of the array at (x,y,z) 
ArrayBaseType * GetPValue(const int x,const int y,const int z) {checkbounds(x,y,z); return & TheArray[z*DimVec[0]*DimVec[1]+y*DimVec[0]+x];}

  /// returns poiter to first element of array 
ArrayBaseType * StartPos(void) {return TheArray;}
  /// returns poiter to last element of array 
ArrayBaseType * EndPos(void) {return LastElem;}

ArrayBaseType * begin(void) {return TheArray;}
ArrayBaseType * end(void) {return LastElem+1;}
size_t size(void) {return end()-begin();}

};




#endif
