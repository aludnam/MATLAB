#ifndef tarrayfft_h  		// -*- C++ -*- 
#define tarrayfft_h

#include <complex>
#include "stlfixes.h"
#include "rawarray.h"

#ifdef THREADS
const int NTHREADS=2;
// #include "srfftw_threads.h"
// #include "sfftw_threads.h"
#include "rfftw_threads.h"
#include "fftw_threads.h"
#else
#ifdef OLDFFTNAME
#include "rfftw.h"
#include "fftw.h"
#else
#include "srfftw.h"
#include "sfftw.h"
#endif
#endif

// The 3-D array is taken to be the most general case for now.

// Lower Dimension-Arrays are derived from it, so other templates work on them too, if TArray is required.

/// to make arbitrary numbers (complex or real) real

// This class adapts basically any class which is derived from tarray in rawarray.h


template<class ArrayBType, class ArrayType = TArray3d<ArrayBType> >    // another array has to be used as base, which will then be fftable
class TFFTArray : public ArrayType { 
public:  
  typedef ArrayBType FArrayBaseType;  // get the basis type of myself and make it accessable to the outside

protected:
  fftwnd_plan FFTPlanFwd;
  fftwnd_plan FFTPlanBwd;
  rfftwnd_plan RFFTPlanFwd;
  rfftwnd_plan RFFTPlanBwd;

  int usedDimensions;

  // typedef typename ArrayType::ArrayBaseType ArrayBaseType;  // get the basis type of myself
  typedef ArrayType MyArrayType;  // make the normal array accessible  

  int * Dims;   // Array for storing the sizes in inverse order (z,y,x)
  bool IsComplex;

  bool TestComplex(complex<float>)
  {
    return true;
  }

  bool TestComplex(float)
  {
    return false;
  }

  bool TestComplex(complex<double>)
  {
    return true;
  }

  bool TestComplex(double)
  {
    return false;
  }


  /*
  template <class BaseT1, class ArrayType2>
  friend void Copy(TFFTArray<BaseT1> * arr1,ArrayType2 * arr2); // these functions need access to member ExtraX

  template <class BaseT1, class BaseT2>
  friend void Copy(TFFTArray<BaseT1> * arr1,TFFTArray<BaseT2> * arr2);

  template <class ArrayType1, class BaseT2>
  friend void CopyRtoC(ArrayType1 *,TFFTArray<BaseT2> *);

  template <class ArrayType1, class BaseT2>
  friend void Copy(ArrayType1 * arr1,TFFTArray<BaseT2> * arr2);

  template <class BaseT1, class BaseT2>
  friend void CopyCtoC(TFFTArray<BaseT1> * arr1,TFFTArray<BaseT2> * arr2);
  template <class ArrayType1, class BaseT2>
  friend void CopyCtoR(ArrayType1 * arr1,TFFTArray<BaseT2> * arr2); 
  */

public:

  int ExtraX;    // number of X positions extra   , would have been greate to keep this private!!! Windows has problems

/// empty constructor
  TFFTArray() : MyArrayType(), FFTPlanFwd(0), FFTPlanBwd(0), RFFTPlanFwd(0), RFFTPlanBwd(0), IsComplex(TestComplex( FArrayBaseType(0))),ExtraX(0)
{
  usedDimensions = ArrayType::TheArray.Dimension();
  Dims = new int[4];
}

/// empty constructor
TFFTArray(int SizeX,int SizeY=1,int SizeZ=1) : MyArrayType(SizeX,SizeY,SizeZ), 
  FFTPlanFwd(0), FFTPlanBwd(0), RFFTPlanFwd(0), RFFTPlanBwd(0), IsComplex(TestComplex( FArrayBaseType(0))),ExtraX(0)
{
  usedDimensions = ArrayType::TheArray.Dimension();
  Dims = new int[4];
  Resize(SizeX,SizeY,SizeZ);  // Now adjust again to needed size
}

/// rezises array to new sizes (will be cleared)
void DoResize(int SizeX,int SizeY=1,int SizeZ=1)
{
  if (usedDimensions == 1)
    {
      Dims[0] = SizeX;
      Dims[1] = 0;
    }
  else if (usedDimensions == 2)
    {
      Dims[0] = SizeY;
      Dims[1] = SizeX;
      Dims[2] = 0;
    }
  else if (usedDimensions == 3)
    {
      Dims[0] = SizeZ;
      Dims[1] = SizeY;
      Dims[2] = SizeX;
      Dims[3] = 0;
    }

  if (! IsComplex)   // If this is a float->complex transformable array
    {
      ExtraX = 2 * (SizeX/2+1) - SizeX;  
      SizeX = 2 * (SizeX/2+1);  // make some space in the last dimension
      cout << "Making extra space in float FFT array \n";
    }


  MyArrayType::Resize(SizeX,SizeY,SizeZ);
}

int GetSize(int dim)
{
	return ArrayType::GetSize(dim);
}

int GetLogicSize (int dim)  // returns the size for real space data (after subtraction of extra allocated positions)
{
  if (dim == 0)
    return GetSize(dim) - ExtraX;
  else
    return GetSize(dim);
}

// Function operates on the range of the real arrays
double MinPosReal(int & mx, int & my, int & mz)
{
  int x,y,z;
  mx=0,my=0,mz=0;
  int sx=GetLogicSize(0);
  int sy=GetLogicSize(1);
  int sz=GetLogicSize(2);
  double minVal = CastToReal(ArrayType::Value(0,0,0));

    for (z=0;z<sz;z++)              // copy into tmpArray
    for (y=0;y<sy;y++)              // copy into tmpArray
    for (x=0;x<sx;x++)              // copy into tmpArray
	{
	if (CastToReal(ArrayType::Value(x,y,z)) < minVal)
			{minVal = CastToReal(ArrayType::Value(x,y,z)); mx=x;my=y;mz=z; }
	}
  return minVal;
}


/// rezises array to new sizes (will be cleared)
void Resize(int SizeX,int SizeY=1,int SizeZ=1)
{
  // cout << "resizing FFTable array \n";

  if (FFTPlanFwd)
    fftwnd_destroy_plan(FFTPlanFwd);
  if (FFTPlanBwd)
    fftwnd_destroy_plan(FFTPlanBwd);
  if (RFFTPlanFwd)
    rfftwnd_destroy_plan(RFFTPlanFwd);
  if (RFFTPlanBwd)
    rfftwnd_destroy_plan(RFFTPlanBwd);

  FFTPlanFwd=0;
  FFTPlanBwd=0;
  RFFTPlanFwd=0;
  RFFTPlanBwd=0;

  DoResize(SizeX,SizeY,SizeZ);   // will pad the last dim by two to make space for fft of a float array

}

/// rezises array to same as *AP (this array will be cleared)
void Resize(TFFTArray * Ap)
{
  Resize(Ap->GetSize(0),Ap->GetSize(1),Ap->GetSize(2));
}

/// rezises array to same as *AP (this array will be cleared)
void Resize(MyArrayType * Ap)
{
  Resize(Ap->GetSize(0),Ap->GetSize(1),Ap->GetSize(2));
}


void AdjustX(int px) // after cyclic operations in float arrays the right half is off by one pixel and the unused line is at position x. This is corrected
{
  if (IsComplex)
    return;

  int DimX=GetSize(0), DimY=GetSize(1), DimZ=GetSize(2);
  px = ArrayType::IntoRange(px,DimX) - ExtraX;

  for (int z=0;z<DimZ;z++)
    for (int y=0;y<DimY;y++)
      for (int x=px;x<DimX-ExtraX;x++)
	ArrayType::SetValue(x,y,z,ArrayType::Value(x+ExtraX,y,z));
}

void rotate(int OfX,int OfY,int OfZ)
{
  MyArrayType::rotate(OfX,OfY,OfZ);
  AdjustX(OfX);
}

 /// (UN) scrambles the fft-data into normal reziprocal space (0 in the middle)
void scramble(int unscramble=0)
{
  int DimX=Dims[usedDimensions-1], DimY=GetSize(1), DimZ=GetSize(2);
  int OfX=DimX/2, 
      OfY=DimY/2,
      OfZ=DimZ/2;

  if (unscramble)
    {
      OfX=-OfX; 
      OfY=-OfY; 
      OfZ=-OfZ; 
    }

  rotate(OfX,OfY,OfZ);
}

/// copies mirrors and conjugates a plane having a distance to the center of CtrDist
/// This has to be done, when a scaling is performed by appending zeroes.
void CopyFFTPlane(int dim,int CtrDist)
{
  int x,y;
  int DimX, DimY, DimZ;
  
  int * px, *py, *pz;  // pointers are used to switch dimensions meaning

  if (dim == 0)
    {
      DimZ=GetSize(0);
      px= &CtrDist;
      DimY=GetSize(1);
      py= &y;
      DimX=GetSize(2);
      pz= &x;
    }

  if (dim == 1) 
    {
      DimZ=GetSize(1);
      py= &CtrDist;
      DimY=GetSize(2);
      pz= &y;
      DimX=GetSize(0);
      px= &x;
    }

  if (dim == 2)
    {
      DimZ=GetSize(2);
      pz= &CtrDist;
      DimY=GetSize(1);
      py= &y;
      DimX=GetSize(0);
      px= &x;
    }

  int CtrX=GetSize(0)/2;
  int CtrY=GetSize(1)/2;
  int CtrZ=GetSize(2)/2;

  int X,Y,Z;
  for (y=-DimY/2;y<DimY-DimY/2;y++)
    {
      for (x=-DimX/2;x<DimX-DimX/2;x++)
	{
	  X=CtrX-*px;
	  Y=CtrY-*py;
	  Z=CtrZ-*pz;
	  if (X >= GetSize(0)) X=0;
	  if (Y >= GetSize(1)) Y=0;
	  if (Z >= GetSize(2)) Z=0;
	  ArrayType::SetValue(CtrX+*px,CtrY+*py,CtrZ+*pz,Conj(ArrayType::Value(X,Y,Z)));
	}
    }
}

/// the same as 2*mkreal()
void AddConj(void)
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;
  FArrayBaseType val;

  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
	{
	  val = ArrayType::Value(x,y,z);
	  ArrayType::SetValue(x,y,z,val+Conj(val));
	}
}

void  PhaseFilter(CalcFloatType wiener)  // works also on a scrambled fft-image
{
  // Why is this not working ??? :
  // transform(TheArray.begin(),TheArray.end(),TheArray.begin(),
  //     bind2nd(Phasor<ArrayBaseType>(),wiener)); // see above
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;
  FArrayBaseType val;

  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
	{
	  val = ArrayType::Value(x,y,z);
	  ArrayType::SetValue(x,y,z,Phasor(val,wiener));
	}
}

void  GaussFilter(CalcFloatType GaussWidth, CalcFloatType GaussWidthZ)  // works only on a scrambled fft-image
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;
  cout << "Applied Gauss weight, width "<< GaussWidth <<" \n";
  ArrayType::ApplyDistWeight(GaussWidth, GaussWidthZ, &ArrayType::GaussWeightor);
}

/// returns the energy content of frequencies lying in a ring (low to high)
CalcFloatType  RingFreqEnergy(CalcFloatType  low,CalcFloatType high)  // works on a scrambled fft-scrambled image !
{
  int x,y,z;
  CalcFloatType x2,y2,z2;
  CalcFloatType r2;
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2);
  CalcFloatType SumSqr=0;
  low *= low; // square it
  high *= high; // square it

  for(z=0;z<DimZ;z++)
    for(y=0;y<DimY;y++)
      for(x=0;x<DimX;x++)
	{
	  x2=x/CalcFloatType(DimX),y2=y/CalcFloatType(DimY),z2=z/CalcFloatType(DimZ);
	  if (x > 1.0/2) x2 -= 1.0;
	  if (y > 1.0/2) y2 -= 1.0;
	  if (z > 1.0/2) z2 -= 1.0;
	  r2 = x2*x2+y2*y2+z2*z2;
	  if ((r2 > low) && (r2 < high))  // compare to square
	    SumSqr += abs(ArrayType::Value(x,y,z))*abs(ArrayType::Value(x,y,z));
	}
  return SumSqr;
}


CalcFloatType  HighFreqEnergy(CalcFloatType  reldistance)  // works on a scrambled fft-scrambled image !
{
  return RingFreqEnergy(reldistance,1.0);
}

int myfftwFwd(unsigned char * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwFwd(char * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwFwd(short int * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwFwd(int * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwFwd(double * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwBwd(unsigned char * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwBwd(char * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwBwd(short int * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwBwd(int * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

int myfftwBwd(double * Arr, bool DoNorm=true)
{ cerr << "FATAL ERROR : Tried to call FFTW with non complex array at "<< Arr<< " ! \n"; exit(-1); }

void myfftwFwd(float * Arr, bool DoNorm=true)  // compressed arrays !! real -> complex
{
  // cout << "Doing float-fft Fwd \n";
  if (! RFFTPlanFwd)
    {

#ifdef FFT_MEASURE
      RFFTPlanFwd=rfftwnd_create_plan(usedDimensions,Dims, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_MEASURE);
#else
      RFFTPlanFwd=rfftwnd_create_plan(usedDimensions,Dims,FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif
	cerr << "Plan 1\n";
#ifdef THREADS
      if(fftw_threads_init() != 0)
	cerr << "Fatal error ! Problems initializing multithreading !\n",exit(-1);
#endif
    }

#ifdef THREADS
  rfftwnd_threads_one_real_to_complex(NTHREADS,RFFTPlanFwd, (fftw_real *) Arr , NULL);  // If arg 2 gives an error : FFTW was compiled without --enable-float
#else
  rfftwnd_one_real_to_complex(RFFTPlanFwd, (fftw_real *) Arr , NULL);  // If arg 2 gives an error : FFTW was compiled without --enable-float
#endif
  if (DoNorm)
    Mul((FArrayBaseType) (1.0/sqrt(double(Dims[usedDimensions-1]*GetSize(1)*GetSize(2)))));
}

void myfftwBwd(float * Arr,bool DoNorm=true)  // compressed arrays !! complex -> real
{ 
  // cout << "Doing float-fft Bwd \n";
  if (! RFFTPlanBwd)
    {
#ifdef FFT_MEASURE
      RFFTPlanBwd=rfftwnd_create_plan(usedDimensions ,Dims, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_MEASURE);
#else
      RFFTPlanBwd=rfftwnd_create_plan(usedDimensions ,Dims, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
#endif

#ifdef THREADS
      if(fftw_threads_init() != 0)
	cerr << "Fatal error ! Problems initializing multithreading !\n",exit(-1);
#endif
    }

#ifdef THREADS
  rfftwnd_threads_one_complex_to_real(NTHREADS,RFFTPlanBwd,(FFTW_COMPLEX *) Arr , NULL);
#else
  rfftwnd_one_complex_to_real(RFFTPlanBwd,(FFTW_COMPLEX *) Arr , NULL);
#endif
  if (DoNorm)
    Mul((FArrayBaseType) (1.0/sqrt(double(Dims[usedDimensions-1]*GetSize(1)*GetSize(2)))));
}

void Magnitude(void)  // calculates the absolute maginitude
{
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2);
  for (int z=0;z<DimZ;z++)
    for (int y=0;y<DimY;y++)
      for (int x=0;x<DimX;x++)
	{
	  ArrayType::SetValue(x,y,z,abs(ArrayType::Value(x,y,z)));  // the real part will be the absolute value
	}
}

void ConvMul(TFFTArray * ArgArr, bool DoConjugate=false)  // This does the inplace multiplication of two array beeing fft of a real array !
{
  /* int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2);
  for (int z=0;z<DimZ;z++)
    for (int y=0;y<DimY;y++)
      for (int x=0;x<DimX;x++)
	{
	  ArrayType::SetValue(x,y,z,ArrayType::Value(x,y,z)*conj(ArgArr->Value(x,y,z)));  // the real part will be the absolute value
	}
  */

  /* For some reason the stuff below gives trouble */
  // typedef complex<FArrayBaseType> * cfp;
  typedef typename CastToComplex<FArrayBaseType>::Type * cfp;
  cfp beg1 = cfp (ArrayType::TheArray.StartPos());
  cfp beg2 = cfp (ArgArr->TheArray.StartPos());
  int allsize=GetSize(2)*GetSize(1)*GetSize(0);
  if (! IsComplex)
    allsize=GetSize(2)*GetSize(1)*GetSize(0)/2;  // in complex
  if (! DoConjugate)
    transform(beg1,beg1+allsize,beg2,beg1, multiplies<typename CastToComplex<FArrayBaseType>::Type> ()); 
  else
    transform(beg1,beg1+allsize,beg2,beg1, mulconj<typename CastToComplex<FArrayBaseType>::Type> ());   // A* conj(B)

}


void myfftwBwd(void * Arr, bool DoNorm=true)   // strange .... Bwd is the other way round here ...
{
  // cout << "Doing wfft Fwd \n";
  if (! FFTPlanFwd)
      {

#ifdef FFT_MEASURE
	FFTPlanFwd=fftwnd_create_plan(usedDimensions ,Dims, FFTW_FORWARD, FFTW_IN_PLACE | FFTW_USE_WISDOM | FFTW_MEASURE);
#else
	FFTPlanFwd=fftwnd_create_plan(usedDimensions ,Dims, FFTW_FORWARD, FFTW_IN_PLACE | FFTW_USE_WISDOM);
#endif
#ifdef THREADS
      if(fftw_threads_init() != 0)
	cerr << "Fatal error ! Problems initializing multithreading !\n",exit(-1);
#endif
      }

#ifdef THREADS
  fftwnd_threads(NTHREADS,FFTPlanFwd,1,(FFTW_COMPLEX *) Arr, 1, 0, 0, 0, 0);
#else
  fftwnd(FFTPlanFwd,1,(FFTW_COMPLEX *) Arr, 1, 0, 0, 0, 0);
#endif
  if (DoNorm)
    Mul((FArrayBaseType) (1.0/sqrt(double(Dims[usedDimensions-1])*GetSize(1)*GetSize(2))));
}


void myfftwFwd(void * Arr, bool DoNorm=true)   // strange .... Bwd is the other way round here ...
{
  // cout << "Doing wfft Bwd \n";
  if (! FFTPlanBwd)
      {

#ifdef FFT_MEASURE
	FFTPlanBwd=fftwnd_create_plan(usedDimensions ,Dims, FFTW_BACKWARD, FFTW_IN_PLACE | FFTW_USE_WISDOM | FFTW_MEASURE);
#else
	FFTPlanBwd=fftwnd_create_plan(usedDimensions ,Dims, FFTW_BACKWARD, FFTW_IN_PLACE | FFTW_USE_WISDOM);
#endif
	cerr << "Plan 2\n";
#ifdef THREADS
      if(fftw_threads_init() != 0)
	cerr << "Fatal error ! Problems initializing multithreading !\n",exit(-1);
#endif
      }

#ifdef THREADS
  fftwnd_threads(NTHREADS,FFTPlanBwd,1,(FFTW_COMPLEX *) Arr, 1, 0, 0, 0, 0);
#else
  fftwnd(FFTPlanBwd,1,(FFTW_COMPLEX *) Arr, 1, 0, 0, 0, 0);
#endif

  if (DoNorm)
    Mul((FArrayBaseType) (1.0/sqrt((double) Dims[usedDimensions-1]*GetSize(1)*GetSize(2))));
}

/// does a forward fft to the array
void FFTfwd(bool DoNorm=true)
{
  myfftwFwd(ArrayType::TheArray.StartPos(),DoNorm);  // 2 : Forward, mixed array
}

/// does a backward fft to the array
void FFTbwd(bool DoNorm=true)
{
  myfftwBwd(ArrayType::TheArray.StartPos(),DoNorm);  // 2 : Forward, mixed array
}


/// this function does a real-space shift in the fourier domain (cyclic). scramble has to be used before !
void FFTshift(CalcFloatType sx,CalcFloatType sy=0, CalcFloatType sz=0)
{
  FArrayBaseType val;
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;
  int MidX=DimX/2,MidY=DimY/2,MidZ=DimZ/2;
  const double M__PI = 2.0*asin(1.0);

  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	{
	  val = (FArrayBaseType) (ArrayType::Value(x,y,z)*
                exp(MulI((FArrayBaseType) (-M__PI*
                (((x-MidX)/(DimX/2.0))*sx+
                 ((y-MidY)/(DimY/2.0))*sy+
                 ((z-MidZ)/(DimZ/2.0))*sz))))); // imaginary exponential
	ArrayType::SetValue(x,y,z,(FArrayBaseType) val);
	}
}

/// this function does not need scramble before !
void FFTshift2(CalcFloatType sx,CalcFloatType sy=0, CalcFloatType sz=0)
{
  FArrayBaseType val;
  int DimX=GetSize(0),DimY=GetSize(1),DimZ=GetSize(2),x,y,z;

  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	{
	  val = (FArrayBaseType) (ArrayType::Value(x,y,z)*exp(MulI((FArrayBaseType) (-M_PI*((x/(DimX/2.0))*sx+(y/(DimY/2.0))*sy+(z/(DimZ/2.0))*sz))))); // imaginary exponential
	ArrayType::SetValue(x,y,z,(FArrayBaseType) val);
	}
}

/// calculates shift between allready FFTed data and applies it to this 
CalcFloatType calcShift(TFFTArray * ArgArr,  // This should allread be FFT transformed !
	       TFFTArray * TmpArr,  // Tmp-Array of right sizes
	       int steps,int range, double phases,
	       double & tx, double & ty, double & tz , 
	       double MaxRelDist=0.0, double MaxRelDistZ=0.0, bool forceplanes=false, // calculates and applies shift
	       double GaussF=0.0, double GaussFZ=0.0)  // If > 0.0 the Images will be low-pass filtered during the determination of the cross-correlation
{
 CalcFloatType max,mx,my,mz,j,norm0=1.0,norm1=1.0;
 tx=0;ty=0;tz=0;

 
 if (phases != 0)
   {
     max = abs(ArgArr->Maximum());
     cout << "phases set to " << phases << " times " << max << "\n";
     phases *= max;
   }

 norm0=abs(ArrayType::Value(0,0,0));
 norm1=abs(ArgArr->Value(0,0,0));

 for (j=0;j< steps;j++)
   {
     TmpArr->Copy(this);
     
     TmpArr->Mul(ArgArr); // convolution

     if (phases != 0)
       TmpArr->PhaseFilter(phases); // enhance HF to count mostly phases.  Is equivalent to filtering individual images

     if (GaussF > 0.0)
       TmpArr->GaussFilter(GaussF,GaussFZ); // enhance HF to count mostly phases.  Is equivalent to filtering individual images

     TmpArr->FFTbwd();

     if (MaxRelDist > 0.0 || MaxRelDistZ > 0.0)
       {
         cout << "Applying MaxRelDist:  " <<  MaxRelDist << "\n";
         TmpArr->ApplyDistWeight(MaxRelDist,MaxRelDistZ);
       }
     // TmpArr->DSave(true,"/tmp/test.kdf");
     
     max = TmpArr->WrappedMaxIntensityPosition(&mx,&my,&mz,range,forceplanes);
     cout << "step nr " << j+1 << " :  maximum at  " << 
       mx << "  " << my << "  " << mz << " is " << max << "\n" << flush;
     tx+=mx;ty+=my;tz+=mz;
	     
     scramble();
     FFTshift(mx,my,mz);  // application of the shift
     scramble(1);
   }

 cout << "  total shift applied : " << tx << "  " << ty << "  " << tz << "\n";
 return max; // / norm1; //  / norm0;
}



int DLoad(int kflag,const char * filename,const char * ReqType,int * SizeX,int * SizeY=0, int * SizeZ=0,int ElementNr=-1)
{
  int res;
  MyArrayType * tmp = new MyArrayType();  // load into a different array, which will be destroyed

  cout << "DLoad called for FFTArray\n";
  res=tmp->DLoad(kflag,filename,ReqType,SizeX,SizeY,SizeZ,ElementNr);
  Resize(* SizeX,* SizeY,* SizeZ);
  LCopy(tmp,this);
  delete tmp;

  return res;
}


/** writes array to a file in raw-format. Stream "to" has to be set to right position before. */
void  Write(ofstream * to)
{
  if (! to)
    cerr << "Fatal Error : Output stream not existent for writing !\n",exit(-1);

  
  for(int z=0;z<GetSize(2);z++)
    for(int y=0;y<GetSize(1);y++)
      {
	to->write((char *) ArrayType::Pointer(0,y,z),sizeof(FArrayBaseType)*(GetSize(0)-ExtraX));
	if (to->bad()) 
	  {
	    cerr << "Error writing file, no space on disk \n";
	    break;
	  }
      }
}


void DHeader(int kflag,ofstream & to,int Elements=1) // write header depending on the k-flag
{
  ArrayType::TheArray.DimVec[0] -= ExtraX;            // Uuggh awful trick !
  MyArrayType::DHeader(kflag,to,Elements);   
  ArrayType::TheArray.DimVec[0] += ExtraX;
}

void DSave(int kflag,const char * filename) // Save data in format dependent on kflag
{
  MyArrayType tmp;  // copy into a different array, which will be destroyed

  cout << "DSave called for FFTArray\n";
  tmp.Resize(GetSize(0)-ExtraX,GetSize(1),GetSize(2));
  LCopy(this,& tmp);
  tmp.DSave(kflag,filename);
}

};   // end of class definition

/// template specialization for copiing between FFT-arrays  and ordinary arrays !
template <class BaseT1, class ArrayType2>
void Copy(TFFTArray<BaseT1> * arr1,ArrayType2 * arr2) // copy part of array centered
{
  // cout << "specialized copy FFTArray-> normal called\n";
  Copy<TArray3d<BaseT1> , ArrayType2>(arr1,arr2,- arr1->ExtraX,0);  // this copy works than with logical sizes
}

/// template specialization for copiing between ordinary arrays and FFT-arrays !
template <class ArrayType1, class BaseT2>
void Copy(ArrayType1 * arr1,TFFTArray<BaseT2> * arr2) // copy part of array centered
{
  // cout << "specialized copy normal->FFTArray called\n";
  Copy<ArrayType1,TArray3d<BaseT2> >(arr1,arr2,0,- arr2->ExtraX);  // this copy works than with logical sizes
}

/// template specialization for copiing between ordinary arrays and FFT-arrays !
template <class BaseT1, class BaseT2>
void CopyRtoC(TFFTArray<BaseT1> * arr1,TFFTArray<BaseT2> * arr2) // copy part of array centered
{
  // cout << "specialized copy normal->FFTArray called\n";
  Copy<TArray3d<BaseT1>,TArray3d<BaseT2> >(arr1,arr2,- arr1->ExtraX,0);  // this copy works than with logical sizes
}

/// template specialization for copiing between FFT-arrays !
template <class BaseT1, class BaseT2>
void Copy(TFFTArray<BaseT1> * arr1,TFFTArray<BaseT2> * arr2) // copy part of array centered
{
  // cout << "specialized copy FFFTArray->FFTArray called\n";
  Copy<TArray3d<BaseT1>,TArray3d<BaseT2> >(arr1,arr2,- (arr1->ExtraX),- (arr2->ExtraX));  // this copy works than with logical sizes
}

/// template specialization for copiing between FFT-arrays !
template <class BaseT1, class BaseT2>
void CopyCtoC(TFFTArray<BaseT1> * arr1,TFFTArray<BaseT2> * arr2) // copy part of array centered
{
  // cout << "specialized copy FFFTArray->FFTArray called\n";
  Copy<TArray3d<BaseT1>,TArray3d<BaseT2> >(arr1,arr2,- (arr1->ExtraX),- (arr2->ExtraX));  // this copy works than with logical sizes
}

template <class BaseT1, class BaseT2>
void CopyCtoR(TFFTArray<BaseT1> * arr1,TFFTArray<BaseT2> * arr2) // copy part of array centered
{
  // cout << "specialized copy CtoR normal->FFTArray called\n";
  CopyCtoR<TArray3d<BaseT1>,TArray3d<BaseT2> >(arr1,arr2,0,- (arr2->ExtraX));  // this copy works than with logical sizes
}

#endif

