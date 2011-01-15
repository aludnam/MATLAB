#ifndef ClipProjector_h		// -*- C++ -*- 
#define ClipProjector_h

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


// #include "point.h"
#include "projector.h"
      
/** The projector and its childs deal with arrays, that know about clipping
 */
template<class ArrayType>
class ClipProjector : public Projector
{
private:
  void UpdatePos2(void);

public:
  /// Length of clipping along direction vector
  CalcFloatType Length;
  /// Width of clipping perpendicular direction vector
  CalcFloatType Width;
  /// Width of clipping perpendicular direction vector
  CalcFloatType Depth;


  ClipProjector(const int Dim) : Projector(Dim) , Length(1.0), Width(1.0), Depth(1.0)
    { };

  /// Copy Constructor has to allocate space for Vectors, but not for Array
  ClipProjector(const ClipProjector & CopyFromPrj) : Projector(CopyFromPrj),Length(CopyFromPrj.Length), Width(CopyFromPrj.Width), Depth(CopyFromPrj.Depth)
    { }    

  /// just copy, without constructing.
  void copy(const ClipProjector * APrj) 
    {
      Projector::copy(APrj);
      Length=APrj->Length;
      Width=APrj->Width;
      Depth=APrj->Depth;
    }

  /// copy without allocating space.
  void copy(const Vector* pItsPos,
	    const Vector* pItsDir,
	    CalcFloatType L, CalcFloatType W,CalcFloatType D) 
    {
        Projector:: copy(pItsPos,pItsDir);
	Length = L;
	Width = W;
	Depth = D;
    }

  void SetClip(CalcFloatType L, CalcFloatType W=0, CalcFloatType D=0) 
    {
	Length = L;
	Width = W;
	Depth = D;
    }

  void GiveClip(CalcFloatType * L, CalcFloatType*  W, CalcFloatType*  D)
    {
      (*L) = Length;
      (*W) = Width;
      (*D) = Depth;
    }


  /** computes intersection of itself with one pixel.
  pPosition is Position in Image (image-coordinates)
  computes intersection of itself with one pixel. 
  Returns value between 0.0 and 1.0*Value
  This Function is to be written by the user
  Return value of 1.0 means full contribution, 0 means none. */

CalcFloatType ProjectYourselfFwd(ArrayType * pArr)	   
{
  return 1.0;
  pArr=0;  // to prevent warning
}

 /// computes intersection of itself with one pixel for the backward Projection.See  ProjectYourselfFwd for description
CalcFloatType ProjectYourselfBwd(ArrayType * pArr)	
{
  return 1.0;
  pArr=0;  // to prevent warning
}

};

#endif

