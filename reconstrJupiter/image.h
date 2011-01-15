#ifndef Image_h		// -*- C++ -*- 
#define Image_h

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


#include "point.h"
#include "projector.h"
#include "comperator.h"

/** This class contains an image. It is used to store the reconstructed image and pre-iterations of it. It is the class for the image to reconstruct.*/
template <class TArray,class TProjector>
class Image	
{

public:
  typedef TArray ArrayType;
  typedef TProjector ProjectorType;

protected:
  typedef typename TArray::ArrayBaseType ArrayBaseType;

typedef ArrayBaseType (* PCompareFkt) (ArrayBaseType,ArrayBaseType);

  /// can be applied to the Array
  PCompareFkt CompareFunc;

  /// The PadValue (normally == 0) is used, when forward-projecting positions outside the data-covered region.
  ArrayBaseType PadValue;

public:
TArray Array;

 /// Constructor
  Image(void) : PadValue(0) // set PadValue to 0
{
  TakeNewComperator(& MulComperator);
}

  /// changes the PadValue of the class instance
void TakeNewPadVal(ArrayBaseType Val)
{
  PadValue = Val;
}

void TakeNewComperator(PCompareFkt aNewCompareFunc)
{
  CompareFunc = aNewCompareFunc;
}

void PrepareArrayForClipping(TProjector * AProjector) 
{
  Array.InitArrayClippingFirst(AProjector->GiveDirection(),AProjector->Length,AProjector->Width,AProjector->Depth);
}

/// Computes only one beam/PSF (projector) through an image
ArrayBaseType ProjProjectorFwd(TProjector * AProjector)
{
  ArrayBaseType sum=0, 
    * APos = Array.InitArrayClippingSecond(AProjector->GivePosition());     // InitArray for clipping with Projectors ClipRanges

  if (! APos)
    return PadValue;

  while (APos)
    {
      // Give the newly clipped pixel with image direction to Prj, for overlap computation;
      sum += (* APos)*AProjector->ProjectYourselfFwd(& Array);
      APos = (Array.GiveNextClipped());  // next Clipped ArrayPosition changes also Positions
    }
  // if (sum != 0) cout << "Nonzero Sum : " << sum << "\n";

  return sum;
}

/// Backprojects one beam/PSF (projector) into the image
void ProjProjectorBkwd(TProjector * AProjector,ArrayBaseType Weight)
{  
  ArrayBaseType * APos= Array.InitArrayClippingSecond(AProjector->GivePosition());      // The Projector takes care for initializing the clipping in the Array

  while (APos)
    {
      (* APos) += Weight * AProjector->ProjectYourselfBwd(& Array);  // LastValid Pixel is in Image coordinates
 
      APos = (Array.GiveNextClipped());  // returns pointer to next clipped pixel, returns 0 if none exists
    }
}

/// Correct Yourself by doing the multiplicativer change with CorrectionImg
void CorrectImg(Image<TArray,TProjector> * CorrectionImg, CalcFloatType OverRelax)
{
 // CorrectionImg->SubOne();                               // Because of the linear Projection (Values are added) this can be done before backprojecting.

  // CorrectionImg->ReplaceVal(-1.0,0.0);                  // if any value is -1.0 it was not influenced by backprojection. Therefore it will be set to 0 (assumed to be right)  
  CorrectionImg->Array.Apply(CompareFunc,& Array);         // Use user definable compare function

  // CorrectionImg->Mul(this);                             // hands itself to be multiplied into correction image
  CorrectionImg->Array.Mul(OverRelax);
  Array.Add(& CorrectionImg->Array);                       // Now really change the image
}


};


#endif
