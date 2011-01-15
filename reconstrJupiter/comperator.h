#ifndef COMPERATORS_H
#define COMPERATORS_H

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

#include "rawarray.h"
#include <math.h>

static int PRJS=1;    // Number of Projections
// This subtraction is part of the ML/EM Algorithm and can be carried out in the projection domain due to lineatity of the backprojection

void ChangeComperatorPrjsNumber(int NewPrjs)
{
  PRJS = NewPrjs;    // Making a Comperator Class is far better ...
}

template <class ArrayBaseType>
ArrayBaseType DivComperator(ArrayBaseType MeasuredValue,ArrayBaseType ProjectedReconstructionImageVal)    
// This is the standard comperator for Projection - Projection Comparisons
// First Argument is Measured Projection, second is ProjectedReconstructionImageValue
{
  if (ProjectedReconstructionImageVal != 0)
    return (MeasuredValue/ProjectedReconstructionImageVal - 1.0)/PRJS;
  else
    if (MeasuredValue == 0)
      return 0.0;
    else
      return MAXARRAYTYPENR;
}


// Argument one us the CorrectionImageValue, agument two is the Imagevalue to correct
template <class ArrayBaseType>
ArrayBaseType MulComperator(ArrayBaseType ImageValue,ArrayBaseType CorrectionValue)    // This is the standard comperator for Image - Image Comparisons
{
      return CorrectionValue*ImageValue;          // This subtraction is part of the ML/EM Algorithm and can be carried out in the projection domain due to lineatity of the backprojection
}

// Argument one is the CorrectionImageValue, agument two is the Imagevalue to correct
template <class ArrayBaseType>
ArrayBaseType SubComperator(ArrayBaseType ImageValue,ArrayBaseType CorrectionValue)    // 
{
      return -CorrectionValue;          //
}

template <class ArrayBaseType>
ArrayBaseType DivOnlyComperator(ArrayBaseType MeasuredValue,ArrayBaseType ProjectedReconstructionImageVal)    
// This is the standard comperator for Projection - Projection Comparisons
// First Argument is Measured Projection, second is ProjectedReconstructionImageValue
{
  if (ProjectedReconstructionImageVal != 0)
    return (MeasuredValue/ProjectedReconstructionImageVal)/PRJS;
  else
    if (MeasuredValue == 0)
      return 0.0;
    else
      return MAXARRAYTYPENR;
}

// Argument one us the CorrectionImageValue, agument two is the Imagevalue to correct
template <class ArrayBaseType>
ArrayBaseType WillOnlyMulComperator(ArrayBaseType ImageValue,ArrayBaseType CorrectionValue)    // 
{
      return ImageValue*CorrectionValue-ImageValue;          // precompensates for F k+1  = F k * Correct
}


/*
template <class ArrayType>
ArrayType CtPrjComperator(ArrayType MeasuredValue,ArrayType ProjectedReconstructionImageVal)  // This is a comperator for CT-Type data (logarithmic)
{
  return  (exp(- ProjectedReconstructionImageVal) - MeasuredValue)/PRJS;
}   // exp(-MeasuredValue) - exp(...)

template <class ArrayType>
ArrayType CtImgComperator(ArrayType ImageValue,ArrayType CorrectionValue)  // This is a comperator for CT-Type data (logarithmic)
{
  return 0.1*CorrectionValue;
}
*/

template <class ArrayBaseType>
ArrayBaseType CtPrjComperator(ArrayBaseType MeasuredValue,ArrayBaseType ProjectedReconstructionImageVal)  // This is a comperator for CT-Type data (logarithmic)
{
  return  (MeasuredValue - ProjectedReconstructionImageVal)/PRJS;         // - log(Measured)
}   // exp(-MeasuredValue) - exp(...)          // -log(Measured / exp(Proj))/PRJS

template <class ArrayBaseType>
ArrayBaseType CtImgComperator(ArrayBaseType ImageValue,ArrayBaseType CorrectionValue)  // This is a comperator for CT-Type data (logarithmic)
{
  return CorrectionValue;
}

#endif
