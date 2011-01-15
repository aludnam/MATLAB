#ifndef Projection_h		// -*- C++ -*- 
#define Projection_h

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


#include "image.h"
#include "debug.h"
#include "comperator.h"
#include "projector.h"
#include "rawarray.h"

/*
/// little helper funktion for eventual conversion to float/double
template <class ResType>
inline ResType CastIfNeeded (ResType ResultDummy, ResType Val) {return Val;}

inline float   CastIfNeeded (float ResultDummy, float_complex Val) {return CastToReal(Val);} 
inline double  CastIfNeeded (double ResultDummy, double_complex Val) {return CastToReal(Val);} */

/** Has ability to project back and forth. It is the class for the projected (measured) data.
   This class represents a conplete projection of an image. It consists of many projectors. 
   Its instances (objects) can either be measured data or projected estimations of the image to reconstruct. */
template <class TArray,class TImage,class TProjectorSelf>  // TImage is the type of the image I am talking to
class Projection : public Image<TArray,TProjectorSelf>     // TArray is my own ArrayType and TProjectorSelf a projector that works on my Array
// Due to an internal compiler error it is impossible to leave  the Image type undefined using function templates :(
{
protected:
  typedef typename Image<TArray,TProjectorSelf>::ArrayBaseType ArrayBaseType;  // Why is that not inherited ??
  /// This projector stores the actual position in space and time temporary and is moved according to its functionality
  typedef typename TImage::ProjectorType TProjector;  // This is the Projectortype of the image to talk to
  // typedef typename TImage::ProjectorType TProjector;  // This is the Projectortype of the image to talk to

public:
  TProjector  CalcPrj;  // WHy can this not be protected ???

  // ArrayBasetype is defined in image

  /** Is needed to find a start, because CalcProjector might change during the iteration. 
  Everytime a new Projections starts, CalcPrj is resetted to this one. */
  TProjector  FirstPrj;  

  ofstream * FwdFile;

/// Initialization by copy-constructors. First Projector is dublicated, and then  used as calculation Projector
Projection(TProjector * AProjector) : Image<TArray,TProjectorSelf> () , CalcPrj(* AProjector), FirstPrj(* AProjector), FwdFile(0)
{
  Image<TArray,TProjectorSelf>::TakeNewComperator(& DivComperator);         // classes own comperator is set to Standard ML/EM Algorithm
}

  virtual ~Projection() {};  // To satisfy the compiler

/// Computes the complete projection of an image
virtual  void ProjectImage(TImage * ImageToProj)  // can be virtual if NormalProjectImage is ever uses directly
{
  NormalProjectImage(ImageToProj);
}

/** Computes the complete projection of an image onto this projection
   The algorithm traverses the projection array (Values) and projects each clipped image point 
   with the appropiate weight thought the projector.
   Optained Values are added to the array. This means IT NORMALLY HAS TO BE CLEARED before. */
void NormalProjectImage(TImage * ImageToProj)
{
  // Use CalcPrj. StartPosition of Projection is the one of the first projector
  ArrayBaseType nVal=0,dummy;
  int movedim;

  CalcPrj.copy(&FirstPrj);   // initialize for a first projection
  /* This function runs through all projectors of the projection and backprojects them into the image */
  ImageToProj->PrepareArrayForClipping(&FirstPrj);

  CalcPrj.InitMove();  // Array has to be prepared before !

  dummy= Image<TArray,TProjectorSelf>::Array.FirstValue();   // Sets the Position of the projector to the first arrayindex, sets RelPos to 0 and linear initializes Array indices.
  
  do {
    // now AddValue instead of "Array.ChangeValue"
    // Calculate the intensity given by the intersection of this projector with the image. Set Value in the array
    Image<TArray,TProjectorSelf>::Array.AddValue(ImageToProj->ProjProjectorFwd(&CalcPrj));
    // cout << "SetVal at "; CalcPrj.show(); cout << "New Val : " << nVal+PVal << "\n";
    movedim = Image<TArray,TProjectorSelf>::Array.NextValue(& nVal);        // Sets (RelPos= Vector to last coordinate) and (nVal=Value in the array)
    if (CalcPrj.Move(movedim) != 0)
       ImageToProj->PrepareArrayForClipping(& CalcPrj);               // The clipping has to be reinitialisized each time direction/cliprange changes
  } while (movedim != -1);
}

  /// Computes a complete projection into an image
virtual void ProjectYourselfIntoImage(TImage * TheImg) // can be virtual, if NormalProjectIntoImage is ever used directly
{
  NormalProjectYourselfIntoImage(TheImg);
}

/** Computes a complete backprojection into an image, adding projected values to it
  The direction of the complete projection shall indicate the direction in which individual projectors lie.
  The direction of CalcPrj is the direction used for the first projector.
  For further  projectors this direction may be changed in this function. */
void NormalProjectYourselfIntoImage(TImage * TheImg)
{
  // Use CalcPrj. StartPosition of Projection is the one of the first projector
  ArrayBaseType nVal=0;
  int movedim;

  CalcPrj.copy(&FirstPrj);   // initialize for a first projection

  TheImg->PrepareArrayForClipping(&FirstPrj);        // Initailizes ArrayClipping in image. (--> InitArrayClippingFirst)
  /* This function runs through all projectors of the projection and backprojects them into the image */

  CalcPrj.InitMove();  // Array has to be prepared before !

  nVal= Image<TArray,TProjectorSelf>::Array.FirstValue();
  do {
    TheImg->ProjProjectorBkwd(&CalcPrj,nVal);  // Project this beam / projector into the image. If necessary convert to The ArrayBaseType used here

    movedim = Image<TArray,TProjectorSelf>::Array.NextValue(& nVal);               // Sets (RelPos= Vector between successive projection-coordinates in image systhem) and (nVal=Value in the projector array)
    if (CalcPrj.Move(movedim) != 0)
       TheImg->PrepareArrayForClipping(& CalcPrj);   // The clipping has to be reinitialisized each time direction/cliprange changes
  } while (movedim != -1);
}


/** Computes itself as an ML correction projection, overwriting old contents. 
 Still has to be backprojected, 1 has to be subtracted in the image and added to orig. */
void CmpMlCorrection(TImage * ReconstructionImage,Projection * MeasuredPrj)
{
  // this is very important, because projecting into an image only adds the values.
  Image<TArray,TProjectorSelf>::Array.Set(0); 
  ProjectImage(ReconstructionImage);    

  if (FwdFile) Image<TArray,TProjectorSelf>::Array.Write(FwdFile); // save forward projection

#ifdef DEBUG_SAVE  
  Image<TArray,TProjectorSelf>::Array.KSave("/tmp/fwd.raw");                  // save Fwd projected image
  MeasuredPrj->Array.KSave("/tmp/msr.raw");                  // save Fwd projected image
#endif

  Image<TArray,TProjectorSelf>::Array.Apply(Image<TArray,TProjectorSelf>::CompareFunc,& MeasuredPrj->Array);   // Use user definable compare function

#ifdef DEBUG_SAVE
  Image<TArray,TProjectorSelf>::Array.KSave("/tmp/fwddiv.raw");               // save Fwd projected image
#endif 
}

/// Computes correction img and adds it into existing correction img
void CorrectToCorrectionImg(TImage * ReconstructionImage,Projection* MeasuredPrj,TImage * CorrectionImg)
{
#ifdef DEBUG_PRINT
  cout << "Projecting forward and calculating ML-correction\n"; 
#endif
  CmpMlCorrection(ReconstructionImage,MeasuredPrj);  // Fwd-Projection and Quotient

#ifdef DEBUG_SAVE
  Image<TArray,TProjectorSelf>::Array.KSave("/tmp/fwddiv2.raw");  
#endif
#ifdef DEBUG_PRINT
  cout << "Projecting back now \n"; 
#endif

  ProjectYourselfIntoImage(CorrectionImg);
#ifdef DEBUG_SAVE
  CorrectionImg->Array.KSave("/tmp/correct3.raw");
#endif
}

/// Uses itself for cumputing FwdPrj and Quotient. Does One Step ! CorrectionImg is changed !
void IterateImg(TImage * ReconstructionImage,Projection* MeasuredPrj,TImage * CorrectionImg, CalcFloatType OverRelax)
{
  CorrectionImg->Array.Clear();
  CorrectToCorrectionImg(ReconstructionImage,MeasuredPrj,CorrectionImg);   // Fills in the correction image
  ReconstructionImage->CorrectImg(CorrectionImg,OverRelax);                // Corrects the reconstruction image
}

/// just print that I am a projection and own a FirstProjector
void show(void)
{
  cout << "Projection \n";
  cout << "  owning Projector : ";
  FirstPrj.show();
}

};

#endif































