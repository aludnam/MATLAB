#ifndef VecProjection_h		// -*- C++ -*- 
#define VecProjection_h

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

#include "projection.h"

/** This projection class projects vectorial represantations of Spots. To make projections fast, a clipping in the projection is incorporated.
This module works only with projectorpos used as projector, because it uses specific features !! */
template <class TArray,class TImage,class TProjectorSelf>
class VecProjection : public Projection<TArray,TImage,TProjectorSelf>
{
private:
  TImage * ClipImg;
  CalcFloatType Background;
  TProjectorSelf SelfProjector;  // does clipping inside this projection

  // TProjector is defined in Projection.h
public:
  typedef typename Projection<TArray,TImage,TProjectorSelf>::ArrayBaseType ArrayBaseType;  // Why is that not inherited ?
  typedef typename Projection<TArray,TImage,TProjectorSelf>::TProjector TProjector;  // Why is that not inherited ?

VecProjection(TProjector * AProjector,TImage * ClipI,TProjectorSelf * SProjector) : Projection<TArray,TImage,TProjectorSelf> (AProjector), Background(0),SelfProjector( * SProjector)
{
  ClipImg=ClipI;           // May be the same as Imgages in ProjectImage or ...
}

void SetBackground(CalcFloatType Val)
{
  Background=Val;
}

void SetClip(CalcFloatType Length,CalcFloatType Width,CalcFloatType Depth)
{
  SelfProjector.SetClip(Length,Width,Depth);
}

void SetPsfFac(CalcFloatType Val)
{
  // SelfProjector.SetPsfFac(Val);
  this.CalcPrj.SetPsfFac(Val);
}

/** Computes the complete projection of a vector-image onto this projection
  The algorithm traverses clipped parts of the projection array and projects each clipped image row thought the projector. */
virtual void ProjectImage(TImage * ImageToProj)
{
  // Use CalcPrj. Position of Projection is the one of the first projector
  ArrayBaseType * APos=0;
  int SpotNr;
  Projection<TArray,TImage,TProjectorSelf>::Array.Set((double) Background);      // Clear the Array
  const Vector Zero(0.0,0.0,0.0),One(1.0,0.0,0.0);


  PrepareArrayForClipping(&SelfProjector);  // prepare clipping in this projection

  // CalcPrj is used for Image, SelfProjector is used for self.
  Projection<TArray,TImage,TProjectorSelf>::CalcPrj.copy(&Zero,& One,5.0,0.1,0.1);
  // Pos,dir,clip,interp --> only one row in image is clipped
  typedef typename TImage::ProjectorType TProjector;

  ImageToProj->PrepareArrayForClipping(& Projection<TArray,TImage,TProjectorSelf>::CalcPrj);

for (SpotNr=0;SpotNr<ImageToProj->Array.GetSize(1);SpotNr++)  // Y-Size beeing number of Spots
  if (Projection<TArray,TImage,TProjectorSelf>::CalcPrj.IsActiveSpot(SpotNr))  // only active spots are projected
    {
      Vector ClipVecCenter(ClipImg->Array.Value(1,SpotNr),
		     ClipImg->Array.Value(2,SpotNr),
		     ClipImg->Array.Value(3,SpotNr));
      Vector CalcPrjPos(0.0,SpotNr,0.0);
	
	// First read middle clippingvector from image (== value,vecx,vecy,..)
	APos = Projection<TArray,TImage,TProjectorSelf>::Array.InitArrayClippingSecond(& ClipVecCenter); 
    // Sets Position of CalcPrj due to values in image. SO FAR ONLY THREE-D !!!

    Projection<TArray,TImage,TProjectorSelf>::CalcPrj.SetCoordinate(& CalcPrjPos);  // change position of projector
   
    while (APos)
      {
	// CalcPrj.SetPrjPos(Array.LastValidPixel);         // change position of projector.  If you get an error here : USE PROJECTORPOS AS PROJECTOR !!!
	Projection<TArray,TImage,TProjectorSelf>::CalcPrj.SetPrjPos(Projection<TArray,TImage,TProjectorSelf>::Array.Clipperator.TXPos(),Projection<TArray,TImage,TProjectorSelf>::Array.Clipperator.TYPos(),Projection<TArray,TImage,TProjectorSelf>::Array.Clipperator.TZPos()); // change position of projector.  If you get an error here : USE PROJECTORPOS AS PROJECTOR !!!
       (* APos) += ImageToProj->ProjProjectorFwd(&Projection<TArray,TImage,TProjectorSelf>::CalcPrj); // LastValid Pixel is in Image coordinates
       APos = (Projection<TArray,TImage,TProjectorSelf>::Array.GiveNextClipped());                    // returns pointer to next clipped pixel, returns 0 if none exists
      }
    }
}

/// Computes a complete backprojection into an image, adding projected values to it
virtual void ProjectYourselfIntoImage(TImage * TheImg)
{
  // Use CalcPrj. Position of Projection is the one of the first projector
  ArrayBaseType * APos=0;
  int SpotNr;
  static const Vector Zero(0.0,0.0,0.0),One(1.0,0.0,0.0);
  
  PrepareArrayForClipping(&SelfProjector);           // prepare clipping in this projection

  // CalcPrj is used for Image, SelfProjector is used for self.
  Projection<TArray,TImage,TProjectorSelf>::CalcPrj.copy(&Zero,& One,5.0,0.0,0.0);      // Pos,dir,clip,interp --> only one row in image is clipped
  TheImg->PrepareArrayForClipping(&this.CalcPrj);    

for (SpotNr=0;SpotNr<TheImg->Array.GetSize(1);SpotNr++)              // Y-Size beeing number of Spots
  if (this.CalcPrj.IsActiveSpot(SpotNr))  // only active spots are backprojected
    {
      Vector ClipVecCenter(ClipImg->Array.Value(1,SpotNr),
		     ClipImg->Array.Value(2,SpotNr),
		     ClipImg->Array.Value(3,SpotNr));
      Vector CalcPrjPos(0.0,SpotNr,0.0);
	
    // First read middle clippingvector from image (== value,vecx,vecy,..)
    APos = Projection<TArray,TImage,TProjectorSelf>::Array.InitArrayClippingSecond(& ClipVecCenter); 
    // Sets Position of CalcPrj due to values in image. SO FAR ONLY THREE-D !!!
 
    Projection<TArray,TImage,TProjectorSelf>::CalcPrj.SetCoordinate(& CalcPrjPos);  // change position of projector
   
    while (APos)
      {
        // change position of projector.  If you get an error here : USE PROJECTORPOS AS PROJECTOR !!!
	Projection<TArray,TImage,TProjectorSelf>::CalcPrj.SetPrjPos(Projection<TArray,TImage,TProjectorSelf>::Array.Clipperator.TXPos(),Projection<TArray,TImage,TProjectorSelf>::Array.Clipperator.TYPos(),Projection<TArray,TImage,TProjectorSelf>::Array.Clipperator.TZPos());
	// CalcPrj.SetPrjPos(Array.LastValidPixel);       // change position of projector
       TheImg->ProjProjectorBkwd(& Projection<TArray,TImage,TProjectorSelf>::CalcPrj,(* APos)); // LastValid Pixel is in Image coordinates
       APos = (Projection<TArray,TImage,TProjectorSelf>::Array.GiveNextClipped());              // returns pointer to next clipped pixel, returns 0 if none exists
      }
    }
}

};

#endif































