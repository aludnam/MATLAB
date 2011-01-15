#ifndef clipperator2d_h  		// -*- C++ -*- 
#define clipperator2d_h

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

#include "clipperator.h"
#include <limits.h>

/// fractional distances are allways measured in SUBUNITS of a pixel.
#ifndef BIGSUBUNITS
const int SUBUNITS=(256*256*256); // better would be local *4, 256*256*256*16 ist too big for confospot program
const int SUBUNITSY=256;        // for measuring Y-distance to center of cliprange
#else
const int SUBUNITS=(8*256*256*256); // better would be local *4, 256*256*256*16 ist too big for confospot program
const int SUBUNITSY=256;        // for measuring Y-distance to center of cliprange
#endif

template<class FType>
int arint(FType val)  // round value to nearest integer value
{
  assert(fabs(val) < INT_MAX);
  if (val > 0)
    return int(val+0.5);
  else
    return -int((-val)+0.5);
}


/// A modified 2D Version of ClipperatorClip
class Clipperator2DClip : public ClipperatorNeighbClip {
public:

  typedef int INT32;
protected:
  IndexType StepYRel;
  /// the coordinates of the last clipped pixel in the transformed unit-system (flipped, but NOT sheared!)
  IndexType PosX,PosY;// Readable from others !
  INT32 SPY;  // SUBPIXELY  position to avoid conversion
  CalcFloatType FPY;  // float position to avoid conversion

  CalcFloatType DirX,DirY;

  /// Distance from the middle of the pixel to the middle of the beam (cliprange) taken in unified coordinates measured in SUBUNITS
  INT32 PixXToBeam;
  INT32 PixXSToBeam;
  /// will be susbtracted for every y-step (measured in SUBUNITS)
  INT32 XDecrement;
  /// Left and right border in SUBUNITS for X
  INT32 MinXUnits,MaxXUnits;
  /// Pixels in Y-Direction for cliprange
  INT32 YUnits;
  /// Border-Pixel cooredinates (X-Start) in unit-system
  /// Borders of the array in the  unitsystem
  IndexType BorderXL,BorderYL;

  /// variables for coordinate transform
  short TX,TY,XCXY;              // FlipX, FlipY, Exchange X and Y

public:



Clipperator2DClip() 
{
  XDecrement = 0;
  TX = TY = 1;
  XCXY = 0;

}
  /// returns Float distance to clipping center in the TRANSFORMED coordinate system !  Overwritten Method !
CalcFloatType XDistToCenter(void) const
  {
    return (CtrX-(CalcFloatType) PosX);
  }

  /// returns Float distance to clipping center in the TRANSFORMED coordinate system !  Overwritten Method !
CalcFloatType YDistToCenter(void)   const
  {
    return (CtrY-(CalcFloatType) PosY);
  }

/// returns Position in transformed coordinate system  Overwritten Method !
IndexType XPos(void) const
  {
    return PosX;
  }
/// returns Position in transformed coordinate system  Overwritten Method !
IndexType YPos(void) const
  {
    return PosY;
  }

INT32 SYPos(void) const  // Y-distance of actual pixel to focus in SUBUNITSY
  {
    return SPY;
  }

CalcFloatType FYPos(void) const  // Y-distance of actual pixel to focus in float
  {
    return FPY;
  }

  /// returns tan of beamangle*SUBUNITS
INT32 BeamXTan(void) const
  {
    return XDecrement;
  }

  /// returns PixXtoBeam in SUBUNITS
INT32 XDistToLine(void) const
  {
    return PixXToBeam;
  }


  /// returns Position in original coordinate system
IndexType TXPos(void) const
  {
    if (XCXY)
      return PosY*TX; 
    else
      return PosX*TX;
  }

  /// returns Position in original coordinate system
IndexType TYPos(void) const
  {
    if (XCXY)
      return PosX*TY; 
    else
      return PosY*TY;
  }

/// calculates the nessesary transformation from the unitsystem to the real array
void CalcXForm(CalcFloatType DirectionX,CalcFloatType DirectionY) 
{
  DirX=DirectionX;
  DirY=DirectionY;

  TX=TY=1;
  XCXY=0;  // exchange X<->Y ...
  // cerr << " CalcXForm TheDir "; TheDirection->show();  cerr << "TheClip";  TheClipRange->show();  cerr << "\n";
  StepX=1; 
  StepY=SizeX;  // Other definition of Y-direction !!

  if (DirX < 0) 
    {
      TX=-1;
      DirX=-DirX;
      StepX = -StepX;
    }

  if (DirY < 0) 
    {
      TY=-1;
      DirY=-DirY;
      StepY = -StepY;
    }

  if (DirX > DirY) 
    {
      XCXY=1;
      SWAP(DirX,DirY);
      SWAP(StepX,StepY);
    }

  //  cerr << "CalcXForm TX " << TX << " TY "<< TY << " XC " << XCXY <<"\n";  TheDirection->show();  TheDirection2->show(); TheClipRange->show();
}

void XFormFwd(CalcFloatType& X,CalcFloatType& Y) const  // from orginal coordinates into Clipping coordinates
{ 
  X*=TX; Y*=TY; 
  if (XCXY != 0) 
      SWAP(Y,X);
}
void XFormFwd(IndexType& X,IndexType& Y)  const // from orginal coordinates into Clipping coordinates
{ 
  X*=TX; Y*=TY;
  if (XCXY != 0) 
      SWAP(Y,X);
}

/// transforms from the transformed coordinate system into normal array system
void XForm(IndexType& X,IndexType& Y) const
{ 
  if (XCXY != 0) SWAP(X,Y);

  X *= TX;  Y *= TY;

  //  cerr << x2 << "  " << y2 << " Tx " << TX << "  " << TY << "\n";

  /* if ((X < 0 ||  X >= SizeX)  // check if pixel is valid
      || (Y < 0 ||  Y >= SizeY) )
      return 0; */
  
      // cerr << "Pixel good , at X,Y " << X << " "<< Y << "\n";

  // Pixel X,Y are now in original array
}

/// Length : Half Length of whole clipping region;  Width : Half width of clipping region
void InitClippingFirst3D(CalcFloatType DirectionX,CalcFloatType DirectionY,CalcFloatType DirectionZ,CalcFloatType Length,CalcFloatType Width,CalcFloatType Depth)
{
  CalcXForm(DirectionX,DirectionY);    // From now on everything can be considered to be in the first octant. 
  StepY = - StepY;  // points into the other direction

  // then compute width of clipping function
  // This model is simplebut fast. Clipping varies only in Y-direction

  // Next calculate decrement for DirectionDistance
  XDecrement=arint(SUBUNITS*DirX/DirY);  // Has to be subtracted in every Y-Step for X-border

  MaxXUnits = (INT32) fabs((Width/DirY)* SUBUNITS);//  + (SUBUNITS + XDecrement)/2;

  MinXUnits = (INT32) - MaxXUnits;

  YUnits = (INT32) fabs(Length * DirY);
}

int InitClippingSecond3D(CalcFloatType CenterX,CalcFloatType CenterY,CalcFloatType CenterZ)
{
  /// The relative Starting ofset
  INT32 BCenterY; // Transformed Center-Point in Y

  CalcFloatType ofx;
  int XB1=(SizeX-1);
  int YB1=(SizeY-1); // imageborders
  int XB2=0;
  int YB2=0;
  // Transform borders and OffSet into new coordinate system and 
  // calculate offset (lower left corner) in transformed coordinates of the central pixel

  CtrX=CenterX;  CtrY=CenterY;

  XFormFwd(XB1,YB1);
  // XFormFwd(XB2,YB2); // Not necessary, because they are 0 !
  XFormFwd(CtrX,CtrY);
  // CtrX, CtrY contain now the transformed center coordinates

  BCenterY = (INT32) arint(CtrY);  // Save Variable
  BorderYL = BCenterY - YUnits;   // Border of lower y-cliprange

  // ofy = CtrY + YUnits;          // Float

  PosY = BCenterY + YUnits;  // is rounded Y offset of lower left edge of clipping range in array (this y is really touched by the cliprange !)

  // Now clip at the upper Y-Border if Necessary !
  if (YB1 < YB2) 
    SWAP(YB1,YB2);   // So now YB1 is guaranteed to be >= YB2 !

  if (YB1 < PosY)   // clipping at upper y-Border
      PosY = YB1;   // Start at upper border. PixToBeams, PosX  will be calculated below !

  if (YB2 > BorderYL)
    BorderYL = YB2;   // Set lower Border according to Image-Border instead of cliprange (above)
  
  if (PosY < BorderYL)
    return 0;  // No pixel in range

  // Positions below define startingposition (X+width,Y+height), == top left corner of cliprange
  ofx = CtrX + (PosY-CtrY)*(CalcFloatType) XDecrement / (CalcFloatType) SUBUNITS;  // is middle of beam at first y-position

  PosX = (IndexType) arint (0.5 + ofx + MinXUnits / (CalcFloatType)SUBUNITS);  // calc first pixel, that is really touched

  // Now correct Pos and set startingpositions if there y-cliprange is too big.
  // Shift measures the shift relative to starting points x-position in the negative direction !

  PixXToBeam = arint ((PosX-ofx) * (CalcFloatType)SUBUNITS);   // will be negative


  if (XB1 < XB2) 
    SWAP(XB1,XB2);   // So now XB1 is guaranteed to be >= XB2 !

  BorderXL = XB2;    // Set lower Border
  if (PosX < BorderXL)    // clipping at lower x-Border
    {
      if (PixXToBeam+(BorderXL-PosX)*CalcFloatType(SUBUNITS) > MaxXUnits)
	return 0;
      PixXToBeam += (BorderXL-PosX)*SUBUNITS;
      PosX = BorderXL;
    }

  BorderXH = XB1;    // Set higher Border
  if (PosX > BorderXH)    // clipping at upper x-Border. A little lower in y it might enter a legal range !
    {
      INT32 tmp;  // how many pixels must I go down ?
      CalcFloatType PixXFloat=PixXToBeam;
      if (XDecrement == 0) // out of range
	return 0;
      PixXFloat -= CalcFloatType(PosX-BorderXH)*SUBUNITS;  // This MAY BE DANGEROUS if long is too little range
      tmp = 1 + int((MinXUnits-PixXFloat) /  CalcFloatType(XDecrement)); // This is hopefully the right formula !
      PosY -= tmp;
      PixXFloat += CalcFloatType(tmp)* XDecrement;
      PixXToBeam=arint(PixXFloat);
      if (PosY < BorderYL) // out of range
	return 0;
      PosX = BorderXH;
    }

  PosXS = PosX;
  PixXSToBeam=PixXToBeam;
  StepYRel=0;  

  if (PixXToBeam > MaxXUnits) // not in beam anymore
    if (GiveNextClipped() == 0)
      return 0;  // no pixel in range

  PositionX=PosX,PositionY=PosY;
  SPY= int(CalcFloatType(PosY-CtrY)*CalcFloatType(SUBUNITSY)+0.5); // Y-distance of pixel to focus  should be positive
  FPY= CalcFloatType(PosY-CtrY); // Y-distance of pixel to focus  should be positive
  XForm(PositionX,PositionY);  // transform back to normal coordinates

  return 1;
}

/// returns next clipped pixel. However : one pixel per row is returned at least
IndexType GiveNextClipped(void)     
{
  IndexType Step;
  // calculate new position and then return it;
  PixXToBeam += SUBUNITS;
  if ((PixXToBeam <= MaxXUnits)
      && (PosX < BorderXH) // else step in z or y
      )
    {
      ++PosX;
      StepYRel-=StepX;
      return StepX;
    }

  Step=StepYRel; // go back to PosXS
  
  do {
    Step+=StepY;  // Go one step downwards
    SPY-=SUBUNITSY;  // A whole pixel downwards
    FPY--;  // A whole pixel downwards
    if (--PosY < BorderYL)
      return 0;  // no more pixels

    PixXSToBeam += XDecrement; // we are closer to the center now
    if (PixXSToBeam - MinXUnits >= SUBUNITS) // one left ?
      {
	if (PosXS > BorderXL) // >= BorderXL+1
	  {
	    PixXSToBeam -= SUBUNITS;
	    Step -= StepX;
	    --PosXS;
	  }
	else
	  if (PixXSToBeam > MaxXUnits) return 0; // End of beam at left (lower x-) side
      }

  StepYRel=0;  // prepare for next Step into -Y direction
  PosX = PosXS;
  PixXToBeam=PixXSToBeam;
  } while (PixXToBeam > MaxXUnits);

  return Step;
}
};

#endif
