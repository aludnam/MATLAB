#ifndef clipperator3d_h  		// -*- C++ -*- 
#define clipperator3d_h

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

#include "clipperator2d.h"


class Clipperator3DClip : public Clipperator2DClip{
protected:

  IndexType StepZRel;

  /// the coordinates of the last clipped pixel in the transformed unit-system (flipped, but NOT sheared!)
  IndexType PosZ;  // Readable from others !

  CalcFloatType DirZ;

  /// Distance from the middle of the pixel to the middle of the beam (cliprange) taken in unified coordinates measured in SUBUNITS
  INT32 PixZToBeam;
  INT32 PixZSToBeam;
  /// will be susbtracted for every y-step (measured in SUBUNITS)
  INT32 ZDecrement;
  /// Left and right border in SUBUNITS for Z
  INT32 MinZUnits,MaxZUnits;
  /// Border-Pixel cooredinates (X-Start, Z-Start) in unit-system
  IndexType PosZS;  // PosYS is not used in CLIP and CYCLIC
  /// Borders of the array in the  unitsystem
  IndexType BorderZL;

  /// variables for coordinate transform
  short TZ,XCYZ;              // FlipX, FlipY, Exchange X and Y

public:


Clipperator3DClip() : Clipperator2DClip()
{
  ZDecrement = 0;
  TZ = 1;
  XCYZ = 0;

}

/// returns Float distance to clipping center in the TRANSFORMED coordinate system !  Overwritten Method !
CalcFloatType ZDistToCenter(void)
  {
    return (CtrZ-(CalcFloatType) PosZ);
  }

/// returns Position in transformed coordinate system
IndexType ZPos(void)
  {
    return PosZ;
  }

  /// returns tan of beamangle*SUBUNITS
INT32 BeamZTan(void)
  {
    return ZDecrement;
  }

  /// returns PixZtoBeam in SUBUNITS
INT32 ZDistToLine(void)
  {
    return PixZToBeam;
  }

  /// returns Position in original coordinate system
IndexType TZPos(void)
  {
    if (XCYZ)
      if (XCXY)
	return PosY*TZ; 
      else
	return PosY*TZ;
    else
      if (XCXY)
	return PosZ*TZ; 
      else
	return PosZ*TZ;
  }
  /// returns Position in original coordinate system
IndexType TXPos(void)
  {
    if (XCYZ)
      if (XCXY)
       return PosZ*TX; 
      else
	return PosX*TX;
    else
       if (XCXY)
	return PosY*TX; 
      else
	return PosX*TX;
  }

  /// returns Position in original coordinate system
IndexType TYPos(void)
  {
    if (XCYZ)
      if (XCXY)
	return PosX*TY; 
      else
	return PosZ*TY;
    else
      if (XCXY)
	return PosX*TY; 
      else
	return PosY*TY;
  }

/// calculates the nessesary transformation from the unitsystem to the real array
void CalcXForm(CalcFloatType DirectionX,CalcFloatType DirectionY, CalcFloatType DirectionZ)
{
  DirX=DirectionX;
  DirY=DirectionY;
  DirZ=DirectionZ;

  TX=TY=TZ=1;
  XCXY=XCYZ=0;  // exchange X<->Y ...
  // cerr << " CalcXForm TheDir "; TheDirection->show();  cerr << "TheClip";  TheClipRange->show();  cerr << "\n";
  StepX=1; 
  StepY=SizeX;  // Other definition of Y-direction !!
  StepZ=SizeX*SizeY;

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

  if (DirZ < 0) 
    {
      TZ=-1;
      DirZ=-DirZ;
      StepZ = -StepZ;
    }

  if (DirX > DirY) 
    {
      XCXY=1;
      SWAP(DirX,DirY);
      SWAP(StepX,StepY);
    }

  if (DirZ > DirY) 
    {
      XCYZ=1;
      SWAP(DirZ,DirY);
      SWAP(StepZ,StepY);
    }

  //  cerr << "CalcXForm TX " << TX << " TY "<< TY << " XC " << XCXY <<"\n";  TheDirection->show();  TheDirection2->show(); TheClipRange->show();
}

void XFormFwd(CalcFloatType& X,CalcFloatType& Y,CalcFloatType& Z)  // from orginal coordinates into Clipping coordinates
{ 
  X*=TX; Y*=TY; Z*=TZ;
  if (XCXY != 0) 
      SWAP(Y,X);
  if (XCYZ != 0) 
      SWAP(Y,Z);
}
void XFormFwd(IndexType& X,IndexType& Y,IndexType& Z)  // from orginal coordinates into Clipping coordinates
{ 
  X*=TX; Y*=TY; Z*=TZ;
  if (XCXY != 0) 
      SWAP(Y,X);
  if (XCYZ != 0) 
      SWAP(Y,Z);
}

/// transforms from the transformed coordinate system into normal array system
void XForm(IndexType& X,IndexType& Y,IndexType& Z)
{ 
  if (XCYZ != 0) SWAP(Z,Y);
  if (XCXY != 0) SWAP(X,Y);

  X *= TX;  Y *= TY;  Z *= TZ; 

  //  cerr << x2 << "  " << y2 << " Tx " << TX << "  " << TY << "\n";

  /* if ((X < 0 ||  X >= SizeX)  // check if pixel is valid
      || (Y < 0 ||  Y >= SizeY)
      || (Z < 0 ||  Z >= SizeZ ))
      return 0; */
  
      // cerr << "Pixel good , at X,Y " << X << " "<< Y << " "<< Z << "\n";

  // Pixel X,Y,Z are now in original array
}


// ************************************ Clipping routines start here ********************

/// Length : Half Length of whole clipping region;  Width : Half width of clipping region
void InitClippingFirst3D(CalcFloatType DirectionX,CalcFloatType DirectionY,CalcFloatType DirectionZ,CalcFloatType Length,CalcFloatType Width, CalcFloatType Depth)
{
  CalcXForm(DirectionX,DirectionY,DirectionZ);    // From now on everything can be considered to be in the first octant. 
  StepY = - StepY;  // points into the other direction

  // then compute width of clipping function
  // This model is simplebut fast. Clipping varies only in Y-direction

  // Next calculate decrement for DirectionDistance
  XDecrement=rint(SUBUNITS*DirX/DirY);  // Has to be subtracted in every Y-Step for X-border
  ZDecrement=rint(SUBUNITS*DirZ/DirY);  // Has to be subtracted in every Y-Step for Z-border

  MaxXUnits = (INT32) fabs((Width/DirY)* SUBUNITS) 
    + (SUBUNITS + XDecrement)/2;    // If removed -> clipping of pixelcenters, not edges

  MinXUnits = (INT32) - MaxXUnits;

  MaxZUnits = (INT32) fabs((Depth/DirY)* SUBUNITS) 
    + (SUBUNITS + ZDecrement)/2;    // If removed -> clipping of pixelcenters, not edges

  MinZUnits = (INT32) - MaxZUnits;

  YUnits = (INT32) fabs(2*Length * DirY+0.5);  // Length in both directions of clipping center
}

int InitClippingSecond3D(CalcFloatType CenterX,CalcFloatType CenterY,CalcFloatType CenterZ)
{
  /// The relative Starting ofset
  INT32 BCenterY; // Transformed Center-Point in Y

 CalcFloatType ofx,ofz;
  int XB1=(SizeX-1);
  int YB1=(SizeY-1); // imageborders
  int ZB1=(SizeZ-1);
  int XB2=0;
  int YB2=0;
  int ZB2=0;
  // Transform borders and OffSet into new coordinate system and 
  // calculate offset (lower left corner) in transformed coordinates of the central pixel

  CtrX=CenterX;  CtrY=CenterY;  CtrZ=CenterZ;

  XFormFwd(XB1,YB1,ZB1);
  // XFormFwd(XB2,YB2,ZB2); // Not necessary, because they are 0 !
  XFormFwd(CtrX,CtrY,CtrZ);
  // CtrX, CtrY and CtrZ contain now the transformed center coordinates

  BCenterY = (INT32) rint(CtrY);  // Save Variable
  BorderYL = BCenterY - YUnits;   // Border of lower y-cliprange

  // ofy = CtrY + YUnits;          // Float

  PosY = BCenterY + YUnits;  // is rounded Y offset of lower left edge of clipping range in array (this y is really touched by the cliprange !)

  // Now clip at the upper Y-Border if Necessary !
  if (YB1 < YB2) 
    SWAP(YB1,YB2);   // So now YB1 is guaranteed to be >= YB2 !

  if (YB1 < PosY)   // clipping at upper y-Border
      PosY = YB1;   // Start at upper border. PixToBeams, PosX and PosZ will be calculated below !

  if (YB2 > BorderYL)
    BorderYL = YB2;   // Set lower Border according to Image-Border instead of cliprange (above)
  
  if (PosY < BorderYL)
    return 0;  // No pixel in range

  // Positions below define startingposition (X+width,Y+height,Z+depth), == top left corner of cliprange
  ofx = CtrX + (PosY-CtrY)*XDecrement / (CalcFloatType) SUBUNITS;  // is middle of beam at first y-position
  ofz = CtrZ + (PosY-CtrY)*ZDecrement / (CalcFloatType) SUBUNITS;  // -is middle of beam at first y-position 

  PosX = (IndexType) rint (0.5 + ofx + MinXUnits / (CalcFloatType)SUBUNITS);  // calc first pixel, that is really touched
  PosZ = (IndexType) rint (0.5 + ofz + MinZUnits / (CalcFloatType)SUBUNITS);  // calc first pixel, that is really touched

  // Now correct Pos and set startingpositions if there y-cliprange is too big.
  // Shift measures the shift relative to starting points x-position in the negative direction !

  PixXToBeam = rint ((PosX-ofx) * (CalcFloatType)SUBUNITS);   // will be negative
  PixZToBeam = rint ((PosZ-ofz) * (CalcFloatType)SUBUNITS);   // will be negative


  if (XB1 < XB2) 
    SWAP(XB1,XB2);   // So now XB1 is guaranteed to be >= XB2 !

  BorderXL = XB2;    // Set lower Border
  if (PosX < BorderXL)    // clipping at lower x-Border
    {
      PixXToBeam += (BorderXL-PosX)*SUBUNITS;
      PosX = BorderXL;
      if (PixXToBeam > MaxXUnits)
	return 0;
    }

  BorderXH = XB1;    // Set higher Border
  if (PosX > BorderXH)    // clipping at upper x-Border. A little lower in y it might enter a legal range !
    {
      INT32 tmp;
      if (XDecrement == 0) // out of range
	return 0;
      PixXToBeam -= (PosX-BorderXH)*SUBUNITS;  // This MAY BE DANGEROUS if long is too little range
      tmp = 1+(MinXUnits-PixXToBeam) / XDecrement; // This is hopefully the right formula !
      PosY -= tmp;
      PixXToBeam += tmp*XDecrement;
      if (PosY < BorderYL) // out of range
	return 0;
      PosX = BorderXH;
    }


  if (ZB1 < ZB2) 
    SWAP(ZB1,ZB2);   // So now ZB1 is guaranteed to be >= ZB2 !

  BorderZL = ZB2;    // Set lower Border
  if (PosZ < BorderZL)    // clipping at lower z-Border
    {
      PixZToBeam += (BorderZL-PosZ)*SUBUNITS;
      PosZ = BorderZL;
      if (PixZToBeam > MaxZUnits)
	return 0;
    }
  
  BorderZH = ZB1;    // Set higher Border
  if (PosZ > BorderZH)    // clipping at upper z-Border
    {
      INT32 tmp;
      if (ZDecrement == 0) // out of range
	return 0;
      PixZToBeam -= (PosZ-BorderZH)*SUBUNITS;  // This MAY BE DANGEROUS if long is too little range
      tmp = 1+(MinZUnits-PixZToBeam) / ZDecrement; // This is hopefully the right formula !
      PosY -= tmp;
      PixZToBeam += tmp*ZDecrement;
      if (PosY < BorderYL) // out of range
	return 0;
      PosZ = BorderZH;
    }
  
  PositionX=PosX,PositionY=PosY,PositionZ=PosZ;

  XForm(PositionX,PositionY,PositionZ);  // transform back to normal coordinates
  PosXS = PosX;
  PixXSToBeam=PixXToBeam;
  PosZS = PosZ;
  PixZSToBeam=PixZToBeam;
  StepYRel=StepY;  // stepY is negative Ystep
  StepZRel=StepZ;
  
  return 1;
}

/// returns next clipped pixel. However : one pixel per row is returned at least
IndexType GiveNextClipped(void)     
{
  IndexType Step;
  // calculate new position and then return it;
  PixXToBeam += SUBUNITS;
  if ((PixXToBeam <= MaxXUnits) 
      && (PosX < BorderXH)) // else step in z or y
    {
      PosX++;
      StepYRel-=StepX;
      StepZRel-=StepX;
      return StepX;
    }

  PixZToBeam += SUBUNITS;
  if ((PixZToBeam <= MaxZUnits) 
      && (PosZ < BorderZH)) // else step in y
    {
      PosZ++;
      PosX = PosXS;
      PixXToBeam=PixXSToBeam;
      Step = StepZRel; //  + StepZ is allready preset in StepZRel !;
      StepZRel = StepZ;  // Set again for next stepz
      StepYRel -= Step;
      return Step;
    }

  if (--PosY < BorderYL)
    return 0;  // no more pixels
  
  Step = StepYRel; //  StepY is allready preset !;
  PixXSToBeam += XDecrement; // we are closer to the center now
  if (PixXSToBeam - MinXUnits >= SUBUNITS) // one left ?
    {
      if (PosXS > BorderXL) // >= BorderXL+1
	{
	  PixXSToBeam -= SUBUNITS;
	  Step -= StepX;
	  PosXS--;
	}
      else
	if (PixXSToBeam > MaxXUnits) return 0; // End of beam at left (lower x-) side
    }
  
  PixZSToBeam += ZDecrement; // we are closer to the center now
  if (PixZSToBeam - MinZUnits >= SUBUNITS) // one in -z direction ?
    {
      if (PosZS > BorderZL) // >= BorderZL+1
	{
	  PixZSToBeam -= SUBUNITS;
	  Step -= StepZ;
	  PosZS--;
	}
      else
	if (PixZSToBeam > MaxZUnits) return 0; // End of beam at lower z-side
    }
  
  StepYRel= StepY,StepZRel= StepZ;
  PosX = PosXS;
  PixXToBeam=PixXSToBeam;
  PosZ = PosZS;
  PixZToBeam=PixZSToBeam;
  return Step;
}
};

#endif
