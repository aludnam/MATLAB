#ifndef clipperator_h  		// -*- C++ -*- 
#define clipperator_h

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


template<class BaseType> void SWAP(BaseType & x,BaseType & y) {BaseType z=x; x=y; y=z;}


#define rint(x) (((x) > 0) ? (int((x)+0.5)) : (- (int((-(x))+0.5))))

/// type for indices in arrays and transformed arrays (-> signed type)
typedef int IndexType;

/** This class knows how to perform clipping in arrays */
/// The Clipperators know only about indices in arrays
class ClipperatorNoClip {  /// == NOCLIP No Clipping
public:
  IndexType SizeX,SizeY,SizeZ;
  /// the coordinates of the last clipped pixel in the unit-system (flipped, but NOT sheared!)
  IndexType PositionX,PositionY,PositionZ;  // Readable from others !
  /// the coordinates of the (transformed) center measured in float Pixels
  CalcFloatType CtrX,CtrY,CtrZ; 

void InitClippingFirst3D(CalcFloatType DirectionX,CalcFloatType DirectionY,CalcFloatType DirectionZ,CalcFloatType Length,CalcFloatType Width,CalcFloatType Depth)
{
  return ;  // nothing to be done here
}

int InitClippingSecond3D(CalcFloatType CenterX,CalcFloatType CenterY,CalcFloatType CenterZ) // returns boolean if cliprange is valid
{
  CtrX=CenterX;  CtrY=CenterY;  CtrZ=CenterZ;
  PositionX = rint(CenterX);   PositionY = rint(CenterY);   PositionZ = rint(CenterZ);
  return 1;
}

IndexType GiveNextClipped(void)
{
  return 0;
}
};

class ClipperatorAllClip : public ClipperatorNoClip {
public:
void InitClippingFirst3D(CalcFloatType DirectionX,CalcFloatType DirectionY,CalcFloatType DirectionZ,CalcFloatType Length,CalcFloatType Width,CalcFloatType Depth)
{
  return ;  // nothing to be done here
}

int InitClippingSecond3D(CalcFloatType CenterX,CalcFloatType CenterY,CalcFloatType CenterZ)
{
  CtrX=CenterX;  CtrY=CenterY;  CtrZ=CenterZ;
  PositionX = 0;   PositionY =0;   PositionZ =0;
  return 1;
}

/// returns all pixels of the array
IndexType GiveNextClipped(void)    
{
    if (++PositionX >= SizeX)
      {
	PositionX=0;
	if (++PositionY >= SizeY)  
	  {
	    PositionY=0;
	    if (++PositionZ >= SizeZ)
	      {
		PositionX = -1;
		PositionZ = 0;
		return 0; // Ready
	      }
	  }
      }
    return 1;  // just advance one index
}
};


/// clips neighbouring pixels
class ClipperatorNeighbClip : public ClipperatorNoClip {
protected:
  /// incremental steps when incrementing in the unitsystem
  IndexType StepX, StepZ, StepY;

  IndexType PosXS,PosYS,PosZS;
  IndexType BorderXH,BorderYH,BorderZH;

public:

  /// returns Float distance to clipping center in the TRANSFORMED coordinate system !
CalcFloatType XDistToCenter(void)
  {
    return (CtrX-(CalcFloatType) PositionX);
  }

  /// returns Float distance to clipping center in the TRANSFORMED coordinate system !
CalcFloatType YDistToCenter(void)  
  {
    return (CtrY-(CalcFloatType) PositionY);
  }

/// returns Position in transformed coordinate system
IndexType XPos(void)
  {
    return PositionX;
  }
/// returns Position in transformed coordinate system
IndexType YPos(void)
  {
    return PositionY;
  }

/// returns Float distance to clipping center in the TRANSFORMED coordinate system !
CalcFloatType ZDistToCenter(void)
  {
    return (CtrZ-(CalcFloatType) PositionZ);
  }

/// returns Position in transformed coordinate system
IndexType ZPos(void)
  {
    return PositionZ;
  }



void InitClippingFirst3D(CalcFloatType DirectionX,CalcFloatType DirectionY,CalcFloatType DirectionZ,CalcFloatType Length,CalcFloatType Width,CalcFloatType Depth)
{
  return ;  // nothing to be done here
}

int InitClippingSecond3D(CalcFloatType CenterX,CalcFloatType CenterY,CalcFloatType CenterZ)
{
  CtrX=CenterX;  CtrY=CenterY;  CtrZ=CenterZ;

  StepX=-1;   // These Steps decribe, how to into next dimension, when one dimension is finished in the normal case (inside the array)
  StepY=-SizeX+1;  // Other definition of Y-direction !!
  StepZ=-SizeX*SizeY+SizeX+1;

  BorderXH = BorderYH = BorderZH = 1;
  PosXS = PosYS = PosZS = 1; // will be used as counter to 0
  
  if (CenterX < -1.0)
    return 0;
  if (CenterY < -1.0)
    return 0;
  if (CenterZ < -1.0)
    return 0;

  PositionX = (int)(CenterX+PosXS);  // first pixel X-pos
  PositionY = (int)(CenterY+PosYS);  // Is allways >= 0 !!
  PositionZ = (int)(CenterZ+PosZS);

  // will all pixels be inside X ?
  if (PositionX >= SizeX) 
    if (PositionX > SizeX)  // outside array
      return 0;
    else
      {    // X-direction is limited to only one pixel
	BorderXH = 0;
	PosXS = 0;
	PositionX --;
	StepY -= 1;  // do not add an x-step
	StepZ -= 1;  // do not add an x-step
      }
  
  if (PositionX <= 0)
    {    // X-direction is limited to only one pixel
      BorderXH = 0;
      PosXS = 0;
      StepY -= 1;  // at turnaround do a step in -x too
      StepZ -= 1;  // at turnaround do a step in -x too
    }
  
  // will all pixels be inside Y ?
  if (PositionY >= SizeY)
    if (PositionY > SizeY)
      return 0;
    else
      {    // Y-direction is limited to only one pixel
	BorderYH = 0;
	PosYS = 0;
	PositionY --;
	StepZ-=SizeX;  // at turnaround do a step in -y too
      }
  
  if (PositionY <= 0)
    {    // Y-direction is limited to only one pixel
      BorderYH = 0;
      PosYS = 0;
      StepZ-=SizeX;  // at turnaround do a step in -y too
    }
  
  // will all pixels be inside Z ?
  if (PositionZ >= SizeZ)
    if (PositionZ > SizeZ)
      return 0;
    else
      {    // Z-direction is limited to only one pixel
	BorderZH = 0;
	PosZS = 0;
	PositionZ --;
      }
  
  if (PositionZ <= 0)
    {    // Z-direction is limited to only one pixel
      BorderZH = 0;
      PosZS = 0;
    }

  
  if (((double) ((int)(CenterX)) == CenterX)  // directly on the pixel :  much faster in this case
     && BorderXH == 1)
    {
      BorderXH = 0;
      PosXS = 0;
      PositionX = (int)(CenterX);  // first pixel X-pos
      StepY -= 1;  // do not add an x-step
      StepZ -= 1;  // do not add an x-step
    }
  if (((double) ((int)(CenterY)) == CenterY)  // directly on the pixel
     && BorderYH == 1)
    {
      BorderYH = 0;
      PosYS = 0;
      PositionY = (int)(CenterY);  // first pixel X-pos
      StepZ-=SizeX;  // at turnaround do a step in -y too
    }
  if (((double) ((int)(CenterZ)) == CenterZ)  // directly on the pixel
     && BorderZH == 1)
    {
      BorderZH = 0;
      PosZS = 0;
      PositionZ = (int)(CenterZ);  // first pixel X-pos
    }

  return 1;
}

IndexType GiveNextClipped(void)     // returns index to be added to actual position
{
  -- PositionX;
  if (--PosXS >= 0)
    {
      return StepX;
    }

  -- PositionY;
  if (-- PosYS >= 0)
    {
      PosXS = BorderXH;
      PositionX += BorderXH+1;
      return StepY;
    }
  
  -- PositionZ;
  if (-- PosZS >= 0)  // then : must be == 0
    { 
      PosYS = BorderYH;  // PosYS max must be 1 ! Otherwise PosZS should be < 0 by mow !
      PosXS = BorderXH;
      PositionY += BorderYH +1;
      PositionX += BorderXH +1;
      return StepZ;
    }
  else
    return 0;
}

};
 

#endif
