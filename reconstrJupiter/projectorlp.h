#ifndef ProjectorLp_h		// -*- C++ -*- 
#define ProjectorLp_h

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
#include "clipprojector.h"


/** represents a multi (up to 3-d) dimensional linear (in each dimension) interpolation between data-values in array */
template <class ArrayType>
class ProjectorLp : public ClipProjector<ArrayType>       // Position is allways Center.
{
protected:
CalcFloatType LinWght(CalcFloatType val)
  {
    if (val > 0.0)
      if (val < 1.0)
	return (1.0-val);
      else
	return 0.0;
    else
      if (val > -1.0)
	return (1.0+val);
      else
	return 0.0;
  }

public:

ProjectorLp(const int Dim,CalcFloatType Length, CalcFloatType Width, CalcFloatType Depth) : ClipProjector<ArrayType> (Dim)
    {ClipProjector<ArrayType>::SetClip(Length,Width,Depth);}

/// Copy Constructor has to allocate space for Vectors, but not for Array
ProjectorLp(const ProjectorLp & CopyFromPrj) : ClipProjector<ArrayType> (CopyFromPrj) {}

/** computes intersection of itself with one pixel.
 pPosition is Position in Image (image-coordinates), p Direction is direction of image
 computes intersection of itself with one pixel. */
CalcFloatType ProjectYourselfFwd(ArrayType * pArr)
{
  // return 1.0;
  return LinWght(pArr->Clipperator.XDistToCenter()) * LinWght(pArr->Clipperator.YDistToCenter()) * LinWght(pArr->Clipperator.ZDistToCenter());
  // return Position.MulLinDist(pArr->LastValidPixel);
}


  /** computes intersection of itself with one pixel.
   interpolation for the backward (into-image) projection, if the raster is finer in the image then in the projection
   computes intersection of itself with one pixel. */
CalcFloatType ProjectYourselfBwd(ArrayType * pArr)
{  
   // return 1.0;
  return LinWght(pArr->Clipperator.XDistToCenter()) * LinWght(pArr->Clipperator.YDistToCenter()) * LinWght(pArr->Clipperator.ZDistToCenter());
  // return Position.MulLinDist(pArr->LastValidPixel);
}

};

template <class ArrayType>
class ProjectorLpZ : public ProjectorLp<ArrayType> // zoomed LinWeight projector
{
private:
  CalcFloatType ScaleX,ScaleY,ScaleZ;
  int OnlyNearest;

protected:

  CalcFloatType LinWght(CalcFloatType dist,CalcFloatType size) // Takes source and dest Pixels as beeing rectangular
    {
      CalcFloatType Max,Min;
      
      Max=dist+size/2.0;
      if (Max > 0.5)
	Max=0.5;

      Min=dist-size/2.0;
      if (Min < -0.5)
	Min=-0.5;

      return Clip(Max-Min);
    }

  CalcFloatType LinWghtT(CalcFloatType dist) // Takes source Pixels as beeing Triangular and destination pixels as deltas
    {                                        // == normal linear interpolation
      CalcFloatType d=fabs(dist);
      if (d> 1.0)
	return 0.0;
      else
	return 1.0-d;
    }
 
  CalcFloatType LinWghtTR(CalcFloatType dist,CalcFloatType size) // Takes source Pixels as beeing Triangular (1.0 wide) and destination pixels as rectangular (size wide)
    {                                        // == normal linear interpolation
      CalcFloatType Max,Min,h1,h2;
      
      Max=dist+size/2.0;
      if (Max > 0.5)
	Max=0.5;

      Min=dist-size/2.0;
      if (Min < -0.5)
	Min=-0.5;


      if (Max <= Min)
	return 0.0;

      h1 = (0.5-fabs(Max));  // if Max < 0 this will yield in a correct substraction
      h2 = (0.5-fabs(Min));  // if Min < 0 this will yield in a correct substraction
	
      return (Max*(1.0+h1)- Min*(1.0+h2))/size;  // / 2.0 * 2.0  // / size is for integral norm
    }
 

public:

ProjectorLpZ(const int Dim,CalcFloatType Length, CalcFloatType Width, CalcFloatType Depth=0) : ProjectorLp<ArrayType> (Dim,Length,Width,Depth), ScaleX(1.0), ScaleY(1.0), ScaleZ(1.0), OnlyNearest(0)
{}

ProjectorLpZ(const ProjectorLp<ArrayType> & CopyFromPrj) : ProjectorLp<ArrayType> (CopyFromPrj), ScaleX(1.0), ScaleY(1.0), ScaleZ(1.0), OnlyNearest(0)
{}

ProjectorLpZ(const ProjectorLpZ & CopyFromPrj) : ProjectorLp<ArrayType> (CopyFromPrj), ScaleX(CopyFromPrj.ScaleX),ScaleY(CopyFromPrj.ScaleY),ScaleZ(CopyFromPrj.ScaleZ), OnlyNearest(CopyFromPrj.OnlyNearest)
    { } 


CalcFloatType Clip(CalcFloatType val)
{
  return (val > 0.0) ? val : 0.0;
}

void copy(const ProjectorLpZ * APrj)
{
  ProjectorLp<ArrayType> ::copy(APrj);
  ScaleX = APrj->ScaleX;
  ScaleY = APrj->ScaleY;
  ScaleZ = APrj->ScaleZ;
  OnlyNearest = APrj->OnlyNearest;
}

void SetZoom(CalcFloatType X, CalcFloatType Y,CalcFloatType Z)
  {
    ScaleX=fabs(X);
    ScaleY=fabs(Y);
    ScaleZ=fabs(Z);
  }

void SetOnlyNearest(int val)   // This is not beautiful, but allows for fast projections, when the rage is zero
{
  OnlyNearest=val;
}

CalcFloatType ProjectYourselfFwd(ArrayType * pArr)
{
  if (! OnlyNearest)
    return LinWghtTR(pArr->Clipperator.XDistToCenter(),ScaleX)*LinWghtTR(pArr->Clipperator.YDistToCenter(),ScaleY)*LinWghtTR(pArr->Clipperator.ZDistToCenter(),ScaleZ);
  else
    return 1.0/(ScaleX*ScaleY*ScaleZ);
}


CalcFloatType ProjectYourselfBwd(ArrayType * pArr)
{  
  if (! OnlyNearest)
    return LinWghtTR(pArr->Clipperator.XDistToCenter(),ScaleX)*LinWghtTR(pArr->Clipperator.YDistToCenter(),ScaleY)*LinWghtTR(pArr->Clipperator.ZDistToCenter(),ScaleZ);
  else
    return 1.0/(ScaleX*ScaleY*ScaleZ);
}

};



#endif

