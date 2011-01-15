#ifndef ProjectorPos_h		// -*- C++ -*- 
#define ProjectorPos_h

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
#include "tarrayclipf.h"  // fast 2 and 3d array templates


/** knows how to project spots at a given position (float !) and of given intensity.
 Projector is asked about the weight it gives a pixel at a given position "Position"
 The zero Vector component is the "element= value, pos1,pos2,pos3, ..."
 The pPosition indicates the spotnumber asked and which component (value, posx,...)
 This strange projector has to know of the image array (SpotList), to read out the positionvalues of the spots
 which dimension in vector stores the element information (VALUE,POSX,POSY,..)
 Handles the weight of generalized image-pixels (positions, values of spots) */
template <class ArrayType,int SpotDim, class ArrayClassType,class TSpotArray>
class ProjectorPos : public ClipProjector<ArrayType> 
{
private:

  Vector * PsfVec;
  /// Vector of Projection (what is normaly Position)
  Vector * PrjPos;
  /// Arrays containing PSF for VALUE, POSX, POSY, POSZ, ...
  ArrayClassType ** PSFArrays;
  /// is used in the backprojection
  CalcFloatType * PSFIntegral;

  /// Factor for enlargement of psf
  CalcFloatType PsfFac;

public:

/// Should be the same as the image to reconstruct
TSpotArray * SpotList;

/// Marks if a spot is active at the moment (0 : not active)
char * ActiveList;

  /// Constructor
ProjectorPos(TSpotArray * SList,char * AList,
	     ArrayClassType ** PSFs,CalcFloatType * pInt,CalcFloatType Length,CalcFloatType Width, CalcFloatType Depth) : ClipProjector<ArrayType> (3), PsfFac(1.0)
  {
    static const Vector YOne(0.0,1.0,0.0);
    SpotList=SList;
    ActiveList=AList;
    PSFArrays=PSFs;              // This is a vector of arrays belonging to the different PSFs of x-coordinates of SpotList
    PsfVec= new Vector(SpotDim); // SpotDim Dimensions
    PrjPos= new Vector(SpotDim);
    ClipProjector<ArrayType>::SetClip(Length,Width,Depth);
    ClipProjector<ArrayType>::SetDir(& YOne);
    PSFIntegral=pInt;
  }

  /// Copy Constructor has to allocate space for Vectors, but not for Array
ProjectorPos(const ProjectorPos & CopyFromPrj) : ClipProjector<ArrayType> (CopyFromPrj)
  {
    SpotList=CopyFromPrj.SpotList;
    ActiveList=CopyFromPrj.ActiveList;
    PSFArrays=CopyFromPrj.PSFArrays;              // This is a vector of arrays belonging to the different PSFs of x-coordinates of SpotList
    PsfVec= new Vector(SpotDim); // SpotDim Dimensions
    PrjPos= new Vector(SpotDim);
    PSFIntegral=CopyFromPrj.PSFIntegral;
    PsfFac=CopyFromPrj.PsfFac;
  }

void SetPrjPos(CalcFloatType X,CalcFloatType Y, CalcFloatType Z)
  {
    PrjPos->copy(X,Y,Z); // what if vector only 2d ??
  }

void SetPrjPos(Vector * pVec)
  {
    PrjPos->copy(pVec);
  }

void SetPsfFac(CalcFloatType factor)
  {
    PsfFac = factor;
  }

int IsActiveSpot(int SpotNr)
  {
    // return 1;
    // return (SpotList->Value(NoSel,0,SpotNr) != 0);
    return (ActiveList[SpotNr] != 0);
  }

/// computes intersection of itself with one pixel.
CalcFloatType ProjectYourselfFwd(ArrayType * pArr)
{
  int i;
  const int ElementNr = pArr->Clipperator.TXPos();  // X-component is Element (PosX,PosY,PosZ,..)
  const int SpotNr    = pArr->Clipperator.TYPos();

  if (SpotList->Value(0,SpotNr) == 0.0) return 0.0;   // Value == 0 -> These Spots are not active yet
  else
  {
  for (i=0;i<SpotDim;i++)
     PsfVec->ChangeComp(i,SpotList->Value(i+1,SpotNr));    // 0-Position is Value of Spot, therefore one has to be added

  if (ElementNr == 0)   // case VALUE : The Positioncomponents of the spot tell, which value (from the psf-array) to return
    return PSFArrays[0]->Value(CenterSel+ClipSel+InterpSel,PrjPos,PsfVec,PsfFac); 
  else  
    return 0.0;
  }
}

/// computes intersection of itself with one pixel.
CalcFloatType ProjectYourselfBwd(ArrayType * pArr)
{
  int i;
  const int ElementNr = pArr->Clipperator.TXPos();  // X-component is Element (PosX,PosY,PosZ,..)
  const int SpotNr    = pArr->Clipperator.TYPos();

  if (SpotList->Value(0,SpotNr) == 0.0) return 0.0;   // Value == 0 -> These Spots are not active yet
  else
  {
  for (i=0;i<SpotDim;i++)
     PsfVec->ChangeComp(i,SpotList->Value(i+1,SpotNr));    // 0-Position is Value of Spot, therefore one has to be added

  if (ElementNr == 0) // case VALUE : The Positioncomponents of the spot tell, which value (from the psf-array) to return
    // the multiplication with the spot-intensity is, because Sub-comperator is used
   return -1.4/(* PSFIntegral)* PsfFac *PsfFac*PsfFac*SpotList->Value(0,SpotNr)*PSFArrays[0]->Value(CenterSel+ClipSel+InterpSel,PrjPos,PsfVec,PsfFac);   
    // This multiplication with the integral Intensity ensures equal convergence
  else
    return 3.2/(* PSFIntegral)* PsfFac*PSFArrays[ElementNr]->Value(CenterSel+ClipSel+InterpSel,PrjPos,PsfVec,PsfFac); //  5/SpotList->Value(ElementNr,SpotNr)*PSF...;    // Return derivation   * PsfFac ?
  // return 3.2/(* PSFIntegral)* PsfFac*PsfFac*PsfFac*PsfFac*PSFArrays[ElementNr]->Value(CenterSel+ClipSel+InterpSel,PrjPos,PsfVec,PsfFac); //  5/SpotList->Value(ElementNr,SpotNr)*PSF...;    // Return derivation   * PsfFac ?
  }
}

};

#endif

