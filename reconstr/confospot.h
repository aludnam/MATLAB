
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

#include <iostream>		// -*- C++ -*- 
#include <math.h>
#include <string>
#include "vec.h"
#include "point.h"
#include "tarrayclipf.h"
#include "parseargs.h"

#include "projector.h"
#include "projectorpos.h"

#include "vecprojection.h"
#include "comperator.h"
#include "clipperator3d.h"

using namespace std;

// This program tries to select spots in an confocal image and determines their position

static const int SpotDimension=4;   // Value, Xpos, YPos, ZPos
int SpotNumber=55;                 // Maximal Number of Spots to detect, might be changed
CalcFloatType BACKGR=10; // // 67;            // A zero is dagerous (div by zero ..), might be changed later

bool verbose=false;   // print verbose information ?

double NewSpotLimit=3.0;      // 45;  100; // 22;   // Above that new Spots are added (value after subtraction of background)
double ThrowAwayLimit=0;   // Spots below are remover for the last iteration
double IdentDistance=5.0;          // New Spots can be identified as frozen spots
double FreezeLimit=20;             // If spot goes below this limit, it will be frozen
double SmallestIntenCorr=0.1; // If correction is below that, a new spot is added
int AddTogether=1;            // Number of Spots that are added together each round
int region=1;                 // region to be integrated when maximum is found for starting pos  
int IterationsBtw=300;        // Maximum number of iterations between adding a new spot
int IterationsBtw2=3;         // Iterations after the limit is exeeded, but still with only some spots
int IterationsAfter=50;       // Iterations after no more spots have been added
double GuessFactor=1.0;       // Intensity is guessed a little higher than maximum in remain
const int GearFak=20;         // Faktor for Overrelaxation when convergence too slow

int PSFSizeX=12; // 14; // The PSF has to be zero at the border of the region !!
int PSFSizeY=12; // 14
int PSFSizeZ=12; // 14

static double InitShiftX=0.0;    // -0.3;  // Shift to spot applied when adding a spot ( estimated maximum of PSF to Center)
static double InitShiftY=0.0;    // -0.3;  // These values are computed from the PSF
static double InitShiftZ=0.0;    // 0.3; 

static int ActualSpotNumber=0;    // Momentary Number of Spots in iteration 

string RemainFile, PrjFile;

// typedefs for global types

typedef float ArrayBType;
typedef float PSFBType;

typedef Dyn3dArray(ArrayBType,Clipperator3DClip)  TPrjArray;  // Clipping is needed 
typedef Dyn2dArray(ArrayBType,Clipperator3DClip)  TImgArray;  // Spots array with clipping 
typedef Dyn3dArray(PSFBType,ClipperatorNeighbClip)  TPSF;  

typedef ProjectorPos<TImgArray,SpotDimension-1,TPSF,TImgArray> TProjector; // without value
typedef ClipProjector<TPrjArray> TProjectorSelf;   
typedef Image<TImgArray,TProjector>                TImage;
typedef VecProjection<TPrjArray,TImage,TProjectorSelf> TProjection;

// following objects are global, because they are too big to be on the main stack

static CalcFloatType MaxInten=0.0;

// Images needed

static   TImage ReconSpotsImg; // SpotWeights // Zero as first guess  // Array, ItsPosition, ItsDirection
static   TImage CorrectionImg;
static   TImgArray OldCorr,SaveArray;

// Needed projections follow. They dont need clipping.

static   TPSF ValuePSFArray,XPSFArray,YPSFArray,ZPSFArray;         // z-Array is empty..
static   TPSF * PSFArrays[SpotDimension]={&ValuePSFArray,&XPSFArray,&YPSFArray,&ZPSFArray}; // for every iteration component there is a PSF

CalcFloatType PSFIntegral=1.0;

static   char ActiveSpots[20000];       // not more spots than this allowed !!
static   TProjector  SpotProjector(& ReconSpotsImg.Array,ActiveSpots,PSFArrays,& PSFIntegral,1.0,1.0,1.0); // can also clip in projection
static   TProjectorSelf SelfProjector(3);

static   TProjection MeasuredPrj(& SpotProjector,& ReconSpotsImg,& SelfProjector); // Projector,pPos and pDir 
static   TProjection ReconPrj(& SpotProjector,& ReconSpotsImg,& SelfProjector);    // Serves as data memory for computation of intermediate projection


void Constrain(TImage * CorrImg,int j ,double Overrelax, double MaxMove=0.5)      // This procedure constrains correction-values to certain maximals
{
  int i;
  CalcFloatType tmp;

      if ((fabs((OldCorr.Value(0,j) - CorrImg->Array.Value(0,j))/CorrImg->Array.Value(0,j)) < 3.0/GearFak)      // Correction did not change ! -> next gear
        && (fabs((OldCorr.Value(1,j) - CorrImg->Array.Value(1,j))/CorrImg->Array.Value(1,j)) < 4.0/GearFak)
        && (fabs((OldCorr.Value(2,j) - CorrImg->Array.Value(2,j))/CorrImg->Array.Value(2,j)) < 4.0/GearFak)
        && (fabs((OldCorr.Value(3,j) - CorrImg->Array.Value(3,j))/CorrImg->Array.Value(3,j)) < 4.0/GearFak)
	  &&  fabs(CorrImg->Array.Value(0,j)) > 0.0001 )          // No truncation errors
          {
#ifdef PRINTALL
	    cout << " Spot " << j +1 <<" Covergence to slow, increasing overrelax !\n";
            cout << flush;
#endif
	    for (i=0;i<4;i++) 
	      CorrImg->Array.SetValue(i,j,CorrImg->Array.Value(i,j) * GearFak);
	  }

      if ((tmp=fabs(CorrImg->Array.Value(1,j))) > MaxMove/Overrelax)    //  No Position change shall be greater than 1/2 Pixel
	{
	  for (i=1;i<4;i++)  // just restrict the movement but allow full intensity correction
	    CorrImg->Array.SetValue(i,j,CorrImg->Array.Value(i,j) * MaxMove/tmp);
	}

      if ((tmp=fabs(CorrImg->Array.Value(2,j))) > MaxMove/Overrelax)
	{
	  for (i=1;i<4;i++)  // just restrict the movement but allow full intensity correction
	    CorrImg->Array.SetValue(i,j,CorrImg->Array.Value(i,j) * MaxMove/tmp);
	}

      if ((tmp=fabs(CorrImg->Array.Value(3,j))) > MaxMove/Overrelax)
	{
#ifdef PRINTALL
	  cout << "constrained z \n";
#endif
	  for (i=1;i<4;i++)  // just restrict the movement but allow full intensity correction
	    CorrImg->Array.SetValue(i,j,CorrImg->Array.Value(i,j) * MaxMove/tmp);
	}

      // if (CorrImg->Array.Value(0,j) > - 0.8*OldCorr.Value(0,j))   // No undamped oszillations
      // for (i=0;i<4;i++) CorrImg->Array.SetValue(i,j,CorrImg->Array.Value(i,j) * (-0.8* OldCorr.Value(0,j)) / CorrImg->Array.Value(0,j));
}

void Constrains(TImage * CorrImg,double Overrelax, double MaxMove=0.5)
{
  int j;
  for (j=0;j<ActualSpotNumber;j++)
    {
      Constrain(CorrImg,j,Overrelax, MaxMove);
    }

  OldCorr.Copy(& CorrImg->Array);
}

double IterationStep(CalcFloatType Overrelax, double MaxMove=0.5)
{    
  int j,maxSpot=-1;
  CalcFloatType maxic=0.0;

      for (j=0;j<ActualSpotNumber;j++)
        {
#ifdef PRINTALL
          cout << "Spot Nr. " << j+1 << "\n";
 	  cout << " Val : "  << ReconSpotsImg.Array.Value(0,j) << "\n";
	  cout << " XPos : " << ReconSpotsImg.Array.Value(1,j) << " X-Dist to ref. : " << ReconSpotsImg.Array.Value(1,j) - ReconSpotsImg.Array.Value(1,0) << "\n";
	  cout << " YPos : " << ReconSpotsImg.Array.Value(2,j) << " Y-Dist to ref. : " << ReconSpotsImg.Array.Value(2,j) - ReconSpotsImg.Array.Value(2,0) << "\n";
	  cout << " ZPos : " << ReconSpotsImg.Array.Value(3,j) << " Z-Dist to ref. : " << ReconSpotsImg.Array.Value(3,j) - ReconSpotsImg.Array.Value(3,0) << "\n";
#endif
	} 


      // ReconPrj.IterateImg(&ReconSpotsImg,& MeasuredPrj,& CorrectionImg,Overrelax);
      CorrectionImg.Array.Clear();
      ReconPrj.CorrectToCorrectionImg(&ReconSpotsImg,& MeasuredPrj,& CorrectionImg);
      Constrains(& CorrectionImg,Overrelax,MaxMove);
      ReconSpotsImg.CorrectImg(& CorrectionImg, Overrelax);


      for (j=0;j<ActualSpotNumber;j++)
	{
#ifdef PRINTALL     
	cout << "Spot " << j+1 << " Correct : Val " << CorrectionImg.Array.Value(NoSel,0,j) << 
	  " X " << CorrectionImg.Array.Value(NoSel,1,j) << 
	  " Y " << CorrectionImg.Array.Value(NoSel,2,j) << 
	  " Z " << CorrectionImg.Array.Value(NoSel,3,j) << "\n"; 
#endif
	if ( fabs(CorrectionImg.Array.Value(NoSel,0,j)) > fabs(maxic)) maxic=CorrectionImg.Array.Value(NoSel,0,j), maxSpot=j;
	}
if (verbose) {        
      cout << "Maximum Value Correction = " << maxic << " at Spot Nr. : " << maxSpot+1 
	   << " dX : " << CorrectionImg.Array.Value(NoSel,1,maxSpot) 
	   << " dY : " << CorrectionImg.Array.Value(NoSel,2,maxSpot) 
	   << " dZ : " << CorrectionImg.Array.Value(NoSel,3,maxSpot) << "\n";
       cout << flush;
       }
       ReconPrj.Array.Clear();
       return fabs(maxic);
}

/// activates a spot and its near neighbours
void Activate(int SpotNr)
{
  int i;
  ActiveSpots[SpotNr] = 1;

  for (i=0;i<ActualSpotNumber;i++)
    {
        if ( ReconSpotsImg.Array.Value(0,i) != 0.0)  // frozen or deleted spots cannot be activated
        if ((fabs(ReconSpotsImg.Array.Value(NoSel,1,i) - ReconSpotsImg.Array.Value(NoSel,1,SpotNr)) < PSFSizeX+1+2)
      && (fabs(ReconSpotsImg.Array.Value(NoSel,2,i) - ReconSpotsImg.Array.Value(NoSel,2,SpotNr)) < PSFSizeY+1+2)
      && (fabs(ReconSpotsImg.Array.Value(NoSel,3,i) - ReconSpotsImg.Array.Value(NoSel,3,SpotNr)) < PSFSizeZ+1+2))
	{
	  ActiveSpots[i] = 1;
                if (verbose) 
                       cout << "activated spot Nr. " << i+1 << "\n";
	}
    }
 
}

void ActivateAll(void)
{
  int i;
  for (i=0;i<ActualSpotNumber;i++) ActiveSpots[i] = 1;
}

void ActivateNone(void)
{
  int i;
  for (i=0;i<ActualSpotNumber;i++) ActiveSpots[i] = 0;
}


void  ReactivateNearestFrozenDist()  // uses the last spot and joines it eventually with a frozen one
{
   double px =  ReconSpotsImg.Array.Value(1,ActualSpotNumber-1);
   double py =  ReconSpotsImg.Array.Value(2,ActualSpotNumber-1);
   double pz =  ReconSpotsImg.Array.Value(3,ActualSpotNumber-1);
   int minnr=0;
   double dx,dy,dz,dist2,mindist2;
   mindist2 = IdentDistance +1.0;               // start with something bigger than IdentDistance to find minimum
   mindist2 *= mindist2;
   // cout << "\n Testing reactivation ... spot "<< ActualSpotNumber -1 << "\n";

   for (int j=0;j<ActualSpotNumber-1;j++)
     {
  //        cout << "\nspot " << j <<" Value " << ReconSpotsImg.Array.Value(0,j) << "\n";
         if (ReconSpotsImg.Array.Value(0,j) == 0)   // was this spot marked as silent ?
         {
                        dx =    ReconSpotsImg.Array.Value(1,j) - px;
                        dy =    ReconSpotsImg.Array.Value(2,j) - py;
                        dz =    ReconSpotsImg.Array.Value(3,j) - pz;
                        dist2 = dx*dx+dy*dy+dz*dz;
  //                      cout << "\nTesting for reactivation: spot " << j <<" dist " << sqrt(dist2) << "\n";
                        if (dist2 < mindist2)
             {
                            mindist2 = dist2;
                            minnr=j;
             }
          }
      }
      if (mindist2 <= IdentDistance*IdentDistance)  // Spots are joined
      {
                  cout << "JOINED spots " << ActualSpotNumber << " and " << minnr << "at distance " << sqrt(mindist2) << "\n";
                  ReconSpotsImg.Array.SetValue(0,minnr, ReconSpotsImg.Array.Value(0,ActualSpotNumber-1));
                  ReconSpotsImg.Array.SetValue(1,minnr,px);
                  ReconSpotsImg.Array.SetValue(2,minnr,py);
                  ReconSpotsImg.Array.SetValue(3,minnr,pz);
                  ReconSpotsImg.Array.SetValue(0,ActualSpotNumber-1,0.0);
                  ReconSpotsImg.Array.SetValue(1,ActualSpotNumber-1,-1000.0);
                  ReconSpotsImg.Array.SetValue(2,ActualSpotNumber-1,-1000.0);
                  ReconSpotsImg.Array.SetValue(3,ActualSpotNumber-1,-1000.0);
                  ActualSpotNumber--;
      }
}


int AddNewSpot(void)
{
CalcFloatType Intensity;
static Vector  MaxVector(0,0,0);

  if (ActualSpotNumber >= SpotNumber)
    {
    cout << "Maximal number of " << SpotNumber << " Spots reached ! Cannot add more spots \n";
    return 0;
    }

  ActivateAll();   // spots will be activated for projection
  ReconPrj.Array.Set(0);
  ReconPrj.ProjectImage(&ReconSpotsImg); // only adds into image

  if (PrjFile != "")
    ReconPrj.Array.DSave(true,PrjFile.c_str());

  ReconPrj.Array.SubsDivSqrtAdd(&MeasuredPrj.Array);

  if (RemainFile != "")
    ReconPrj.Array.DSave(true,RemainFile.c_str());
 
  ReconPrj.Array.MaxIntensityVec(& MaxVector, & MaxInten, region);   // find the vector to the maximum Intensity

  double px=MaxVector.comp(0),
    py=MaxVector.comp(1),
    pz=MaxVector.comp(2);

  MaxInten=ReconPrj.Array.MeanRegionIntensity(rint(px),rint(py),rint(pz),region,region,region);
   Intensity=MeasuredPrj.Array.MeanRegionIntensity(rint(px),rint(py),rint(pz),region,region,region) - BACKGR;

  // Intensity=MeasuredPrj.Array.MaxIntensityNeighbour( px,py,pz,region)-BACKGR;  // What was this?
  
  if (MaxInten >= NewSpotLimit)
    {
      ReconSpotsImg.Array.SetValue(0,ActualSpotNumber,GuessFactor*Intensity);
      // ReconSpotsImg.Array.SetValue(1,ActualSpotNumber,MaxVector.comp(0)-InitShiftX);
      // ReconSpotsImg.Array.SetValue(2,ActualSpotNumber,MaxVector.comp(1)-InitShiftY);
      // ReconSpotsImg.Array.SetValue(3,ActualSpotNumber,MaxVector.comp(2)-InitShiftZ);
      ReconSpotsImg.Array.SetValue(1,ActualSpotNumber,px-InitShiftX);
      ReconSpotsImg.Array.SetValue(2,ActualSpotNumber,py-InitShiftY);
      ReconSpotsImg.Array.SetValue(3,ActualSpotNumber,pz-InitShiftZ);

      cout << "Added Spot # " << ActualSpotNumber << " Remain : "<< MaxInten << " Intensity " << GuessFactor*Intensity << " starting at ";
      MaxVector.show();
      cout << flush;

      ActualSpotNumber++;

      // now activate only this spot and surroundings
      ActivateNone();  
      Activate(ActualSpotNumber-1);

      IterationStep(0.8);                  // correct others
      ReconSpotsImg.Array.SetValue(0,ActualSpotNumber-1,GuessFactor*Intensity);   // reset intensity of this spot
      IterationStep(0.8);                  // correct others
      ReconSpotsImg.Array.SetValue(0,ActualSpotNumber-1,GuessFactor*Intensity);
      IterationStep(0.8);                  // correct others
      ReconSpotsImg.Array.SetValue(0,ActualSpotNumber-1,GuessFactor*Intensity);

      ReactivateNearestFrozenDist();  // identifies an old spot and eventually joines the two

      return 1;
    }
  else
    {
      cout << "Spot of intensity " << MaxInten << " is lower than " << NewSpotLimit << " . No more Spots to be found !\n";
      return 0;
    }

}
int AddNewSpots(void)
{
  // int i;

  //  for (i=0;i<AddTogether;i++)
    //  {
      if ((AddNewSpot()) == 0) return 0;
    //  }

  OldCorr.Set(0.0);
  return 1;      
}

void ClearSpot(int i)
{
   ReconSpotsImg.Array.SetValue(0,i,0.0);              // delete all entries
   ReconSpotsImg.Array.SetValue(1,i,-1000.0);          // clipping will make it faster :)
   ReconSpotsImg.Array.SetValue(2,i,-1000.0);
   ReconSpotsImg.Array.SetValue(3,i,-1000.0);
}

void FreezeSpot(int i)
{
   ReconSpotsImg.Array.SetValue(0,i,0.0);              // delete intensity but keep spot position
   if (SaveArray.GetSize(1) >  i)
   if (SaveArray.Value(0,i) != 0)
   {
           ReconSpotsImg.Array.SetValue(1,i,SaveArray.Value(1,i));     
           ReconSpotsImg.Array.SetValue(2,i,SaveArray.Value(2,i));     // copy the saved values into this frozen spot
           ReconSpotsImg.Array.SetValue(3,i,SaveArray.Value(3,i));     // this corrects for possible drifts, that have occured
   }
   // cout << "Froze spot Nr. " << i << "\n\n";
}

void FreezeDimSpots(void)
{
   // cout << "Trying to freeze\n";

   for (int j=0;j<ActualSpotNumber;j++)
     {
         // cout << "testing ... spot " << j << " val " << ReconSpotsImg.Array.Value(0,j)  << " pos (" <<
         //        ReconSpotsImg.Array.Value(1,j) << ", " << ReconSpotsImg.Array.Value(2,j) << ", " << ReconSpotsImg.Array.Value(3,j) <<")\n";
       if ( ReconSpotsImg.Array.Value(0,j) <= FreezeLimit)
	   FreezeSpot(j);
     }
}

void RemoveSmallSpots(void)
{
   for (int j=0;j<ActualSpotNumber;j++)
     {
       if ( ReconSpotsImg.Array.Value(0,j) <= ThrowAwayLimit)
	   ClearSpot(j);
     }
}


void ComputePSFCenter(TPSF * PSF)         // Computes the Center of Mass of the PSF for initial guesses
{
  Vector COM(3);
  PSF->ClippedCOM(0.0,& COM);  // COM of positive PSF
  InitShiftX= COM.comp(0) - PSFSizeX/2.0;
  InitShiftY= COM.comp(1) - PSFSizeY/2.0;
  InitShiftZ= COM.comp(2) - PSFSizeZ/2.0;
  cout << " PSF ShiftX : " << InitShiftX << " Y : " << InitShiftY << " Z : " << InitShiftZ << "\n";
}
