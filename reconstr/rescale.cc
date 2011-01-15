// rescale : this program reads in a 3-D datafile and rescales (translates, rotates and scales) it into the outputfile

/*   This file is part of a software package written by 
     Rainer Heintzmann
     Institute of Applied Optics and Information Processing
     Albert Ueberle Strasse 3-5
     69120 Heidelberg
     Tel.: ++49 (0) 6221  549264
     Current Address : Max Planck Inst. for Biophysical Chemistry, Am Fassberg 11, 37077 Goettingen, Germany
     Tel.: ++49 (0) 551 201 1029, e-mail: rheintz@gwdg.de  or rainer@heintzmann.de
     No garantee, whatsoever is given about functionallaty and  savety.
     No warranty is taken for any kind of damage it may cause.
     No support for it is provided !

     THIS IS NOT FREE SOFTWARE. All rights are reserved, and the software is only
     given to a limited number of persons for evaluation purpose !
     Please do not modify and/or redistribute this file or any other files, libraries
     object files or executables of this software package !
*/

#include <iostream>
#include <string>
#include "vec.h"
#include "projection.h"
#include "projectorlp.h"
#include "tarrayclipf.h"
#include "parseargs.h"
#include "clipperator3d.h"
#include "fftarray.h"


int MKEven(int number) // returns even number, that is lower than number
{
  return ((int) (number/2))*2;
}

Vector Dir(0.0,1.0,0.0);

typedef float ArrayBType;
typedef Dyn3dArray(ArrayBType,Clipperator3DClip)  TImgArray;
typedef Dyn3dArray(ArrayBType,Clipperator3DClip)  TPrjArray;   // Larger Clipping is needed
typedef Dyn3dArray(ArrayBType,Clipperator3DClip)  TPrjArray2;
// typedef Dyn3dArray(ArrayBType,ClipperatorNoClip)  TFastArray; 
typedef complex<float> ArrayBCType;

typedef TFFTArray<ArrayBCType>  TInterArray; // Only for FFT methods
typedef ProjectorLp<TImgArray>                     TProjector;
typedef ProjectorLpZ<TImgArray>                    TProjector2;
// typedef Projector                    TFastProjector;

typedef Image<TPrjArray,TProjector2>               TImage;       // TProjector2 ?
typedef Image<TPrjArray2,TProjector2>              TImage2;       // TProjector2 ?
typedef Projection<TPrjArray,TImage2,TProjector2>  TProjection;       // TProjector2 ?
typedef Projection<TPrjArray2,TImage,TProjector2> TProjection2;       // TProjector2 ?   SetZoom must be activated !!
// typedef Projection<TFastArray,TImage,TFastProjector> TFastProjection;    

static const ArrayBType FillValue=0.0;        

static TProjector2    SmallPrj(3,0.5,0.5,0.5);
// static TFastProjector    FastPrj(3);
static TProjector2 SmallPrjBwd(3,0,0,0);
static TProjection2  InputImg(& SmallPrjBwd);
static TProjection2 * pOutputPrj;
// static TFastProjection * pOutputPrjFast;

static TImgArray Vectors;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-sX X ...] [-i inputfile] [-o outputfile] [-fX X ...] [-mX X ...] [-b]\n" << flush;
  cerr <<  "[-sX -sY -sZ]  : sizes of destination image\n" << flush;
  cerr <<  "[-iX -iY -iZ]  : sizes of source image\n" << flush;
  cerr <<  "[-fX -fY -fZ]  : factors for multiplication in directi X,Y and Z\n" << flush;
  cerr <<  "[-mX -mY -mZ]  : Midpoint of destination image in koordinates of the source image \n" << flush;
  cerr <<  "[-k] :           use KDF format for in- and output \n" << flush;
  cerr <<  "[-b] :           integrate in source image instead in destination image (better for shrinkage)\n" << flush;
  cerr <<  "[-method nr] :    1,2: trilinear, 3: fft base method, 4: trilinear (triangle in source)\n" << flush;
  exit(-1);
}


int main(int argc, char *argv[])
{ 

  int OUTPUTSizeX=0;  // determine automatically
  int OUTPUTSizeY=0; 
  int OUTPUTSizeZ=0; 
  int INPUTSizeX=32;  // These is the standart size, if raw data is used
  int INPUTSizeY=32;  
  int INPUTSizeZ=32; 
  int VECSizeX=7;  // These is the standart size, if raw data is used
  int VECSizeY=1;  
  int VECSizeZ=1; 

int kflag=0,method=1,Elem,Elements,smx=0,smy=0,smz=0;

string INPUTFileName,OUTPUTFileName,VFileName;

double FacX=1.0;       // Scaling Factor for X
double FacY=1.0;       // Y
double FacZ=1.0;       // .0/3.09204;       // Z

double MidX=OUTPUTSizeX/2.0, MidY=OUTPUTSizeY/2.0, MidZ=OUTPUTSizeZ/2.0;

float  Fill=0.0;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-b",parg)) {method=2;continue;}  // integrate in source image
   if (readArg("-f",parg)) {method=3;continue;}  
   if (readArg("-method",& method,parg)) continue;  
   if (readArg("-i", INPUTFileName, parg)) continue;
   if (readArg("-vecin", VFileName, parg)) continue;
   if (readArg("-o", OUTPUTFileName, parg)) continue;
   if (readArg("-iX",& INPUTSizeX,parg)) continue; // input size X
   if (readArg("-iY",& INPUTSizeY,parg)) continue;
   if (readArg("-iZ",& INPUTSizeZ,parg)) continue;
   if (readArg("-sX",& OUTPUTSizeX,parg)) continue; // output size X
   if (readArg("-sY",& OUTPUTSizeY,parg)) continue;
   if (readArg("-sZ",& OUTPUTSizeZ,parg)) continue;
   if (readArg("-Pad",& Fill,parg)) continue;  // Fill outside with

   if (readArg("-fX",& FacX, parg)) continue; // scale factors
   if (readArg("-fY",& FacY, parg)) continue; 
   if (readArg("-fZ",& FacZ, parg)) continue; 

   if (readArg("-mX",& MidX, parg)) {smx=1;continue;} // center of output will be coord mid in input
   if (readArg("-mY",& MidY, parg)) {smy=1;continue;}
   if (readArg("-mZ",& MidZ, parg)) {smz=1;continue;} 

    usage(argv[0]);
  }

  if (VFileName!="")
    {
      Vectors.DLoad(kflag,VFileName.c_str(),"Float",& VECSizeX,& VECSizeY, & VECSizeZ,0);
      MidX=Vectors.Value(1,0,0);
      MidY=Vectors.Value(2,0,0);
      MidZ=Vectors.Value(3,0,0);
      smx=smy=smz=1;  // Do not calculate again !
    }


  Elements=InputImg.Array.DLoad(kflag,INPUTFileName.c_str(),"Float",& INPUTSizeX,& INPUTSizeY, & INPUTSizeZ,-1);

  if (OUTPUTSizeX <= 0)
    {
    OUTPUTSizeX = int(INPUTSizeX*FacX);
    if (OUTPUTSizeX <= 0) OUTPUTSizeX=1;
    cerr << "calculated OUTPUTSizeX to be " << OUTPUTSizeX << "\n";
    }
  if (OUTPUTSizeY <= 0)
    {
    OUTPUTSizeY = int(INPUTSizeY*FacY);
    if (OUTPUTSizeY <= 0) OUTPUTSizeY=1;
    cerr << "calculated OUTPUTSizeY to be " << OUTPUTSizeY << "\n";
    }
  if (OUTPUTSizeZ <= 0)
    {
    OUTPUTSizeZ = int(INPUTSizeZ*FacZ);
    if (OUTPUTSizeZ <= 0) OUTPUTSizeZ=1;
    cerr << "calculated OUTPUTSizeZ to be " << OUTPUTSizeZ << "\n";
    }

  if (!smx)  // If no midpoint given determine middle of Inputimage
    {
      MidX= INPUTSizeX/2.0;
      cerr << "calculated MidPointX to be " << MidX << "\n";
    }
  if (!smy)
    {
      MidY= INPUTSizeY/2.0;
      cerr << "calculated MidPointY to be " << MidY << "\n";
    }
  if (!smz)
    {
      MidZ= INPUTSizeZ/2.0;
      cerr << "calculated MidPointZ to be " << MidZ << "\n";
    }

Vector PosZero(MidX,MidY,MidZ), // Position in Image, that will be at the corner of the projection
       PosMid(OUTPUTSizeX/2.0/FacX,OUTPUTSizeY/2.0/FacY,OUTPUTSizeZ/2.0/FacZ),    // Middle of Projection
       PixVecX(1.0/FacX,0.0,0.0),PixVecY(0.0,1.0/FacY,0.0),PixVecZ(0.0,0.0,1.0/FacZ),
       PosZeroBwd(OUTPUTSizeX/2.0,OUTPUTSizeY/2.0,OUTPUTSizeZ/2.0), // Position in Image, that will be at the corner of the projection
       PosMidBwd(FacX*MidX,FacY*MidY,FacZ*MidZ),    // Middle of Projection
       PixVecXBwd(FacX,0.0,0.0),PixVecYBwd(0.0,FacY,0.0),PixVecZBwd(0.0,0.0,FacZ);

 PosZero.Sub(& PosMid);
 PosZeroBwd.Sub(& PosMidBwd);

 ofstream to(OUTPUTFileName.c_str());

for (Elem=0;Elem<Elements;Elem++)
  {
    cout << "Loading Element " << Elem << "\n";
    Elements=InputImg.Array.DLoad(kflag,INPUTFileName.c_str(),"Float",& INPUTSizeX,& INPUTSizeY, & INPUTSizeZ,Elem);
  if (VFileName!="")
    {
      Vectors.DLoad(kflag,VFileName.c_str(),"Float",& VECSizeX,& VECSizeY, & VECSizeZ,Elem);
      MidX=Vectors.Value(1,0,0);
      MidY=Vectors.Value(2,0,0);
      MidZ=Vectors.Value(3,0,0);
      smx=smy=smz=1;  // Do not calculate again !
    }

    SmallPrj.SetPos(& PosZero);
    SmallPrj.SetDir(& Dir);
    SmallPrj.TakeNewChangeVec(0,&PixVecX);
    SmallPrj.TakeNewChangeVec(1,&PixVecY);
    SmallPrj.TakeNewChangeVec(2,&PixVecZ);
    // SmallPrj.SetClip(& Vector(FacX/2.0+0.5,FacY/2.0+0.5,FacZ/2.0+0.5));
    SmallPrj.SetZoom(FacX,FacY,FacZ);

    InputImg.FirstPrj.SetPos(& PosZeroBwd);
    InputImg.FirstPrj.SetDir(& Dir);
    InputImg.FirstPrj.TakeNewChangeVec(0,&PixVecXBwd);
    InputImg.FirstPrj.TakeNewChangeVec(1,&PixVecYBwd);
    InputImg.FirstPrj.TakeNewChangeVec(2,&PixVecZBwd);
    double Length=FacY/2.0,Width=FacX/2.0, Depth=FacZ/2.0;  // the different order is due to the Dir-Vector
    if (Length < 0.5) Length=0.5;
    if (Width < 0.5) Width=0.5;
    if (Depth < 0.5) Depth=0.5;
    InputImg.FirstPrj.SetClip(Length,Width,Depth);  // Length, Width, Depth

    InputImg.FirstPrj.SetZoom(FacX,FacY,FacZ);

    pOutputPrj= new TProjection2(& SmallPrj);
    pOutputPrj->Array.Resize(OUTPUTSizeX,OUTPUTSizeY,OUTPUTSizeZ);
	
    pOutputPrj->Array.Set(0);
    InputImg.TakeNewPadVal(Fill);
    pOutputPrj->TakeNewPadVal(Fill);

    if (method == 3) // fourier method. For speed enhancement some a preselection in the source may be very useful !
      {  // at the moment zoom is only sizes
	/* int ZoomedX=(int) (INPUTSizeX*FacX +0.5);
	   int ZoomedY=(int) (INPUTSizeY*FacY +0.5);
	   int ZoomedZ=(int) (INPUTSizeZ*FacZ +0.5); */
	
	int EINPUTSizeX = INPUTSizeX;
	int EINPUTSizeY = INPUTSizeY;
	int EINPUTSizeZ = INPUTSizeZ;
	// ehemals: int EINPUTSizeZ = MKEven(INPUTSizeZ+1);
	
	// int ZoomedX=MKEven((int) (EINPUTSizeX*FacX+0.5));  // has to be even !!
	int ZoomedX=(int) (EINPUTSizeX*FacX+0.5);
	if (ZoomedX <= 0) ZoomedX=1;
	int ZoomedY=(int) (EINPUTSizeY*FacY+0.5);
	if (ZoomedY <= 0) ZoomedY=1;
	int ZoomedZ=(int) (EINPUTSizeZ*FacZ+0.5);
	if (ZoomedZ <= 0) ZoomedZ=1;
	
	double RFacX=double(ZoomedX) / double(EINPUTSizeX);
	double RFacY=double(ZoomedY) / double(EINPUTSizeY);
	double RFacZ=double(ZoomedZ) / double(EINPUTSizeZ);
	
        double ShiftX,ShiftY,ShiftZ;  // are applied to intermediate image
	
	cout << "FFT Magnifications are : " << RFacX << ", " << RFacY << ", " << RFacZ << "\n";
	cout << "MidX : " << MidX << ", MidY : " << MidY << ", MidZ : " << MidZ << "\n";
	cout << "Zoomed Sizes  : " << ZoomedX << ", " << ZoomedY << ", " << ZoomedZ << "\n";
	
	TInterArray Intermediate,Zoomed;
	Intermediate.Resize(EINPUTSizeX,EINPUTSizeY,EINPUTSizeZ); // will also set contens to zero
	Zoomed.Resize(ZoomedX,ZoomedY,ZoomedZ); // will also set contens to zero
	
	Copy(& InputImg.Array,& Intermediate); // template function

	Intermediate.FFTfwd();
	Intermediate.scramble();

	ShiftX = MidX-EINPUTSizeX/2.0;
	if (even(OUTPUTSizeX) && uneven(ZoomedX))
	  ShiftX = MidX-EINPUTSizeX/2.0-0.5/RFacX;
	if(uneven(OUTPUTSizeX) && even(ZoomedX))
	  if (RFacX == 1.0)
	    ShiftX = MidX-EINPUTSizeX/2.0-0.5/RFacX;
	  else
	    ShiftX = MidX-EINPUTSizeX/2.0+0.5/RFacX;
	if(uneven(OUTPUTSizeX) && uneven(ZoomedX))
	  ShiftX = MidX-EINPUTSizeX/2.0+0.0/RFacX;

	ShiftY = MidY-EINPUTSizeY/2.0;
	if (even(OUTPUTSizeY) && uneven(ZoomedY))
	  ShiftY = MidY-EINPUTSizeY/2.0-0.5/RFacY;
	if(uneven(OUTPUTSizeY) && even(ZoomedY))
	  if (RFacY == 1.0)
	    ShiftY = MidY-EINPUTSizeY/2.0-0.5/RFacY;
	  else
	    ShiftY = MidY-EINPUTSizeY/2.0+0.5/RFacY;
	if(uneven(OUTPUTSizeY) && uneven(ZoomedY))
	  ShiftY = MidY-EINPUTSizeY/2.0+0.0/RFacY;

	ShiftZ = MidZ-EINPUTSizeZ/2.0;
	if (even(OUTPUTSizeZ) && uneven(ZoomedZ))
	  ShiftZ = MidZ-EINPUTSizeZ/2.0-0.5/RFacZ;
	if (uneven(OUTPUTSizeZ) && even(ZoomedZ))
	  if (RFacZ == 1.0)
	    ShiftZ = MidZ-EINPUTSizeZ/2.0-0.5/RFacZ;
	  else
	    ShiftZ = MidZ-EINPUTSizeZ/2.0+0.5/RFacZ;
	if (uneven(OUTPUTSizeZ) && uneven(ZoomedZ))
	  ShiftZ = MidZ-EINPUTSizeZ/2.0+0.0/RFacZ;

	Intermediate.FFTshift(ShiftX, ShiftY, ShiftZ);
	
	Copy(& Intermediate,& Zoomed); // centered

	if ((ZoomedX > EINPUTSizeX) && even(EINPUTSizeX))
	  Zoomed.CopyFFTPlane(0,EINPUTSizeX/2);
	  // duplicates the highest frequency, because it doesnt exist in the other direction
	if ((ZoomedY > EINPUTSizeY) && even(EINPUTSizeY))
	  Zoomed.CopyFFTPlane(1,EINPUTSizeY/2);
	  // duplicates the highest frequency, because it doesnt exist in the other direction
	if ((ZoomedZ > EINPUTSizeZ) && even(EINPUTSizeZ))
	  Zoomed.CopyFFTPlane(2,EINPUTSizeZ/2);
	  // duplicates the highest frequency, because it doesnt exist in the other direction

	// Zoomed.KSave("/tmp/xxx1");
	Zoomed.scramble(1);  // unscramble
	// Zoomed.KSave("/tmp/xxx2");
	Zoomed.FFTbwd();
	// Zoomed.KSave("/tmp/xxx3");
	CopyCtoR(& Zoomed,& pOutputPrj->Array);  // centered at 0
      }
    else
      {
	if (method == 1)
	  {
	    if (1.0/FacX != int(1.0/FacX) || 1.0/FacY != int(1.0/FacY) || 1.0/FacZ != int(1.0/FacZ))
	      cerr << "WARNING ! Non-integer fractional Scaling factors result into severe artifacts !\n";

	    if (FacX > 1.0 || FacY > 1.0 || FacZ > 1.0)
	      cerr << "WARNING ! Scaling factors > 1.0 result into missing lines in the image ! !\n";
	      
	    InputImg.FirstPrj.SetClip(0.0,0.0,0.0);  // Length, Width, Depth
	    InputImg.FirstPrj.SetOnlyNearest(1);
	  }

	pOutputPrj->Array.Set(Fill);
	InputImg.Array.Add(-Fill);
	if (method == 4 ) // backwards
	  {
	    pOutputPrj->FirstPrj.SetClip(0.5,0.5,0.5);  // Length, Width, Depth
	    pOutputPrj->FirstPrj.SetZoom(1.0,1.0,1.0);

	    pOutputPrj->ProjectImage(& InputImg);
	    // InputImg.ProjectYourselfIntoImage(pOutputPrj);
	  }
	else // method == 2
	  InputImg.ProjectYourselfIntoImage(pOutputPrj);
      }
    
    if (Elem == 0)
      pOutputPrj->Array.DHeader(kflag,to,Elements);
    
    pOutputPrj->Array.Write(& to);
  }
 if (to) 
  to.close(); 
}
