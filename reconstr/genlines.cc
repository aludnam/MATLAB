// genlines : draws lines into images. Lines are read from a vector file

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
#include "tarrayclipf.h"    // in this case really needed
#include "parseargs.h"
#include "clipperator3d.h"

Vector Dir(0.0,1.0,0.0);

typedef float ArrayBType;
typedef Dyn3dArray(ArrayBType,Clipperator3DClip)  TImgArray;  // ClipperatorNeighbClip can generated segmentation faults, when the positions of the line are outside the array
typedef Dyn3dArray(ArrayBType,Clipperator3DClip)  TPrjArray;   // Clipping is needed, Clipperator3DClip is also possible
typedef ProjectorLpZ<TImgArray>                       TProjector;
typedef ProjectorLpZ<TPrjArray>                       TProjector2;
typedef Image<TImgArray,TProjector>                   TImage;       // TProjector2 ?
typedef Projection<TPrjArray,TImage,TProjector>       TProjection;       // TProjector2 ?   SetZoom must be activated !!

static const ArrayBType FillValue=0.0;        

static TProjector    SmallPrj(3,0.5,0.5,0.5);
static TProjector    SmallPrjBwd(3,0,0,0);
static TProjection   InputImg(& SmallPrjBwd);
static TProjection * pOutputPrj=0;
TArray3d<float>       InputVecs, OutputImg;


void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-i inputfile] [-o outputfile]\n" << flush;
  cerr <<  "[-iX -iY -iZ]  : sizes of source image\n" << flush;
  cerr <<  "[-sX -sY -sZ]  : start point in 3D\n" << flush;
  cerr <<  "[-eX -eY -eZ]  : end point in 3D\n" << flush;
  cerr <<  "[-k] :           use KDF format for in- and output \n" << flush;
  exit(-1);
}


int main(int argc, char *argv[])
{ 

int INPUTSizeX=32;  // These is the standart size, if raw data is used
int INPUTSizeY=32;  
int INPUTSizeZ=32; 

int kflag=0,Elem,Elements;

string IFileName, OFileName, IVFileName;

double value=255.0;
bool add=false;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-add",parg)) {add=true;continue;}
   if (readArg("-i", IFileName, parg)) continue;     // input image to which to draw to
   if (readArg("-iv", IVFileName, parg)) continue;   // vector file with the coordinates
   if (readArg("-o", OFileName, parg)) continue;    // Image with the lines drawn
   if (readArg("-iX",& INPUTSizeX,parg)) continue; // input size X
   if (readArg("-iY",& INPUTSizeY,parg)) continue;
   if (readArg("-iZ",& INPUTSizeZ,parg)) continue;
   if (readArg("-v",& value, parg)) continue;

   usage(argv[0]);
  }


Vector NullVec(0,0,0);

 ofstream to(OFileName.c_str());

 Elements=1;
 int VecSizeX,VecSizeY,VecSizeZ,VElements=1,vy;

 VElements=InputVecs.DLoad(kflag,IVFileName.c_str(),"Float",& VecSizeX,& VecSizeY, & VecSizeZ,0);
 if (VElements != 1)
 {
   cerr << "ERROR: Vector file contains more than one Element! Stack the time in depth direction in the vector file!\n";
   exit(-1);
 }

float DX=0,PX=0,PXO=0,DY=0,PY=0,PYO=0,DZ=0,PZ=0,PZO=0,MSize;

for (Elem=0;Elem<Elements;Elem++)
  {
    cout << "Loading Element " << Elem << "\n";
    Elements=OutputImg.DLoad(kflag,IFileName.c_str(),"Float",& INPUTSizeX,& INPUTSizeY, & INPUTSizeZ,Elem);

    if (Elem == 0)
    {
      InputImg.Array.Resize(& OutputImg);
      OutputImg.DHeader(kflag,to,Elements);
    }

    if (!add)
     InputImg.Array.Clear();
     
    if (Elem >= 1)
    for (vy=0;vy<VecSizeY;vy++)
    {
    /*if (InputVecs(0,vy,Elem) <= 0)   // This spot is valid
    {
      if (InputVecs(0,vy,Elem-1) > 0)  // It was also valid in the previous time step
        {
           InputVecs(0,vy,Elem) = InputVecs(0,vy,Elem-1);
           InputVecs(1,vy,Elem) = InputVecs(1,vy,Elem-1);
           if (VecSizeX > 2)
             InputVecs(2,vy,Elem) = InputVecs(2,vy,Elem-1);
           if (VecSizeX > 3)
             InputVecs(3,vy,Elem) = InputVecs(3,vy,Elem-1);
        }
     }
     else*/
    if (InputVecs(0,vy,Elem) > 0)   // This spot is valid
      if (InputVecs(0,vy,Elem-1) > 0)  // It was also valid in the previous time step
        {
          PX = InputVecs(1,vy,Elem);
          PXO = InputVecs(1,vy,Elem-1);
          DX = PX-PXO;
          if (VecSizeX > 2)
          {
            PY = InputVecs(2,vy,Elem);
            PYO = InputVecs(2,vy,Elem-1);
            DY = PY-PYO;
          }
          if (VecSizeX > 3)
          {
            PZ = InputVecs(3,vy,Elem);
            PZO = InputVecs(3,vy,Elem-1);
            DZ = PZ-PZO;
          }
          MSize = sqrt(DX*DX+DY*DY+DZ*DZ);
          Vector PixVecX(DX/MSize,DY/MSize,DZ/MSize);
          Vector PosStart(PXO,PYO,PZO), PosEnd(PX,PY,PZ),DistVec(DX,DY,DZ);

          cout << "From " << PXO << ", " << PYO << ", " << PZO << " to " << PX <<", " << PY << ", " << PZ << "\n";
          SmallPrj.SetPos(& PosStart);   // When projecting forward ...
          SmallPrj.SetDir(& Dir);
          SmallPrj.TakeNewChangeVec(0,&PixVecX);
          SmallPrj.TakeNewChangeVec(1,&NullVec);
          SmallPrj.TakeNewChangeVec(2,&NullVec);
          SmallPrj.SetClip(1.0,1.0,1.0);  // Length, Width, Depth
          SmallPrj.SetZoom(1.0,1.0,1.0);  // ????

          if (pOutputPrj) delete(pOutputPrj);
          pOutputPrj= new TProjection(& SmallPrj);
          pOutputPrj->Array.Resize(int(MSize+1.5),1,1);  // project pixels

          pOutputPrj->Array.Set(value);         // color

          pOutputPrj->ProjectYourselfIntoImage(&InputImg);
        }
    }
    // if (add)
     OutputImg.Add(&InputImg.Array);
    /* else
     OutputImg.Copy(&InputImg.Array);*/

   OutputImg.Write(& to);
 }
 
 if (to) 
  to.close(); 
}
