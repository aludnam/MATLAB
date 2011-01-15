// scanline : scans a line in 3D through a datastack using linear interpolation with a specified width

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
static TProjection * pOutputPrj;


void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-i inputfile] [-o outputfile] [-s samples per x-pixel-size] \n" << flush;
  cerr <<  "[-iX -iY -iZ]  : sizes of source image\n" << flush;
  cerr <<  "[-sX -sY -sZ]  : start point in 3D\n" << flush;
  cerr <<  "[-eX -eY -eZ]  : end point in 3D\n" << flush;
  cerr <<  "[-w scanwidth] : width of interpolation\n" << flush;
  cerr <<  "[-k] :           use KDF format for in- and output \n" << flush;
  exit(-1);
}


int main(int argc, char *argv[])
{ 

int OUTPUTSize=32;
int INPUTSizeX=32;  // These is the standart size, if raw data is used
int INPUTSizeY=32;  
int INPUTSizeZ=32; 

int kflag=0,Elem,Elements;

string IFileName, OFileName, AFileName, LFileName;

double SamplesPP=2.0,StartX=0.0,StartY=0.0,StartZ=0.0,EndX=0.0,EndY=0.0,EndZ=0.0,scanwidth=0.5,xscaling=1.0,LVal=255.0;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-ao", AFileName, parg)) continue;
   if (readArg("-lo", LFileName, parg)) continue;
   if (readArg("-lv", & LVal, parg)) continue;
   if (readArg("-iX",& INPUTSizeX,parg)) continue; // input size X
   if (readArg("-iY",& INPUTSizeY,parg)) continue;
   if (readArg("-iZ",& INPUTSizeZ,parg)) continue;
   if (readArg("-sX",& StartX,parg)) continue; // output size X
   if (readArg("-sY",& StartY,parg)) continue;
   if (readArg("-sZ",& StartZ,parg)) continue;

   if (readArg("-eX",& EndX, parg)) continue; // scale factors
   if (readArg("-eY",& EndY, parg)) continue; 
   if (readArg("-eZ",& EndZ, parg)) continue; 
   if (readArg("-s",& SamplesPP, parg)) continue; 
   if (readArg("-w",& scanwidth, parg)) continue; 
   if (readArg("-axs",& xscaling, parg)) continue; 

   usage(argv[0]);
  }


Vector PosStart(StartX,StartY,StartZ), // Position in Image, that will be at the corner of the projection
       PosEnd(EndX,EndY,EndZ), DistVec(EndX-StartX,EndY-StartY,EndZ-StartZ);

double dist = DistVec.Norm();
double Samples = dist*SamplesPP;

Vector NullVec(0,0,0),PixVecX((EndX-StartX)/Samples,(EndY-StartY)/Samples,(EndZ-StartZ)/Samples);
Vector PosZeroBwd((EndX-StartX)/Samples,(EndY-StartY)/Samples,(EndZ-StartZ)/Samples);

double sampledist = PixVecX.Norm();

 DistVec.Mul(1.0/dist);  // normalize

 ofstream to(OFileName.c_str());

 ofstream * ato=0,* to3d=0;

 if (AFileName != "") ato = new ofstream(AFileName.c_str());
 if (LFileName != "") to3d = new ofstream(LFileName.c_str());

 Elements=1;

for (Elem=0;Elem<Elements;Elem++)
  {
    cout << "Loading Element " << Elem << "\n";
    Elements=InputImg.Array.DLoad(kflag,IFileName.c_str(),"Float",& INPUTSizeX,& INPUTSizeY, & INPUTSizeZ,Elem);

    double Length=0.5,Width=scanwidth/2.0, Depth=scanwidth/2.0;  // the different order is due to the Dir-Vector
    if (Length < 0.5) Length=0.5;
    if (Width < 0.5) Width=0.5;
    if (Depth < 0.5) Depth=0.5;

    SmallPrj.SetPos(& PosStart);   // When projecting forward ...
    SmallPrj.SetDir(& Dir);
    SmallPrj.TakeNewChangeVec(0,&PixVecX);
    SmallPrj.TakeNewChangeVec(1,&NullVec);
    SmallPrj.TakeNewChangeVec(2,&NullVec);
    SmallPrj.SetClip(1.0,1.0,1.0);  // Length, Width, Depth
    SmallPrj.SetZoom(1.0,1.0,1.0);  // ????

    pOutputPrj= new TProjection(& SmallPrj);
    OUTPUTSize = int(dist/sampledist);
    pOutputPrj->Array.Resize(OUTPUTSize,1,1);

    pOutputPrj->Array.Set(0);

    // InputImg.ProjectYourselfIntoImage(pOutputPrj);
    pOutputPrj->ProjectImage(&InputImg);
    
    if (ato)
      {
	int x;
	(* ato) << "\n";
	for (x=0;x < OUTPUTSize;x++)
	  (* ato) << x*sampledist*xscaling << "     " << pOutputPrj->Array.Value(x,0,0) << "\n";
      }
    if (Elem == 0)
      {
	pOutputPrj->Array.DHeader(kflag,to,Elements);
	if (to3d)
	  InputImg.Array.DHeader(kflag,(* to3d),Elements);
      }

    pOutputPrj->Array.Write(& to);
    if (to3d)
      {
        pOutputPrj->Array.Set(LVal);
        pOutputPrj->ProjectYourselfIntoImage(&InputImg);
        InputImg.Array.Write(to3d);
      }
  }

 if (to) 
  to.close(); 
}
