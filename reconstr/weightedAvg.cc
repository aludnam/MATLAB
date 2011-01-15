// weightedAvg.cc : computes the weighted average of two files.
// The noise levels, the OTFs and the goal function has to be given

/*   This file is part of a software package written by 
     Rainer Heintzmann
     MPI bpc
     Am Fassberg 11
     37077 Goettingen
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
#include "rawarray.h"
#include "fftarray.h"
#include "parseargs.h"

typedef float ArrayBType;
typedef complex<float> ArrayBCType;

typedef TArray3d<ArrayBType>    TImgArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TFFTArray<ArrayBCType>  TCImgArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

static TImgArray  OTF1,OTF2,GoalOTF;
static TCImgArray Img1,Img2;


void usage(char * filename)
{
  cerr <<  "usage: " << filename << " -i1 inputfile -i2 input -otf1 Transfer1 -otf2 Transfer2 -sigma1 StdDev1 -sigma2 StdDev2 -o outputfile [-sX SizeX] [-sY SizeY] [-sZ SizeZ] [-integral] [-val value]\n" << flush;
  cerr << "computes the weigthed average of i1 and i2 \n";
  exit(-1);
}

int main(int argc, char *argv[])
{ 

static int SizeX=64,SizeY=64,SizeZ=64,method=0;
bool kflag=false, computeGoal=false;

string I1FileName,I2FileName,OTF1FileName,OTF2FileName,OFileName,GOALFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning
double alpha1,alpha2, Sigma1, Sigma2,DX=1.0,DY=1.0,DZ=1.0;
 ArrayBCType val;

 while (* parg)
  {
    if (readArg("-sX", & SizeX, parg)) continue; 
    if (readArg("-sY", & SizeY, parg)) continue;
    if (readArg("-sZ", & SizeZ, parg)) continue;
    if (readArg("-i1", I1FileName, parg)) continue;
    if (readArg("-i2", I2FileName, parg)) continue;
    if (readArg("-sigma1", & Sigma1, parg)) continue;
    if (readArg("-sigma2", & Sigma2, parg)) continue;
    if (readArg("-otf1", OTF1FileName, parg)) continue;
    if (readArg("-otf2", OTF2FileName, parg)) continue;
    if (readArg("-DX", & DX, parg)) continue;  // Borders in relative units
    if (readArg("-DY", & DY, parg)) continue;
    if (readArg("-DZ", & DZ, parg)) continue;
    if (readArg("-method", & method,parg)) continue;  // use linear decaying function
    if (readArg("-goal", GOALFileName, parg)) continue;  // if not given, an artificial OTF with boarder freq.s will be used
    if (readArg("-o", OFileName, parg)) continue;
    if (readArg("-k", parg)) {kflag=true;continue;}
    usage(argv[0]);
  }

  if (OFileName=="") 
    usage(argv[0]);

   ofstream to(OFileName.c_str());

   int Elements = 1;

   double Sigma1Sqr = Sigma1*Sigma1;
   double Sigma2Sqr = Sigma2*Sigma2;
   double Goal = 1.0,rrel=1.0;

   double DXSqr,DYSqr,DZSqr;

   if (GOALFileName=="") computeGoal = true;
   else computeGoal = false;

 for (int Elem=0;Elem<Elements;Elem++)
   {
     cout << "Loading Element " << Elem << "\n";
     Elements= Img1.DLoad(kflag,I1FileName.c_str(),"Complex",& SizeX, &SizeY, & SizeZ,Elem);
     Elements= Img2.DLoad(kflag,I2FileName.c_str(),"Complex",& SizeX, &SizeY, & SizeZ,Elem);
     Elements= OTF1.DLoad(kflag,OTF1FileName.c_str(),"Float",& SizeX, &SizeY, & SizeZ,Elem);
     Elements= OTF2.DLoad(kflag,OTF2FileName.c_str(),"Float",& SizeX, &SizeY, & SizeZ,Elem);
     if (! computeGoal)
       Elements= GoalOTF.DLoad(kflag,GOALFileName.c_str(),"Float",& SizeX, &SizeY, & SizeZ,Elem);

     int DimX=Img1.GetSize(0),DimY=Img1.GetSize(1),DimZ=Img1.GetSize(2);
     int MidX=DimX/2,MidY=DimY/2,MidZ=DimZ/2;

     DXSqr= DX * (DimX/2.0);
     DYSqr = DY * (DimY/2.0);
     DZSqr = DZ * (DimZ/2.0);
     DXSqr *= DXSqr;
     DYSqr *= DYSqr;
     DZSqr *= DZSqr;

     if (! Img1.SizesEqual(& Img2) ||
	 ! Img1.SizesEqual(& OTF1) ||
	 ! Img1.SizesEqual(& OTF2) ||
	 ((!computeGoal) && ! Img1.SizesEqual(& GoalOTF)) )
       cerr << "All sizes of all images have to be identical ! Bailing out \n", exit (-1);


     for(int z=0;z<DimZ;z++)
       for(int y=0;y<DimY;y++)
	 for(int x=0;x<DimX;x++)
	   {
	     if (computeGoal)
	       {
		 rrel = (x-MidX)*(x-MidX)/DXSqr+(y-MidY)*(y-MidY)/DYSqr+(z-MidZ)*(z-MidZ)/DZSqr;
		 if (rrel < 1.0)
		   {
		     switch (method)
		       {
		       case 1: Goal = 1.0 - sqrt(rrel); break;
		       case 2: Goal = 1.0 - rrel;break;
		       case 3: Goal = (1.0+cos(M_PI*sqrt(rrel)))/2.0;break;
		       default: Goal = (1.0+cos(M_PI*sqrt(rrel)))/2.0; Goal *= Goal;
		       }
		   }
		 else
		   Goal = 0.0;
	       }
	     else
	       {
		 Goal = GoalOTF.Value(x,y,z);
	       }
	     alpha1 = Goal / OTF1.Value(x,y,z);
	     alpha2 = Goal / OTF2.Value(x,y,z);
	     if (Goal != 0.0)
	       {
		 val = Img1.Value(x,y,z)/ArrayBCType(alpha1*Sigma1Sqr) + Img2.Value(x,y,z)/ArrayBCType(alpha2*Sigma2Sqr);
		 val /= 1.0 / (Sigma1Sqr*alpha1*alpha1) +  1.0 / (Sigma2Sqr*alpha2*alpha2) ;
		 Img1.SetValue(x,y,z, val); 
	       }
	     else
	       Img1.SetValue(x,y,z, 0); 
	   }


     if (Elem == 0)
       {
	 Img1.DHeader(kflag,to,Elements);
       }
	 
     Img1.Write(& to);
   }

}
