// ringtrans :  transformes either sum or squared energy into a ring

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
#include "vec.h"
// #include "tarrayclipf.h"
#include "rawarray.h"
#include "parseargs.h"

typedef float ArrayBType;

typedef TArray3d<ArrayBType> TArray;
typedef TArray3d<ArrayBType> TOutArray;
// typedef Dyn3dArray(ArrayBType,ClipperatorNoClip)  TArray;
// typedef Dyn1dArray(ArrayBType,ClipperatorNoClip)  TOutArray;

static TArray Img;
static TOutArray Result;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-iX X ...] [-bins Bins]\n" << flush;
  cerr <<  "[-k] :           use KDF format for in- and output \n" << flush;
  exit(-1);
}


int main(int argc, char *argv[])
{ 

int INPUTSizeX=32;  // These is the standart size, if raw data is used
int INPUTSizeY=32;  
int INPUTSizeZ=32; 

int kflag=0,bins=256;

char IFileName[100];
char OFileName[100];
IFileName[0]=0,OFileName[0]=0;


char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-i", & IFileName, parg)) continue;
   if (readArg("-o", & OFileName, parg)) continue;
   if (readArg("-iX",& INPUTSizeX,parg)) continue; // input size X
   if (readArg("-iY",& INPUTSizeY,parg)) continue;
   if (readArg("-iZ",& INPUTSizeZ,parg)) continue;
   if (readArg("-bins",& bins,parg)) continue;  // Fill outside with

    usage(argv[0]);
  }

 if (! IFileName[0]) cerr << "Error, missing input file. Input is required !\n",exit(-1);
 if (! OFileName[0]) cerr << "Error, missing output file name, required argument !\n",exit(-1);

 int Elem, Elements=1;

 ofstream to(OFileName);
 
 double xscale=1.0;

 for (Elem=0;Elem<Elements;Elem++)
   {
     cout << "Loading Element " << Elem << "\n";
     Elements=Img.DLoad(kflag,IFileName,"Float",&INPUTSizeX,&INPUTSizeY,&INPUTSizeZ,Elem);

     Result.Resize(bins,1,1);

     xscale=Img.NormedRingSum(& Result);
     cout << "1 histo-bin (along X) corresponds to " << xscale << " voxels\n";
 
     if (Elem == 0)
       Result.DHeader(kflag,to,Elements);
	 
     Result.Write(& to);
   }
}
