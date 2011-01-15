// pamgenpattern.cc : shifts a pattern pointwise to generate a series of scanning patterns

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
#include "parseargs.h"
#include "rawarray.h"
#include "fftarray.h"


typedef float ArrayBType;

typedef TArray3d<ArrayBType> TImgArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

TImgArray InputImg,OutputImg;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-iX X ...] [-i inputimage] [-o outputfile] [-stepX num] [-stepY num] [-numX num] [-numY] \n" << flush;
  cerr <<  "[-iX -iY -iZ]  : sizes of images\n" << flush;
  cerr <<  "[-k] :      use KDF format for in- and output \n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 
  
  int SizeX=32,SizeY=32, SizeZ=32, NumX=5, NumY=5, StepX=1, StepY=1;
  
  int kflag=0;
  
  string IFileName, OFileName;
  
  char ** parg= & argv[1];
  argc=0;  // to prevent warning
  
  while (* parg)
    {
      if (readArg("-k",parg)) {kflag=1;continue;}
      if (readArg("-i", IFileName, parg)) continue;
      if (readArg("-o", OFileName, parg)) continue;
      if (readArg("-numX",& NumX,parg)) continue; // number of scans in X direction
      if (readArg("-numY",& NumY,parg)) continue; // number of scans in Y direction
      if (readArg("-stepX",& StepX,parg)) continue; // Steps to move
      if (readArg("-stepY",& StepY,parg)) continue; // Steps to move
      if (readArg("-iX",& SizeX,parg)) continue; // input size X
      if (readArg("-iY",& SizeY,parg)) continue;
      if (readArg("-iZ",& SizeZ,parg)) continue;
      usage(argv[0]);
    }
  

  ofstream to(OFileName.c_str());
  if (! to) cerr << "Fatal Error : Could not open outputfile " << OFileName << "\n", exit(-1);
 
  InputImg.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY, & SizeZ); // Load Image
  if (SizeZ != 1)
    cerr << "Fatal error : planes-file should have size 1 in depth direction. Store different planes as elements!\n",exit(-1);
  OutputImg.Resize(SizeX,SizeY,SizeZ);
  OutputImg.DHeader(kflag,to,NumX*NumY);  
 
  for (int y=0;y< NumY*StepY;y+=StepY)
    for (int x=0;x< NumX*StepX;x+=StepX)
      {
	CyclicCopy(InputImg,OutputImg,x,y,0,project1st<ArrayBType,ArrayBType>());  // copy and shift
	OutputImg.Write(& to);
      }
  if (to)
    to.close();
   
}
