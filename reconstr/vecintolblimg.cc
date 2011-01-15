
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

// wrapexpand.cc : expands one image to the sizes of a second one by wrapping

#include <iostream>
#include <string>
#include <math.h>
#include "rawarray.h"
// #include "tarrayclipf.h"
#include "parseargs.h"
#include "stlfixes.h"

typedef float ArrayBType;
typedef TArray3d<ArrayBType>  TImgArray;  // Spots array with clipping 

static TImgArray   Array1,Array2,Array3;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] -i1 image1 -i2 image2 -o outputfile\n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

static int kflag=0,el,x,y,z;
static int I1SizeX=0,I1SizeY=0,I1SizeZ=0,elements1=0;
static int I2SizeX=0,I2SizeY=0,I2SizeZ=0,elements2=0;

string OFileName,I1FileName,I2FileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-i1", I1FileName, parg)) continue; // label image
   if (readArg("-i2", I2FileName, parg)) continue; // vector image
    usage(argv[0]);
  }

  if (OFileName=="" || I1FileName=="" || I2FileName=="") 
    usage(argv[0]);

  elements1=Array1.DLoad(kflag,I1FileName.c_str(),"Float",& I1SizeX,& I1SizeY,& I1SizeZ,0);
  elements2=Array2.DLoad(kflag,I2FileName.c_str(),"Float",& I2SizeX,& I2SizeY,& I2SizeZ,0); // for sizes

  Array3.Resize(I1SizeX,I1SizeY,I2SizeX);

  ofstream of(OFileName.c_str());
  if (! of)
    cerr << "Error ! Couldnt open " << OFileName << " for writing !! \n",exit(-1);
  
  float val;
  int ival;
  for (el=0;el< elements2;el++)
    {

 	for (z=0;z<Array3.GetSize(2);z++)
 	  for (y=0;y<Array3.GetSize(1);y++)
 	    for (x=0;x<Array3.GetSize(0);x++)
		{
		val=Array1.Value(x,y,z%I1SizeZ);
	   	ival = int(val);
		if (val > 0)
			{
			if (ival < I2SizeY)
				{
				Array3.SetValue(x,y,z,Array2.Value(z,ival,z%I1SizeZ)); 
				}
			}
		}

      if (el == 0)
	Array3.DHeader(kflag,of,elements2);
      Array3.Write(& of);
    }

  of.close();
}
