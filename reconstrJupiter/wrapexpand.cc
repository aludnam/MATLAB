
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

static TImgArray   Array1;
static TImgArray   Array2;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] -i1 image1 -i2 image2 -o outputfile\n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

static int kflag=0,el,blend=0;
static int I1SizeX=0,I1SizeY=0,I1SizeZ=0,elements1=0;
static int I2SizeX=0,I2SizeY=0,I2SizeZ=0,elements2=0;
bool keepx=false,keepy=false,keepz=false,keepe=false;

int offe=0,offx=0,offy=0,offz=0;

string OFileName,I1FileName,I2FileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-keepx",parg)) {keepx=true;continue;}
   if (readArg("-keepy",parg)) {keepy=true;continue;}
   if (readArg("-keepz",parg)) {keepz=true;continue;}
   if (readArg("-keepe",parg)) {keepe=true;continue;}
   if (readArg("-blend",parg)) {blend=1;continue;}
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-i1", I1FileName, parg)) continue; // image to expand
   if (readArg("-i2", I2FileName, parg)) continue; 
   if (readArg("-oe", & offe, parg)) continue; 
   if (readArg("-ox", & offx, parg)) continue; 
   if (readArg("-oy", & offy, parg)) continue; 
   if (readArg("-oz", & offz, parg)) continue; 
    usage(argv[0]);
  }

  if (OFileName=="" || I1FileName=="" || I2FileName=="") 
    usage(argv[0]);

  elements1=Array1.DLoad(kflag,I1FileName.c_str(),"Float",& I1SizeX,& I1SizeY,& I1SizeZ,0);
  elements2=Array2.DLoad(kflag,I2FileName.c_str(),"Float",& I2SizeX,& I2SizeY,& I2SizeZ,0); // for sizes

  if (keepx) I2SizeX = I1SizeX;
  if (keepy) I2SizeY = I1SizeY;
  if (keepz) I2SizeZ = I1SizeZ;
  if (keepe) elements2 = elements1;
  
  Array2.Resize(I2SizeX,I2SizeY,I2SizeZ);
  
  ofstream of(OFileName.c_str());
  if (! of)
    cerr << "Error ! Couldnt open " << OFileName << " for writing !! \n",exit(-1);
  
  if (offe < 0) offe = elements1-((-offe)%elements1);
  for (el=0;el< elements2;el++)
    {
      int nelem=((el-offe) % elements1);
      if (nelem < 0) nelem = elements1 + nelem;
      cout << "processing element " << el << " being original element" << nelem<< "\n";
      Array1.DLoad(kflag,I1FileName.c_str(),"Float",& I1SizeX,& I1SizeY,& I1SizeZ,nelem);

      if (! blend)
	  CyclicCopy(Array1, Array2,offx,offy,offz,project1st<ArrayBType,ArrayBType>());
	  // CyclicCopy(Array1, Array2,offx,offy,offz);
      else BlendCopy(& Array1, & Array2);

      if (el == 0)
	Array2.DHeader(kflag,of,elements2);
      Array2.Write(& of);
    }

  of.close();

}
