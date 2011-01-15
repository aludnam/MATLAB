// subbackgr : Subtracts the background of an image (measured as mean value of lower x% of pixels)

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
#include "parseargs.h"
#include "khoros.h"

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-i inputstack] -d dimension \n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

int SizeX=32;  // These is the standart size, if raw data is used
int SizeY=32;  
int SizeZ=32; 
int Elements=1;
int Times=1;

int  dim=-1,ptype=0;

string IFileName;
char TypeString[100];

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-d",& dim,parg)) continue;
   if (readArg("-t",& ptype,parg)) continue;  // Print the datatype as a number
    usage(argv[0]);
  }


 if (IFileName=="")
   cerr << "Missing Input Filename\n",usage(argv[0]);

ifstream file(IFileName.c_str());

if (!file)
  cerr << "ERROR : could not open file " << IFileName << " for reading\n"<< flush,exit(-1);

  Times=ReadKhorosHeader(& file,TypeString,SizeX,SizeY,SizeZ,Elements);
  
  if (ptype!=0) cout << TypeString;

  if (dim==0)
    cout << SizeX;

  if (dim==1)
    cout << SizeY;

  if (dim==2)
    cout << SizeZ;

  if (dim==3)
    cout << Elements;

  if (dim==4)
    cout << Times;

file.close();
}
