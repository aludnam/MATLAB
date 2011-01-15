
/*   This file is part of a software package written by 
     Rainer Heintzmann
     Max Planck Inst. for Biophysical Chemistry, Am Fassberg 11, 37077 Goettingen, Germany
     Tel.: ++49 (0) 551 201 1029, e-mail: rheintz@gwdg.de  or rainer@heintzmann.de
     No garantee, whatsoever is given about functionallaty and  savety.
     No warranty is taken for any kind of damage it may cause.
     No support for it is provided !

     THIS IS NOT FREE SOFTWARE. All rights are reserved, and the software is only
     given to a limited number of persons for evaluation purpose !
     Please do not modify and/or redistribute this file or any other files, libraries
     object files or executables of this software package !
*/

// delvecs.cc : deletes entries in a file containing vectors
#include <iostream>
#include <math.h>
#include <string>
#include "rawarray.h"
#include "parseargs.h"
#include "projector.h"
#include "projectorpos.h"

#include "vecprojection.h"
#include "comperator.h"
#include "poidev.h"
#include "clipperator3d.h"


typedef float ArrayBType;
typedef TArray3d<ArrayBType>  TArray;  // Spots array 

static TArray   Array1,Array2;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] -i1 image1 -i2 image2 -o outputfile -d distancesFile\n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

static int kflag=0,x,y,greater=0;
static int SizeX=0,SizeY=0,SizeZ=0,elements=0,valid=0;
double thresh=0;
int column=0;

string OFileName, IFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-greater",parg)) {greater=1;continue;}
   if (readArg("-thresh", &thresh, parg)) continue;
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-i",  IFileName, parg)) continue;
   if (readArg("-column", & column, parg)) continue;
    usage(argv[0]);
  }

  if (OFileName=="" || IFileName=="") 
    usage(argv[0]);

  elements=Array1.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY,& SizeZ,0);
  if (SizeZ != 1)
    cerr << "Error : Datastack is 3D ! use Elements for multiple processing!\n",exit(-1);

  Array2.Resize(SizeX,SizeY,1);

ofstream * oi=0;

 if (OFileName!="")
   {
     oi=new ofstream(OFileName.c_str());
     if (! oi)
       cerr << "Error ! Couldnt open " << OFileName << " for writing !! \n",exit(-1);
   }

if (elements > 1)
	{
	  cerr << "WARNING! Datafile has " << elements << " elements. Datasize will not be changed !\n";
	}

 for (int  Elem=0;Elem<elements;Elem++)
  {
      if (Elem > 0)
        elements=Array1.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY,& SizeZ,Elem);

  if (column >= SizeX)
    cerr << "Error ! Datafile has only " << SizeX << " columns. Tried to select the " << column <<"s column !\n",exit(-1);
  

  for (y=0;y<SizeY;y++)
      if (greater)
	{
	if (Array1.Value(column,y,0) > thresh) valid++;
	}
      else
	{
	  if (Array1.Value(column,y,0) < thresh) valid++;
	}


    if (elements <= 1)
    {
             Array2.Resize(SizeX,valid);
    }
    if (Elem == 0)
          Array2.DHeader(kflag,*oi,elements);  // determines its sizes on its own

  valid=0;
  for (y=0;y<SizeY;y++)
      if (greater)
	{
	if (Array1.Value(column,y,0) > thresh) 
	  {
	  for (x=0;x<SizeX;x++)
	    Array2.SetValue(x,valid,0,Array1.Value(x,y));
	  valid ++;
	  }
	}
      else
	{
	  if (Array1.Value(column,y,0) < thresh)
	    {
	    for (x=0;x<SizeX;x++)
	      Array2.SetValue(x,valid,0,Array1.Value(x,y)) ;
	    valid ++;
	    }
	}
    if (elements > 1)
      for (y=valid;y<SizeY;y++)
	    for (x=0;x<SizeX;x++)
	      Array2.SetValue(x,valid,0,0.0) ;

  cout << "Deleted " << SizeY - valid << " vectors out of " << SizeY << "\n";

   Array2.Write(oi);
  }

if (oi)  oi->close();

  // Array2.DSave(kflag,OFileName);
}
