
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

// cmpvecs.cc : compares vector images.
// result is a difference image of nearest vectors in second datastack to first datastack
// vectors are assumend to have the following format :
// Intensity,XPos,YPos,ZPos, something else
#include <iostream>
#include <math.h>
#include "tarrayclipf.h"
#include "parseargs.h"
#include "projector.h"
#include "projectorpos.h"

#include "vecprojection.h"
#include "comperator.h"
#include "poidev.h"
#include "clipperator3d.h"

using namespace std;

typedef float ArrayBType;
typedef Dyn2dArray(ArrayBType,ClipperatorNoClip)  TArray;  // Spots array with clipping 

static TArray   Array1;
static TArray   Array2;
static TArray   Distances;
static TArray   ArrayOut;
int IgnoreEq=0;

double Thedistance(int i,int j)
{
  double dx,dy,dz;
  dx= Array1.Value(1,i)- Array2.Value(1,j);  // same scaling is assumed
  dy= Array1.Value(2,i)- Array2.Value(2,j);
  dz= Array1.Value(3,i)- Array2.Value(3,j);
  return sqrt(dx*dx+dy*dy+dz*dz);
}

/// searches the nearest in Array2 to nr in Array1 using distance(nr,index2)
int SearchNearest(int nr,double & mindist)
{
  int minnr=0,j;
  double dist;
  mindist=1e4;

  for (j=0;j<Array2.GetSize(1);j++)
    if ((dist=Thedistance(nr,j)) < mindist)
      if (IgnoreEq)
	{
	  if(dist != 0.0)
	  mindist=dist,minnr=j;
	}
      else
	mindist=dist,minnr=j;

  return minnr;
}


void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] -i1 image1 -i2 image2 -o outputfile -d distancesFile\n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

static int kflag=0,el,spot,i,j,k;
static int I1SizeX=0,I1SizeY=0,I1SizeZ=0,elements1=0;
static int I2SizeX=0,I2SizeY=0,I2SizeZ=0,elements2=0;
double mindist;

string OFileName, DFileName, I1FileName, I2FileName;

bool noSubtract=false, useprev=false;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-up",parg)) {useprev=true;continue;}
   if (readArg("-nosubtract",parg)) {noSubtract=true;continue;}
   if (readArg("-ie",parg)) {IgnoreEq=1;continue;}
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-d", DFileName, parg)) continue;
   if (readArg("-i1", I1FileName, parg)) continue;
   if (readArg("-i2", I2FileName, parg)) continue;  // every vector of this array is compared to all vectors in array1 and the nearest is found.
    usage(argv[0]);
  }

  if (OFileName=="" || I1FileName=="" ||  I2FileName=="" || DFileName=="") 
    usage(argv[0]);

  elements1=Array1.DLoad(kflag,I1FileName.c_str(),"Float",& I1SizeX,& I1SizeY,& I1SizeZ,0);
  elements2=Array2.DLoad(kflag,I2FileName.c_str(),"Float",& I2SizeX,& I2SizeY,& I2SizeZ,0);
  if (I1SizeZ != 1 || I2SizeZ != 1)
    cerr << "Error : Datastack is 3D ! use Elements for multiple processing!\n",exit(-1);

  ArrayOut.Resize(I1SizeX,I1SizeY);
  Distances.Resize(1,I1SizeY);

  cout << "Reading vectors in Format : Intensity,XPos,YPos,ZPos, ...\n";

  ofstream of(OFileName.c_str());
  if (! of)
    cerr << "Error ! Couldnt open " << OFileName << " for writing !! \n",exit(-1);
  
  ofstream odf(DFileName.c_str());
  if (! odf)
    cerr << "Error ! Couldnt open " << DFileName << " for writing !! \n",exit(-1);
  
  int spots=0;

  int elems= elements1;
  if (elems < elements2) elems = elements2;   // Use the bigger one
  if (elements1 != elements2) cerr << "WARNING: Element size are not equal! Using a total of " << elems << " elements and wrapping the smaller size\n";
  
  for (el=0;el< elems;el++)
    {
      if (el > 0)
	{
                  if (useprev)
              {
                              cout << "Using previous result for comparison of element " << el << "\n";
                              Array1.Set(0);
                              Array1.Copy(& ArrayOut);
              }
                  else
         {
        	        Array1.DLoad(kflag,I1FileName.c_str(),"Float",& I1SizeX,& I1SizeY,& I1SizeZ,el % elements1);
          }
	  cout << "loading element " << el << "\n";
	  Array2.DLoad(kflag,I2FileName.c_str(),"Float",& I2SizeX,& I2SizeY,& I2SizeZ,el % elements2);
	}
      else // el == 0
	{
	  for (i=0;i<Array1.GetSize(1);i++)
	    if (Array1.Value(0,i) != 0)
		spots++;  // count spots
	  ArrayOut.Resize(I1SizeX,spots);
	  Distances.Resize(1,spots);
	  ArrayOut.DHeader(kflag,of,elems);
	  Distances.DHeader(kflag,odf,elems);
	}

      ArrayOut.Set(0);
      Distances.Set(0);
      spot=0;
      for (i=0;i<Array1.GetSize(1);i++)   // For every spot in array1 do
	if ((Array1.Value(0,i) != 0) && (spot < spots))
	  {
	    j = SearchNearest(i,mindist);   // search for nearest spot in array2
	    Distances.SetValue(0,spot,mindist);

	    for (k=0;k < Array1.GetSize(0);k++)
                      if (noSubtract)
        	        ArrayOut.SetValue(k,spot,Array2.Value(k,j));
                      else
	        ArrayOut.SetValue(k,spot,Array1.Value(k,i) - Array2.Value(k,j));
	    spot++;
	  }
       ArrayOut.Write(& of);
       Distances.Write(& odf);
    }

   if (of)
        of.close();
   if (odf)
      odf.close();

}
