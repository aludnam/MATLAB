// genwave.cc : Generates a complex spherical wave

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
#include "rawarray.h"
#include "parseargs.h"
#include <complex>
#include "clipperator3d.h"

typedef complex<float> ArrayBType;
typedef TArray3d<ArrayBType>  TImgArray;  

static TImgArray     Img;

using namespace std;

complex<float> CalcRefrAmpl(float sinBeta,float n1,float n2,float A,float kr)
{
  float cosAlpha,cosBeta,Alpha,Beta;

  if ((sinBeta*n2/n1) > 1.0)  // total reflection
    return ArrayBType(0,0);

  Beta=asin(sinBeta);
  cosBeta=cos(Beta);
  Alpha=asin(sinBeta*n2/n1);
  cosAlpha=cos(Alpha);

  // return transmission coefficient
  return ArrayBType(2.0/(cosAlpha/cosBeta+n1/n2)) *  // assuming parallel Polarization everywhere : E||
    exp(ArrayBType(0,-kr*(n2*cos(Alpha-Beta)-n1)*A/cosAlpha));      // Phase correction A: Distance of object to surface
}

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] -o outputfile [-sX SizeX] [-sY SizeY] [-sZ SizeZ]\n" << flush;
  exit(-1);
}

float sqr(float x) { return x*x;}

int main(int argc, char *argv[])
{ 

static int kflag=0,sflag=0,Aberation=0,SizeX=64,SizeY=64,SizeZ=64,x,y,z;
static float Wavelength=2.0,Freq=2.0,r,kr,Scale=1.0,dr,n1=1.0,n2=1.0,A=0.0,ScaleX=1.0,ScaleY=1.0,ScaleZ=1.0;
complex<float> Val(1.0,0.0);

string OFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-shell",parg)) {sflag=1;continue;}
   if (readArg("-wl", & Wavelength, parg)) continue; // 1.0/(Vakuum wavelength in Pixels)
   if (readArg("-s", & Scale, parg)) continue; 
   if (readArg("-sX", & SizeX, parg)) continue; 
   if (readArg("-sY", & SizeY, parg)) continue;
   if (readArg("-sZ", & SizeZ, parg)) continue;
   if (readArg("-scaleX", & ScaleX, parg)) continue; // scales the shell in that direction
   if (readArg("-scaleY", & ScaleY, parg)) continue;
   if (readArg("-scaleZ", & ScaleZ, parg)) continue;
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-n1",& n1,parg)) continue;    // Refraction index of first medium
   if (readArg("-ab",parg)) {Aberation=1;continue;}
   if (readArg("-abn2",& n2,parg)) continue;
   if (readArg("-abA",& A,parg)) continue;  // Distance Object->Surface
    usage(argv[0]);
  }

  if (OFileName=="") 
    usage(argv[0]);

  Img.Resize(SizeX,SizeY,SizeZ);

  Freq = n1*2*M_PI/(Wavelength);

  if (! sflag)
    {
      for (z=0;z<SizeZ;z++)
	for (y=0;y<SizeY;y++)
	  for (x=0;x<SizeX;x++)
	    {
	      r = Scale*sqrt(sqr(x-SizeX/2.0)+sqr(y-SizeY/2.0)+sqr(z-SizeZ/2.0));
	      kr = r*Freq;
	      if (r == 0) r=1.0*Scale;
	      // if (x-SizeX/2.0 > 0) 
	      Img.SetValue(x,y,z,exp(ArrayBType(0,-kr))/r);
	    }
      }
  else
    {  // draw a shell
      if (Aberation)
	Freq = n2*2*M_PI/(Wavelength); // reduce k-vector in outgoing wave

      for (z=0;z<SizeZ;z++)
	for (y=0;y<SizeY;y++)
	  for (x=0;x<SizeX;x++)
	    {
	      r = Scale*sqrt(sqr((x-SizeX/2.0)/ScaleX)+sqr((y-SizeY/2.0)/ScaleY)+sqr((z-SizeZ/2.0)/ScaleZ));

	      dr = abs(r-SizeX*Freq/2/M_PI);
	      if (dr < 1.0)
		{
		  Val=complex<float> (1.0,0.0);

		  if (Aberation)
		    {
		      Val = CalcRefrAmpl(sqrt(sqr((y-SizeY/2.0)/ScaleY)+sqr((z-SizeZ/2.0)/ScaleZ))/(r/Scale),n1,n2,A,Freq);
		      // Phase has to be calculated here !!
		    }

		  Img.SetValue(x,y,z,Val*ArrayBType((1.0-dr))*
                      exp(ArrayBType(0,M_PI*(x+y+z)))); // to shift to mid-image
		}
	    }
    }

    Img.DSave(kflag,OFileName.c_str());
}
