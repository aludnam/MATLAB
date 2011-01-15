// genrecmat : generates a matrix to be used for saturated patterned illumination reconstruction
/*   This file is part of a software package written by 
     Rainer Heintzmann
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

// the delta peaks in the fft will be counted in the following way:
// example : 0 order and +- 1 order
//      0   1   2
//      3   4   5
//      6   7   8

// Thus the zero order is position 4
// The matrix will describe the generation of images

#include <iostream>
#include <string>
#include "rawarray.h"
#include "parseargs.h"

typedef complex<float> ArrayBType;
typedef TArray3d<ArrayBType> TImgArray;

static TImgArray Matrix;
const double Pi= 2.0* asin(1.0);


void usage(char * filename)
{
  // cout << "PI "  << Pi << "\n";
  cerr <<  "usage: " << filename << " [-k] [-sx PattternStepsX] [-sY PatternStepsY] [-nx orders along x] [-ny orders along y] [-o outputmatrix] \n" << flush;
  exit(-1);
}


int main(int argc, char *argv[])
{ 
bool diagscan=false,noZeroOrder=false;
int PatternStepsX=3, PatternStepsY=3, ordersx=1,ordersy=1, totalordersX=3, totalordersY=3;
int kflag=0,i,j,x,y;
float val,kx,ky;

string OFile("/tmp/matrix.raw");

string OUTPUTFileName;
string OFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-nx", & ordersx, parg)) continue; // orders to be used (+-)
   if (readArg("-ny", & ordersy, parg)) continue; // orders to be used (+-)
   if (readArg("-o", & OFileName, parg)) continue;
   if (readArg("-sx",& PatternStepsX, parg)) continue;  // steps assume a rectangular grid
   if (readArg("-sy",& PatternStepsY,parg)) continue;
   if (readArg("-diagscan", parg))  {diagscan=true; continue;}// the order vectors will be rotated by this angle 
   if (readArg("-nozeroorder", parg))  {noZeroOrder=true; continue;}// the order vectors will be rotated by this angle 
    usage(argv[0]);
  }

  if (OFileName=="") OUTPUTFileName=OFile;
    else OUTPUTFileName=OFileName;

  kx = 1.0 / PatternStepsX ;  
  ky = 1.0 / PatternStepsY ;  // frequency for a full repetition


  totalordersX=2*ordersx+1;
  totalordersY=2*ordersy+1;
  cout << "total orders in each direction is " << totalordersX*totalordersY << "\n";
  Matrix.Resize(totalordersX*totalordersY,PatternStepsX*PatternStepsY,1);

  cout << "Matrix  = exp (2 pi i * \n";
  for (j=-ordersy;j<=ordersy;j++)   // starting top left order
    for (i=-ordersx;i<=ordersx;i++)
      {
	for (y=0;y<PatternStepsY;y++)   // counting the sub-shifts of the illumination pattern
	  {
	  for (x=0;x<PatternStepsX;x++)
	    {
	      // if (diagscan)  // even in this case one cycle is completed for the first order, but also a cycle along order y
		// val = (x+y) * float(i) * kx + (y-x) * float(j) * ky;
	      // else
		// val = x * float(i) * kx + y * float(j) * ky;
	      if (! noZeroOrder)
		      val = x * float(i) * kx + y * float(j) * ky;
		else
		      val = 2 * x * float(i) * kx + 2 * y * float(j) * ky;
		      
	      if (diagscan)  // even in this case one cycle is completed for the first order, but also a cycle along order y
		val += y * float(i) * ky -x * float(j) * kx;
		    
	      // cout << " orderx " << i << " ordery " << j << " val " << val << "\n";
	      Matrix.SetValue((j+ordersy)*totalordersX+(i+ordersx),y*PatternStepsX+x,0, 
			      exp(ArrayBType(0.0f,1.0f) * float(2.0 * Pi * val)));
	      printf("%5.2f  ",val);
	      // cout << val << "\t";
	    }
      }
	  cout << val << "\n";
      }
  cout << ")\n";
  Matrix.DSave(kflag,OUTPUTFileName.c_str());
}
