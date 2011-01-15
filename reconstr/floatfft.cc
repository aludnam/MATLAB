// floatfft : computes the fft of arbitrary sized float images

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
typedef complex<float> ArrayBCType;

typedef TArray3d<ArrayBType> TImgArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TFFTArray<ArrayBCType>  TImgCArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TFFTArray<ArrayBCType>  TSliceArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

TImgArray InputImg;
TImgCArray OutputImg;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-iX X ...] [-i inputimage] [-o outputfile] [-us]\n" << flush;
  cerr <<  "[-iX -iY -iZ]  : sizes of images\n" << flush;
  cerr <<  "[-us]:      unscramble\n" << flush;
  cerr <<  "[-k] :      use KDF format for in- and output \n" << flush;
  cerr <<  "[-complex] :      do full complex In/Output \n" << flush;
  cerr <<  "[-inverse] :    do inverse fft \n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

int SizeX=32;  // These is the standart size, if raw data is used
int SizeY=32;  
int SizeZ=32; 
int Elements=1;

int kflag=0,unscramble=0,Elem, inverse=0, compx=0,ndepth=0,onlydepth=0,center0=0;

string IFileName, OFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-complex",parg)) {compx=1;continue;}  // treat both input/output as complex
   if (readArg("-nodepth",parg)) {ndepth=1;continue;}  // perform FFT slicewise
   if (readArg("-onlydepth",parg)) {onlydepth=1;continue;}  // perform FFT slicewise
   if (readArg("-us",parg)) {unscramble=1;continue;}  // Fourier-space data is scrambled
   if (readArg("-inverse",parg)) {inverse=1;continue;}  // performe inverse FFT
   if (readArg("-center0",parg)) {center0=1;continue;}  // real space data should have the zero in the center int(Size /2)
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-o",  OFileName, parg)) continue;
   if (readArg("-iX",& SizeX,parg)) continue; // input size X
   if (readArg("-iY",& SizeY,parg)) continue;
   if (readArg("-iZ",& SizeZ,parg)) continue;
    usage(argv[0]);
  }



 Elements=1;

 ofstream to(OFileName.c_str());

 for (Elem=0;Elem<Elements;Elem++)
   {
     cout << "Loading Element " << Elem << "\n";
     if (inverse)
       {
	 Elements=OutputImg.DLoad(kflag,IFileName.c_str(),"Complex",& SizeX,& SizeY, & SizeZ,Elem);
	 if (Elem==0)
	   InputImg.Resize(SizeX,SizeY,SizeZ);
	 
	 if (onlydepth)
	   {
	     TImgCArray Slice;
	     Slice.Resize(SizeZ,1,1);
	     int x,y,z;

	     for (y=0;y < SizeY;y++)
	       for (x=0;x < SizeX;x++)
		 {
		   for (z=0;z < SizeZ;z++)
		       Slice.SetValue(z,0,0,OutputImg.Value(x,y,z));  // not very fast ...
		   if (unscramble)
		     Slice.scramble(1);
		   Slice.FFTbwd();
		   if (center0)
		     Slice.scramble(1);
		   for (z=0;z < SizeZ;z++)
		     OutputImg.SetValue(x,y,z,Slice.Value(z,0,0));
		 }
	   }
	 else
	 if (! ndepth) // full 3D fft
	   {
	     if (unscramble)
	       OutputImg.scramble(1);
	     OutputImg.FFTbwd();
	     if (center0)
	       OutputImg.scramble(1);
	   }
	 else  // slicewise
	   {
	     TImgCArray Slice;
	     Slice.Resize(SizeX,SizeY,1);
	     int z;

	     for (z=0;z < SizeZ;z++)
       {
		 Slice.Copy(OutputImg.Pointer(0,0,z));  // not very fast ...
		 if (unscramble)
		   Slice.scramble(1);
		 Slice.FFTbwd();
		 if (center0)
		   Slice.scramble(1);
		 Slice.CopyTo(OutputImg.Pointer(0,0,z));
       }
	   }

	 if (! compx)
	   {
	     CopyCtoR(& OutputImg,& InputImg);
	 
	     if (Elem == 0)
	       InputImg.DHeader(kflag,to,Elements);
	 
	     InputImg.Write(& to);
	   }
	 else
	   {
	     if (Elem == 0)
	       OutputImg.DHeader(kflag,to,Elements);
	 
	     OutputImg.Write(& to);
	   }
       }
     else  //forward fft
       {
	 if (! compx)
	   {
	     Elements=InputImg.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY, & SizeZ,Elem);
	     if (Elem==0)
	       OutputImg.Resize(SizeX,SizeY,SizeZ);
	 
	     Copy(& InputImg,& OutputImg);
	   }
	 else
	   {
	     Elements=OutputImg.DLoad(kflag,IFileName.c_str(),"Complex",& SizeX,& SizeY, & SizeZ,Elem);
	   }

	 if (onlydepth)
	   {
	     TImgCArray Slice;
	     Slice.Resize(SizeZ,1,1);
	     int x,y,z;

	     for (y=0;y < SizeY;y++)
	       for (x=0;x < SizeX;x++)
		 {
		   for (z=0;z < SizeZ;z++)
		       Slice.SetValue(z,0,0,OutputImg.Value(x,y,z));  // not very fast ...
		   if (center0)
		     Slice.scramble();
		   Slice.FFTfwd();
		   if (unscramble)
		     Slice.scramble();
		   for (z=0;z < SizeZ;z++)
		     OutputImg.SetValue(x,y,z,Slice.Value(z,0,0));
		 }
	   }
	 else
	 if (! ndepth)
	   {
	     if (center0)
	       OutputImg.scramble();
	     OutputImg.FFTfwd();
	     if (unscramble)
	       OutputImg.scramble();
	   }
	 else
	   {
	     TImgCArray Slice;
	     Slice.Resize(SizeX,SizeY,1);
	     int z;

	     for (z=0;z < SizeZ;z++)
	       {
		 Slice.Copy(OutputImg.Pointer(0,0,z));  // not very fast ...
		 if (center0)
		   Slice.scramble();
		 Slice.FFTfwd();
		 if (unscramble)
		   Slice.scramble();
		 Slice.CopyTo(OutputImg.Pointer(0,0,z));
	       }
	   }

	 if (Elem == 0)
	   OutputImg.DHeader(kflag,to,Elements);
	 
	 OutputImg.Write(& to);
       }
   }

 if (to) 
   to.close(); 
}
