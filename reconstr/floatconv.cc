// floatconvolve : computes the fft of arbitrary sized float images

/*   This file is part of a software package written by 
     Rainer Heintzmann
     Institute of Applied Optics and Information Processing
     Albert Ueberle Strasse 3-5
     69120 Heidelberg
     Tel.: ++49 (0) 6221  549264
     Current Address : Max Planck Inst. for Biophysical Chemistry, Am Fassberg 11, 37077 Goettingen, Germany
     Tel.: ++49 (0) 551 201 1029, e-mail: rheintz@gwdg.de  or rainer@heintzmann.de
     No garantee, whatsoever is given about functionallity and  savety.
     No warranty is taken for any kind of damage it may cause.
     No support for it is provided !

     THIS IS NOT FREE SOFTWARE. All rights are reserved, and the software is only
     given to a limited number of persons for evaluation purpose !
     Please do not modify and/or redistribute this file or any other files, libraries
     object files or executables of this software package !
*/

#include <iostream>
#include "parseargs.h"
#include "rawarray.h"
#include "fftarray.h"

typedef float ArrayBType;

typedef TFFTArray<ArrayBType>  TImgCArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TArray3d<ArrayBType> TImgArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

TImgArray InputImg2;
TImgCArray FImg2,InputImg1;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-iX X ...] [-i inputimage] [-o outputfile] [-us]\n" << flush;
  cerr <<  "[-iX -iY -iZ]  : sizes of images\n" << flush;
  cerr <<  "[-cpsf]:      unscramble\n" << flush;
  cerr <<  "[-k] :      use KDF format for in- and output \n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

int SizeX=32;  // These is the standart size, if raw data is used
int SizeY=32;  
int SizeZ=32; 
int SizePX=32; 
int SizePY=32;  
int SizePZ=32; 
int Elements=1,Elements2=1;

int kflag=0,cpsf=0,Elem;

string I1FileName, I2FileName,OFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-cpsf",parg)) {cpsf=1;continue;}  // the psf has be be thought of as centered in the middle of image
   if (readArg("-i1",  I1FileName, parg)) continue;
   if (readArg("-i2",  I2FileName, parg)) continue;
   if (readArg("-o",  OFileName, parg)) continue;
   if (readArg("-iX",& SizeX,parg)) continue; // input size X
   if (readArg("-iY",& SizeY,parg)) continue;
   if (readArg("-iZ",& SizeZ,parg)) continue;
    usage(argv[0]);
  }

 Elements=1;

 ofstream to(OFileName.c_str());

 if (! to) cerr << "Fatal Error : Could not open outputfile " << OFileName << "\n", exit(-1);

 for (Elem=0;Elem<Elements;Elem++)
   {
     cout << "Loading Element " << Elem << "\n";
     Elements=InputImg1.DLoad(kflag,I1FileName.c_str(),"Float",& SizeX,& SizeY, & SizeZ,Elem);

     if (Elem==0)
       Elements2=InputImg2.DLoad(kflag,I2FileName.c_str(),"Float",& SizePX,& SizePY, & SizePZ,Elem);
     else
       if (Elem >= Elements2)
	 {
	   cerr << "WARNING : To few elements in PSF-File ! Taking element # " << Elements2-1 << " for convolution !\n";
	   Elements2=InputImg2.DLoad(kflag,I2FileName.c_str(),"Float",& SizePX,& SizePY, & SizePZ,Elements2-1);
	 }
       else
	 Elements2=InputImg2.DLoad(kflag,I2FileName.c_str(),"Float",& SizePX,& SizePY, & SizePZ,Elem);


     FImg2.Resize(SizeX,SizeY,SizeZ);
     Copy(& InputImg2,& FImg2); // Centered

     // FImg2.DSave(1,"/tmp/0.kdf");
     if (cpsf)
       FImg2.rotate(-SizeX/2,-SizeY/2,-SizeZ/2);  // move center to edge, will automatically adjust the extra size
     // FImg2.DSave(1,"/tmp/1.kdf");
     
     FImg2.Mul(1.0/(SizeX*SizeY*SizeZ));      // special normalization, because FFT is used in fast, unnormalized mode
     // normalized fft(true) whould need FImg2.Mul(sqrt(SizeX*SizeY*SizeZ));
  
     // cout << "FFT1\n";
     InputImg1.FFTfwd(false);
     // cout << "FFT2\n";
     FImg2.FFTfwd(false);
     
     // cout << "Multipl.\n";
     InputImg1.ConvMul(& FImg2);
     
     // cout << "FFTBwd\n";
     InputImg1.FFTbwd(false);
     
     if (Elem == 0)
       InputImg1.DHeader(kflag,to,Elements);
     
     InputImg1.Write(& to);
   }

 if (to) 
   to.close(); 
}
