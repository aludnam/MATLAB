// waveletspots : a trous - wvelet based recongnition of spots (J.C Olivio-Marin, Pattern Recognition 35 (2002) 1989-1996)

/*   This file is part of a software package written by 
     Rainer Heintzmann
     Institute of Applied Optics and Information Processing
     Albert Ueberle Strasse 3-5
     69120 Heidelberg
     Tel.: ++49 (0) 6221  549264
     Current Address : Max Planck Inst. for biophysical Chemistry, Am Fassberg 11, 37077 Goettingen, Germany
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
#include "tarrayclipf.h"
#include "parseargs.h"

typedef float ArrayBType;
typedef TArray3d<ArrayBType>   TImgArray;

static TImgArray InputImg, Processed, Wavelet, Product;

int mwrap(int pos,int dsize)   // returns position wrapped by mirroring at size
{
  if (pos < 0) pos = -pos;
  int dpos=(pos%(2*dsize));
  return (dpos >= dsize) ? (dsize - (dpos - dsize +1)) : dpos;
}

ArrayBType computeSigma( TImgArray & wavelet, TImgArray & help)
{
  int dimx = wavelet.GetSize(0);
  int dimy = wavelet.GetSize(1);
  int dimz = wavelet.GetSize(2);


  ArrayBType mean=wavelet.Mean();

  cout << " compute sigma: mean : " << mean << "\n";

  for (int z=0;z<dimz;z++)
    for (int y=0;y<dimy;y++)
      for (int x=0;x<dimx;x++)
	{
	  help.SetValue(x,y,z,fabs(wavelet.Value(x,y,z) - mean));  
	}
  help.Sort();

  ArrayBType * pp = help.Pointer(0,0,0);

  return pp[(dimx*dimy*dimz)/2 ] / 0.67;   // The MAD estimate (Sadler et al. IEEE Trans. Inform. Theory 45 (3), 1999, 1048-1051)
  
}

void applyFilter(TImgArray & img, TImgArray & img2, TImgArray & wavelet, int dir, int step)
{
  int dimx = img.GetSize(0);
  int dimy = img.GetSize(1);
  int dimz = img.GetSize(2);

  int dx=0,dy=0,dz=0,x,y,z;
  ArrayBType val;

  cout << "applying filter, dir " << dir << " step " << step << "\n";

  if (dir == 0)
    dx = step;
  if (dir == 1)
    dy = step;
  if (dir == 2)
    dz = step;

  for (z=0;z<dimz;z++)
    for (y=0;y<dimy;y++)
      for (x=0;x<dimx;x++)
	{
	  val=3.0/8.0* img.Value(x,y,z); 
	  val+=1.0/4.0 * img.Value(mwrap(x+dx,dimx),mwrap(y+dy,dimy),mwrap(z+dz,dimz)); 
	  val+=1.0/4.0 * img.Value(mwrap(x-dx,dimx),mwrap(y-dy,dimy),mwrap(z-dz,dimz)); 
	  val+=1.0/16.0 * img.Value(mwrap(x+2*dx,dimx),mwrap(y+2*dy,dimy),mwrap(z+2*dz,dimz)); 
	  val+=1.0/16.0 * img.Value(mwrap(x-2*dx,dimx),mwrap(y-2*dy,dimy),mwrap(z-2*dz,dimz)); 
	  wavelet.SetValue(x,y,z,img2.Value(x,y,z) - val);  // works only in 2D ! Since img2 contains the old data
	  img2.SetValue(x,y,z,val);
	}
}

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-o outputfile]\n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 
int Elements=1,i,x,y,z,wavelets=3,minwavelet=1;
static int INPUTSizeX=32;  // These is the standart size, if raw data is used
static int INPUTSizeY=32;  // 256
static int INPUTSizeZ=32;  // 22

 ArrayBType k=3.0;
int kflag=0,krelgiven=false;

string IFileName,O1FileName, O2FileName,O3FileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-i",  IFileName, parg)) continue;
   if (readArg("-nw", & wavelets, parg)) continue;
   if (readArg("-mw", & minwavelet, parg)) continue;  // number of the first wavelet (starting at zero)
   if (readArg("-krel", & k, parg)) {krelgiven=1;continue;}
   if (readArg("-o1",  O1FileName, parg)) continue;
   if (readArg("-o2",  O2FileName, parg)) continue;
   if (readArg("-o3",  O3FileName, parg)) continue;
   if (readArg("-IX",& INPUTSizeX, parg)) continue;
   if (readArg("-IY",& INPUTSizeY,parg)) continue;
   if (readArg("-IZ",& INPUTSizeZ,parg)) continue;
    usage(argv[0]);
  }

 if ( IFileName=="") usage(argv[0]);

 ofstream to1(O1FileName.c_str());
 ofstream to2(O2FileName.c_str());
 ofstream to3(O3FileName.c_str());

 if (O1FileName!="")
     if (! to1 )
       {
	 cerr << "Couldn't open file " << O1FileName << " for writing !!\n" << flush;
	 exit(-1);
       }

 if (O2FileName!="")
     if (! to2 )
       {
	 cerr << "Couldn't open file " << O2FileName << " for writing !!\n" << flush;
	 exit(-1);
       }

 if (O3FileName!="")
     if (! to3 )
       {
	 cerr << "Couldn't open file " << O3FileName << " for writing !!\n" << flush;
	 exit(-1);
       }

  for (i=0;i<Elements;i++)
  {

    Elements=InputImg.DLoad(kflag,IFileName.c_str(),"Float", & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,i);

    if (i==0)
      {
	Wavelet.Resize(INPUTSizeX,INPUTSizeY,INPUTSizeZ);
	Processed.Resize(INPUTSizeX,INPUTSizeY,INPUTSizeZ);
	Product.Resize(INPUTSizeX,INPUTSizeY,INPUTSizeZ);
	if (kflag)
	  {
	    if (to1 && kflag) WriteKhorosHeader(& to1,"Generated by waveletspots 2003","Float",INPUTSizeX,INPUTSizeY,INPUTSizeZ,wavelets*Elements);
	    cerr << "writing file " << O1FileName << " \n" << flush;
	    if (to2 && kflag) WriteKhorosHeader(& to2,"Generated by waveletspots 2003","Float",INPUTSizeX,INPUTSizeY,INPUTSizeZ,wavelets*Elements);
	    cerr << "writing file " << O1FileName << " \n" << flush;
	    if (to3 && kflag) WriteKhorosHeader(& to3,"Generated by waveletspots 2003","Float",INPUTSizeX,INPUTSizeY,INPUTSizeZ,Elements);
	    cerr << "writing file " << O1FileName << " \n" << flush;
	  }
      }

    int plvl = 1;
    ArrayBType sigma;

    for (int lvl = 0 ; lvl < wavelets;lvl++)
      {
	applyFilter(InputImg,Processed, Wavelet, 0 , plvl);  // steps is allways by 1 bigger than number of zeroes in between
	applyFilter(Processed,InputImg, Wavelet, 1 , plvl);
	// applyFilter(InputImg,Processed, Wavelet, 2 , plvl -1);   // For now, do not apply filter in Z-direction
	// InputImg.Copy(& Processed);  
	plvl *= 2;
	if (to1) InputImg.Write(& to1);

	sigma=computeSigma(Wavelet,Processed);
	cout << "sigma : " << sigma << " \n";
	if (krelgiven)
	  Wavelet.ClipAt(k * sigma,0.0);

	if (to2) Wavelet.Write(& to2);

	if (lvl==minwavelet)
	  Product.Copy(&Wavelet);
	else if (lvl > minwavelet)
	  Product.Mul(&Wavelet);
      }

    if (to3) Product.Write(& to3);

  }
  if (to1) to1.close();
  if (to2) to2.close();
  if (to3) to3.close();

}
