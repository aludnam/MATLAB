// window : this program reads in a 3-D datafile and applies a windowing

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
#include <time.h>

#include "fftarray.h"
#include "parseargs.h"
#include "window.h"

typedef float ArrayBType;
typedef complex<float> ArrayBCType;

typedef TArray3d<ArrayBType>    TImgArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TFFTArray<ArrayBCType>  TCImgArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

static TImgArray     InputImg;
static TCImgArray     CInputImg;

bool rectborder=false,ring=false,center=true;
float gaussstretch=1.0,OutX=1.0,InX=0.0,OutY=1.0,InY=0.0,OutZ=1.0,InZ=0.0;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-w windowtype] [-rect] [-i inputfile] [-o outputfile] [-relo OuterBorder] [-reli InnerBorder] \n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 
int Elements=1,i;
static int INPUTSizeX=32;  // These is the standart size, if raw data is used
static int INPUTSizeY=32;  // 256
static int INPUTSizeZ=32;  // 22
bool kflag=false,cdata=false,mirrored=false;
 int window=4;  // Window code

string INPUTFileName,OUTPUTFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=true;continue;}
   if (readArg("-c",parg)) {cdata=true;continue;}
   if (readArg("-nocenter",parg)) {center=false;continue;}
   if (readArg("-i",  INPUTFileName, parg)) continue;
   if (readArg("-o", OUTPUTFileName, parg)) continue;
   if (readArg("-sX",& INPUTSizeX, parg)) continue;
   if (readArg("-sY",& INPUTSizeY,parg)) continue;
   if (readArg("-sZ",& INPUTSizeZ,parg)) continue;
   if (readArg("-window", & window, parg)) continue; 
   if (readArg("-mirror", parg)) {mirrored=true;continue;}
   if (readArg("-rect", parg)) {rectborder=true;continue;}
   if (readArg("-ring", parg)) {ring=true;continue;}
   if (readArg("-relXo", & OutX, parg)) continue; 
   if (readArg("-relXi", & InX, parg)) continue; 
   if (readArg("-relYo", & OutY, parg)) continue; 
   if (readArg("-relYi", & InY, parg)) continue; 
   if (readArg("-relZo", & OutZ, parg)) continue; 
   if (readArg("-relZi", & InZ, parg)) continue; 
   if (readArg("-gs", & gaussstretch, parg)) continue;
    usage(argv[0]);
  }


  ofstream to(OUTPUTFileName.c_str());
      
  if (! to )
    {
      cerr << "Couldn't open file " << OUTPUTFileName << " for writing !!\n" << flush;
      exit(-1);
    }
	  
  if (kflag)
    {
      ifstream dummy(INPUTFileName.c_str());
      char name[100];
      if (! dummy)
	cerr << "Couldn't open inputfile " << OUTPUTFileName << "\n" << flush, exit(-1);

      ReadKhorosHeader(& dummy,name,INPUTSizeX,INPUTSizeY,INPUTSizeZ,Elements);
      if ((strcmp(name,"Complex")==0) || (strcmp(name,"Double Complex")==0))
	cdata=1;  // Force to complex data !
      dummy.close();
    }

  for (i=0;i<Elements;i++)
  {
    if (! cdata)
      Elements=InputImg.DLoad(kflag,INPUTFileName.c_str(),"Float",
			      & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,i);
    else
      Elements=CInputImg.DLoad(kflag,INPUTFileName.c_str(),"Complex",
			      & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,i);

      if (INPUTSizeX <= 1)
      {
        InX = 10000.0; OutX = 20000.0;
      }
      if (INPUTSizeY <= 1)
      {
        InY = 10000.0; OutY = 20000.0;
      }
      if (INPUTSizeZ <= 1)
      {
        InZ = 10000.0; OutZ = 20000.0;
      }

    if (! cdata)
      switch (window)
	{
	case 1:
	  {
	  GaussWindow<TImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.gaussstretch=gaussstretch;
	  Win.ApplyWindow(&InputImg,center,mirrored);
	  }
	  break;
	case 2:
	  {
	  EdgeWindow<TImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.ApplyWindow(&InputImg,center,mirrored);
	  }
	  break;
	case 3:
	  {
	  LinWindow<TImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.ApplyWindow(&InputImg,center,mirrored);
	  }
	  break;
	case 4:
	  {
	  SinWindow<TImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.ApplyWindow(&InputImg,center,mirrored);
	  }
	  break;
	case 5:
	  {
	  NearGaussWindow<TImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.ApplyWindow(&InputImg,center,mirrored);
	  }
	  break;
	default:
	  usage(argv[0]);
	}
    else
      switch (window)
	{
	case 1:
	  {
	  GaussWindow<TCImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.gaussstretch=gaussstretch;
	  Win.ApplyWindow(&CInputImg,center,mirrored);
	  }
	  break;
	case 2:
	  {
	  EdgeWindow<TCImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.ApplyWindow(&CInputImg,center,mirrored);
	  }
	  break;
	case 3:
	  {
	  LinWindow<TCImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.ApplyWindow(&CInputImg,center,mirrored);
	  }
	  break;
	case 4:
	  {
	  SinWindow<TCImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.ApplyWindow(&CInputImg,center,mirrored);
	  }
	  break;
	case 5:
	  {
	  NearGaussWindow<TCImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.ApplyWindow(&CInputImg,center,mirrored);
	  }
	  break;
	case 6:
	  {
	  TunnelWindow<TCImgArray> Win(InX,OutX,InY,OutY,InZ,OutZ,rectborder,ring);
	  Win.gaussstretch=gaussstretch;
	  Win.ApplyWindow(&CInputImg,center,mirrored);
	  }
	  break;
	default:
	  usage(argv[0]);
	}

    if (i==0)
      {
	if (kflag)
	  {
	    if (! cdata)
	      WriteKhorosHeader(& to,"Generated by window.cc 2000","Float",INPUTSizeX,INPUTSizeY,INPUTSizeZ,Elements);
	    else
	      WriteKhorosHeader(& to,"Generated by window.cc 2000","Complex",INPUTSizeX,INPUTSizeY,INPUTSizeZ,Elements);
	    cerr << "writing file " << OUTPUTFileName << " \n" << flush;
	  }
	// else nothing  -> only data is written
	
      }

    if (! cdata)
      InputImg.Write(& to);
    else
      CInputImg.Write(& to);
  }
  to.close();

}
