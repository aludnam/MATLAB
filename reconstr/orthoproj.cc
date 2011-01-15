// orthoproj : This program produces 3 othogonal projections from a datastack

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

#include <time.h>
#include "poidev.h"

long rval=-5;       // initvalue of random-generator, can be changed

typedef float ArrayBType;
typedef TArray3d<ArrayBType> TProjArray;
typedef TArray3d<ArrayBType>  TImgArray;

static TImgArray InputImg;
static TProjArray O1Img,O2Img,O3Img;


void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-av] [-mip1] [-mip2] [-o1 outputfile] [-o2 outputfile] [-o3 outputfile] \n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 
int Elements=1,i,method=1,x,y,z, maxp,minp;
ArrayBType sum,max,min,val;
static int INPUTSizeX=32;  // These is the standart size, if raw data is used
static int INPUTSizeY=32;  // 256
static int INPUTSizeZ=32;  // 22

int kflag=0, average=0;

string  IFileName,O1FileName,O2FileName,O3FileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-i",  IFileName, parg)) continue;
   if (readArg("-o1",  O1FileName, parg)) continue;
   if (readArg("-o2",  O2FileName, parg)) continue;
   if (readArg("-o3",  O3FileName, parg)) continue;
   if (readArg("-IX",& INPUTSizeX, parg)) continue;
   if (readArg("-IY",& INPUTSizeY,parg)) continue;
   if (readArg("-IZ",& INPUTSizeZ,parg)) continue;
   if (readArg("-method", & method, parg)) continue;  // 1 : sum , 2 : avg, 3 : mip, 4 : max- min, 5: max+min-2avg
   if (readArg("-av",parg)) {average=1;continue;}
    usage(argv[0]);
  }

 if (average == 1) method = 2; // force the method to average for compatibility reasons

 switch (method) {
 case 1:
   cout << "Projection Method is SUM\n"; break;
 case 2:
   cout << "Projection Method is AVERAGE\n"; break;
 case 3:
   cout << "Projection Method is MAX\n"; break;
 case 4:
   cout << "Projection Method is MAX-MIN\n"; break;
 case 5:
   cout << "Projection Method is MAX+MIN-2AVG\n"; break;
 case 6:
   cout << "Projection Method is MAX Position\n"; break;
 case 7:
   cout << "Projection Method is MIN Position\n"; break;
 }

 if (IFileName=="") usage(argv[0]);

 ofstream to1(O1FileName.c_str());
 if (O1FileName!="")
     if (! to1 )
       {
	 cerr << "Couldn't open file " << O1FileName << " for writing !!\n" << flush;
	 exit(-1);
       }
 ofstream to2(O2FileName.c_str());
 if (O2FileName!="")
     if (! to2 )
       {
	 cerr << "Couldn't open file " << O2FileName << " for writing !!\n" << flush;
	 exit(-1);
       }
 ofstream to3(O3FileName.c_str());
 if (O2FileName!="")
     if (! to3 )
       {
	 cerr << "Couldn't open file " << O3FileName << " for writing !!\n" << flush;
	 exit(-1);
       }

 float divisor = 1.0;
      
  for (i=0;i<Elements;i++)
  {

    Elements=InputImg.DLoad(kflag,IFileName.c_str(),"Float",
			  & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,i);

    if (i==0)
      {
	O1Img.Resize(INPUTSizeX,INPUTSizeY,1);
	O2Img.Resize(INPUTSizeX,INPUTSizeZ,1);
	O3Img.Resize(INPUTSizeY,INPUTSizeZ,1);
      }

    divisor=1.0/INPUTSizeZ;

    for (y=0;y<INPUTSizeY;y++)
      for (x=0;x<INPUTSizeX;x++)
	{
	  val=InputImg.Value(x,y,0); 
	  sum=val;
	  max=val;
	  min=val;
	  maxp=0,minp=0;
	  for (z=1;z<INPUTSizeZ;z++)
	    {
	      val = InputImg.Value(x,y,z);
	      sum+=val;
	      if (val > max) max = val,maxp = z;
	      if (val < min) min = val,minp = z;
	    }
	  switch (method) {
	  case 1:
	    O1Img.SetValue(x,y,0,sum); break;
	  case 2:
	    O1Img.SetValue(x,y,0,sum*divisor); break;
	  case 3:
	    O1Img.SetValue(x,y,0,max); break;
	  case 4:
	    O1Img.SetValue(x,y,0,max-min); break;
	  case 5:
	    O1Img.SetValue(x,y,0,max+min-2*sum*divisor); break;
	  case 6:
	    O1Img.SetValue(x,y,0,maxp); break;
	  case 7:
	    O1Img.SetValue(x,y,0,minp); break;
	  }
	}

    divisor=1.0/INPUTSizeY;

    for (z=0;z<INPUTSizeZ;z++)
      for (x=0;x<INPUTSizeX;x++)
	{
	  val=InputImg.Value(x,0,z); 
	  sum=val;
	  max=val;
	  min=val;
	  maxp=0,minp=0;
	  for (y=1;y<INPUTSizeY;y++)
	    {
	      val = InputImg.Value(x,y,z);
	      sum+=val;
	      if (val > max) max = val, maxp = y;
	      if (val < min) min = val, minp = y;
	    }
	  switch (method) {
	  case 1:
	    O2Img.SetValue(x,z,0,sum); break;
	  case 2:
	    O2Img.SetValue(x,z,0,sum*divisor); break;
	  case 3:
	    O2Img.SetValue(x,z,0,max); break;
	  case 4:
	    O2Img.SetValue(x,z,0,max-min); break;
	  case 5:
	    O2Img.SetValue(x,z,0,max+min-2*sum*divisor); break;
	  case 6:
	    O2Img.SetValue(x,z,0,maxp); break;
	  case 7:
	    O2Img.SetValue(x,z,0,minp); break;
	  }
	}

    divisor=1.0/INPUTSizeX;

    for (z=0;z<INPUTSizeZ;z++)
      for (y=0;y<INPUTSizeY;y++)
	{
	  val=InputImg.Value(0,y,z); 
	  sum=val;
	  max=val;
	  min=val;
	  maxp=0,minp=0;
	  for (x=1;x<INPUTSizeX;x++)
	    {
	      val = InputImg.Value(x,y,z);
	      sum+=val;
	      if (val > max) max = val, maxp = x;
	      if (val < min) min = val, minp = x;
	    }
	  switch (method) {
	  case 1:
	    O3Img.SetValue(y,z,0,sum); break;
	  case 2:
	    O3Img.SetValue(y,z,0,sum*divisor); break;
	  case 3:
	    O3Img.SetValue(y,z,0,max); break;
	  case 4:
	    O3Img.SetValue(y,z,0,max-min); break;
	  case 5:
	    O3Img.SetValue(y,z,0,max+min-2*sum*divisor); break;
	  case 6:
	    O3Img.SetValue(y,z,0,maxp); break;
	  case 7:
	    O3Img.SetValue(y,z,0,minp); break;
	  }
	}

    if (i==0)
      {
	if (kflag)
	  {
	    if (to1 && kflag) WriteKhorosHeader(& to1,"Generated by poissit 1998","Float",INPUTSizeX,INPUTSizeY,1,Elements);
	    cerr << "writing file " << O1FileName << " \n" << flush;
	    if (to2 && kflag) WriteKhorosHeader(& to2,"Generated by poissit 1998","Float",INPUTSizeX,INPUTSizeZ,1,Elements);
	    cerr << "writing file " << O1FileName << " \n" << flush;
	    if (to3 && kflag) WriteKhorosHeader(& to3,"Generated by poissit 1998","Float",INPUTSizeY,INPUTSizeZ,1,Elements);
	    cerr << "writing file " << O1FileName << " \n" << flush;
	  }
	// else nothing
      }

    if (to1) O1Img.Write(& to1);
    if (to2) O2Img.Write(& to2);
    if (to3) O3Img.Write(& to3);
  }
  if (to1) to1.close();
  if (to2) to2.close();
  if (to2) to2.close();

}
