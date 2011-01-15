// imprawser : Impoort a series of raw data

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
#include "tarrayclipf.h"
#include "parseargs.h"

#include <time.h>

typedef unsigned short ArrayBType;
typedef TArray3d<ArrayBType> TImgArray;  // Clipping is needed
typedef unsigned char ArrayBBType;
typedef TArray3d<ArrayBBType> TImgBArray;  // Clipping is needed

static TImgArray InputImg;
static TImgBArray InputBImg;

static TImgArray OutputImg;
static TImgBArray OutputBImg;


void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-o outputfile] \n" << flush;
  exit(-1);
}


int main(int argc, char *argv[])
{ 
int i;
bool ByteType=false;
static int INPUTSizeX=256;  // These is the standart size, if raw data is used
static int INPUTSizeY=256;  // 256
static int INPUTSizeZ=1;  // 256
static int INPUTSizeT=1;  // 256
static int INPUTSizeE=1;  // 256
 static int RESAMPLEX=1;
 static int RESAMPLEY=1;
 static int RESAMPLEZ=1;

// static int INPUTSizeZ=1;  // 22

int offset=0;

string  IFileName,O1FileName,LoadFile, LoadFormat, Ext, Format("%.2d");

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
    if (readArg("-k",parg)) {continue;}  // ignore
   if (readArg("-bt",parg)) {ByteType=true;continue;}
   // if (readArg("-us",parg)) {ushort=true;continue;}
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-o", O1FileName, parg)) continue;
   if (readArg("-sizex",& INPUTSizeX, parg)) continue;
   if (readArg("-sizey",& INPUTSizeY,parg)) continue;
   if (readArg("-sizez",& INPUTSizeZ,parg)) continue;
   if (readArg("-sizet",& INPUTSizeT,parg)) continue;
   if (readArg("-sizee",& INPUTSizeE,parg)) continue;
   if (readArg("-binx",& RESAMPLEX, parg)) continue;
   if (readArg("-biny",& RESAMPLEY,parg)) continue;
   if (readArg("-binz",& RESAMPLEZ,parg)) continue;
   if (readArg("-offset",& offset,parg)) continue;
    usage(argv[0]);
  }
 cout << "Original Data Sizes: " << INPUTSizeX << "x" << INPUTSizeY << "x" << INPUTSizeZ << "\n" << flush;

 if (IFileName=="") usage(argv[0]);

 ofstream to1(O1FileName.c_str());
 if (O1FileName!="")
     if (! to1 )
       {
           cerr << "Couldn't open file " << O1FileName << " for writing !!\n" << flush;
           exit(-1);
       }

 if (ByteType)
   OutputBImg.Resize(INPUTSizeX/RESAMPLEX,INPUTSizeY/RESAMPLEY,INPUTSizeZ/RESAMPLEZ); // (last-first)/step+1);
 else       
   OutputImg.Resize(INPUTSizeX/RESAMPLEX,INPUTSizeY/RESAMPLEY,INPUTSizeZ/RESAMPLEZ); // (last-first)/step+1);

 if (ByteType)
   OutputBImg.DHeader(true,to1,INPUTSizeE,INPUTSizeT);  // always save khoros format !
 else
   OutputImg.DHeader(true,to1,INPUTSizeE,INPUTSizeT);  // always save khoros format !

 int j=0;
 for (i=0;i<INPUTSizeE*INPUTSizeT;i++)
   {
       if (ByteType)
       {
            InputBImg.DLoad(false,IFileName.c_str(),& INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,0,offset); 
	    OutputBImg.ResampleFrom(&InputBImg,RESAMPLEX,RESAMPLEY,RESAMPLEZ);
            OutputBImg.Write(&to1);
            offset += INPUTSizeX*INPUTSizeY*INPUTSizeZ;
       }
       else
       {
          InputImg.DLoad(false,IFileName.c_str(),& INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,0,offset); 
	  OutputImg.ResampleFrom(&InputImg,RESAMPLEX,RESAMPLEY,RESAMPLEZ);
          OutputImg.Write(&to1);
          offset += 2*INPUTSizeX*INPUTSizeY*INPUTSizeZ;
       }
   }

   if (ByteType)
     OutputBImg.Write(&to1);
   else
     OutputImg.Write(&to1);

 if (to1)
   to1.close();
}
