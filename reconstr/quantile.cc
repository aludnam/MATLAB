
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

// quantile.cc : determines the quantile. 0.5 : median, 0.9 : start of 90% quantile
#include <iostream>
#include <math.h>
#include <string>
#include "rawarray.h"
#include "parseargs.h"


typedef float ArrayBType;
typedef TArray3d<ArrayBType>  TArray;  // Spots array 

static TArray   Array1, Gate, Array2;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] -i1 image1 -i2 image2 -o outputfile -d distancesFile\n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

static int kflag=0,x,y,greater=0;
bool DoClip=false;
static int SizeX=0,SizeY=0,SizeZ=0,elements=0,valid=0;
double quantile=0.5,qmax=1.0;
int column=0;

string OFileName="", IFileName="", GFileName="", TFileName="";

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-clip",parg)) {DoClip=true;continue;}
   if (readArg("-quantile", &quantile, parg)) continue;
   if (readArg("-qmax", &qmax, parg)) continue;
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-g", GFileName, parg)) continue;
   if (readArg("-t", TFileName, parg)) continue;  // for a thresholded output
   if (readArg("-i",  IFileName, parg)) continue;
    usage(argv[0]);
  }

  if (OFileName=="" || IFileName=="") 
    usage(argv[0]);

  elements=Array1.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY,& SizeZ,0);
  Array2.Resize(1,1,1);

ofstream * oi=0;

 if (OFileName!="")
   {
     oi=new ofstream(OFileName.c_str());
     if (! oi)
       cerr << "Error ! Couldnt open " << OFileName << " for writing !! \n",exit(-1);
     Array2.DHeader(kflag,*oi,elements);  // determines its sizes on its own
   }

ofstream * ot=0;

 if (TFileName!="")
   {
     ot=new ofstream(TFileName.c_str());
     if (! ot)
       cerr << "Error ! Couldnt open " << TFileName << " for writing !! \n",exit(-1);
     Array1.DHeader(kflag,*ot,elements);  // determines its sizes on its own
   }

 int InvalPixels=0;
 ArrayBType Res=0.0,ResMax=0.0;
 int qpos,qposx,qposy,qposz;

 for (int  Elem=0;Elem<elements;Elem++)
  {
   if (Elem > 0)
       elements=Array1.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY,& SizeZ,Elem);
   if (GFileName != "")
	{
	Gate.DLoad(kflag,GFileName.c_str(),"Float",& SizeX,& SizeY,& SizeZ,Elem);
   	if (! Array1.SizesEqual(& Gate))
           cerr << "Error ! Gate has to be equal in size to data! \n",exit(-1);
	InvalPixels=Gate.CntBelowToMax(0.5,& Array1,HUGE_VAL);
	}
   Array1.Sort();  
   qpos = (int) (quantile*(SizeX*SizeY*SizeZ-InvalPixels)+0.5);
   if (qpos >= SizeX*SizeY*SizeZ) qpos = SizeX*SizeY*SizeZ-1;
   qposz=qpos/(SizeX*SizeY);
   qposy=(qpos-(SizeX*SizeY)*qposz)/SizeX;
   qposx=qpos % SizeX;
   Res = Array1.Value(qposx,qposy,qposz);
   cout << "Quantile position: " << qposx << ", " << qposy << ", " << qposz << ", total : " << qpos << " invalid : " << InvalPixels << "\n";
   cout << "Quantile " << 100*quantile <<"\% starts at value " << Res << "\n";

   qpos = (int) (qmax*(SizeX*SizeY*SizeZ-InvalPixels)+0.5);
   if (qpos >= SizeX*SizeY*SizeZ) qpos = SizeX*SizeY*SizeZ-1;
   qposz=qpos/(SizeX*SizeY);
   qposy=(qpos-(SizeX*SizeY)*qposz)/SizeX;
   qposx=qpos % SizeX;
   ResMax = Array1.Value(qposx,qposy,qposz);
   cout << "MaxQuantile position: " << qposx << ", " << qposy << ", " << qposz << ", total : " << qpos << " invalid : " << InvalPixels << "\n";
   cout << "MaxQuantile " << 100*quantile <<"\% starts at value " << Res << "\n";

   Array2.SetValue(0,0,0,Res);

   Array2.Write(oi);

   if (ot != 0)
	{
         elements=Array1.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY,& SizeZ,Elem);
	 if (DoClip)
		if (qmax < 1.0)
	 	Array1.ClipBetween(Res,ResMax,0);	
		else
	 	Array1.ClipAt(Res,0);	
	 else
		if (qmax < 1.0)
	 	Array1.ThreshBetween(Res,ResMax,0,1);	
		else
	 	Array1.ThreshAt(Res,0,1);	
   	if (GFileName != "")
		Array1.Mul(&Gate);
         Array1.Write(ot);
	}
  }

if (oi)  oi->close();
if (ot)  ot->close();

}
