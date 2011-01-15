// normaize.cc : Normalizes float image by different methods

/*   This file is part of a software package written by 
     Rainer Heintzmann
     MPI bpc
     Am Fassberg 11
     37077 Goettingen, Germany
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

typedef float ArrayType;
typedef TArray3d<ArrayType>  TImgArray; 
typedef TArray3d<char>  BImgArray; 

static TImgArray    Img, Scale;
static BImgArray    BImg;
static char dirstring[]="xyz????";

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " -i inputfile -o outputfile [-sX SizeX] [-sY SizeY] [-sZ SizeZ] [-sE Elements] [-integral] [-val value] [-byte] [-series] [-slicedir directionnumber] [-min]\n" << flush;
  cerr << "The zero will always stay where it is ! If -int is not selected it will be   normalized to maximum, else to integral\n";
  cerr << "-byte converts result into byte format\n-series will save a series of slices\n-slicedir number  indicates which direction to slice (0-2), elements will allways be sliced. Naming conventions are compatible with Leica TCS\n";
  exit(-1);
}

int main(int argc, char *argv[])
{ 

static int SizeX=64,SizeY=64,SizeZ=64;
static bool integral=false,kflag=false,indiv=false,nonormalize=false,color=false,appletoutput=false,byte=false,series=false,min=false;
ArrayType val=1;
 int Elem=0,Elements=1,exelem=0,slicedir=2;
 double scale=1.0;

string IFileName,OFileName,ScaleFileName;

char sliceName[5000];
char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
    if (readArg("-byte", parg)) {byte=true;continue;}
    if (readArg("-series", parg)) {series=true;continue;}  // saves the data slice by slice. Compatible with the Leica conventions
    if (readArg("-slicedir", &slicedir, parg)) continue;  // Direction to separate into slices
    if (readArg("-integral", parg)) {integral=true;continue;}
    if (readArg("-min", parg)) {min=true;continue;}
    if (readArg("-nonormalize", parg)) {nonormalize=true;continue;}  // if selected, no normalization will be performed, only possible element expansion
    if (readArg("-color", parg)) {color=true;continue;}
    if (readArg("-expand", & exelem, parg)) continue; 
    if (readArg("-indiv", parg)) {indiv=true;continue;}  // individually
    if (readArg("-val", & val, parg)) continue; 
    if (readArg("-sX", & SizeX, parg)) continue; 
    if (readArg("-sY", & SizeY, parg)) continue;
    if (readArg("-sZ", & SizeZ, parg)) continue;
    if (readArg("-sE", & Elements, parg)) continue;
    if (readArg("-i",  IFileName, parg)) continue;
    if (readArg("-o",  OFileName, parg)) continue;
    if (readArg("-sf",  ScaleFileName, parg)) continue;
    if (readArg("-k", parg)) {kflag=true;continue;}
    if (readArg("-a", parg)) {appletoutput=true;continue;}
    usage(argv[0]);
  }

  if (slicedir > 2 || slicedir<0) slicedir=2;

  if (OFileName=="") 
    usage(argv[0]);

   ofstream to(OFileName.c_str());

   ofstream * sfile= 0;
   if (ScaleFileName!="") sfile = new ofstream(ScaleFileName.c_str());

   double MaxVal = -1e30;
   double MinVal = 1e30;

   if (! indiv)
     {
       double MaxV,MinV;
       for (Elem=0;Elem<Elements;Elem++)
	 {
	   if (! appletoutput)
	     cout << "Loading Element " << Elem << "\n";
	 Elements= Img.DLoad(kflag,IFileName.c_str(),"Float",& SizeX, &SizeY, & SizeZ,Elem);
	 MinV =  Img.Minimum();
	 if(MinV < MinVal) MinVal=MinV;
	 MaxV =  Img.Maximum();
	 if(MaxV > MaxVal) MaxVal=MaxV;
	 }
       if (! appletoutput)
	{
	 cout << "Maximum Value is " << MaxVal << "\n";
	 cout << "Minimum Value is " << MinVal << "\n";
	}
     }

 for (Elem=0;Elem<Elements;Elem++)
   {
     if (! appletoutput)
       cout << "Loading Element " << Elem << "\n";
     Elements= Img.DLoad(kflag,IFileName.c_str(),"Float",& SizeX, &SizeY, & SizeZ,Elem);
     if (indiv)
         MinVal =  Img.Minimum();
     if (min) Img.Sub(MinVal);

     if (!nonormalize)
       {
     if (integral)
       {
	 if (indiv)
	   {
	     scale = val / Img.Integral();
	     if (! appletoutput)
	       cout << "Elements will be scaled by : "<<scale<<"\n";
	     // Img.Normalize(val);
	   }
	 else if (Elem==0)
	   {
	     scale = val / Img.Integral();
	     if (! appletoutput)
	       cout << "Common scale from Elemnt 0 : "<<scale<<"\n";
	   }
       }
     else
       {
	 if (indiv)
	   {
	     MaxVal =  Img.Maximum();
	     if (MaxVal == 0)
		 {
		   if (! appletoutput)
		     cerr << "WARNING, Element " << Elem+1 << " could not be normalized. Max == 0 !\n";
		   scale = 1.0;
		 }
	     else
	       scale = val / MaxVal; // Img.Maximum();
	   }
	 else
	   scale = val / MaxVal; // Img.Maximum();
       }
     if (! appletoutput)
       cout << "Elements will be scaled by : "<<scale<<"\n";
     else
       {
	 cerr << "<param name=scalev"<<Elem+1<<" value=" << 1.0/scale <<">\n";
	 cerr << "<param name=offsetv"<<Elem+1<<" value=" << 0.0 <<">\n";
	 cout << "<param name=scalev"<<Elem+1<<" value=" << 1.0/scale <<">\n";
	 cout << "<param name=offsetv"<<Elem+1<<" value=" << 0.0 <<">\n";
       }

     Img.Mul(scale);
       }


     if (Elem == 0)
       {
	 Scale.Resize(1,1,1);
	 if (exelem > 0 && Elements < exelem)
	   if (byte)
	    BImg.DHeader(kflag,to,exelem);
	   else
	    Img.DHeader(kflag,to,exelem);
	 else
	   if (byte)
	    BImg.DHeader(kflag,to,Elements);
	   else
	    Img.DHeader(kflag,to,Elements);

	 if (sfile)
	   if (exelem > 0 && Elements < exelem)
	     Scale.DHeader(kflag,* sfile,exelem);
	   else
	     Scale.DHeader(kflag,* sfile,Elements);
       }
     Scale.SetValue(0,0,0,scale);

    if (byte)
    {
	if (series)
	{
	if (Elem == 0 && slicedir == 0)
		BImg.Resize(1,Img.GetSize(1),Img.GetSize(2));
	if (Elem == 0 && slicedir == 1)
		BImg.Resize(Img.GetSize(0),1,Img.GetSize(2));
	if (Elem == 0 && slicedir == 2)
		BImg.Resize(Img.GetSize(0),Img.GetSize(1),1);
	for (int slicepos=0;slicepos<Img.GetSize(slicedir);slicepos++)
		{
		BImg.extract(slicedir,slicepos,& Img); // copy each element into the Byte Array
		//sprintf(sliceName,"%s%c%.3d_ch%.2d.raw",OFileName.c_str(),dirstring[slicedir],slicepos,Elem);
		sprintf(sliceName,OFileName.c_str(),slicepos,Elem);
		BImg.DSave(kflag,sliceName);
		}
	}
	else 
	{
	if (Elem == 0)
		BImg.Resize(Img.GetSize(0),Img.GetSize(1),Img.GetSize(2));
	Copy(& Img,& BImg); // copy each element into the Byte Array
      	BImg.Write(& to);
	}
    }
    else
      Img.Write(& to);
     if (sfile)
       Scale.Write(sfile);
   }

 if (exelem > 0 && Elements < exelem)
   {
     if (color) 
       if (byte)
         BImg.Clear();
       else
         Img.Clear();
   for (Elem=0;Elem<exelem-Elements;Elem++)
     {
       if (byte)
         BImg.Write(& to);
       else
         Img.Write(& to);
     }
   }
 if (sfile)
   sfile->close();
 to.close();
}
