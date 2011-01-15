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

static TImgArray OutputImg;
static TImgBArray OutputBImg;


void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-o outputfile] \n" << flush;
  exit(-1);
}

ifstream * NewIFStream(string Front, string Format, int num, string Ext)
{
 char  numbers[1000];
 sprintf(numbers,Format.c_str(),num);
 return new ifstream((Front+numbers+Ext).c_str());
 cout << "testing: " << num <<"#\n" << flush;
}

void parseFileName(string & FileName, string & Format, string & Ext, int & first)
{
  char  numbers[1000];
  int endnum=FileName.find_last_of("0123456789");
  string Rest(FileName,0,endnum+1);
  int begnum=Rest.find_last_not_of("0123456789");
  sprintf(numbers,"%d",endnum-begnum);
  Format=string("%.")+numbers+"d";
  Ext= FileName;
  Ext.erase(0,endnum+1);
  string Num=Rest;
  Num.erase(0,begnum+1);
  first=atoi(Num.c_str());
  FileName.erase(begnum+1,FileName.length()-(begnum+1));
  // cout << "Parameters: " << FileName <<"#\n";
  // cout << "Rest:" << Rest <<"#\n";
  // cout << "Ext:" << Ext <<"#\n";
  // cout << "Num:" << Num <<"#\n";
  // cout << "first:" << first <<"#\n";
  // cout << "endnum:" << endnum <<"#\n";
  // cout << "begnum:" << begnum <<"#\n";
}

void checkparameters(string & FileName, string & Format, string & Ext,int & first,int & last) // fills this values in automatically if possible
{
 ifstream * in;
 if (FileName=="") return;
 if (Format != ""  && last-first >= 0)
 {
   in = NewIFStream(FileName,Format,first,Ext);
   if (* in)
 	{
	in->close();
	return;
	}
 }
 in = new ifstream(FileName.c_str()); 
 if (* in)
	{
	parseFileName(FileName,Format,Ext,first);
	int i=first;
	do	{
		in->close();
		i++;
		in=NewIFStream(FileName,Format,i,Ext);
		}
	while (* in);
	last=i-1;
	}
}

int main(int argc, char *argv[])
{ 
int i;
bool ByteType=false;
static int INPUTSizeX=256;  // These is the standart size, if raw data is used
static int INPUTSizeY=256;  // 256
static int INPUTSizeZ=1;  // 256
static int INPUTSizeT=-1;  // 256
static int INPUTSizeE=-1;  // 256

// static int INPUTSizeZ=1;  // 22

int first=0,last=10,step=1,offset=0;

string  IFileName,O1FileName,LoadFile, LoadFormat, Ext, Format("%.2d");

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-bt",parg)) {ByteType=true;continue;}
   // if (readArg("-us",parg)) {ushort=true;continue;}
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-o",  O1FileName, parg)) continue;
   if (readArg("-sizex",& INPUTSizeX, parg)) continue;
   if (readArg("-sizey",& INPUTSizeY,parg)) continue;
   if (readArg("-sizez",& INPUTSizeZ,parg)) continue;
   if (readArg("-sizet",& INPUTSizeT,parg)) continue;
   if (readArg("-sizee",& INPUTSizeE,parg)) continue;
   if (readArg("-offset",& offset,parg)) continue;
   if (readArg("-format", Format,parg)) continue;
   if (readArg("-ext",Ext,parg)) continue;
   // if (readArg("-sizez",& INPUTSizeZ,parg)) continue;
   if (readArg("-first",& first,parg)) continue;
   if (readArg("-last",& last,parg)) continue;
   if (readArg("-step",& step,parg)) continue;
    usage(argv[0]);
  }

 if (IFileName=="") usage(argv[0]);

 ofstream to1(O1FileName.c_str());
 if (O1FileName!="")
     if (! to1 )
       {
           cerr << "Couldn't open file " << O1FileName << " for writing !!\n" << flush;
           exit(-1);
       }


 checkparameters(IFileName,Format,Ext,first,last);  // fills this values in automatically if possible

 int DataSizeZ=INPUTSizeZ;
 int SizeE;
 int SizeT=1;

 if (INPUTSizeE>0)
   SizeE=INPUTSizeE;
 else
   SizeE=(last-first)/step+1;
 
 if (INPUTSizeT>0)
   SizeT=INPUTSizeT;

 if (INPUTSizeZ==1) 
	{DataSizeZ=SizeE,SizeE=1;}

 if (INPUTSizeZ<1) 
        {INPUTSizeZ=1;DataSizeZ=1;}

 if (ByteType)
   OutputBImg.Resize(INPUTSizeX,INPUTSizeY,DataSizeZ); // (last-first)/step+1);
 else       
   OutputImg.Resize(INPUTSizeX,INPUTSizeY,DataSizeZ); // (last-first)/step+1);

 if (ByteType)
   OutputBImg.DHeader(true,to1,SizeE,SizeT);  // always save khoros format !
 else
   OutputImg.DHeader(true,to1,SizeE,SizeT);  // always save khoros format !

 int j=0;
 for (i=first;i<=last;i+=step)
   {
       char  numbers[1000];
       sprintf(numbers,Format.c_str(),i);
       // cout << "Numbers: " << numbers <<"\n";

       LoadFile=IFileName+numbers+Ext;

       if (ByteType)
       {
	 if (INPUTSizeZ==1 && DataSizeZ != 1)
          OutputBImg.LoadSlice(LoadFile.c_str(),j,offset);
	 else
	  {
            OutputBImg.DLoad(false,LoadFile.c_str(),& INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,0,offset); 
            OutputBImg.Write(&to1);
	  }
       }
       else
       {
	 if (INPUTSizeZ==1 && DataSizeZ != 1)
          OutputImg.LoadSlice(LoadFile.c_str(),j,offset);
	 else
	  {
          OutputImg.DLoad(false,LoadFile.c_str(),& INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,0,offset); 
          OutputImg.Write(&to1);
	  }
       }
       j++;
       if (j >= DataSizeZ) j=0;
   }
 if (INPUTSizeZ==1 && DataSizeZ != 1)
   if (ByteType)
     OutputBImg.Write(&to1);
   else
     OutputImg.Write(&to1);

 if (to1)
   to1.close();
}
