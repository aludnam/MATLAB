// reads and writes Andreas Schoenle's data format

// andreasio.cc : expands one image to the sizes of a second one by wrapping

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

#include <stdlib.h>
#include <alloca.h>
#include <iostream>
#include <string>
#include <math.h>
#include "rawarray.h"
#include "parseargs.h"

typedef float ArrayBType;  // has to be double for Andreas' programs to work
typedef TArray3d<ArrayBType>  TImgArray;  // Spots array with clipping 

static TImgArray   Array1;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-out] -i image1  -o outputfile\n" << flush;
  exit(-1);
}


int _flip()
{
 int test = 1;
 return (0 ==  (* ((unsigned char *) & test)));
}

int omas_dbl_write(const char * fnam,   // Filename
		   int * res, // resolution int [4]
		   float * len, // length of dataset for scaling
		   double * data, // data to be written
		   int bFlt) // true -> save as float
{
 unsigned char * buf, * bbf, * ddt;
 FILE * fp;

 double * dpt;
 float  * fpt;

 unsigned char header[128];

 size_t flt = sizeof(float);
 size_t dbl = bFlt ? flt : sizeof(double);

 size_t siz = res[0] * res[1] * res[2] * res[3] * dbl;

 size_t cnt = siz / BUFSIZ + 1;
 size_t rst = siz % BUFSIZ;
 size_t num = BUFSIZ / dbl;

 long idx = cnt;
 int  flp = _flip();

 int jdx, kdx;

 if (! (fp=fopen(fnam,"wb")))
   return 0;

 if (bFlt)
 {
  fpt = (float  *) malloc(siz);
  dpt = (double *) data;
  ddt = (unsigned char *) fpt;


  idx = siz / flt;

  while(idx--) *fpt++ = (float) *dpt++;

 } else {

  ddt = (unsigned char *) data;
 }

 header[0] = res[2] / 256;
 header[1] = res[2] % 256;

 header[2] = res[1] / 256;
 header[3] = res[1] % 256;

 header[4] = res[0] / 256;
 header[5] = res[0] % 256;

 header[6] = res[3] / 256;
 header[7] = res[3] % 256;

 *((float *) (header + 8 + 0 * flt)) = len[2];
 *((float *) (header + 8 + 1 * flt)) = len[1];
 *((float *) (header + 8 + 2 * flt)) = len[0];
 *((float *) (header + 8 + 3 * flt)) = len[3];


 if (flp)
 {
  ddt = header + 8; bbf = header + 8 + 4*flt;

  jdx = 4; while(jdx--)
  {
   ddt += flt;
   kdx = flt; while(kdx--) *bbf++ = *(--ddt);
   ddt += flt;
  }

  memcpy(header + 8, header + 8 + 4*flt, 4*flt);
 }

 header[24] = bFlt ? 1 : 0;
 if (fwrite(header,1,128,fp) != 128) 
   cerr << "Error writing file, no space on disk \n";

 if (flp)
 {
  if (flp) buf = (unsigned char *) alloca(BUFSIZ);

  while(cnt--)
  {

  if (! cnt) num = rst / dbl;

  bbf = buf; jdx = num; while(jdx--)
  {
   ddt += dbl;
   kdx = dbl; while(kdx--) *bbf++ =  *(--ddt);
   ddt += dbl;
  }

  }

  if (! cnt) 
    {
      if (fwrite(buf,1,rst,fp) != rst) 
	cerr << "Error writing file, no space on disk \n";
      cout << " RSTSize : " << rst << "\n";
    }
  else {
    if (fwrite(buf,1,BUFSIZ,fp) != BUFSIZ) 
      cerr << "Error writing file, no space on disk \n";
    cout << " BufSize : " << BUFSIZ << "\n";
  }

 } else 
   {
     if (fwrite(ddt,1,siz,fp) != siz) 
      cerr << "Error writing file, no space on disk \n";
     cout << " size : " << siz << "\n";
   }


 if (bFlt) free(ddt);

 fclose(fp);
 return 1;
}

int omas_dbl_read (const char * fnam, int * res, float * len, double ** data)
{
 unsigned char * buf, * bbf, * ddt;

 double * dpt;
 float  * fpt;

 FILE * fp;

 unsigned char header[128];

 size_t dbl = sizeof(double);
 size_t flt = sizeof(float);

 size_t siz, cnt, rst, num;
 long idx;

 int jdx, kdx;
 int bFlt;

 int  flp = _flip();

 if ((fp = fopen(fnam,"rb")) == NULL) 
   return 0;

 fread( header, 128,1,fp);
 bFlt = header[24];

 if (bFlt) dbl = flt;

 res[2] = header[0]*256 + header[1];
 res[1] = header[2]*256 + header[3];
 res[0] = header[4]*256 + header[5];
 res[3] = header[6]*256 + header[7];

 siz = res[0] * res[1] * res[2] * res[3] * dbl;

 cnt = siz / BUFSIZ + 1;
 rst = siz % BUFSIZ;
 num = BUFSIZ / dbl;

 idx = cnt;

 if (flp)
 {
  ddt = header + 8; bbf = header + 8 + 4*flt;

  jdx = 4; while(jdx--)
  {
   ddt += flt;
   kdx = flt; while(kdx--) *bbf++ = *(--ddt);
   ddt += flt;
  }

  memcpy(header + 8, header + 8 + 4*flt, 4*flt);
 }

 len[2] = *((float *) (header + 8 + 0 * flt));
 len[1] = *((float *) (header + 8 + 1 * flt));
 len[0] = *((float *) (header + 8 + 2 * flt));
 len[3] = *((float *) (header + 8 + 3 * flt));

 *data = (double *) malloc(siz * (bFlt ? 2 : 1));
 ddt = (unsigned char *) (* data);

 if (flp)
 {
  buf = (unsigned char *) alloca(BUFSIZ);

  while(cnt--)
  {

  if (! cnt)
  {
   num = rst / dbl;
   fread( buf, rst,1,fp);

  } else 
   fread( buf, BUFSIZ,1,fp);

  bbf = buf; jdx = num; while(jdx--)
  {
   ddt += dbl;
   kdx = dbl; while(kdx--) *(--ddt) = *bbf++;
   ddt += dbl;
  }

  }

 } else 
   fread( ddt, siz,1,fp);

 fclose(fp);

 if (bFlt)
 {
  siz /= dbl;

  fpt = ((float *)  *data) + siz - 1;
  dpt = ((double *) *data) + siz - 1;

  while(siz--) *dpt-- = *fpt--;
 }

 return 1;
}



int main(int argc, char *argv[])
{ 

static int kflag=0, in=0;
static int I1SizeX=0,I1SizeY=0,I1SizeZ=0,elements=0,el;

string OFileName,I1FileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 int res[4];
 float len[4];

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-in",parg)) {in=1;continue;}
   if (readArg("-o", & OFileName, parg)) continue;
   if (readArg("-i", & I1FileName, parg)) continue; // image to expand
    usage(argv[0]);
  }

  if (OFileName=="" || I1FileName=="") 
    usage(argv[0]);

  if (!in)
    {
      elements=Array1.DLoad(kflag,I1FileName.c_str(),"Float",& I1SizeX,& I1SizeY,& I1SizeZ,0);
      double * data= (double *) malloc(I1SizeX*I1SizeY*I1SizeZ*elements*sizeof(double));
      double * dat=data;
      cout << "Data read\n";
      for (el=0;el< elements;el++)
	{
	  Array1.DLoad(kflag,I1FileName.c_str(),"Float",& I1SizeX,& I1SizeY,& I1SizeZ,el);
	  for (int z=0;z< I1SizeZ;z++)
	    for (int y=0;y< I1SizeY;y++)
	      for (int x=0;x< I1SizeX;x++)
		(* dat++) = Array1.Value(x,y,z);
	}
      res[0] = I1SizeX;
      res[1] = I1SizeY;
      res[2] = I1SizeZ;
      res[3] = elements;

      len[0] = I1SizeX;
      len[1] = I1SizeY;
      len[2] = I1SizeZ;
      len[3] = elements;

      cerr << "writing ... \n";
      if (! omas_dbl_write(OFileName.c_str(), res, len, data, 1)) // true -> save as float
      cerr << "Error writing file " << OFileName << "\n",exit(-1);
    }
  else  // read data
    {
      double * data;
      cerr << "reading ... \n";
      if (! omas_dbl_read (I1FileName.c_str(), res, len, & data))
	cerr << "Error reading file " << I1FileName << "\n",exit(-1);
      cerr << "Data read\n";
      ofstream of(OFileName.c_str());
      Array1.Resize(res[0],res[1],res[2]);
      double * dat = data;
      cout << "read : " << res[3] << " elements \n";
	for (el=0;el< res[3];el++)
	{
	  for (int z=0;z< res[2];z++)
	    for (int y=0;y< res[1];y++)
	      for (int x=0;x< res[0];x++)
		Array1.SetValue(x,y,z,(* dat++));
	  if (el == 0)
	    Array1.DHeader(kflag,of,res[3]);
	  Array1.Write(& of);
	  cout << "Wrote file : " << res[0] << " x " << res[1] << " x " << res[2] << "\n";
	}
	of.close();
    }

}
