		// -*- C++ -*- 
#ifndef asi_h
#define asi_h

/*   This file is part of a software package written by 
     Rainer Heintzmann
     Institute of Applied Optics and Information Processing
     Albert Ueberle Strasse 3-5
     69120 Heidelberg
     Tel.: ++49 (0) 6221  549264
     No garantee, whatsoever is given about functionallaty and  savety.
     No warranty is taken for any kind of damage it may cause.
     No support for it is provided !

     THIS IS NOT FREE SOFTWARE. All rights are reserved, and the software is only
     given to a limited number of persons for evaluation purpose !
     Please do not modify and/or redistribute this file or any other files, libraries
     object files or executables of this software package !
*/


#define AREAD(stream,Var) {len+= sizeof(Var);stream->read((char *) & (Var),sizeof(Var));}
// should this be "unsigned char *" above?? It used to be

void ReadASIHeader(ifstream * from, int & ValueDimX,int & ValueDimY,int & ValueDimZ)  // ASI spectral raw file
{
  unsigned long headersize=1024,i,len=0;
  unsigned char bdummy;
  unsigned short dummy;
  unsigned long cx,cy;
  unsigned short nch,nx,ny;

  AREAD(from,headersize);
  cerr << "ASI spectral file read, header size: ("<< headersize << ")\n";
  if (headersize != 1024)
    {
      cerr << "Unknown file-format, headersize wrong\n";
      exit(0);
    }
  for (i=0;i<26;i++) AREAD(from,dummy);
  AREAD(from,cx);
  AREAD(from,cy);
  cerr << "ASI file read, CCD size: ("<< cx << ", " << cy << ")\n";
  for (i=0;i<147-32;i++) AREAD(from,dummy);
  AREAD(from,nx);
  AREAD(from,ny);
  cerr << "len: "<< len << ", " << len/2 << "\n";
  while (len/2 < 156) AREAD(from,dummy);
  AREAD(from,nch);

  ValueDimX = nx; ValueDimY = ny; ValueDimZ = nch;

  cerr << "ASI file read, sizes: ("<< ValueDimX << ", " << ValueDimY << ", " << ValueDimZ << ")\n";

  for (i=len;i < headersize;i++)
    AREAD(from,bdummy);

}

unsigned short GetValue(ifstream * from, int pos)
{
  unsigned short val;
      int len=0;
  from->seekg(pos*sizeof(unsigned short));
    AREAD(from,val);
	return val;
}


void ReadSDIHeader(ifstream * from, int & ValueDimX,int & ValueDimY,int & ValueDimZ, int mode=0, int posx=-1, int posy=-1, int posz=-1, int headersize=-1)   // streaming SDI file (raw data)
{
  unsigned short nch,nx,ny;

  if (mode <= 0)
    {
      if (posx < 0) posx = 326;
      if (posy < 0) posy = 327;
      if (posz < 0) posz = 325;
      if (headersize < 0) headersize = 2*334;
    }
  else
    {
      if (posx < 0) posx = 224;
      if (posy < 0) posy = 225;
      if (posz < 0) posz = 231;
      if (headersize < 0) headersize = 2*325+1;
    }

  nx = GetValue(from, posx);
  ny = GetValue(from, posy);
  nch = GetValue(from, posz);
  from->seekg(headersize);

  ValueDimX = nx; ValueDimY = ny; ValueDimZ = nch;

  cerr << "SDI interferometer file read, sizes: ("<< ValueDimX << ", " << ValueDimY << ", " << ValueDimZ << ")\n";

}

void ReadRAWHeader(ifstream * from, int & ValueDimX,int & ValueDimY,int & ValueDimZ)
{
  int len=0,i;
  char bdummy;
  short dummy;
  short nch,nch_1,nx,ny,OrderMode;
  long coeff_wlen;
  short dlt_ch,LamOffset;
  float fcoeff_wlen,fdlt_ch,fLamOffset;

  AREAD(from,nch);
  AREAD(from,nch_1);
  AREAD(from,dummy);
  AREAD(from,nx);
  AREAD(from,dummy);
  AREAD(from,dummy);
  AREAD(from,ny);
  AREAD(from,dummy);
  AREAD(from,dummy);
  AREAD(from,dummy);
  AREAD(from,OrderMode);
  if (OrderMode != 0)
    { ValueDimX = nx; ValueDimY = ny; ValueDimZ = nch;}
  else
    { ValueDimX = nch; ValueDimY = nx; ValueDimZ = ny;}

  cerr << "SDM interferometer file read, sizes: ("<< ValueDimX << ", " << ValueDimY << ", " << ValueDimZ << ")\n";

  AREAD(from,coeff_wlen);
  for (i=0;i<4;i++) AREAD(from,dummy);
  AREAD(from,dlt_ch);
  AREAD(from,dummy);
  AREAD(from,LamOffset);

  AREAD(from,fcoeff_wlen);
  AREAD(from,fdlt_ch);
  AREAD(from,fLamOffset);

  cerr << "spectral conversion (nm) :"<< fcoeff_wlen << "/ ( channel + " << fdlt_ch << ") - " << fLamOffset << "\n";
  for (i=len;i < 256;i++)
    AREAD(from,bdummy);
  return;  // from is now pointing to first segment data !
}



#endif
