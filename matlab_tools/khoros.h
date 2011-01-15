		// -*- C++ -*- 
#ifndef khoros_h
#define khoros_h

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


#include <iostream>
#include <fstream>
#include <string>
#include <complex>

using namespace std;

typedef complex<float> float_complex;
typedef complex<double> double_complex;

const char * TypeString(unsigned char arg)
{
 static const char rt[] = "Unsigned Byte";
 arg=0;
 return rt;
}

const char * TypeString(int arg)
{
 static const char rt[] = "Integer";
 arg=0;
 return rt;
}

const char * TypeString(long arg)
{
 static const char rt[] = "Long";
 arg=0;
 return rt;
}

const char * TypeString(float arg)
{
 static const char rt[] = "Float";
 arg=0;
 return rt;
}

const char * TypeString(double arg)
{
 static const char rt[] = "Double";
 arg=0;
 return rt;
}

const char * TypeString(complex<float> arg)
{
 static const char rt[] = "Complex";
 arg=0;
 return rt;
}

const char * TypeString(unsigned short arg)
{
 static const char rt[] = "Unsigned Short";
 arg=0;
 return rt;
}

const char * TypeString(complex<double> arg)
{
 static const char rt[] = "Double Complex";
 arg=0;
 return rt;
}

#define WRITE(stream,Var) stream->write((char *) & (Var),sizeof(Var))
// used to be unsigned char *

void WRITEINT(ofstream * to,int Var)
{
     char dummy;
     dummy =  char(Var & 255);
     WRITE(to,dummy);
     dummy =  char((Var>>8) & 255);
     WRITE(to,dummy);
     dummy =  char((Var>>16) & 255);
     WRITE(to,dummy);
     dummy =  char((Var>>24) & 255);
     WRITE(to,dummy);
}

void WriteKhorosHeader(ofstream * to, const char * Date, const char * TypeString, 
		       int ValueDimX, int ValueDimY, int ValueDimZ, int ValueDimE, int ValueDimT = 1)
{

  if (ValueDimX <= 0 ) cerr << "Warning ! Datasize X was set to zero ! Setting to one \n", ValueDimX=1;
  if (ValueDimY <= 0 ) cerr << "Warning ! Datasize Y was set to zero ! Setting to one \n", ValueDimY=1;
  if (ValueDimZ <= 0 ) cerr << "Warning ! Datasize Z was set to zero ! Setting to one \n", ValueDimZ=1;
  if (ValueDimE <= 0 ) cerr << "Warning ! Datasize E was set to zero ! Setting to one \n", ValueDimE=1;

  const char MagicNumVer[6]={1,3,19,94,0,2};
#ifndef SGI
  const char MachineType=65;   // This is PC format
#else
  const char MachineType=33;  // This is SGI format
#endif
  const int NumSets=1;
  const int NumBlocks=2;

  const char ObjAttrName[] ="";
  const int   ObjSegNr=2;

  const char SegAttrName[]="date";
  const int SegAttrNum=1;

  const int DateDim=1;
  const char DateType[] ="String";
  const char EOA[] ="<>";            // end of attribute tag

  const char Seg2AttrName[]="locationGrid";
  const int Seg2AttrNum=1;

  const int LocationDim=1;
  const char LocationType[] ="Integer";
  const int Location=0;
  const char EOA2[] ="<>";            // end of attribute tag

  const char ValueAttrName[] ="value";
  const int ValueAttrNum=0;

  const int ValueDim=5;
  // const int  ValueDimT = 1;
  const int  ValueOrder[] ={2,3,4,5,6};
  const int FixedDim=-1,FixedIndex=-1;

  WRITE(to,MagicNumVer);
  WRITE(to,MachineType);
  WRITE(to,NumSets);
  WRITE(to,NumBlocks);

  WRITE(to,ObjAttrName);
  WRITE(to,ObjSegNr);

  WRITE(to,SegAttrName);
  WRITE(to,SegAttrNum);
  WRITE(to,DateDim);
  WRITE(to,DateType);
  to->write((char *) Date,strlen(Date)+1);
  WRITE(to,EOA);

  WRITE(to,Seg2AttrName);
  WRITE(to,Seg2AttrNum);
  WRITE(to,LocationDim);
  WRITE(to,LocationType);
  WRITE(to,Location);
  WRITE(to,EOA2);

  WRITE(to,ValueAttrName);
  WRITE(to,ValueAttrNum);
  WRITE(to,ValueDim);
  to->write((char *) TypeString,strlen(TypeString)+1);
  WRITE(to,ValueDimX);
  WRITE(to,ValueDimY);
  WRITE(to,ValueDimZ);
  WRITE(to,ValueDimT);
  WRITE(to,ValueDimE);
  for (int i=0; i <  5; i++)
    WRITE(to,ValueOrder[i]);

  WRITE(to,FixedDim);
  WRITE(to,FixedIndex);
}

#define READ(stream,Var) stream->read((char *) & (Var),sizeof(Var))
// used to be unsigned char *

#define READSTR(stream,Str) {int ip=0; char x=1; for(ip=0;x != 0;ip++) stream->read(& x,sizeof(char)), Str[ip]=x; }

int ReadKhorosHeader(ifstream * from,char * TypeString, // will be filled in
                      int & ValueDimX,int & ValueDimY,int & ValueDimZ,int & ValueDimE)
{
  int i,j,k, ValueDimT=1;
  char MagicNumVer[6]={1,3,19,94,0,2};
#ifndef SGI
  const char MachineType=65;   // This is PC format
#else
  const char MachineType=33;  // This is SGI format
#endif
  int NumSets=1;
  int NumBlocks=2;

  char ObjAttrName = 0;
  int   ObjSegNr=2;

  char SegAttrName[100];
  int SegAttrNum=1,count=1;

  int SegDim=1;
  char SegType[100];
  char EOA[] ="<>";            // end of attribute tag

  char AString[100];
  int  AnInt;
  unsigned long  AULong;
  long  ALong;
  float  AFloat;
  double  ADouble;

  int  SegDims[100];
  int  SegOrders[100];

  int FixedDim=-1,FixedIndex=-1;

  READ(from,MagicNumVer);
  if (MagicNumVer[0] != 1 || MagicNumVer[1] != 3) cerr << "Error: No Khoros file format\n", myexit(-1);
  char MachTyp;
  READ(from,MachTyp);
  if (MachTyp != MachineType)
   { cerr << "This file is from a different machine (type " << int(MachTyp) <<" expecting " << int(MachineType) << ") ! Cannot read. Use native khoros glyph (e.g. extract) first\n";
           myexit (-1);
   }
  READ(from,NumSets);
  if (NumSets != 1) cerr << "Error: Can only read 1-Set Kdf data !\n", myexit(-1);
  READ(from,NumBlocks);

// So called "header" is done
// Now read Object Attribute Block

  READ(from,ObjAttrName);
  if (ObjAttrName != 0) cerr << "Error : Unexpected char as ObjAttrName in Kdf !\n", myexit(-1);
  READ(from,ObjSegNr);

// Now Read Atrribute Block and Attributes for each Segment

  for (i=0;i<ObjSegNr+1;i++)
    {
      READSTR(from,SegAttrName);
#ifdef K_PH
      cout << "KHeader : Attr. " << SegAttrName << " found\n";
#endif
      READ(from,SegAttrNum);
#ifdef K_PH
      cout << "KHeader : SegAttrNum " << SegAttrNum << " found\n";
#endif
       
      READ(from,SegDim);
      READSTR(from,SegType);
#ifdef K_PH
      cout << "KHeader :   Dim " << SegDim << " Type  " << SegType << " found\n";
#endif
      if (SegAttrNum==0)
	count = SegAttrNum+1;
      else
	count = SegAttrNum;

      for (k=0;k<count;k++) // What happenes with more than one attribute ??
	{
	  if (SegDim>1)   // Read Sizes and Orders
	    {
	      for (j=0;j<SegDim;j++)
		{
		  READ(from,SegDims[j]);
#ifdef K_PH
		  cout << "KHeader :   Dim " << j << " Size  " << SegDims[j] << "\n";
#endif
		}
	      for (j=0;j<SegDim;j++)
		{
		  READ(from,SegOrders[j]);
#ifdef K_PH
		  cout << "KHeader :   Dim " << j << " Order  " << SegOrders[j] << "\n";
#endif
		  if (SegOrders[j] != j+2)
		    cerr << "Warning : Orders of segment nr. " << j << " is anormal !\n";
		}
	      READ(from,FixedDim);
	      READ(from,FixedIndex);


	      if (strcmp(SegAttrName,"value")==0)
		{
		  ValueDimX=SegDims[0];
		  ValueDimY=SegDims[1];
		  ValueDimZ=SegDims[2];
		  ValueDimE=SegDims[4];
		  ValueDimT=SegDims[3];
		  strcpy(TypeString,SegType);
#ifdef K_PH
		  cout << "values copied.\n";
#endif
		}
       	      continue;
	    }
	  else
	  if (SegDim==1) {
	      if (strcmp(SegType,"String")==0)
		{
		  READSTR(from,AString);
#ifdef K_PH
		  cout << "KHeader :   String :  " << AString << "\n";  
#endif
		  continue;
		}
	      if (strcmp(SegType,"Integer")==0)
		{
		  READ(from,AnInt);
#ifdef K_PH
		  cout << "KHeader :   Int :  " << AnInt << "\n";  
#endif
		  continue;
		}
	      if (strcmp(SegType,"Long")==0)
		{
		  READ(from,ALong);
#ifdef K_PH
		  cout << "KHeader :   Long :  " << ALong << "\n";  
#endif
		  continue;
		}
	      if (strcmp(SegType,"Float")==0)
		{
		  READ(from,AFloat);
#ifdef K_PH
		  cout << "KHeader :   Float :  " << AFloat << "\n";  
#endif
		  continue;
		}
	      if (strcmp(SegType,"Double")==0)
		{
		  READ(from,ADouble);
#ifdef K_PH
		  cout << "KHeader :   Double :  " << ADouble << "\n";  
#endif
		  continue;
		}
	      if (strcmp(SegType,"Unsigned Long")==0)
		{
		  READ(from,AULong);
#ifdef K_PH
		  cout << "KHeader :   Unsigned Long :  " << AULong << "\n";  
#endif
		  continue;
		}
	      cerr << "Unknown 1 Dim Segmenttype " << SegType << "!\n",myexit(-1);
	    }
	}
  if (strcmp(SegAttrName,"value")!=0)
    {
      READSTR(from,EOA);
      if (strcmp(EOA,"<>") != 0) cerr << "Error : EOA not found in KDF\n";
    }
  else
    {
      float AFloat;
      int i;
      if (SegAttrNum != 0) 
	for (i=0;i<9;i++) READ(from,AFloat);  // Dirty Method for scipping something unknown, which causes trouble
    }
  }

return ValueDimT;  // from is now pointing to first segment data !
}



#endif
