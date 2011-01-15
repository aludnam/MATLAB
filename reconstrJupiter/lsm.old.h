		// -*- C++ -*- 
#ifndef lsm_h
#define lsm_h

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

#include <tiffio.h>

template <class vclass>
void VREAD(ifstream * stream,vclass & Var) {stream->read((char *) & (Var),sizeof(Var));cout << "read " << ": " << Var << "\n";}

//#define VREAD(stream,Var) {stream->read((char *) & (Var),sizeof(Var)); cout << "read "##Var << ": " << Var << "\n";}

const int MAXELEM = 40;
const unsigned short TIF_IMAGEWIDTH=256,
    TIF_IMAGELENGHT=257,
    TIF_BITSPERSAMPLE=258,
    TIF_COMPRESSION=259,
    TIF_STRIPOFFSETS=273,
    TIF_CZLSMINFO=34412;

class LSMParser {
  public:
    ifstream * from;
    unsigned short sdummy;
    unsigned long offset;
    unsigned short numtags;

    // Directory Tag:
    unsigned short dummy,dtag,dtype,compress;
    unsigned long dlng,dval;
    unsigned long width,height,bits,dataoffsets[MAXELEM],infooffset;

    // Infor Structure
    unsigned long magic;
    long isize,SizeX,SizeY,SizeZ,SizeE,SizeT,DatType,ThumX,ThumY,datsize;
    double VoxX,VoxY,VoxZ;
    
  LSMParser(ifstream * fr) : from(fr)
  {
  }

  void ParseHeader() {
    VREAD(from,sdummy);
    VREAD(from,sdummy);
    VREAD(from,offset);
    ReadDir();
  }

  void ParseInfo() {
    cout << "Parsing Info\n"<<flush;
  	from->seekg(infooffset,ios::beg);
    VREAD(from,magic);
    VREAD(from,isize);
    VREAD(from,SizeX);
    VREAD(from,SizeY);
    VREAD(from,SizeZ);
    VREAD(from,SizeE);
    VREAD(from,SizeT);
    VREAD(from,DatType);
    VREAD(from,ThumX);
    VREAD(from,ThumY);
    VREAD(from,VoxX);
    VREAD(from,VoxY);
    VREAD(from,VoxZ);
    datsize = SizeX*SizeY*SizeZ*DatType;
    cout << "Sizes: " << SizeX << ", " << SizeY << ", " << SizeZ << ", " << SizeE << "\n" << flush;
    }

  void ScanTable(unsigned long tableoffset, unsigned long tablesize) {
    cout << "Scanning Offset-Table\n"<<flush;
  	from->seekg(tableoffset,ios::beg);
    for (int i=0;i < tablesize;i++)
    {
      VREAD(from,dataoffsets[i]);
      }
    }
    
  bool ReadDir() {
    unsigned long tableoffset=0;
    unsigned long tablesize=0;
    cout << "Read Dir\n"<<flush;
    if (offset == 0) return false;
  	from->seekg(offset,ios::beg);
    VREAD(from,numtags);

    infooffset=0;
    for (int i=0;i < numtags;i++)
    {
    cout << "New Tag\n"<<flush;
      VREAD(from,dtag);
      VREAD(from,dtype);
      VREAD(from,dlng);
      if (dtag == TIF_COMPRESSION)
      {
      VREAD(from,compress);
      VREAD(from,dummy);
        }
      else
      VREAD(from,dval);
      switch(dtag) {
        /*case TIF_IMAGEWIDTH:
        width = dval;
        break;
        case TIF_IMAGELENGHT:
        height = dval;
        break;
        case TIF_BITSPERSAMPLE:
        bits = dval;
        break;*/
        case TIF_STRIPOFFSETS:
        tableoffset = dval;
        tablesize = dlng;
        break;
        case TIF_CZLSMINFO:
        infooffset = dval;
        default:
        break;        
      }
    }
    
    VREAD(from,offset);
    if (infooffset != 0) ParseInfo();
    if (tableoffset != 0) ScanTable(tableoffset, tablesize);
    if (compress != 1)
      {cerr << "FATAL Error! File is compressed, cannot uncompress!\n";
      // exit (-1);
      }
    return true;  
  }

  bool ReadSlice(void * data, int element=0) {   // The offset for the next data is allready stored
    cout << "Reading Slice Nr " << element << ", offset is " << dataoffsets[element] << ", size: " << datsize << "\n"<<flush;
    if (dataoffsets[element] != 0)
      from->seekg(dataoffsets[element],ios::beg);
    else
      {
      cerr << "ERROR: No more data found\n";
      exit(-1);
      }

    from->read((char *) data,datsize);
 

    return ReadDir();  // position to next directory
  }
};




#endif
