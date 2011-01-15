// applies any function to a number of datasets (0 to 6 datasets)

// THIS PROGRAM CAN BE COMPILED WITH OR WITHOUT THE RUN FLAG  !
// If RUN flag is present, it is meant for excecution otherwise for compilation


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

typedef float ArrayBType;
typedef TArray3d<ArrayBType> TImgArray;

#define MAXIMG   5
#define MAXIMGS "5"
static TImgArray InputImg[MAXIMG], MaskImg, OutputImg;

/*
export int INPUTSizeX[MAXIMG];  // These is the standart size, if raw data is used
export int INPUTSizeY[MAXIMG];  // These sizes can be used in the C-program
export int INPUTSizeZ[MAXIMG];
export int INPUTSizeE[MAXIMG];
export int MaskSizeX;
export int MaskSizeY;
export int MaskSizeZ;
export int MaskSizeE;
export bool NoMask;
*/
int INPUTSizeX[MAXIMG];  // These is the standart size, if raw data is used
int INPUTSizeY[MAXIMG];  // These sizes can be used in the C-program
int INPUTSizeZ[MAXIMG];
int INPUTSizeE[MAXIMG];
int MaskSizeX;
int MaskSizeY;
int MaskSizeZ;
int MaskSizeE;
bool NoMask;

#ifdef RUN
void ProcessVoxel(int x, int y, int z, int c, float w, float h, float d, float e, float & out, float & i1, float & i2, float & i3, float & i4, float & i5);
#endif

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-i1  ... -i6 inputfiles] [-o outputfile] \n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 
static int OUTPUTSizeX=0;  // These is the standart size, if raw data is used
static int OUTPUTSizeY=0;
static int OUTPUTSizeZ=0;
static int OUTPUTSizeE=0;

INPUTSizeX[0] = 256;
INPUTSizeY[0] = 256;
INPUTSizeZ[0] = 16;
INPUTSizeE[0] = 1;
MaskSizeE = 1;
NoMask = true;

bool kflag=false, iter1=false, ClearZ=false,ClearE=false,ClearY=false;

string INPUTFileName[MAXIMG], MaskFileName, OUTPUTFileName,
          OUTPUTProgramTxt,OUTPUTProgramObj,OUTPUTProgramExe, Program,OldProgram;

char ** parg= & argv[1];
printCall(argc,argv);

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=true;continue;}
   if (readArg("-iter1",parg)) {iter1=true;continue;}
   if (readArg("-clearE",parg)) {ClearE=true;continue;}
   if (readArg("-clearZ",parg)) {ClearZ=true;continue;}
   if (readArg("-clearY",parg)) {ClearY=true;continue;}
   if (readArg("-i1", INPUTFileName[0], parg)) continue;
   if (readArg("-i2", INPUTFileName[1], parg)) continue;
   if (readArg("-i3", INPUTFileName[2], parg)) continue;
   if (readArg("-i4", INPUTFileName[3], parg)) continue;
   if (readArg("-i5", INPUTFileName[4], parg)) continue;
   if (readArg("-mask", MaskFileName, parg)) {NoMask = false;continue;}
   if (readArg("-prog", Program, parg)) continue;  // this contains the little program to be executed for every pixel
   if (readArg("-o", OUTPUTFileName, parg)) continue;
   if (readArg("-pout", OUTPUTProgramTxt, parg)) continue;  // Filename of output program
   if (readArg("-pobj", OUTPUTProgramObj, parg)) continue;  // Filename of object file to link with
   if (readArg("-pexe", OUTPUTProgramExe, parg)) continue;  // Filename of output program
   if (readArg("-sX",& OUTPUTSizeX, parg)) continue;
   if (readArg("-sY",& OUTPUTSizeY,parg)) continue;
   if (readArg("-sZ",& OUTPUTSizeZ,parg)) continue;
   if (readArg("-se",& OUTPUTSizeE,parg)) continue;
    usage(argv[0]);
  }

// cout << INPUTFileName+" usw"+" und noch mehr" << " ... worked\n";
#ifndef RUN  
// create a subrouting by adding a header and footer:
 Program = "#include <math.h>\nextern int INPUTSizeX["+string(MAXIMGS)+
                     "];\nextern int INPUTSizeY["+string(MAXIMGS)+"];\nextern int INPUTSizeZ["+string(MAXIMGS)+
                     "];\nextern int INPUTSizeE["+string(MAXIMGS)+
                     "];\nextern bool NoMask;\n"+
                     "void ProcessVoxel(int x, int y, int z, int e, float w, float h, float d, float re, float & out, float & i1, float & i2, float & i3, float & i4, float & i5)\n{\n"+Program;
 Program += "\nreturn;\n}\n";  

  ifstream oldp(OUTPUTProgramTxt.c_str());
  string txt;
  if (oldp)
    {
              while (! oldp.eof())
       {
                    getline(oldp,txt);
                   OldProgram += txt;
                   if (! oldp.eof()) OldProgram += "\n";
       }
              oldp.close();
    }

    // cout << "Previous Program:\n" << OldProgram;
    cout << "New Program:\n" << Program;

  int sucess = 0;
  if (OldProgram == Program)
  {
        cout << "Programs are identical, no need to recompile\n";
  }
  else
  {
        cout << "Programs was changed, recompiling ...\n";
          ofstream to(OUTPUTProgramTxt.c_str());
          to << Program;  // wrties the program to the file
          to.close();
        
          system(("echo \"cp "+OUTPUTProgramTxt+" "+OUTPUTProgramTxt+".cc;g++ -c "  +  OUTPUTProgramTxt+".cc\"").c_str());
          system(("echo \"rm "+  OUTPUTProgramTxt+"\".cc").c_str());
          system(("echo \"g++ -o "+OUTPUTProgramExe+"  " + OUTPUTProgramTxt + ".o "+ OUTPUTProgramObj+"\"").c_str());
          system(("cp "+OUTPUTProgramTxt+" "+OUTPUTProgramTxt+".cc;\n").c_str());
          sucess = system(("g++ -o "+ OUTPUTProgramTxt +".o -c "  +  OUTPUTProgramTxt+".cc").c_str());
          // cout << "Exit code: " << sucess << "\n" << flush;
          if (sucess != 0)
            {
              system(("echo \"ERROR\" >> "+OUTPUTProgramTxt+"\n").c_str());
              cerr << "Compilation did not succeed!\n"; exit(-1);}
          system(("rm "+  OUTPUTProgramTxt+".cc").c_str());
          sucess=system(("g++ -o "+OUTPUTProgramExe+"  " + OUTPUTProgramTxt + ".o "+ OUTPUTProgramObj).c_str());
          // cout << "Exit code: " << sucess << "\n" << flush;
          if (sucess != 0)
            { cerr << "Linking did not succeed!\n"; exit(-1);}
  }


  string callstring=OUTPUTProgramExe;
  for (int i=1;i< argc;i++)
  {
        if (string(argv[i])== string("-prog"))            // scip the program argument
      {i++; continue;}
       callstring += string(" ") + argv[i];
  }
    
    system(("echo \"running program: "+callstring+"\"").c_str());    // Run the program
    return system(callstring.c_str());    // Run the program   

#else    // not we have to analyse the data

static int IterSizeX=0;  // These is the standart size, if raw data is used
static int IterSizeY=0;
static int IterSizeZ=0;
static int IterSizeE=0;
  
  ofstream to(OUTPUTFileName.c_str());
      
  if (! to )
    {
      cerr << "Couldn't open file " << OUTPUTFileName << " for writing !!\n" << flush;
      exit(-1);
    }

      int x,y,z,j;
      float w,h,d,e;
      float val[MAXIMG];
      for (int j=0;j< MAXIMG;j++)
          val[j]=0.0;

  IterSizeE=1;
  for (j=0;j< MAXIMG;j++)
    INPUTSizeE[j] = 1;

  for (int i=0;i<IterSizeE;i++)
  {
    if (INPUTFileName[0] != "")
                INPUTSizeE[0]=InputImg[0].DLoad(kflag,INPUTFileName[0].c_str(),"Float",
			  & INPUTSizeX[0],& INPUTSizeY[0],& INPUTSizeZ[0],i % INPUTSizeE[0]);
    if (MaskFileName != "")
                MaskSizeE=MaskImg.DLoad(kflag,MaskFileName.c_str(),"Float",
			  & MaskSizeX,& MaskSizeY,& MaskSizeZ,i % MaskSizeE);

    if (i == 0)
        {
           if (OUTPUTSizeX == 0)
                 OUTPUTSizeX = INPUTSizeX[0];
           if (OUTPUTSizeY == 0)
                OUTPUTSizeY = INPUTSizeY[0];
           if (OUTPUTSizeZ == 0)
               OUTPUTSizeZ = INPUTSizeZ[0];
           if (OUTPUTSizeE == 0)
               OUTPUTSizeE = INPUTSizeE[0];

           if (iter1 && INPUTFileName[0] != "")
           {IterSizeX=INPUTSizeX[0];IterSizeY=INPUTSizeY[0];IterSizeZ=INPUTSizeZ[0];IterSizeE=INPUTSizeE[0];}
              else
           {IterSizeX=OUTPUTSizeX;IterSizeY=OUTPUTSizeY;IterSizeZ=OUTPUTSizeZ;IterSizeE=OUTPUTSizeE;}

           OutputImg.Resize(OUTPUTSizeX,OUTPUTSizeY,OUTPUTSizeZ);
           OutputImg.DHeader(kflag,to,OUTPUTSizeE);
         }

         
    // cout << "Iterating over " << IterSizeE << " Elements \n";
               
    // cout << "Loaded file nr. 0\n";
    for (j=1;j<MAXIMG;j++)
      {
        // cout << "Loading file nr. " << j <<" : " << INPUTFileName[j] << "\n";
         if (INPUTFileName[j] != "")
                INPUTSizeE[j]=InputImg[j].DLoad(kflag,INPUTFileName[j].c_str(),"Float",
			  & INPUTSizeX[j],& INPUTSizeY[j],& INPUTSizeZ[j],i % INPUTSizeE[j]);
        // cout << "Loaded file nr. " << j <<" : " << INPUTFileName[j] << "\n" << flush;
     }

      float CenterX = OUTPUTSizeX/2;  // Integer divisions !
      float CenterY = OUTPUTSizeY/2;
      float CenterZ = OUTPUTSizeZ/2;
      float CenterE = OUTPUTSizeE/2;
      float WX = OUTPUTSizeX/2;  // Integer divisions !
      float WY = OUTPUTSizeY/2;
      float WZ = OUTPUTSizeZ/2;
      float WE = OUTPUTSizeE/2;
      if (WX < 1) WX = 1;
      if (WY < 1) WY = 1;
      if (WZ < 1) WZ = 1;
      if (WE < 1) WE = 1;

      if (ClearE)
        OutputImg.Clear();

      for (z=0;z<IterSizeZ;z++)
      {
        if (ClearZ)
          OutputImg.Clear();

         for (y=0;y<IterSizeY;y++)
         {
           if (ClearY)
             OutputImg.Clear();
           for (x=0;x<IterSizeX;x++)
            if (NoMask || MaskImg.Value(x % MaskSizeX,y % MaskSizeY,z % MaskSizeZ) > 0)
             {
                    w = (x-CenterX) / WX;
                    h = (y-CenterY) / WY;
                    d = (z-CenterZ) / WZ;
                    e = (i-CenterE) / WE;
                    // .SetValue(x,y,z,out);
                    for (j=0;j<MAXIMG;j++)
                      if (INPUTFileName[j] != "")
                            val[j] = InputImg[j].Value(x % INPUTSizeX[j],y % INPUTSizeY[j],z % INPUTSizeZ[j]);

                    ProcessVoxel(x, y,  z, i, w,  h,  d,  e, OutputImg(x%OUTPUTSizeX,y%OUTPUTSizeY,z%OUTPUTSizeZ),  val[0],  val[1],  val[2],  val[3], val[4]);
               }
           }
      }
    if (i >= (IterSizeE-OUTPUTSizeE))
      OutputImg.Write(& to);
  }

  to.close();

#endif
}
