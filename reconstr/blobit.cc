// this program reads in a 3-D datafile and puts out a pov-ray blob by applying a threshold to the data

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
#include "vec.h"
#include "parseargs.h"
#include "rawarray.h"
#include "stat.h"

static const char * filename = "/tmp/start.raw";

typedef float ArrayBType;

typedef TArray3d<ArrayBType>    TArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
TArray3d<char>    LArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
TArray3d<int>    NumArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

static int genVRML=0;
static int nospix=0; // Do not interpolate positions
static double ZExpand=1.0,Zoom=1.0,XCenter=0.5,YCenter=0.5,ZCenter=0.5,blobs=1,polycol=0;
static double ImageDX=1.0,ImageDY=1.0,ImageDZ=1.0;
double smoothfac=2.5;
TArray Measured;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-t threshhold] [-z zscale] [-tf texturefilesname] -i inputfile -ao ascoutput \n" << flush;
  cerr << "[- bf begin-file] [-tf texture-files] [-ef end-file] \n";
  cerr <<  "Texturefiles will be extended by #Elem.txb for begin-file and #Elem.txe for end-file\n" << flush;
  exit(-1);
}
void CalcDataSize(double & XC, double & YC, double & ZC, double & XS, double & YS, double & ZS) // in case only vectors are given
{
   int y;
   ArrayBType I;
   Stat<ArrayBType> X,Y,Z,V;

   int DimY=Measured.GetSize(1);

   for (y=0;y<DimY;y++)
       {
	 I=Measured.Value(ClipSel,0,y,0);
	 if (I != 0)
	   {
	     X.Register(Measured.Value(ClipSel,1,y,0));
	     Y.Register(Measured.Value(ClipSel,2,y,0));
	     Z.Register(Measured.Value(ClipSel,3,y,0));
	     V.Register(Measured.Value(ClipSel,4,y,0));
	   }
       }
   XS=X.Max()-X.Min();
   YS=Y.Max()-Y.Min();
   ZS=Z.Max()-Z.Min();
   XC=X.Min()+XS/2.0;
   YC=Y.Min()+YS/2.0;
   ZC=Z.Min()+ZS/2.0;
   
   cerr << "Sizes were : " << XS << " x " << YS << " x " << ZS << "\n";
   cerr << "Centers were : " << XC << " x " << YC << " x " << ZC << "\n";
}

/// generates the spheres from vectorial input
void VecSpheres(ostream * ostrm,double XS,double XC,double YS,double YC,double ZS,double ZC,int elem)
{
   int y;
   ArrayBType I,V;
   CalcFloatType px,py,pz;
   int DimY=Measured.GetSize(1);
   double cr,cg,cb,alpha=0,betha=0,posx,posy,posz;

   for (y=0;y<DimY;y++)
       {
	 I=Measured.Value(ClipSel,0,y,0);
	 px=Measured.Value(ClipSel,1,y,0);
	 py=Measured.Value(ClipSel,2,y,0);
	 pz=Measured.Value(ClipSel,3,y,0);  // changed to -Measured ...
	 V=Measured.Value(ClipSel,4,y,0);
	 if (I != 0)
	   {
	     alpha += M_PI/2/2.862;    // computer cyclic colors
	     betha += M_PI/2/3.8182;
	     cb = fabs(cos(betha));
	     cr = fabs(sin(betha)*sin(alpha));
	     cg = fabs(sin(betha)*cos(alpha));

	     posx = -Zoom*((double(px)-ImageDX*XCenter)/(ImageDX));
	     // "," << -Zoom*ZExpand*((double(z)-dz1+dz2-DimZ*ZCenter)/DimX)  << 
	     posy = -Zoom*ZExpand*((double(pz+1.0)-ImageDZ*ZCenter)/(ImageDX)) ;
	     posz = Zoom*((double(py)-ImageDY*YCenter)/(ImageDX));

	     if (genVRML)
	       {
		 (* ostrm) <<  "Transform { translation " << posx << "  " << posy << "  " << posz << "  children [ USE MyObject"<<elem<<" ] },\n";
	       } 
	     else
	       {
		 (* ostrm) << "sphere { <" << posx << "," << posy  << "," << posz << ">,SR";

		 /* (* ostrm) << "sphere { <" << Zoom*(X-XC-XS*(XCenter-0.5))/XS << "," <<  
		    Zoom*ZExpand*(Z-ZC-ZS*(ZCenter-0.5))/YS  << "," << 
		    Zoom*(Y-YC-YS*(YCenter-0.5))/ZS<< ">,SR"; */
		 if (blobs)
		   {
		     (* ostrm) << ", 1 "; 
		   }
		 if (polycol)
		   (* ostrm) <<  "  pigment{rgb <" <<cr << ","<< cg << "," << cb <<">}";

		 (* ostrm) << " } \n"; // make z point upward  DimX is used for scaling
	       }
	   }
       }
}

int inRange(int x, int y, int z)
{
  int DimX=Measured.GetSize(0),
      DimY=Measured.GetSize(1),
      DimZ=Measured.GetSize(2);

    if (x < 0) return 0;
    if (y < 0) return 0;
    if (z < 0) return 0;
    if (x >= DimX) return 0;
    if (y >= DimY) return 0;
    if (z >= DimZ) return 0;
    return 1;
}

int GenPlanes (ostream * ostrm,double Threshhold) // generates triangularization
{
  int DimX=Measured.GetSize(0),
    DimY=Measured.GetSize(1),
    DimZ=Measured.GetSize(2),
    x,y,z,dx,dy,dz,dx2,dy2,dz2,count=0,num=1;
  double px,py,pz, posx,posy,posz;

  LArray.Resize(DimX,DimY,DimZ);
  NumArray.Resize(DimX,DimY,DimZ);

  (* ostrm) << "geometry IndexedFaceSet {\n ";
  (* ostrm) << "coord Coordinate { point [ ";

  Measured.ClipAt(Threshhold,Threshhold/2.0); // Clip and replace by 0 -> smother surfaces !
  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	if (Measured.Check6OuterThreshVoxel(Threshhold,x,y,z))
	  {
	    LArray.SetValue(x,y,z,1);  // mark voxels
	    NumArray.SetValue(x,y,z,num++); // give each corner a number
	    px=x,py=y,pz=z;
	    if (! nospix)
	      Measured.RangedCOM(px,py,pz,1,1,1);

	    
	    px = x+smoothfac*(px-x); py = y+smoothfac*(py-y); pz = z+smoothfac*(pz-z);
	    // px = (px+x)/2.0; py = (py+y)/2.0; pz = (pz+z)/2.0;

	    posx = -Zoom*((double(px)-DimX*XCenter)/(DimX));
	    // "," << -Zoom*ZExpand*((double(z)-dz1+dz2-DimZ*ZCenter)/DimX)  << 
	    posy = -Zoom*ZExpand*((double(pz+1.0)-DimZ*ZCenter)/(DimX));
	    posz =  Zoom*((double(py)-DimY*YCenter)/(DimX));

	    (* ostrm) << -posx << "  " << posy << "  " << posz << ",\n" ;
	  }
  (* ostrm) << "] }\n solid FALSE\n creaseAngle 0.5\n coordIndex [";

  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	if (LArray.Value(x,y,z))
	  {
	    LArray.SetValue(x,y,z,0);  // delete voxel, will be completely processed
	    for (dz=-1;dz<=1;dz++)
	      for (dy=-1;dy<=1;dy++)
		for (dx=-1;dx<=1;dx++)
		  {
		    if (dx==0)
		      if (dy==0)
			if (dz==0)
			  continue;

		    if (inRange(x+dx,y+dy,z+dz))
			LArray.SetValue(x+dx,y+dy,z+dz,LArray.Value(x+dx,y+dy,z+dz)*2);  // mark neigbours
		  }

	    for (dz=-1;dz<=1;dz++)
	      for (dy=-1;dy<=1;dy++)
		for (dx=-1;dx<=1;dx++)
		  if (inRange(x+dx,y+dy,z+dz))  // All 3rd members of this 2nd voxel are processed
		    if (LArray.Value(ClipSel,x+dx,y+dy,z+dz) == 2)  // is a marked neigbour
		      {
			if (dx==0)
			  if (dy==0)
			    if (dz==0)
			      continue;

			for (dz2=-1;dz2<=1;dz2++)
			  for (dy2=-1;dy2<=1;dy2++)
			    for (dx2=-1;dx2<=1;dx2++)
			      {
				if (dx2==0)
				  if (dy2==0)
				    if (dz2==0)
				      continue;
				if (LArray.Value(ClipSel,x+dx+dx2,y+dy+dy2,z+dz+dz2) == 2)  // found a complete triangle
				  {
				    (* ostrm) << NumArray.Value(x,y,z)-1 << ", " << 
				      NumArray.Value(x+dx,y+dy,z+dz)-1 << ", " << 
				      NumArray.Value(x+dx+dx2,y+dy+dy2,z+dz+dz2)-1 << ", -1,\n" ;
				    count++;
				  }
			      }
			// inRange is checked above !
			LArray.SetValue(x+dx,y+dy,z+dz,LArray.Value(x+dx,y+dy,z+dz)/2);  // mark neigbours
		      }

	  }
  (* ostrm) << "] }\n ";
  return count;

}

int GenSpheres(ostream * ostrm,double Threshhold,int elem) // writes spheres to ostream return #speres
{
  int x,y,z;
  int count=0;
  CalcFloatType px,py,pz,posx,posy,posz;
  int DimX=Measured.GetSize(0),
      DimY=Measured.GetSize(1),
      DimZ=Measured.GetSize(2);

  for (z=0;z<DimZ;z++)
    for (y=0;y<DimY;y++)
      for (x=0;x<DimX;x++)
	if (Measured.Check6OuterThreshVoxel(Threshhold,x,y,z))
	  {
	    /* 
	    dz1=dz2=0.0;
	    myval=Measured.Value(ClipSel,x,y,z);
	     if ((val=Measured.Value(ClipSel,x,y,z-1)) < Threshhold)
	       dz1=1.0/( (Threshhold-val)/(myval-Threshhold) +1.0);   // linear intersection

	     if ((val=Measured.Value(ClipSel,x,y,z+1)) < Threshhold)
	       dz2=1.0/( (Threshhold-val)/(myval-Threshhold) +1.0);   // linear intersection
	    */
	    px=x,py=y,pz=z;

	    if (! nospix)
	      Measured.RangedCOM(px,py,pz,1,1,1);

	     posx = -Zoom*((double(px)-DimX*XCenter)/(DimX));
	     // "," << -Zoom*ZExpand*((double(z)-dz1+dz2-DimZ*ZCenter)/DimX)  << 
	     posy = -Zoom*ZExpand*((double(pz+1.0)-DimZ*ZCenter)/(DimX));
	     posz =  Zoom*((double(py)-DimY*YCenter)/(DimX));

	     if (genVRML)
	       {
		 (* ostrm) <<  "Transform { translation " << posx << "  " << posy << "  " << posz << "  children [ USE MyObject"<<elem<<" ] },\n";
	       }
	     else
	       {
		 (* ostrm) << "sphere { <" << posx << "," << posy  << "," << posz << ">,SR";
		 if (blobs)
		   (* ostrm) << ", 1 "; 
		 
		 (* ostrm) << "} \n"; 
	       }
	     count++;
	   }
  return count;
}

int main(int argc, char *argv[])
{  
  int khoros=0,count=0,Elem=0,i,VectorInput=0,autocenter=0,faces=0,ScaleBar=1000,relative=0;
  Vector COM(0.0,0.0,0.0);
  double XS,XC,YS,YC,ZS,ZC;  // only needed for vector input

static int MEASUREDSizeX=256; // 64; 
static int MEASUREDSizeY=256; // 64;
static int MEASUREDSizeZ=32;  // 32;

static double Threshhold=1.0,SphereRadius=2.2,PixelSizeXY=78.0,PixelSizeZ=162.0; 

string InFileName, AscFileName, TextureFiles,TextureFile,
  BeginFile, EndFile,Unit, TitleX,TitleY,TitleZ;


Unit="1 µm";
TitleX="X-axis";
TitleY="Y-axis";
TitleZ="Z-axis";

BeginFile="/usr/local/KhorosInstall/goodies/data/part1.pov";
EndFile="/usr/local/KhorosInstall/goodies/data/part2.pov";

char ** parg= & argv[1];

while (* parg)
  {
    if (readArg("-t", & Threshhold, parg)) continue;
    if (readArg("-rel",parg)) {relative=1;continue;}
    if (readArg("-k",parg)) {khoros=1;continue;}
    if (readArg("-noblob",parg)) {blobs=0;continue;}
    if (readArg("-polycol",parg)) {polycol=1;continue;}
    if (readArg("-z",& ZExpand,parg)) continue;
    if (readArg("-pixXY",& PixelSizeXY,parg)) continue;
    if (readArg("-pixZ",& PixelSizeZ,parg)) continue;
    if (readArg("-ac",parg)) {autocenter=1;continue;}
    if (readArg("-zoom",& Zoom,parg)) continue;
    if (readArg("-dx",& ImageDX,parg)) continue;  // Dimensions of virtual image, when vectors are used
    if (readArg("-dy",& ImageDY,parg)) continue;
    if (readArg("-dz",& ImageDZ,parg)) continue;
    if (readArg("-mx",& XCenter,parg)) continue;
    if (readArg("-my",& YCenter,parg)) continue;
    if (readArg("-mz",& ZCenter,parg)) continue;
    if (readArg("-sb",& ScaleBar,parg)) continue; // basic units to form the length of the scale-bar
    if (readArg("-unit",Unit,parg)) continue; // text to be written at the scale bar
    if (readArg("-i",InFileName, parg)) continue;
    if (readArg("-ao",AscFileName, parg)) continue;
    if (readArg("-tf",TextureFiles, parg)) continue;
    if (readArg("-titX", TitleX, parg)) continue;
    if (readArg("-titY",TitleY, parg)) continue;
    if (readArg("-titZ",TitleZ, parg)) continue;
    if (readArg("-bf",BeginFile, parg)) continue;
    if (readArg("-ef",EndFile, parg)) continue;
    if (readArg("-sr",& SphereRadius, parg)) continue;
    if (readArg("-smoothfac",& smoothfac, parg)) continue;
    if (readArg("-v",parg)) {VectorInput=1;continue;}
    if (readArg("-vrml",parg)) {genVRML=1;continue;}
    if (readArg("-faces",parg)) {faces=1;continue;}
    if (readArg("-nospix",parg)) {nospix=1;continue;}
 
    usage(argv[0]);
  }

  if (genVRML)
  {XCenter=0.0;YCenter=0.0;ZCenter=0.0;}  // ignore the image centers
  
if (InFileName=="" || AscFileName=="") usage(argv[0]);

  cerr << "// Threshhold = " << Threshhold << ", now reading file : " << filename << "\n" << flush;

  ofstream * ostrm= new ofstream(AscFileName.c_str());


  if (genVRML)
    (* ostrm) << "#VRML V2.0 utf8\n Group { children [ \n";

  Elem=1;
  for (i=0;i<Elem;i++)
    {
      Elem=Measured.DLoad(khoros,InFileName.c_str(),"Float",& MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ,i);

      if (relative)
	Measured.NormalizeToMax(100.0);

      if (!genVRML)
	(* ostrm) << "// Element Nr. : " << i << "\n";
      if (i==0)
	{
	  if (! VectorInput)
	    {
	      ImageDX=MEASUREDSizeX;
	      ImageDY=MEASUREDSizeY;
	      ImageDZ=MEASUREDSizeZ;
	    }
	  if (!genVRML)                                  
	    {
          int mpos = Unit.find("\\m");
          if (mpos != string::npos) Unit.replace(mpos,2,"µ");

	      (* ostrm) << "#declare ScaleText=\""<< Unit<<"\" ;\n"  ;  // ????
	      (* ostrm) << "#declare Scale=<-1,"<< PixelSizeZ/PixelSizeXY<<",-1> ;\n"  ;  // ????
	      (* ostrm) << "#declare ScaleImage=<1,"<< ImageDZ<<"*"<<PixelSizeZ<<"/("<<PixelSizeXY<<"*"<<ImageDX<<"),1> ;\n"  ;
	      (* ostrm) << "#declare ScaleBar="<<1.0/(PixelSizeXY/float(ScaleBar))/ImageDX << ";\n"  ;  // ????
	      
	      (* ostrm) << "// blobit file : " << filename << "\n";
	      (* ostrm) << "#include \"" << BeginFile <<"\"\n";


	      (* ostrm) << "// Size : " << 
		ImageDX << " x " << ImageDY << " x " << ImageDZ << 
		",  Threshhold : " << Threshhold << "\n";
	    }
	  else
	    {
              double width = SphereRadius * double(fabs(Zoom)) /double(ImageDX);
	      (* ostrm) << "  DEF MYScene Group { children [\n";
	      (* ostrm) << "      DEF MyCoordX Shape { appearance DEF White Appearance { material Material {diffuseColor 1 1 1} }\n";
              (* ostrm) <<"geometry Box { size " << 1.0 << "  " << width << "  " << width << " }\n},\n";
	      (* ostrm) << "      DEF MyCoordY Shape { appearance DEF White Appearance { material Material {diffuseColor 1 1 1} }\n";
              (* ostrm) <<"geometry Box { size " << width << "  " << 1.0 << "  " << width << " }\n},\n";
	      (* ostrm) << "      DEF MyCoordZ Shape { appearance DEF White Appearance { material Material {diffuseColor 1 1 1} }\n";
              (* ostrm) <<"geometry Box { size " << width << "  " << width << "  " << 1.0 << " }\n},\n";
		for (int j=0;j < Elem;j++)
		{
		double cr=(j%3==0),cg=(j%3==1),cb=(j%3==2);
	      (* ostrm) << "      DEF MyObject"<<j<<" Shape { appearance DEF White Appearance { material Material {diffuseColor "<<cr<<" "<<cg<<" "<<cb<<"} }\n"; 
	      if (faces)
		{
		  count = GenPlanes(ostrm,Threshhold);
		  cerr << "// number of planes generated : " << count << "\n" << flush;
		}
	      else
		{
		  (* ostrm) <<"geometry Box { size " << width << "  " << width << "  " << width << " }\n";
		}

              (* ostrm) <<"}, \n"; 
		}
                
	      (* ostrm) <<  "  ] },\n";
	      (* ostrm) << "Viewpoint { position    0 0 2 orientation 0 0 1 0 description \"DefaultView\" }, \n";

	      (* ostrm) << "Transform { scale  -1 " << PixelSizeZ/PixelSizeXY << " -1  children [ \n"  ;  // ????

              (* ostrm) << "Separator {\n Translation { translation -1.5 0 0 } \n Material { emissiveColor 0.8 0.5 0.5 }\n";
              (* ostrm) << "Transform {scale -1 1 1 children [\n";
              (* ostrm) << "FontStyle {size 0.2 family SANS style BOLD } AsciiText { string \""<< TitleX <<"\" } ] }\n},\n";
  	      (* ostrm) << "Transform { translation  -2 0 0 children [ USE MyCoordX ]},\n"  ;  // ????

              (* ostrm) << "Separator {\n Translation { translation 0 -1.5 0 }\n Rotation { rotation 0 0 1 1.571 }\n Material { emissiveColor 0.5 0.8 0.5 }\n";
              (* ostrm) << "Transform {scale -1 1 1 children [\n";
              (* ostrm) << "FontStyle {size 0.2 family SANS style BOLD } AsciiText { string \""<< TitleY <<"\" } ] }\n},\n";
              (* ostrm) << "Transform { translation  0 -2 0 children [ USE MyCoordY ]},\n"  ;  // ????

              (* ostrm) << "Separator {\n Translation { translation 0 0 1.5 }\n Rotation { rotation 0 0 1 1.571 }\n Rotation { rotation 0 1 0 1.571 }\nRotation { rotation 1 0 0 1.571 }\nMaterial { emissiveColor 0.5 0.5 0.8 }\n";
              (* ostrm) << "Transform {scale -1 1 1 children [\n";
              (* ostrm) << "FontStyle {size 0.2 family SANS style BOLD } AsciiText { string \"" << TitleZ << "\" } ] }\n},\n";
	      (* ostrm) << "Transform { translation  0 0 2 children [ USE MyCoordZ ]},\n"  ;  // ????

	    } 
	}

	  if (!genVRML)
	    {
	      (* ostrm) << "#declare SR="<< SphereRadius<<" * " << double(fabs(Zoom)) ;
	      // if (!VectorInput)
	      (* ostrm) << " / " << double(ImageDX);

	      (* ostrm) << ";   // Sphere Radius\n";

	      if (TextureFiles != "")
		{
                     char dum[10];
		  sprintf(dum,"%d.txb",i);
		  TextureFile=TextureFiles+dum;
		  (* ostrm) << "#include \"" << TextureFile <<"\"\n";
		}
	      else
		blobs=0;
	    }

      if (VectorInput)  // place spheres with radius r
	{
	  // XC = 0,YC=0,ZC=0;
	  // XS =-1.0,YS=-1.0,ZS=1.0;
	  // if ((i == 0) && autocenter)
	  //   CalcDataSize(XC,YC,ZC,XS,YS,ZS); 

	  VecSpheres(ostrm,XS,XC,YS,YC,ZS,ZC,i);
	  count=MEASUREDSizeY;
	}
      else // do your own labelling
	{
	  Measured.ClippedCOM(Threshhold,&COM);
	  if (!genVRML)
	    (* ostrm) << "// COM at : " << COM.comp(0) << "  " << COM.comp(1) << "  " << COM.comp(2) << "\n";
	  
	  if (!genVRML || (genVRML && ! faces))
	    {
	      count = GenSpheres(ostrm,Threshhold,i);
	      cerr << "// number of sphere generated : " << count << "\n" << flush;
	    }
	}


      if (TextureFiles != "")
	{
            char dum[10];
            sprintf(dum,"%d.txe",i);
	   TextureFile=TextureFiles+dum;
	 (* ostrm) << "#include \"" << TextureFile <<"\"\n";
	}
    }

  if (genVRML)
    (* ostrm) << "] }\n" <<  
      "   ] } ";
  else
    (* ostrm) << "#include \"" << EndFile <<"\"\n";

    ostrm->close();
}

