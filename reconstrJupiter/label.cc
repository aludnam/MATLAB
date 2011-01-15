// label : this program reads in a 3-D datafile and labels individual regions and determines their COM

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
#include "vec.h"
#include "parseargs.h"
#include "tarrayclipf.h"
#include "label.h"

using namespace std;

typedef float ArrayBType;

typedef Dyn3dArray(ArrayBType,ClipperatorNoClip)  TArray;  // Clipping is needed
typedef Dyn3dArray(int,ClipperatorNoClip)  TArrayL;  // Clipping is needed
typedef Dyn2dArray(ArrayBType,ClipperatorNoClip)  TArrayVec;  // Clipping is needed


static TArray Measured;
static TArrayL Label;

int MinVol=0;      // minimal volume of spots to be couted and displayed
float MinInt=0.0;  // minimal Intensity to be counted
float MinBorderX=0.0,MaxBorderX=0.0,MinBorderY=0.0,MaxBorderY=0.0,MinBorderZ=0.0,MaxBorderZ=0.0;

int NoSubClip=0;   // subtract before computing center of mass !

const int MAXSPOT=10000;   // maximum amount of spots

label<MAXSPOT> * data1, * data2;  // will be allocated later on

#define usedLabel (data1->farr)
#define storedLabel (data2->farr)
#define ActLabelNr (data1->spots)
#define StoredLabels (data2->spots)

// int ActLabelNr=0;
// static float tarr usedLabel;         // saves # Pixels per label in [0] and some more data COMS ...
// int StoredLabels=0;                  // Number of labels stored in storedLabels
// static float tarr storedLabel;       // stores the information of Element=0 for comparison

/// marks labels to be joined to others
int join(int a, int b)
{
  int i,n,j,m;
  
  i=b;
  while((n=(int) usedLabel[i-1].Join) != 0)  // find lowest of b and store in i
    i=n;
  j=a;
  while((m=(int) usedLabel[j-1].Join) != 0)  // find lowest of a and store in j
    j=m;

  if (j > i)  // (a >= b) the most propable case
    {
      usedLabel[j-1].Join = i;
      return i;
    }
  else if (i > j)
    {  
      usedLabel[i-1].Join = j;
      return j;
    }
  else  // i==j -> No joining needed, because regions are allready connected
    return i;
}

/// calculates all important data for sets of labelled voxels
void calcIt(bool applyConstraints)
{
  int i,Vol;
  float inten,Vxx,Vyy,Vxy,tmp1,tmp2,Elow,Ehigh;

  for (i=0;i<ActLabelNr;i++)
    {
      Vol=usedLabel[i].Vol;
      inten = usedLabel[i].SumI;
      if ((! applyConstraints) || ((Vol >= MinVol) && (inten >= MinInt)))
	{
	  // COMS and VARS still contain sums
	  usedLabel[i].COMX /= inten;
	  usedLabel[i].COMY /= inten;
	  usedLabel[i].COMZ /= inten;
	  usedLabel[i].VARX = (usedLabel[i].VARX/inten)-usedLabel[i].COMX*usedLabel[i].COMX; // inten+Vol*ClipVal
	  usedLabel[i].VARY = (usedLabel[i].VARY/inten)-usedLabel[i].COMY*usedLabel[i].COMY;
	  usedLabel[i].VARZ = (usedLabel[i].VARZ/inten)-usedLabel[i].COMZ*usedLabel[i].COMZ;
	  usedLabel[i].VARXY = (usedLabel[i].VARXY/inten)-usedLabel[i].COMX*usedLabel[i].COMY; // inten+Vol*ClipVal
	  // usedLabel[i].MeanI = (inten+ClipVal)/usedLabel[i].Vol;
	  usedLabel[i].MeanI = inten/usedLabel[i].Vol;

                 // perform the normalized diagonalization
                Vxx =     usedLabel[i].VARX;
                Vyy =     usedLabel[i].VARY;
                Vxy =    -  usedLabel[i].VARXY;
                tmp1 = 4*Vxy*Vxy +(Vxx-Vyy)*(Vxx-Vyy);
                if (tmp1 <  0.0)
                        tmp1 = 0.0;
                else
                        tmp1 = sqrt(tmp1);
                
                tmp2 = 4*Vxy*Vxy - (Vxx-Vyy)*(Vyy-Vxx + tmp1);
                if (tmp2 <  0.0)
                        tmp2 = 0.0;
                else      tmp2 = sqrt(tmp2);
                
                Elow =   0.5*(Vxx+Vyy-tmp1);
                Ehigh =  0.5*(Vxx+Vyy+tmp1);
                usedLabel[i].Ellipt = (Ehigh-Elow)/(Ehigh+Elow);
                usedLabel[i].EigenLow = Ehigh;
                usedLabel[i].EigenHigh = Elow;

                // cout << "Eigenvalues: " << Elow <<", "<< Ehigh<<"\n";
                // cout <<"Vxx: " << Vxx << ", Vyy: " << Vyy << ", Vxy: " << Vxy <<", Tmp1 : " << tmp1 << ", tmp2: "<<tmp2<<"\n";
                if (Vxx + Vyy <= 0.0)
            {
                        usedLabel[i].ExcentX = 0.0;
                        usedLabel[i].ExcentY = 0.0;
            }
                 else
            {
                          if (tmp1 == 0.0)  // this means there is no direction (Entartet)
                {
                                usedLabel[i].ExcentX = 0.0;            
                                usedLabel[i].ExcentY = 0.0;
                }
                        else
            {
                         if (tmp2 == 0.0)
              {
                            usedLabel[i].ExcentY = 0.0;
                            usedLabel[i].ExcentX = 1.0;
              }
                          else
              {
                            usedLabel[i].ExcentY = - tmp2 / sqrt(2.0)/(Vxx+Vyy);
                            usedLabel[i].ExcentX = Vxy*sqrt(2.0)*tmp1/(Vxx+Vyy)/tmp2;
              }
            }}

	  if ((usedLabel[i].COMX < MinBorderX) ||
	      (usedLabel[i].COMX > MaxBorderX) ||
	      (usedLabel[i].COMY < MinBorderY) ||
	      (usedLabel[i].COMY > MaxBorderY) ||
	      (usedLabel[i].COMZ < MinBorderZ) ||
	      (usedLabel[i].COMZ > MaxBorderZ) )
	    {
	      usedLabel[i].Nr =0;
	      usedLabel[i].Vol=0;  // delete this label
	    }
	}
      else
	{
	  usedLabel[i].Nr =0;
	  usedLabel[i].Vol=0;  // delete this label
	}
    }
    
}

void CrunchLabels(void)  // crunches all the labels down to a contiguous set and reassigns them in the image.
{
   int nextFree=0,i;
   for (i=0;i<ActLabelNr;i++)
   {
      if (usedLabel[i].Vol >  0)
     {
            usedLabel[i].Join = nextFree + 1;
            nextFree++;
     }
       else   usedLabel[i].Join = 0;
   }
          
  int x,y,z,val,newlabel;
  int MEASUREDSizeX=Label.GetSize(0);
  int MEASUREDSizeY=Label.GetSize(1);
  int MEASUREDSizeZ=Label.GetSize(2);

  for (z=0;z<MEASUREDSizeZ;z++)
    for (y=0;y<MEASUREDSizeY;y++)
       for (x=0;x<MEASUREDSizeX;x++)
     {
             val =  (int) Label.Value(x,y,z);
             if(val > 0)
             {
                                newlabel =   usedLabel[val-1].Join; 
                                if (newlabel == 0)
                                        newlabel = -1;
                                 Label.SetValue(x,y,z, newlabel);
             }
     }
   nextFree=0;
   for (i=0;i<ActLabelNr;i++)
   {
      if (usedLabel[i].Vol >  0)
     {
           usedLabel[nextFree] = usedLabel[i];
           if (nextFree +1 != usedLabel[nextFree].Join)
                cout << "internal ERROR while reassigning labels, Nr. " << nextFree + 1 << " is " << usedLabel[nextFree].Join << "\n";
           usedLabel[nextFree].Nr = nextFree + 1;
           usedLabel[nextFree].Join = 0;
           nextFree++;
     }
   }
     for (i=nextFree;i<ActLabelNr;i++)
        data1->Clear(i);   // clear this usedLabel[i]

     ActLabelNr = nextFree;
}

void computeExcent()  // computes the 2D  (!!) excentricity vector
{
// Equations:
// Vx  = Sum ( [(Px_i^2 -Py_i^2)/sqrt(Px_i^2 + Py_i^2) + 1.0]Intensity_i)  / Sum(Intensity_i)
// Vy  = Sum ( [1.0 - (Px_i^2 -Py_i^2)/sqrt(Px_i^2 + Py_i^2)] Intensity_i)  / Sum(Intensity_i)
// px_i and py_i being difference vectors to center of mass

  int MEASUREDSizeX=Label.GetSize(0);
  int MEASUREDSizeY=Label.GetSize(1);
  int MEASUREDSizeZ=Label.GetSize(2);
  float inten,px,py,dist;
  int val;
  
  for (int z=0;z<MEASUREDSizeZ;z++)
    for (int y=0;y<MEASUREDSizeY;y++)
       for (int x=0;x<MEASUREDSizeX;x++)
	   {
	     if((val=(int) Label.Value(x,y,z)) > 0)
             {
                              inten = Measured.Value(x,y,z);            // This is suboptimal but shall be tested
                              px = x - usedLabel[val-1].COMX;      // This requires tha the COMs are allready calculated
                              py = y - usedLabel[val-1].COMY;
                              dist = sqrt(px*px+py*py);
                              if ( dist >  0.0)
                {
                                usedLabel[val-1].ExcentX += inten*(dist+(px*px-py*py)/dist);
                                usedLabel[val-1].ExcentY += inten*(dist-(px*px-py*py)/dist);        
                }
             }
           }

 for (int i=0;i<ActLabelNr;i++)
    {
          usedLabel[i].ExcentX /= usedLabel[i].SumI;
          usedLabel[i].ExcentY /= usedLabel[i].SumI;
    }
}

/// joins the marked voxels together
void dojoinCOM(float ClipVal)
{
  int x,y,z,val,n,nr=1,surf=0;
  int MEASUREDSizeX=Label.GetSize(0);
  int MEASUREDSizeY=Label.GetSize(1);
  int MEASUREDSizeZ=Label.GetSize(2);
  float inten;

  for (z=0;z<MEASUREDSizeZ;z++)
    for (y=0;y<MEASUREDSizeY;y++)
       for (x=0;x<MEASUREDSizeX;x++)
	   {
	     if((val=(int) Label.Value(x,y,z)) > 0)
	       {
		 while((n=(int) usedLabel[val-1].Join) != 0)  // find lowest
		   val=n;
		 // if (NoSubClip)
		   inten = Measured.Value(x,y,z);            // This is suboptimal but shall be tested
		 // else
		 // inten = Measured.Value(x,y,z) - ClipVal;  // should be bigger than 0
		 
		 usedLabel[val-1].COMX += inten*x;
		 usedLabel[val-1].COMY += inten*y;
		 usedLabel[val-1].COMZ += inten*z;
		 usedLabel[val-1].SumI += inten;
		 usedLabel[val-1].Vol ++;
		 surf=0;
		 if ((Label.Value(ClipSel,x-1,y,z)) == 0)
		   usedLabel[val-1].SurfaceX ++,surf=1;
		 if ((Label.Value(ClipSel,x+1,y,z)) == 0)
		   usedLabel[val-1].SurfaceX ++,surf=1;
		 if ((Label.Value(ClipSel,x,y-1,z)) == 0)
		   usedLabel[val-1].SurfaceY ++,surf=1;
		 if ((Label.Value(ClipSel,x,y+1,z)) == 0)
		   usedLabel[val-1].SurfaceY ++,surf=1;
		 if ((Label.Value(ClipSel,x,y,z-1)) == 0)
		   usedLabel[val-1].SurfaceZ ++,surf=1;
		 if ((Label.Value(ClipSel,x,y,z+1)) == 0)
		   usedLabel[val-1].SurfaceZ ++,surf=1;
		 if (surf)
		   usedLabel[val-1].Surface++;
		 usedLabel[val-1].VARX += ((inten)*x*x);  // inten + ClipVal
		 usedLabel[val-1].VARY += ((inten)*y*y);
		 usedLabel[val-1].VARZ += ((inten)*z*z);
		 usedLabel[val-1].VARXY += ((inten)*x*y);
		 if (usedLabel[val-1].Nr == 0)
		   usedLabel[val-1].Nr = nr++;

		 Label.SetValue(x,y,z,val); // usedLabel[val-1].Nr);
	       }
	   }
}


int NextLabel(void)
{
  ++ActLabelNr;
  if (ActLabelNr >= MAXSPOT)
    cerr << "Fatal error ! Too many (>" << MAXSPOT<<") spots for array ! Try to increase threshhold ! \n",exit(-1);
  return ActLabelNr;
}

void COMS(ostream * ostrm)   // prints the label information
{
  int i;
  const float fac=2.35482; // == sqrt(8*ln(2));

  Vector COM(3);

  for (i=0;i< ActLabelNr;i++)
    if (usedLabel[i].Vol > 0)
      {
	cout << usedLabel[i].Nr << "   ";
	if (ostrm)  (* ostrm) << usedLabel[i].Nr << "   ";

	// to have FWHM : VAR * sqrt(8*ln(2)) == VAR * 2.35482
	
	cout << usedLabel[i].COMX << "  " << usedLabel[i].COMY <<"  " << usedLabel[i].COMZ << 
	  "    " << usedLabel[i].Vol << "   " << usedLabel[i].MeanI << "    " << usedLabel[i].SumI << 
	  "    " << fac*sqrt(usedLabel[i].VARX) << "  " << fac*sqrt(usedLabel[i].VARY) <<"  " << fac*sqrt(usedLabel[i].VARZ) << "    " << usedLabel[i].VARXY << "    "  << usedLabel[i].ExcentX << "  " << usedLabel[i].ExcentY << "  " <<  usedLabel[i].Ellipt << "\n";

	if (ostrm)  (* ostrm)  << usedLabel[i].COMX << "  " << usedLabel[i].COMY <<"  " << usedLabel[i].COMZ << 
	  "    " << usedLabel[i].Vol << "   " << usedLabel[i].MeanI << "    " << usedLabel[i].SumI << 
	  "    " << fac*sqrt(usedLabel[i].VARX) << "  " << fac*sqrt(usedLabel[i].VARY) <<"  " << fac*sqrt(usedLabel[i].VARZ) << "    " << usedLabel[i].VARXY << "    " << usedLabel[i].ExcentX << "  " << usedLabel[i].ExcentY << "  " <<  usedLabel[i].Ellipt << "\n";

      }
  cout << "D First-Last = " << data1->FirstLastDist() << "\n";
}



void usage(char * filename)
{
  cerr << "usage: " << filename << " -i inputfile [-k] [-m] [-t threshhold] [-lo labelimage] [-ao ascii-output] [-SX VoxSizeX] [-SY VoxSizeY] [-SZ VoxSizeZ] [-MD MinDist] [-MV MinVol] [-MB MinBorderDist]\n" << flush;
  cerr << "-k : use Khoros Data format for in- and output\n";
  cerr << "-m : subtract minimum value in array before further evaluation\n";
  cerr << "threshold : minimum value for a pixel to be labeled\n";
  cerr << "labelimage : name of labeled output image (if given)\n";
  cerr << "ascii-output : labelresults are written into this file (if given) and to stdout\n";
  cerr << "[-SX VoxelSizeX] : Size of a voxel in X\n";
  exit(-1);
}

int main(int argc, char *argv[])
{  
  int x,y,z,labelwith,newlabelwith,hasathresh=0;
  int Count=0;
  Vector COM(0.0,0.0,0.0);

  int khoros=0,submin=0;
  string InFileName, LblFileName,AscFileName,BinFileName;

static int MEASUREDSizeX=256; // 256; 
static int MEASUREDSizeY=256; // 256;
static int MEASUREDSizeZ=22;  // 32;

bool labelFirstElement=false;
    
static double Threshhold,RelThreshhold=0.5;  // 8.5

char ** parg= & argv[1];

float PixX=1,PixY=1,PixZ=1,MaxDist=0.0,MinBorder=0.0;

argc=0;  // prevent warning

while (* parg)
  {
    if (readArg("-t", & RelThreshhold, parg)) continue;    // relative to max and min
    if (readArg("-at",& Threshhold, parg)) {hasathresh=1;continue;}
    if (readArg("-labelFirst",parg)) {labelFirstElement=true;continue;}
    if (readArg("-k",parg)) {khoros=1;continue;}
    if (readArg("-m",parg)) {submin=1;continue;} // subtract actual minimum
    if (readArg("-nsc",parg)) {NoSubClip=1;continue;} // do not subtract any clipping values before COM evaluation
    if (readArg("-i",  InFileName, parg)) continue;
    if (readArg("-lo", LblFileName, parg)) continue;
    if (readArg("-ao", AscFileName, parg)) continue;
    if (readArg("-vo", BinFileName, parg)) continue;
    if (readArg("-PX", & PixX, parg)) continue;
    if (readArg("-PY", & PixY, parg)) continue;
    if (readArg("-PZ", & PixZ, parg)) continue;
    if (readArg("-MD", & MaxDist, parg)) continue;
    if (readArg("-MB", & MinBorder, parg)) continue;
    if (readArg("-MV", & MinVol, parg)) continue;
    if (readArg("-MI", & MinInt, parg)) continue;
    if (readArg("-MBX", & MinBorderX, parg)) continue;
    if (readArg("-MBY", & MinBorderY, parg)) continue;
    if (readArg("-MBZ", & MinBorderZ, parg)) continue;

    usage(argv[0]);
  }


if (InFileName=="") usage(argv[0]);

  cerr << "// RelThreshhold = " << RelThreshhold << ", now reading file : " << InFileName << "\n" << flush;
  cerr << "// Pix-X : " << PixX <<" Pix-Y : " << PixY <<" Pix-Z : " << PixZ << "\n" << flush;

double max,min;

int Elem=0,Elems=1;


ofstream * ostrm=0;
if (AscFileName!="") ostrm= new ofstream(AscFileName.c_str());

ofstream * oi=0;
TArrayVec *VecArray=0;

 if (BinFileName!="")
   {
     oi=new ofstream(BinFileName.c_str());
     if (! oi)
       cerr << "Error ! Couldnt open " << BinFileName << " for writing !! \n",exit(-1);
   }

ofstream * loi=0;
 if (LblFileName!="")
   {
     loi=new ofstream(LblFileName.c_str());
     if (! loi)
       cerr << "Error ! Couldnt open " << LblFileName << " for writing !! \n",exit(-1);
   }

 
for (Elem=0;Elem<Elems;Elem++)
  {
  if (Elem==1 && ! labelFirstElement)
    data1->TransferData(data2);

  if (Elem > 0 && ! labelFirstElement)
    ActLabelNr=0;

  Elems=Measured.DLoad(khoros,InFileName.c_str(),"Float",& MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ,Elem);

  MaxBorderX=MEASUREDSizeX-MinBorderX;
  MaxBorderY=MEASUREDSizeY-MinBorderY;
  MaxBorderZ=MEASUREDSizeZ-MinBorderZ;

  if (Elem==0 || ! labelFirstElement)
     Label.Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);  // will only be rezised when necessary
  
  max = Measured.Maximum();
  min = Measured.Minimum();

  if (submin)
      Measured.Sub(min);

  if (Elem==0)
    {
      data1 = new label<MAXSPOT> (PixX,PixY,PixZ,MaxDist,MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ,MinBorder,MinVol,MinInt,Threshhold);
      data2 = new label<MAXSPOT> (PixX,PixY,PixZ,MaxDist,MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ,MinBorder,MinVol,MinInt,Threshhold); // Threshhold not CORRECT !!!!

      if (ostrm)  (* ostrm) << "// label file : " << InFileName << " , Size : " << 
		  MEASUREDSizeX  << " x " << MEASUREDSizeY << " x " << MEASUREDSizeZ << "\n";
      cout << "// label file : " << InFileName << " , Size : " << 
		  MEASUREDSizeX  << " x " << MEASUREDSizeY << " x " << MEASUREDSizeZ << "\n";
    }
      
  if (! hasathresh)
    {
      Threshhold = RelThreshhold*(max-min)+min;  // Value in Image
      // Threshhold = RelThreshhold*max;
      cout << "using relative Threshhold : " << RelThreshhold << "\n";
      cout << "Max : " << max << " Min : " << min << " Threshhold : " << Threshhold << " in original data \n";
      if (submin) Threshhold -= min;   // NEW !! In this way the threshold matches the original data !
    }

  Measured.ClippedCOM(Threshhold,&COM);

  if (ostrm)
    {
     (* ostrm) <<"Element : " << Elem << ",  Threshhold : " << Threshhold << 
		" Max : " << max << " Min : " << min << "\n";
     (* ostrm) << "// COM at : " << COM.comp(0) << "  " << COM.comp(1) << "  " << COM.comp(2) << "\n";
    }

  cout <<"Element : " << Elem << ",  Threshhold : " << Threshhold << 
    " Max : " << max << " Min : " << min << "\n";
  cout << "// COM at : " << COM.comp(0) << "  " << COM.comp(1) << "  " << COM.comp(2) << "\n" << flush;

  if (Elem == 0 || ! labelFirstElement)
  {
    data1->Clear();
    Label.Set(0);

  for (z=0;z<MEASUREDSizeZ;z++)
    for (y=0;y<MEASUREDSizeY;y++)
       for (x=0;x<MEASUREDSizeX;x++)
	   {
	     if(Measured.Value(x,y,z)>= Threshhold)
	       {
		 labelwith=0;
		 if (x>0)
		   {
		     labelwith = (int) Label.Value(x-1,y,z);
		   }
		 if (y>0)
		   {
		     if((newlabelwith= (int) Label.Value(x,y-1,z)) != 0)
		       if (labelwith == 0)
			 labelwith=newlabelwith;
		       else
			 if (newlabelwith!=labelwith) labelwith=join(labelwith,newlabelwith);
		   }
		 if (z>0)
		   {
		     if((newlabelwith= (int) Label.Value(x,y,z-1)) != 0)
		       if (labelwith == 0)
			 labelwith=newlabelwith;
		       else
			 if (newlabelwith!=labelwith) labelwith=join(labelwith,newlabelwith);
		   }

		 if (labelwith) 
		   {
		     Label.SetValue(x,y,z,labelwith);
		     // usedLabel[labelwith-1].Vol +=1;       // count Pixels for label
		   }
		 else 
		   {
		     Label.SetValue(x,y,z,NextLabel());  // usedLabel is set to 1
		     // cout << "New Label : " << ActLabelNr << " at " << x << ", " << y << ", " << z << "\n";
		   }
	       }
	     // else { // Label.SetValue(x,y,z,0);} // is allready done

	   }

  // cout << "MaxLabels found : " << ActLabelNr << " \n";
  cout << "// Label Nr, Cliped COM X (Pixel) , Y (Pixel) , Z (Pixel) ,  Volume [Voxel], mean, itegrated,  FWHM X [Gauss,Pixel] ,  FWHM Y ,FWHM Z, Cross Variance XY, Excentricity X, Excentricity Y, Excentricity total\n";

  dojoinCOM(Threshhold);
  calcIt(true);  // calculates all important label information, and invalidates spots outside the criteria
  CrunchLabels();  // crunches all the labels down to a contiguous set and reassigns them in the image.
 }
 else
 {
     data1->ClearStat();
     dojoinCOM(0);
     calcIt(false);  // calculates all important label information
 }

  //  This old version is not correct: computeExcent();  // requires COMS to be evaluated allready
  COMS(ostrm);                  // just prints the information on the screen

  if (Elem > 0 && ! labelFirstElement)
    {
      cout << " Distances to Element 0 : \n";
      data1->calcshifts(data2);
    }

  if (Elem == 0 &&  loi != 0)
  {
        Label.DHeader(khoros,*loi,Elems);
  }

  if (LblFileName!="") 
       Label.Write(loi);

  if (Elem == 0 && oi != 0)
    {
        Count = 0;
        for (int i=0;i<ActLabelNr;i++)
	if (usedLabel[i].Nr != 0)
	  Count++;
        // VecArray = new TArrayVec(9,ActLabelNr);
        VecArray = new TArrayVec(18,Count);
        VecArray->DHeader(khoros,* oi,Elems);
    }

  if (oi != 0)
    {
      int i,j=0;
      const float fac=2.35482; // == sqrt(8*ln(2));

      for (i=0;i<ActLabelNr;i++)
	if ((usedLabel[i].Nr != 0) && j < Count)
	  {
	    VecArray->SetValue(0,j,usedLabel[i].MeanI);
	    VecArray->SetValue(1,j,usedLabel[i].COMX);
	    VecArray->SetValue(2,j,usedLabel[i].COMY);
	    VecArray->SetValue(3,j,usedLabel[i].COMZ);
	    VecArray->SetValue(4,j,usedLabel[i].Vol);
	    VecArray->SetValue(5,j,usedLabel[i].Surface);
	    VecArray->SetValue(6,j,usedLabel[i].SurfaceX);
	    VecArray->SetValue(7,j,usedLabel[i].SurfaceY);
	    VecArray->SetValue(8,j,usedLabel[i].SurfaceZ);
	    VecArray->SetValue(9,j,sqrt(usedLabel[i].VARX) * fac);  // converted into FWHM
	    VecArray->SetValue(10,j,sqrt(usedLabel[i].VARY) * fac);
	    VecArray->SetValue(11,j,sqrt(usedLabel[i].VARZ) * fac);
	    VecArray->SetValue(12,j,usedLabel[i].Ellipt);
	    VecArray->SetValue(13,j,usedLabel[i].EigenLow);
	    VecArray->SetValue(14,j,usedLabel[i].EigenHigh);
	    VecArray->SetValue(15,j,usedLabel[i].ExcentX);
	    VecArray->SetValue(16,j,usedLabel[i].ExcentY);
	    VecArray->SetValue(17,j,usedLabel[i].Nr); // just the identification code
	    j++;
	  }
      VecArray->Write(oi);
    }
  }

  if (oi)  oi->close();
  if (loi)  loi->close();
  if (ostrm) ostrm->close();
  cout << "Vector OUTPUT:\nwidth position, information\n0   Mean Intensity\n1   Center of Mass X\n2   Center of Mass Y\n3   Center of Mass Z\n4   Volume in Voxels\n5   Number of Voxels having a free surface\n6   Surfaces X\n7   Surfaces Y\n8   Surfaces Z\n9   FWHM from Variance X\n10  FWHM from Variance Y\n11  FWHM from Variance Z\n12  Excentricity (total)\n13  Moment of Inertia (low)\n14  Moment of Inertia (high)\n15  Excentricity Orientation Vector X\n16  Excentricity Orientation Vector Y\n17  Identification Number\n";
}

