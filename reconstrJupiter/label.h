#ifndef label_h
#define label_h

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

#include "stat.h"

// object for easy handling of evaluation of chromatic shifts from confocal data

const char threshtag[]="Threshhold : ";


template <const int MAXSPOT>
class label {

  typedef struct {
    int Nr;
    int Join;
    int Vol;
    float COMX;
    float COMY;
    float COMZ;
    float MeanI;
    float SumI;
    float VARX;
    float VARY;
    float VARZ;
    float VARXY;  // needed for XY-direction excentricity
    int Surface;  // Voxels that lie outside
    int SurfaceX; // Surfaces of Voxels in X-direction
    int SurfaceY; // Surfaces of Voxels in Y-direction
    int SurfaceZ; // Surfaces of Voxels in Z-direction
    float ExcentX; // Surfaces of Voxels in Z-direction
    float ExcentY; // Surfaces of Voxels in Z-direction
    float Ellipt;  // ellipticity
    float EigenLow;  
    float EigenHigh; 
  } sinfo;

typedef sinfo tarr[MAXSPOT];

 protected:
 float ScaleX,ScaleY,ScaleZ;  // may be changed
 float MaxX,MaxY,MaxZ;        // Borders of Imagestack
 float MinBorder;             // minimum distance from image border of spot
 int MinVol;                  // minimum pixel-volume of spot
 float MaxDist,MinInt;               // maximal distance between spots to be evaluated

 public:
 int spots;
 double firstlast,flsoll;


 double thresh;

 public: 
 tarr farr;

label(float SX,float SY,float SZ,float MD,float MX,float MY,float MZ,float MinB,int MinV,float MinI,float thr) : ScaleX(SX),ScaleY(SY),ScaleZ(SZ),
  MaxX(MX),MaxY(MY),MaxZ(MZ),
  MinBorder(MinB),MinVol(MinV),MaxDist(MD),MinInt(MinI),
  spots(0),thresh(thr) 
  {}

  void ClearStat(int i)
  {
      farr[i].Vol=0;
      farr[i].Surface=0;
      farr[i].SurfaceX=0;
      farr[i].SurfaceY=0;
      farr[i].SurfaceZ=0;
      farr[i].COMX=0;
      farr[i].COMY=0;
      farr[i].COMZ=0;
      farr[i].MeanI=0;
      farr[i].SumI=0;
      farr[i].VARX=0;
      farr[i].VARY=0;
      farr[i].VARZ=0;
      farr[i].VARXY=0;
      farr[i].ExcentX=0;
      farr[i].ExcentY=0;
      farr[i].Ellipt=0;
  }

  void Clear(int i)
  {
      farr[i].Nr=0;
      farr[i].Join=0;
      ClearStat(i);
  }
  
void Clear(void)  // clears first array only
{
  int i;
  for (i=0;i<MAXSPOT;i++)
      Clear(i);
}

void ClearStat(void)  // clears first array only
{
  int i;
  for (i=0;i<MAXSPOT;i++)
      ClearStat(i);
}


double distance(label<MAXSPOT> * farr2, int i,int j)
{
  double dx,dy,dz;
  dx= ScaleX*(farr[i].COMX- farr2->farr[j].COMX);  // same scaling is assumed
  dy= ScaleY*(farr[i].COMY- farr2->farr[j].COMY);
  dz= ScaleZ*(farr[i].COMZ- farr2->farr[j].COMZ);
  return sqrt(dx*dx+dy*dy+dz*dz);
}

int SearchNearest(label<MAXSPOT> * farr2,int nr)
{
  int minnr=0,j;
  double dist,mindist=1e4;

  for (j=0;j<farr2->spots;j++)
    if ((dist=distance(farr2,nr,j)) < mindist)
      mindist=dist,minnr=j;

  return minnr;
}

int counts(int i)
{
  // cout << farr[i].Vol << "  " << farr[i].COMX << "   " << farr[i].COMY << "   " << farr[i].COMZ << " \n";
  // cout << MinVol << "  " << MinBorder << "   " << MaxX << "   " << MaxY << "   " << MaxZ << " \n";
  if (farr[i].Vol < MinVol) return 0;
  if (farr[i].SumI < MinInt) return 0;
  if (farr[i].COMX < MinBorder/ScaleX) return 0;
  if (farr[i].COMY < MinBorder/ScaleY) return 0;
  if (farr[i].COMZ < MinBorder/ScaleZ) return 0;
  if (farr[i].COMX > MaxX- MinBorder/ScaleX) return 0;
  if (farr[i].COMY > MaxY- MinBorder/ScaleY) return 0;
  if (farr[i].COMZ > MaxZ- MinBorder/ScaleZ) return 0;
  return 1;
}

void TransferData(label<MAXSPOT> * farr2)  // transfers data to array2
{
  int i;
  farr2->spots=spots;

  for(i=0;i< spots;i++)
     farr2->farr[i] = farr[i];

}

void calcshifts(label<MAXSPOT> * farr2)
{
int i,j;
double DX,DY,DZ,I1,I2;
Stat<double> XStat,YStat,ZStat,I1Stat,I2Stat;

cout << "SpotNr, DX(scaled), DY(scaled), DZ(scaled), ZPOS, ISum1(with thresh), ISum2(with thresh) \n";

for (i=0;i<spots;i++)
 if (farr[i].Vol > 0)
  {
    j = SearchNearest(farr2,i);  // return nearest spotnr in farr2
    DX = ScaleX*(farr[i].COMX-farr2->farr[j].COMX);
    DY = ScaleY*(farr[i].COMY-farr2->farr[j].COMY);
    DZ = ScaleZ*(farr[i].COMZ-farr2->farr[j].COMZ);
    I1 = thresh*farr[i].Vol+farr[i].SumI;
    I2 = farr2->thresh*farr2->farr[i].Vol+farr2->farr[j].SumI;

    // cout << "distance : " << distance(farr2,i,j) << "\n";

    if ((distance(farr2,i,j) < MaxDist) &&
        counts(i) &&
        farr2->counts(j))
      {
	cout << farr[i].Nr << "   " << DX << "   " << DY << "  " << DZ <<
	  "   " << ScaleZ*farr[i].COMZ << "     " << I1 << "   " << I2 << "\n";
	XStat.Register(DX);
	YStat.Register(DY);
	ZStat.Register(DZ);
	I1Stat.Register(I1);
	I2Stat.Register(I2);
      }
    else
      cout << farr[i].Nr << "XXX  " << DX << "   " << DY << "  " << DZ <<
	"   " << ScaleZ*farr[i].COMZ << "\n";
  }

cout << " DX = " << XStat.Mean() << " +/-  " << XStat.StdDev() << "\n";
cout << " DY = " << YStat.Mean() << " +/-  " << YStat.StdDev() << "\n";
cout << " DZ = " << ZStat.Mean() << " +/-  " << ZStat.StdDev() << "\n";
cout << " I1 = " << I1Stat.Mean() << " +/-  " << I1Stat.StdDev() << "\n";
cout << " I2 = " << I2Stat.Mean() << " +/-  " << I2Stat.StdDev() << "\n";
// if (flsoll != 0.0)
//   SumZ *= (flsoll/firstlast), SumSZ *= (flsoll/firstlast)*(flsoll/firstlast);

}

float FirstLastDist(void)
{
  int last;
  for(last=spots-1;((last>0) &&(farr[last].Vol <= 0));) last--;

  cerr << "spot " << farr[last].Nr << " minus spot " << farr[0].Nr << "\n";
  return ScaleZ*(farr[last].COMZ-farr[0].COMZ);
}

/// opens a datafile and parses its contents
int ParseArray(char * FileName)
{
static char line[1024], * ppos;
int spot=0;

ifstream file1(FileName);
  if (! file1) {
    fprintf(stderr, "couldnt open datafile %s \n",FileName);
    exit(-1);
    }

cout << "\n// read file : " << FileName << "\n";

file1.getline(line,1024,'\n');
while (! file1.eof())
  {
  if (line[0] != '/')
   {
    if (! sscanf(line,"%d%f%f%f%f%f\n",& farr[spot].Vol,& farr[spot].COMX,
		 & farr[spot].COMY,& farr[spot].COMZ,& farr[spot].MeanI,& farr[spot].SumI))
       {
	 cerr << "Error in Fileformat !!\n";
	 cerr << "Line reads : " << line;
	 exit(-1);
       }
    // printf("%d I : %g\n",spot,farr[spot][5]);
    spot++;
   }
  else
    {
      cout << line << "\n";
      if ((ppos=strstr(line,threshtag)))
	{
	  sscanf(ppos+strlen(threshtag),"%lf",& thresh);
	  // cout << "\n thresh : " << thresh << "\n";
	}
    }	
 
  file1.getline(line,1024,'\n');
  }
return spot;
}

};
    
#endif
