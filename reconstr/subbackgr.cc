// subbackgr : Subtracts the background of an image (measured as mean value of lower x% of pixels)

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
#include "rawarray.h"
#include "parseargs.h"

typedef float ArrayBType;
typedef TArray3d<ArrayBType> TImgArray;
// typedef Dyn3dArr(ArrayBType)  TImgArray;

TImgArray Img,Histo,ASlice;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-iX X ...] [-i inputimage] [-o outputfile] \n" << flush;
  cerr <<  "[-iX -iY -iZ]  : sizes of images\n" << flush;
  cerr <<  "[-p]:      percentage between min and max to be counted\n" << flush;
  cerr <<  "[-n] :     iterations to count percentage\n" << flush;
  cerr <<  "[-f] :     factor of background to subtract\n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

int SizeX=32;  // These is the standart size, if raw data is used
int SizeY=32;  
int SizeZ=32; 
int Elements=1;

double percentage=10,factor=1.0,scanhp=0.3;
int    iterations=2,sortit=0,histogram=0,bins=256,usebins=5;

int kflag=0,Elem;

string IFileName, OFileName, HOFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-s",parg)) {sortit=1;continue;}
   if (readArg("-histo",parg)) {histogram=1;continue;}   // process data in Histogram mode
   if (readArg("-bins",& bins,parg)) continue; // number of bins for histogram
   if (readArg("-scanhp",& scanhp,parg)) continue; // looks for the maximum in the first ... % of the histogram
   if (readArg("-ub",& usebins,parg)) continue; // use +/- ... bins for determination of the center -> background 
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-oh", HOFileName, parg)) continue;    // Histogram Output file
   if (readArg("-iX",& SizeX,parg)) continue; // input size X
   if (readArg("-iY",& SizeY,parg)) continue;
   if (readArg("-iZ",& SizeZ,parg)) continue;
   if (readArg("-p",& percentage,parg)) continue;
   if (readArg("-n",& iterations,parg)) continue;
   if (readArg("-f",& factor,parg)) continue;
    usage(argv[0]);
  }

 if (IFileName=="")
   cerr << "Missing Input Filename\n",usage(argv[0]);

 if (OFileName=="")
   cerr << "Missing Output Filename\n",usage(argv[0]);


 Elements=1;

 if (! histogram) 
   HOFileName[0] = 0;
 ofstream to(OFileName.c_str());
 ofstream hto(HOFileName.c_str());


 double min,mean,clip,max,py,pz,hmaxp,hmax,binsize; 
 int i;

 if (histogram)
   {
     Histo.Resize(bins,1,1); 
   }
 
 for (Elem=0;Elem<Elements;Elem++)
   {
     cout << "Loading Element " << Elem << "\n";
     Elements=Img.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY, & SizeZ,Elem);

     min = Img.Minimum();
     max = Img.Maximum();
     mean = Img.Mean();

     if (histogram)
       {
         binsize=Img.CreateHisto(& Histo,min,max);  // creates the histogram 
         cout << "Histogram : Binsize = " << binsize<<" starting at " << min << "\n";

         hmax = Histo.Value(0,0,0);
         hmaxp =0;
         for (i=1;i<(bins*scanhp);i++)
           if (Histo.Value(i,0,0) > hmax) 
           {
             hmax=Histo.Value(i,0,0);
             hmaxp=i;
           }
         if (hmaxp+usebins>= bins*scanhp) 
           cerr << "WARNING! Maximum at "<<hmaxp<<" is very close to upper boarder searched !! \n";

         if (hmaxp<usebins) 
         {
           usebins=int(hmaxp+0.5);
           cerr << "WARNING! Number of bins had to be reduced to "<<usebins<<" to accomodate for the maximum found at " << hmaxp << "\n";
         }
         if (hmaxp+usebins>=bins) 
         {
           usebins=int(bins-hmaxp-1.0+0.5);
           cerr << "WARNING! Number of bins had to be reduced to "<<usebins<<" to accomodate for the maximum found at " << hmaxp << "\n";
         }

         py=0;pz=0;
         Histo.RangedCOM(hmaxp,py,pz,usebins,0,0);
         cout << "Maximum in histogram at " << hmaxp << "\n";
         clip=hmaxp*binsize + min;
         cout << "Histogram method is subtracting " << clip << "\n";
       }   // .. if histogram
     else if (! sortit)
       {
         clip = (max-min)*percentage/100.0+min;

         for (i=0;i<iterations;i++)
         {
           cout << "clip is " << clip << "\n";
           clip = Img.Mean(min,clip); // calculate mean of voxels between min and clip
         }

         clip = ((1.0-factor)*clip+factor*mean);
         cout << "subtracting " << clip << "\n";
       }
     else
       {
         Img.Sort();
         clip=Img.Mean(int(percentage/100.0*SizeX*SizeY*SizeZ));
	 // load data again
         Elements=Img.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY, & SizeZ,Elem);
         cout << "subtracting " << clip << "\n";
       }
     Img.Sub(clip);

     if (Elem == 0)
       {
       Img.DHeader(kflag,to,Elements);
       if (hto && histogram)
         Histo.DHeader(kflag,hto,Elements);
       }
	 
     Img.Write(& to);
     if (hto && histogram)
       Histo.Write(& hto);
   }

 if (to) 
   to.close(); 
 if (hto && histogram)
   hto.close(); 
}
