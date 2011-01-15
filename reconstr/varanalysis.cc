// varanalysis : Analyses the mean and variance in a series of datasets

/*   This file is part of a software package written by 
     Rainer Heintzmann
     Institute of Applied Optics and Information Processing
     Albert Ueberle Strasse 3-5
     69120 Heidelberg
     Tel.: ++49 (0) 6221  549264
     Current Address : Max Planck Inst. for biophysical Chemistry, Am Fassberg 11, 37077 Goettingen, Germany
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

TImgArray Img,Mean,Var,Histo,VarHisto, SumHisto, LastImg, MeanDiff;

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


int kflag=0,Elem,differential=0,bins=256;

string IFileName,OFileName,HOFileName,HAOFileName,VFileName;
double binsize=0.0;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-d",parg)) {differential=1;continue;}
   if (readArg("-bins",& bins, parg)) continue;
   if (readArg("-bsize",& binsize, parg)) continue;
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-oh", HOFileName, parg)) continue;    // Histogram Output file
   if (readArg("-ah", HAOFileName, parg)) continue;    // Histogram Axis Output file
   if (readArg("-ov", VFileName, parg)) continue;    // Variance File
   if (readArg("-iX",& SizeX,parg)) continue; // input size X
   if (readArg("-iY",& SizeY,parg)) continue;
   if (readArg("-iZ",& SizeZ,parg)) continue;
    usage(argv[0]);
  }

 if (IFileName=="")
   cerr << "Missing Input Filename\n",usage(argv[0]);

 if (OFileName=="")
   cerr << "Missing Output Filename\n",usage(argv[0]);


 Elements=1;

 ofstream to(OFileName.c_str());
 ofstream hto(HOFileName.c_str());
 ofstream vto(VFileName.c_str());


 double max; 

 Histo.Resize(bins,1,1); 
 VarHisto.Resize(bins,1,1); 
 SumHisto.Resize(bins,1,1); 
 
 if (! differential)
   {
     for (Elem=0;Elem<Elements;Elem++)
       {
	 cout << "Loading Element " << Elem << "\n";
	 Elements=Img.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY, & SizeZ,Elem);
	 
	 if (Elem == 0)
	   {
	     Mean.Resize(& Img);
	     Var.Resize(& Img);
	   }
	 
	 Mean.Add(& Img);
	 Var.SqrAdd(& Img);  // now contains all the sum of squares
       }
     max = Mean.Maximum() / Elements;
     binsize=Mean.CreateVarHisto(& VarHisto, & Histo, &SumHisto, 0,max, & Var, & Mean, Elements, Elements, binsize);  // creates the histogram 
     Var.CalcVarMean(& Mean,Elements); // (SumSqr-sum*sum/Num)/Num
   }
 else // differential
   {
     Elements=LastImg.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY, & SizeZ,0);
     cout << "Elements " << Elements << "\n";
     Mean.Resize(& LastImg);
     MeanDiff.Resize(& LastImg);
     Var.Resize(& LastImg);
     Mean.Add(& LastImg);
	 
     for (Elem=1;Elem<Elements;Elem++)
       {
	 cout << "Loading Element " << Elem << "\n";
	 Elements=Img.DLoad(kflag,IFileName.c_str(),"Float",& SizeX,& SizeY, & SizeZ,Elem);
	 LastImg.Sub(& Img);
	 
	 Mean.Add(& Img);
	 MeanDiff.Add(& LastImg);
	 Var.SqrAdd(& LastImg);
	 LastImg.Copy(& Img);
       }
     max = Mean.Maximum() / Elements;
     binsize=Mean.CreateVarHisto(& VarHisto, & Histo, &SumHisto, 0,max, & Var,& MeanDiff, Elements, Elements -1, binsize);  // creates the histogram 
     Var.CalcVarMean(& MeanDiff,Elements-1); // (SumSqr-sum*sum/Num)/Num
     Mean.Mul(1.0/double(Elements));
   }


 // binsize=Mean.CreateHisto(& Histo,0,max);  // creates the histogram 
 // binsize=Mean.CreateVarHisto(& VarHisto,0,max, & Var);  // creates the histogram 

 /* for (int x=0;x < VarHisto.GetSize(0);x++)
   {
     double num = Histo.Value(x,0,0);
     if (num != 0)
       VarHisto.SetValue(x,0,0, VarHisto.Value(x,0,0)/num);
     else
       VarHisto.SetValue(x,0,0, 0);
       }*/
 cout << " Binsize of Histogram : " << binsize << "\n";

 if (! differential)
   {
     Mean.DHeader(kflag,to,1);
     Mean.Write(& to);
   }
 else
   {
     MeanDiff.DHeader(kflag,to,1);
     MeanDiff.Write(& to);
   }
 Var.DHeader(kflag,vto,1);
 Var.Write(& vto);
 VarHisto.DHeader(kflag,hto,1);
 VarHisto.Write(& hto);

 if (to) 
   to.close(); 
 if (hto)
   hto.close(); 
 if (vto)
   vto.close(); 

 ofstream hato(HAOFileName.c_str());  // Now write the histogram axis
 for (int x=0;x < Histo.GetSize(0);x++)
   {
     Histo.SetValue(x,0,0, 0 + binsize * x);
   }
 Histo.DHeader(kflag,hato,1);
 Histo.Write(& hato);
 if (hato)
   hato.close(); 

}
