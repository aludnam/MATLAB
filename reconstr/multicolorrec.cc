// #define DEBUG_SAVE
// This version for ML-deconvolution expects a psf for every element (color) given for reconstructing a series of 
// (partially segmented) images.
// The program is meant for estimating the most probable chromosome in every reconstructed pixel.

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
#include <time.h>
#include "fftarray.h"
#include "khoros.h"

#include "parseargs.h"
#include "window.h"

static clock_t lasttime;
static double totaltime;

static const int MaxPrjs=1000;

// following objects are global, because they are too big to be on the main stack

typedef float          ArrayBType;

typedef TFFTArray<ArrayBType>  TAllArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TArray3d<ArrayBType> TPSFArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

// Images needed
static   TAllArray * reconImg[MaxPrjs];          // Zero as first guess  // Arra, ItsPos, ItsDir
static   TAllArray * correctionImg[MaxPrjs];     // will collect the correction information before the update is done
static   TAllArray ConstrainArr;      // If given, this image is multiplicated in every step with reconstruction image
static   TAllArray colorMatrix;      // defines the colors of every chromosome
static   TAllArray RecColorMatrix;      // defines the colors of every chromosome
  
static   TAllArray * measuredPrjs[MaxPrjs]; // actually holds the measured data
static   TAllArray * OTFs[MaxPrjs];         // the OTFs

static   TAllArray reconPrj;          // Serves as data memory for computation of intermediate projection


void usage(char * filename)
{
  cerr << "usage: " << filename << " [-k] [-mX X ...] [-i inputfile] [-p psffile] [-o outputfile] [-c corrimage] [-fX X ...] [-n iterations]\n" << flush;
  cerr << "-mX X -mY Y -mZ Z :  X,Y,Z == sizes of measured projection\n";
  cerr << "-pX X -pY Y -pZ Z :  X,Y,Z == sizes of PSF\n";
  cerr << "-k : khoros data format flag\n";
  cerr << "-s file : file for starting iteration with\n";
  cerr << "-f file : file to save forward projection in\n";
  cerr << "-rW file : Width (Pixels) of edge blending window\n";
  cerr << "-mf relative maximal frequency (if <1.0 image will be restrained to this each step)\n";
  exit(-1);
}

void EntropyPenalty(int Signatures, double weight)
{
  int x,y,z,l;
  int sizeX=correctionImg[0]->GetSize(0);
  int sizeY=correctionImg[0]->GetSize(1);
  int sizeZ=correctionImg[0]->GetSize(2);
  double sum,val,penal;

  for (z=0;z<sizeZ;z++)
  for (y=0;y<sizeY;y++)
  for (x=0;x<sizeX;x++)
    {
      sum=0;
      for (l=0;l<Signatures;l++)
	{
	  sum += reconImg[l]->Value(x,y,z);
	}
      for (l=0;l<Signatures;l++)
	{
	  val = reconImg[l]->Value(x,y,z) / sum;
	  penal = weight*(1+log(val));   // This stems from the derivative of the shannon entropy : d/dp_k    ( sum(p_k * ln(p_k))
	  val = correctionImg[l]->Value(x,y,z);
	  correctionImg[l]->SetValue(x,y,z, val + penal);
	}
    }
}

void ConstrainDAPIBrightness(int Signatures)  // assumes the first image (DAPI channel) to define the brightness and normalizes image
{  
  int x,y,z,l;
  int sizeX=correctionImg[0]->GetSize(0);
  int sizeY=correctionImg[0]->GetSize(1);
  int sizeZ=correctionImg[0]->GetSize(2);
  double sum,val;

  for (z=0;z<sizeZ;z++)
  for (y=0;y<sizeY;y++)
  for (x=0;x<sizeX;x++)
    {
      sum=0;
      for (l=0;l<Signatures;l++)
	{
	  sum += reconImg[l]->Value(x,y,z);
	}
      for (l=0;l<Signatures;l++)
	{
	  val = reconImg[l]->Value(x,y,z) / sum * measuredPrjs[0]->Value(x,y,z);
	  reconImg[l]->SetValue(x,y,z, val);
	}
    }
}

void EntropyRegionPenalty(int Signatures, double weight, int regionsize)
{
  int x,y,z,l,dx,dy,dz;
  int sizeX=correctionImg[0]->GetSize(0);
  int sizeY=correctionImg[0]->GetSize(1);
  int sizeZ=correctionImg[0]->GetSize(2);
  double sum,val,penal;

  for (l=0;l<Signatures;l++)
  for (z=0;z<sizeZ;z++)
  for (y=regionsize;y<sizeY-regionsize;y++)
  for (x=regionsize;x<sizeX-regionsize;x++)
    {
      sum=0;
      for (dy=-regionsize;dy <= regionsize;dy++)
      for (dx=-regionsize;dx <= regionsize;dx++)
	{
	  sum += reconImg[l]->Value(x+dx,y+dy,z);
	}

      val = reconImg[l]->Value(x,y,z) / sum;
      penal = weight*(1+log(val));   // This stems from the derivative of the shannon entropy : d/dp_k    ( sum(p_k * ln(p_k))
      val = correctionImg[l]->Value(x,y,z);
      correctionImg[l]->SetValue(x,y,z, val + penal);
    }
}

int main(int argc, char *argv[])
{
  int i,k,l,Elements=1,Signatures=23, FirstConvolution=0;
  double overrelax=1.0,MaxFreq=1.0,MaxZFreq=1.0,HFPercent=0.1,penalty=0.03,regionpenalty=0.03;  // Measure Frequencies > 10% of max
  double WindowPixels=0.0;

  int kflag=0,IterNum=60,applywindow=0,noPSFnorm=0,noDeconv=0,noColorEstimate=0;
  bool HolmesLiu=false;

  string IFileName, OFileName, CMFileName, CFileName, PFileName, OvFileName,  // Overralax-table
  HFFileName, FFileName, SFileName, ConstrainFileName,CMatrixOutFile;

  int MEASUREDSizeX=256, MEASUREDSizeY=256, MEASUREDSizeZ=32;
  int SizeX=23, SizeY=7, SizeZ=1;
  int PSFSizeX=256, PSFSizeY=256, PSFSizeZ=32;

printCall(argc,argv);

char ** parg= & argv[1];
argc = 0; // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   // if (readArg("-norichardson",parg)) {norichardson=1;continue;}   // this is automatically set
   if (readArg("-noPSFnorm",parg)) {noPSFnorm=1;continue;}  // Do not normalize PSF
   if (readArg("-noDeconv",parg)) {noDeconv=1;continue;}  // Do not normalize PSF
   if (readArg("-noColorEstimate",parg)) {noColorEstimate=1;continue;}  // Do not normalize PSF
   if (readArg("-firstConv", & FirstConvolution, parg)) continue;
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-cmatrix", CMFileName, parg)) continue;
   if (readArg("-cmatout", CMatrixOutFile, parg)) continue;
   if (readArg("-c", CFileName, parg)) continue;
   if (readArg("-p", PFileName, parg)) continue;
   if (readArg("-constrain", ConstrainFileName, parg)) continue;   // Image that will be multiplied with reconstruction
   if (readArg("-hf",HFFileName, parg)) continue;   // File for High-Frequency Output
   if (readArg("-hfp",& HFPercent, parg)) continue;   // Percentage where high frequency measure starts
   if (readArg("-penalty",& penalty, parg)) continue;   // Percentage where high frequency measure starts
   if (readArg("-rpenalty",& regionpenalty, parg)) continue;   // Percentage where high frequency measure starts
   if (readArg("-mX",& MEASUREDSizeX, parg)) continue;
   if (readArg("-mY",& MEASUREDSizeY,parg)) continue;
   if (readArg("-mZ",& MEASUREDSizeZ,parg)) continue;
   if (readArg("-pX",& PSFSizeX, parg)) continue;
   if (readArg("-pY",& PSFSizeY,parg)) continue;
   if (readArg("-pZ",& PSFSizeZ,parg)) continue;
   if (readArg("-rW",& WindowPixels,parg)) continue;
   if (readArg("-Over", OvFileName, parg)) continue;  // relaxationtable (ascii-file with overralaxation values)
   if (readArg("-HolmesLiu",parg)) {HolmesLiu=true;continue;}
   if (readArg("-s", SFileName, parg)) continue;  // start-file
   if (readArg("-f", FFileName, parg)) continue;  // forward-file

   if (readArg("-mf",& MaxFreq, parg)) { applywindow=1;continue; }
   if (readArg("-mzf",& MaxZFreq, parg)) { applywindow=1;continue; }

   if (readArg("-n",& IterNum, parg)) continue;

   usage(argv[0]);
  }

 if (IFileName=="") cerr << "Error : No input file given !\n",exit(-1);
 if (OFileName=="") cerr << "Error : No output file given !\n",exit(-1);
 if (PFileName=="") cerr << "Error : No psf file given !\n",exit(-1);
 if (CMFileName=="") cerr << "Error : No color matrix file given! !\n",exit(-1);

 lasttime = clock();
 totaltime = 0.0;
 cout << " start CPU time is " << lasttime/(double) CLOCKS_PER_SEC << "  \n" << flush;
 
 int NPsf=1;
 Elements=1;


 for (i=0;i<Elements;i++)
   {
     if (i >= MaxPrjs) cerr << "Error : Maximal projections number reached\n",exit(-1);
     measuredPrjs[i]=new TAllArray;
     Elements=measuredPrjs[i]->DLoad(kflag,IFileName.c_str(),"Float",& MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ,i);
     measuredPrjs[i]->ClipAt(0.0);  // avoid negative values in measured image !
   }

 reconPrj.Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);
 
 if (ConstrainFileName!="") // will be multiplied with every reconstructed image
   {
     ConstrainArr.DLoad(kflag,ConstrainFileName.c_str(),"Float",& MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ);
     
     if ( ! measuredPrjs[0]->SizesEqual(& ConstrainArr))
       {
	 cerr << "Error : ConstrainImage size is different from InputSize !\n";
	 exit(-1);
       }
   }
 
 double innerborder=1.0-(2*WindowPixels/MEASUREDSizeX);
 GaussWindow<TAllArray> BorderWin(innerborder,1.0,innerborder,1.0,innerborder,1.0,
				  true,false);  // rectangular Gaussian window, no ring shape
 BorderWin.gaussstretch=0.7;
 
 if (WindowPixels > 0)
   for (i=0;i<Elements;i++)
     {
       BorderWin.ApplyWindow(measuredPrjs[i]);
     }
 
 SinWindow<TAllArray> FreqWin(MaxFreq-0.15,MaxFreq+0.15,
			   MaxFreq-0.15,MaxFreq+0.15,
			   MaxZFreq-0.15,MaxZFreq+0.15,false,false);  // elliptical Sinusoidal window
 
 
 TPSFArray * Psf = new TPSFArray();

 NPsf=1;
 if (! noDeconv)
   {
 for (i=0;i< NPsf;i++)
   {
     if (i >= MaxPrjs) cerr << "Error : Maximal psf number reached\n",exit(-1);
     if (i >= Elements)
       {
	 cerr << "Warning ! Too many Psfs given using only Psf up to #" << i-1 <<" \n"; break;
       }
     
     NPsf=Psf->DLoad(kflag,PFileName.c_str(),"Float",& PSFSizeX,& PSFSizeY,& PSFSizeZ,i);
     
     if ( ! measuredPrjs[0]->SizesEqual(Psf))
       {
	 cerr << "Warning : PSF size is different from InputSize !\n";
	 cerr << "Zero-Padding PSF !\n";
       }
     
     if (Psf->Minimum() < 0.0)
       {
	 cerr << "Warning : PSF contains values < 0, clipping them to 0 !\n";
	 Psf->ClipAt(0.0);  // avoid PSF having negative values
       }

     OTFs[i]= new TAllArray(MEASUREDSizeX,
			    MEASUREDSizeY,
			    MEASUREDSizeZ);
     
     OTFs[i] -> Clear();
     
     Copy(Psf,OTFs[i]);       // centered copy algorithm
     
     cout << "Computing OTF #"<<i<<"\n";
     
     // To get the integral intensity correct use following normalization : (0-freq == 1.0 +0*i)
     if (! noPSFnorm)
       OTFs[i]->Normalize(1.0/(MEASUREDSizeX*MEASUREDSizeY*MEASUREDSizeZ));  // Does contain psf at the moment
     else
       OTFs[i]->Mul(1.0/(MEASUREDSizeX*MEASUREDSizeY*MEASUREDSizeZ));  // Assume, that one PSF is allready normalized
     // Normalization is special, since it accounts for the optimized factors during FFTs
     // doing the fft with "true" needs Mul(sqrt(X*Y*Z));

     OTFs[i]->scramble(1);       // Move center to the edge
     OTFs[i]->FFTfwd(false);          // Generate OTF
     
     if (applywindow)  // restrict PSF to low frequencies
       {
	 FreqWin.ApplyWindow(OTFs[i],false);  // may generate small negative numbers in PSF !
       }
   }
 
 if (NPsf < Elements)
   cerr << "Warning ! Fewer Psfs than measured datasets ! using Psf #0 for missing Psfs ! \n";
 for (i=NPsf;i<Elements;i++)
   OTFs[i]= OTFs[0];
 
 delete(Psf); // not needed any more
   }

 colorMatrix.DLoad(kflag,CMFileName.c_str(),"Float",& SizeX,& SizeY,& SizeZ);
 RecColorMatrix.DLoad(kflag,CMFileName.c_str(),"Float",& SizeX,& SizeY,& SizeZ);  // just a copy used during reconstruction

 //colorMatrix.FirstYNormalize(1.0);  // Normalize each color to a DAPI Intensity of one
 //RecColorMatrix.FirstYNormalize(1.0);  // Normalize each color to a DAPI Intensity of one

 if ((SizeY != Elements) || (SizeZ != 1))
   cerr << "Error: Color matrix file is in diagreement with number of colors or has SizeZ != 1!",exit(-1);
 Signatures = SizeX;
 
 for (i=0;i<Signatures;i++)
   {
     reconImg[i]=new TAllArray;
     reconImg[i]->Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);
     reconImg[i]->Set(1.0); // Only accessible measured parts are important
     correctionImg[i]=new TAllArray;
     correctionImg[i]->Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);  
   }

 cout << "Starting Reconstruction .. \n";
 
 // Generate the Forwardfile
 ofstream to,to2,to3;
 
 if (SFileName!="")
   {
     int sx,sy,sz;
     int Sig2=1;
     for (i=0;i<Sig2;i++)
       {
	 Sig2=reconImg[i]->DLoad(kflag,SFileName.c_str(),"Float",&sx,&sy,&sz,i);
     
	 if (! measuredPrjs[0]->SizesEqual(reconImg[i])) 
	   cerr << "Error : StartImage has wrong sizes ! \n", exit(-1);
	 if (Sig2 != Signatures) 
	   cerr << "Error: Color matrix file is in diagreement with number of signatures!",exit(-1);
       }
   }
 
 FILE * overfile=0,* hffile=0;
 
 if (OvFileName!="")
   {
     if (! (overfile=fopen(OvFileName.c_str(),"r"))) cerr << "Error opening File " << OvFileName << "\n",exit(-1);
   }
 
 if (HFFileName!="")
   {
     if (! (hffile=fopen(HFFileName.c_str(),"w"))) cerr << "Error opening File " << HFFileName << "\n",exit(-1);
   }
 
 overrelax = 1.0;
 cout << "Starting first iteration \n";
 for (i=1;i < IterNum+1;i++)
   {
     if (FFileName!="")   // write a new fwd-projectionfile
       {
	 to.open(FFileName.c_str());
	 if (! to )
	   {
	     cerr << "Couldn't open file " << FFileName << " for writing !!\n" << flush;
	     exit(-1);
	   }

	 ArrayBType dummy;
	 
	 if (kflag) 
	   {
	   WriteKhorosHeader(& to,"Generated by Reconstruction Set 1997",TypeString(dummy),MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ,Elements);
	   cerr << "writing file " << FFileName << " \n" << flush;
	   }
	 
       } 

     if (OFileName!="")   // write a new fwd-projectionfile
       {
	 to2.open(OFileName.c_str());
	 if (! to2)
	   {
	     cerr << "Couldn't open file " << OFileName << " for writing !!\n" << flush;
	     exit(-1);
	   }
	 
	 ArrayBType dummy;
	 
	 if (kflag) 
	   {
	   WriteKhorosHeader(& to2,"Generated by Reconstruction Set 2002",TypeString(dummy),MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ,Signatures);
	   cerr << "writing file " << OFileName << " \n" << flush;
	   }
	 
       } 
     if (CFileName!="")   // write a new fwd-projectionfile
       {
	 to3.open(CFileName.c_str());
	 if (! to3)
	   {
	     cerr << "Couldn't open file " << CFileName << " for writing !!\n" << flush;
	     exit(-1);
	   }
	 
	 ArrayBType dummy;
	 
	 if (kflag) 
	   {
	   WriteKhorosHeader(& to3,"Generated by Reconstruction Set 2002",TypeString(dummy),MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ,Signatures);
	   cerr << "writing file " << CFileName << " \n" << flush;
	   }
	 
       } 
     
     
     if (! HolmesLiu)
       {
	 if (! overfile)
	   {
	     overrelax=1.0;
	     
	     if (i > 100) overrelax=1.0; 
	   }
	 else
	   fscanf(overfile,"%lf",&overrelax);
       }
     
     
     
     for (l=0;l<Signatures;l++)  // Clear correction Image
       {
	 correctionImg[l]->Clear();
       }
     
     double likelihoodsum = 0.0,likelihood=0.0;
     if (hffile)
       fprintf(hffile,"%d\t%g\t",i,overrelax);

     for (k=0;k<Elements;k++)  // iterates over different colors
       {
	 reconPrj.Clear();
     // apply the color matrix to the object emission prior to convolution:
	 for (l=0;l<Signatures;l++)  // Now compute the contribution of all chromosomes to this color in every pixel.
	   {
	     reconPrj.AddMul(reconImg[l],RecColorMatrix.Value(l,k,0));
	   }

	 if (! noDeconv && k >= FirstConvolution)
	   {
	     reconPrj.FFTfwd(false);        // Now do first half of convolution with PSF
	     reconPrj.ConvMul(OTFs[k]);
	     reconPrj.FFTbwd(false);        
	   }
	 
	 if (FFileName!="")   
	   reconPrj.Write(& to);

	 if (hffile)
	   {
	     likelihood=reconPrj.ArgDivSelfM1LogLikelihood(measuredPrjs[k]);
	     likelihoodsum += likelihood;
	     fprintf(hffile,"%g\t",likelihood);
	   }
	 else
	   reconPrj.ArgDivSelfM1(measuredPrjs[k]);


	 if (! noDeconv && k >= FirstConvolution)
	   {
	     reconPrj.FFTfwd(false);        // Now convolve with PSF, roominversion is done below
	     reconPrj.ConvMul(OTFs[k],true);  // complex conjugate OTF, to invert room dimensions of PSF
	     reconPrj.FFTbwd(false);        // should be done now for every color
	   }
     //if (! noColorEstimate) // The stuff below was for updating the estimate of the color matrix

     // apply the transpose of the color matrix to the correction image:
	 for (l=0;l<Signatures;l++)  // iterate over chromosomes
	   {
	     correctionImg[l]->AddMul(& reconPrj,RecColorMatrix.Value(l,k,0));  // this is now done in real space
	   }
       }

     if (hffile)
       {
	 fprintf(hffile,"%g\n",likelihoodsum);
       }
     
     // pulling this out of the loop makes it faster for multiple elements
     if (applywindow)   // because of multiplicative algorithm HF can emerge !!
       {
	 cout << " restraining frequencies ...\n";
	 for (l=0;l<Signatures;l++)  // Now compute the contribution of all chromosomes to this color in every pixel.
	   {
	     FreqWin.ApplyWindow(correctionImg[l],false);
	     // correctionImg[l]->DSave(kflag,"/tmp/1.kdf");
	   }
       }

     if (penalty != 0.0)
       EntropyPenalty(Signatures,penalty);  // will be applied to correction Image [l]
     if (regionpenalty != 0.0)
       EntropyRegionPenalty(Signatures,regionpenalty,1);  // will be applied to correction Image [l]

     for (l=0;l<Signatures;l++)  // iterate over chromosomes
       {
	 if (HolmesLiu)
	 {
	   double maxsum=-1e30;
	   overrelax=0.5;  // will be set to 1.0 in the first step
	   likelihoodsum = -1e29;
	   while (likelihoodsum > maxsum)
	     {
	       maxsum = likelihoodsum;
	       overrelax *= 2.0;
	       likelihoodsum = 0.0;
	       for (k=0;k<Elements;k++)
		 {
		   reconPrj.Copy(reconImg[l]);  // just transfer it into the tmp
		   // reconPrj.ClipAt(0.000001);    // Constraint to prevent impossible overrelaxiations : Now done below !
		   reconPrj.MulSelfAdd(correctionImg[l],float(overrelax) / float(Elements));        // Apply correction
		   if (! noDeconv && k >= FirstConvolution)
		     {
		       reconPrj.FFTfwd(false);        // Now do first half of convolution with PSF
		       reconPrj.ConvMul(OTFs[k]);
		       reconPrj.FFTbwd(false);        
		     }
		   
		   likelihoodsum+=reconPrj.ArgDivSelfM1LogLikelihood(measuredPrjs[k]);
		 }
	       printf("overrelax : %g\t likelihoodsum: %g\n",overrelax,likelihoodsum);
	     }
	   overrelax /= 2.0;  // because last step was unsucessful
	 }
     
	 cout << "Iteration Nr. " << i <<", Signature " << l << ", Overrelaxation Factor is " << overrelax << "  \n";

	 //reconImg[l]->MulSelfAdd(correctionImg[l],float(overrelax) / float(Elements));        // Apply correction
	 reconImg[l]->MulSelfAdd(correctionImg[l],float(overrelax));        // Apply correction
	 // also clips at 0.000001
	 
	 if (ConstrainFileName!="")
	   {
	     cout << " constraining image ...\n";
	     reconImg[l]->Mul(&ConstrainArr);
	   }
	 // reconImg[l]->ClipAt(0.000001);    // Constraint to prevent impossible overrelaxiations


	 
	 if (CFileName!="")
	   correctionImg[l]->Write(& to3);
       }

     // if (i == IterNum)  // Only at the final iteration !
     //   ConstrainDAPIBrightness(Signatures);  // Constrains the reconImg data to a sum corresponding to the DAPI image (channel 0)
     // accordingly the colors were normalized to a DAPI intensity of one

     for (l=0;l<Signatures;l++)  // iterate over chromosomes
       {
	 reconImg[l]->Write(& to2);
       }
     
     to2.close();
     to3.close();

     if (FFileName!="")
       to.close();
     
     
     totaltime += (clock()-lasttime)/double(CLOCKS_PER_SEC);
     lasttime = clock();
     cout << " used CPU time is " << totaltime << "  \n" << flush;
   } 
 if (hffile) fclose(hffile);
 if (CMatrixOutFile!="")
   RecColorMatrix.DSave(kflag,CMatrixOutFile.c_str());
 return 0;
}


