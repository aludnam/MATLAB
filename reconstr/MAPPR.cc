// This version for ML-deconvolution expects a psf for every element given for reconstructing one image
// The measured projections have allready to be turned into the right postionens. The psfs have to be turned the same angles but centered
// The program is base on : Image Restoration based on Good's Roughness Penalty with Application to Fluorescence Microscopy
// P. J. Verveer and T. M. Jovin, J. Opt. Soc. Am. A , 15, 1077-1083, 1998
// maximum a posteriori likelihood poisson noise, good's roughness reconstruction  

/*   This file is part of a software package written by 
     Rainer Heintzmann
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
#include <time.h>
#include "fftarray.h"
#include "khoros.h"
#include "MAPPR.h"  // contains the important formulas

#include "parseargs.h"
#include "window.h"

static const int MaxPrjs=100;  // Do not change below 3 (at least 3 parameters are read in)

// following objects are global, because they are too big to be on the main stack

typedef float          ArrayBType;

typedef MAPPRArray<ArrayBType>  TAllArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TArray3d<ArrayBType> TPSFArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

// Images needed
static   TAllArray reconImg;          // One as first guess , this is x_k
static   TAllArray correctionImg;     // will collect the correction information before the update is done
static   TAllArray Hx2[MaxPrjs];               // Holds the convolved x^2, copied from reconPrj
static   TAllArray Hxd[MaxPrjs];               // Holds the convolved x*d
// static   TAllArray NablaPsi;          // is necessary to save the old NablaPsi in the Polak-Ribiere method
static   TAllArray & NablaPsi = Hxd[0];          // reuse Hxd array for NablaPsi
// static   TAllArray Hd2;               // Holds the convolved d^2,  reconPrj is used instead
static   TAllArray ConstrainArr;      // If given, this image is multiplicated in every step with reconstruction image
  
static   TAllArray measuredPrjs[MaxPrjs],  // holds the measured data
                   PatternArr[MaxPrjs];    // excitation patterns
static   TAllArray * OTFs[MaxPrjs];         // the OTFs

static   TAllArray reconPrj[MaxPrjs];       // Serves as data memory for computation of intermediate projection and stores Hd2

double backgr[MaxPrjs];   // C++ should guarantee, that these are cleared to zero !!

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


int main(int argc, char *argv[])
{
  int i,k,Elements=1;
  double alpha=1.0,minalpha,maxalpha,  // limiting alpha before the laplacian has to be recomputed
    sumsqr,sumsqr0,gamma=0.0,gammakorr=0.0, gammakorr2=1.0,  // gammakorr is for correcting gamma
    betha, MaxFreq=1.0,MaxZFreq=1.0,HFPercent=0.1;  // Measure Frequencies > 10% of max
  double weights[MaxPrjs];
  for (i=0;i< MaxPrjs;i++)
    weights[i]=1.0; 

  double WindowPixels=0.0;
  double tmp, dLx, dLd;  // For Polak-Ribiere: , rrold;

  int kflag=0,IterNum=60,applywindow=0,noPSFnorm=0;
  bool Steepest=false, PolakRibiere=false, FirstGuessMeasured=false;

  string IFileName, OFileName, CFileName, PFileName, OvFileName;  // Overralax-table
  string HFFileName, FFileName, SFileName, ConstrainFileName, PatternFileName;

  int MEASUREDSizeX=256, MEASUREDSizeY=256, MEASUREDSizeZ=32;
  int PSFSizeX=256, PSFSizeY=256, PSFSizeZ=32;

printCall(argc,argv);

char ** parg= & argv[1];
argc = 0; // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   // if (readArg("-norichardson",parg)) {norichardson=1;continue;}   // this is automatically set
   if (readArg("-noPSFnorm",parg)) {noPSFnorm=1;continue;}  // Do not normalize PSF
   if (readArg("-Steepest",parg)) {Steepest=true;continue;}  // use steepest descent instead of conjugate gradient
   if (readArg("-PR",parg)) {PolakRibiere=true;continue;}  // use steepest descent instead of conjugate gradient
   if (readArg("-FGM",parg)) {FirstGuessMeasured=true;continue;}  // use steepest descent instead of conjugate gradient
   if (readArg("-i", & IFileName, parg)) continue;
   if (readArg("-o", & OFileName, parg)) continue;
   if (readArg("-c", & CFileName, parg)) continue;
   if (readArg("-p", & PFileName, parg)) continue;
   if (readArg("-constrain", & ConstrainFileName, parg)) continue;   // Image that will be multiplied with reconstruction
   if (readArg("-patterns", & PatternFileName, parg)) continue;   // Accounts for spatially variing excitation patterns
   if (readArg("-hf",& HFFileName, parg)) continue;   // File for High-Frequency Output
   if (readArg("-hfp",& HFPercent, parg)) continue;   // Percentage where high frequency measure starts
   if (readArg("-mX",& MEASUREDSizeX, parg)) continue;
   if (readArg("-mY",& MEASUREDSizeY,parg)) continue;
   if (readArg("-mZ",& MEASUREDSizeZ,parg)) continue;
   if (readArg("-w1",& weights[0],parg)) continue;   // Attach a different weight to log likelihood of element 2
   if (readArg("-w2",& weights[1],parg)) continue;   // Attach a different weight to log likelihood of element 2
   if (readArg("-w3",& weights[2],parg)) continue;
   if (readArg("-pX",& PSFSizeX, parg)) continue;
   if (readArg("-pY",& PSFSizeY,parg)) continue;
   if (readArg("-pZ",& PSFSizeZ,parg)) continue;
   if (readArg("-rW",& WindowPixels,parg)) continue;
   if (readArg("-Gamma", & gamma, parg)) continue;  // regularization parameter
   if (readArg("-Backgr", & backgr[0], parg)) continue;  // regularization parameter
   if (readArg("-Backgr1", & backgr[1], parg)) continue;  // regularization parameter
   if (readArg("-Backgr2", & backgr[2], parg)) continue;  // regularization parameter
   if (readArg("-alpha0", & alpha, parg)) continue;  // starting value for iterating alpha
   if (readArg("-s", & SFileName, parg)) continue;  // start-file
   if (readArg("-f", & FFileName, parg)) continue;  // forward-file

   if (readArg("-mf",& MaxFreq, parg)) { applywindow=1;continue; }
   if (readArg("-mzf",& MaxZFreq, parg)) { applywindow=1;continue; }

   if (readArg("-n",& IterNum, parg)) continue;

   usage(argv[0]);
  }

 if (IFileName=="") cerr << "Error : No input file given !\n",exit(-1);
 if (OFileName=="") cerr << "Error : No output file given !\n",exit(-1);
 if (PFileName=="") cerr << "Error : No psf file given !\n",exit(-1);

 
 int NPsf=1;
 Elements=1;
 int PElements=1; 

 for (i=0;i<Elements;i++)
   {
     if (i >= MaxPrjs) cerr << "Error : Maximal projections number reached\n",exit(-1);
     // measuredPrjs[i]=new TAllArray;
     Elements=measuredPrjs[i].DLoad(kflag,IFileName.c_str(),"Float",& MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ,i);
     if (measuredPrjs[i].Minimum() < 0.0)
       {
	 cerr << "Warning : Measured data element " << i << " contains values < 0, clipping them to 0 !\n";
	 measuredPrjs[i].ClipAt(0.0);  // avoid negative values in measured image !
       }

     if (PatternFileName!="") // accounts for patterned excitation
       {
	 if (i < PElements)
	   PElements=PatternArr[i].DLoad(kflag,PatternFileName.c_str(),"Float",
					   & MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ,i);
	 
	 if ( ! measuredPrjs[i].SizesEqual(& PatternArr[i]))
	   {
	     cerr << "Error : PatternImage size is different from InputSize !\n";
	     cerr << "Pattern Sizes: " << MEASUREDSizeX << " x " << MEASUREDSizeY << " x " << MEASUREDSizeZ << " x " << PElements<< "\n";
	     exit(-1);
	   }
	 PatternArr[i].Normalize(MEASUREDSizeX*MEASUREDSizeY*MEASUREDSizeZ);  // Necessary for stability
       }
   }
 
 if (ConstrainFileName!="") // will be multiplied with every reconstructed image
   {
     ConstrainArr.DLoad(kflag,ConstrainFileName.c_str(),"Float",& MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ);
     
     if ( ! measuredPrjs[0].SizesEqual(& ConstrainArr))
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
       BorderWin.ApplyWindow(& measuredPrjs[i]);
     }
 
 SinWindow<TAllArray> FreqWin(MaxFreq-0.15,MaxFreq+0.15,
			   MaxFreq-0.15,MaxFreq+0.15,
			   MaxZFreq-0.15,MaxZFreq+0.15,false,false);  // elliptical Sinusoidal window
 
 TPSFArray * Psf = new TPSFArray();

 NPsf=1;
 for (i=0;i< NPsf;i++)
   {
     if (i >= MaxPrjs) cerr << "Error : Maximal psf number reached\n",exit(-1);
     if (i >= Elements)
       {
	 cerr << "Warning ! Too many Psfs given using only Psf up to #" << i-1 <<" \n"; break;
       }
     
     NPsf=Psf->DLoad(kflag,PFileName.c_str(),"Float",& PSFSizeX,& PSFSizeY,& PSFSizeZ,i);
     
     if ( ! measuredPrjs[0].SizesEqual(Psf))
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
       {
	 gammakorr += 1.0;   // after collecting the integrals gamma has to be modified, to not cause a huge change in smoothness
	 OTFs[i]->Normalize(sqrt(double(MEASUREDSizeX*MEASUREDSizeY*MEASUREDSizeZ)));  // Does contain psf at the moment
       }
     else
       {
	 // if (i == 0)
	 gammakorr += OTFs[i]->Integral();   // gamma has to be modified, to not cause a huge change in smoothness

	 OTFs[i]->Mul(sqrt(double(MEASUREDSizeX*MEASUREDSizeY*MEASUREDSizeZ)));  // Assume, that one PSF is allready normalized
       }

     OTFs[i]->scramble(1);       // Move center to the edge
     OTFs[i]->FFTfwd();          // Generate OTF
     
     if (applywindow)  // restrict PSF to low frequencies
       {
	 FreqWin.ApplyWindow(OTFs[i],false);  // may generate small negative numbers in PSF !
       }
   }

 // gammakorr *= gammakorr2;
 // Now do the correction of gamma for multiple PSFs and for non-scaled PSFs
 gamma *= gammakorr;
 if (gammakorr != 1.0)
   cout << "Gamma was corrected due to unnormalized PSFs or multiple PSFs by " << gammakorr << " to " << gamma << " \n";
 
 if (NPsf < Elements)
   cerr << "Warning ! Fewer Psfs than measured datasets ! using Psf #0 for missing Psfs ! \n";
 for (i=NPsf;i<Elements;i++)
   OTFs[i]= OTFs[0];
 
 delete(Psf); // not needed any more
 
 reconImg.Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);
 reconImg.Clear();  // clear only at the beginning. Then is needed for the iterations
 correctionImg.Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);  
 correctionImg.Clear();  // clear only at the beginning. Then is needed for the iterations

 for (i=0;i<Elements;i++)
   {
     reconPrj[i].Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);
     reconPrj[i].Clear();  // clear only at the beginning. Then is needed for the iterations
     Hx2[i].Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);
     Hx2[i].Clear();  // clear only at the beginning. Then is needed for the iterations
     Hxd[i].Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);
     Hxd[i].Clear();  // clear only at the beginning. Then is needed for the iterations
     // Hd2.Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);
   }

 // NablaPsi.Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);  
 // NablaPsi.Clear();  // will be needed clean for the Polak-Ribiere method

 cout << "Starting Reconstruction .. \n";
 
 // Generate the Forwardfile
 ofstream to;
 
 if (SFileName!="")
   {
     int sx,sy,sz;
     reconImg.DLoad(kflag,SFileName.c_str(),"Float",&sx,&sy,&sz);
     reconImg.Sqrt(); // take the squareroot everywhere
     
     if (! measuredPrjs[0].SizesEqual(& reconImg)) 
       cerr << "Error : StartImage has wrong sizes ! \n", exit(-1);
   }
 else
   {
     if (FirstGuessMeasured)
       {
	 reconImg.Copy(& measuredPrjs[0]);
	 reconImg.Sqrt(); // take the squareroot everywhere
       }
     else
       {
	 double mean=measuredPrjs[0].Integral() / MEASUREDSizeX /MEASUREDSizeY /  MEASUREDSizeZ;
	 reconImg.Set(sqrt(mean));
       }
     reconImg.Mul(1.0/gammakorr);  // korrect for non-normalized PSFs
   }

 // for (int y=0;y<MEASUREDSizeY;y++)
 //   for (int x=0;x<MEASUREDSizeX;x++)
 //     if ((x+y)%2 >= 1)
 //        reconImg.SetValue(x,y,0,-1.0*reconImg.Value(x,y,0)); 
 
 FILE * hffile=0;
 
 if (HFFileName!="")
   {
     if (! (hffile=fopen(HFFileName.c_str(),"w"))) cerr << "Error opening File " << HFFileName.c_str() << "\n",exit(-1);
   }
 
 
 sumsqr=0;
 sumsqr0=0;

 cout << "Starting first iteration \n";
 for (i=1;i < IterNum+1;i++)
   {
     if (FFileName!="")   // write a new fwd-projectionfile
       {
	 to.open(FFileName.c_str());
	 
	 if (! to )
	   {
	     cerr << "Couldn't open file " << FFileName.c_str() << " for writing !!\n" << flush;
	     exit(-1);
	   }
	 
	 if (kflag) 
	   WriteKhorosHeader(& to,"Generated by Reconstruction Set 1997",TypeString(ArrayBType(0)),MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ,Elements);
	 cerr << "writing file " << FFileName.c_str() << " \n" << flush;
	 
       } 
     
     betha=0;

     for (k=0;k<Elements;k++)
       {
	 // compute Hx^2
	 if (i <= 1)  // only in this case the Forward computation has to be done, otherwise  Hx^2 = Hx2 + 2*alpha*Hxd + alpha^2*Hd2  can be computed
	   {
	     reconPrj[0].Mul(reconImg,reconImg);  // x*x
	     if (PatternFileName!="") // account for patterned excitation
	       reconPrj[0].Mul(&PatternArr[k % PElements]);

	     reconPrj[0].FFTfwd();        // Now do first half of convolution with PSF
	     reconPrj[0].ConvMul(OTFs[k]);
	     reconPrj[0].FFTbwd();          // yield H x*x
	     // reconPrj[0].DSave(true,"/tmp/yyy.kdf");
	   }
	 else  // take the alpha from the last iteration and use it for Hx2
	   reconPrj[0].ComputeHx2(alpha,Hx2[k],Hxd[k],reconPrj[k]);

	 Hx2[k].Copy(& reconPrj[0]);    // save into Hx2, all Hx2 are needed later on for the determination of the optimal alpha
	 
	 if (FFileName!="")   
	   reconPrj[0].Write(& to);

	 // reconPrj[0].DSave(true,"/tmp/tmp1.kdf");
	 // cout << "background elem " << k << " is " << backgr[k] << "\n";
	 reconPrj[0].ComputeOneMinDiv(measuredPrjs[k],backgr[k]);  // reconPrj = 1 - g / (Hx2 + b)
	 reconPrj[0].FFTfwd();        // Now convolve with PSF, roominversion is done below
	 reconPrj[0].ConvMul(OTFs[k],true);  // complex conjugate OTF, to invert room dimensions of PSF
	 reconPrj[0].FFTbwd(); 
	 if (PatternFileName!="") // account for patterned excitation
	    reconPrj[0].Mul(&PatternArr[k % PElements]);
	 // reconPrj.DSave(true,"/tmp/tmp2.kdf");
	 // reconPrj contains Ht * (1- (g/Hx2+b))
	 if (k == 0) // first computation ? -> Hxd[0] has to be cleared
	   Hxd[0].Clear();  // this will temporarily be used for storing the sum of DeltaPsi for all elements

	 if (k < Elements-1)  // not last element ?
	   Hxd[0].ComputeDPsiSumSqrNoLaplace(reconPrj[0],reconImg,weights[k]);  // NablaPsi will be added to in Hxd
	 else  // The code below does not write any longer into Hxd
	   tmp = NablaPsi.ComputeDPsiSumSqr(Hxd[0],reconPrj[0],reconImg,gamma,weights[k]);  // NablaPsi will be constructed from prior additions in Hxd, new SumSqr is returned and r_i+1 * r_i is written
	 // Polak Rebiere: ,gamma,rrold,weights[k]);  // ...  and r_i+1 * r_i is written
       }
     // NablaPsi should contain : nambla psi of all the individual elements
     if (i == 1) 
       {
	 sumsqr = tmp;   // old Nambla psi was the same for first iteration, betha will be one, but d(k-1) is zero
	 sumsqr0 = sumsqr;
       }

     if (PolakRibiere)
       {
	 cerr << "Polak-Ribiere is not supported any longer since it is slower and uses more memory!\n";
	 exit(-1);
	 // betha = (tmp - rrold)/ sumsqr;    // the two signs in rrold cancel
	 // if (betha < 0)
	 //   {
	 //     cout << "Warning Polak-Ribiere had to be restarted \n";
	 //     betha = 0;   // to make it stable
	 //   }
       }
     else // use Fletcher-Reeves
       betha = tmp / sumsqr;  
     sumsqr = tmp;  // store for next iteration in i
	 
     cout << " betha : " << betha << " sumsqr : " << sumsqr << "\n";
     cout << " sumsqr / sumsqr0: " << sumsqr/sumsqr0 << "\n";

     // correctionImg contains d(k-1)
     if (Steepest)
       correctionImg.ComputeSteepest(NablaPsi);  // dk = - DeltaPsi ,   this is also correct for multiple elements 
     else
       correctionImg.ComputeDK(betha,NablaPsi);  // dk = bk * d(k-1) - DeltaPsi ,   this is also correct for multiple elements 

     // correctionImg.Mul(backgr[2]-backgr[1]); // just to fool the algorihtm

     correctionImg.Compute_dLx_dLd(reconImg, dLx, dLd, minalpha, maxalpha);  // will compute the two Laplacian sums, needed later on for the computation of alpha
     
     // since it is changing dk it should be done before alphas are iterated
     if (applywindow)   // because of multiplicative algorithm HF can emerge !!
       {
	 cout << " restraining frequencies ...\n";
	 FreqWin.ApplyWindow(& correctionImg,false);
	 // correctionImg.DSave(kflag,"/tmp/1.kdf");
       }
     
     // below 2x half of the convolution could be saved
     for (k=0;k<Elements;k++)
       {
	 // Now Hxd and Hd2 have to be calculated, Hx2 is already stored 
	 // reconPrj contains backconvolved values, which are not needed anymore. It will store Hd^2 now.
	 reconPrj[k].Mul(correctionImg,correctionImg);  // d*d, overwriting reconPrj
	 reconPrj[k].FFTfwd();          // Now do first half of convolution with PSF
	 reconPrj[k].ConvMul(OTFs[k]);
	 reconPrj[k].FFTbwd();          // yields Hd2 for each element
     	 if (PatternFileName!="") // account for patterned excitation
	    reconPrj[k].Mul(&PatternArr[k % PElements]);

	 Hxd[k].Mul(correctionImg,reconImg);  // Hxd=d*x
	 Hxd[k].FFTfwd();              // Now do convolution with PSF
	 Hxd[k].ConvMul(OTFs[k]);
	 Hxd[k].FFTbwd();              // yields Hxd
      	 if (PatternFileName!="") // account for patterned excitation
	    Hxd[k].Mul(&PatternArr[k % PElements]);
      }
	 
     // Time to start the Newton Iteration, last paramer is the number of iterations
     // cout << "max d : " << correctionImg.Maximum() << " min : "<< correctionImg.Minimum() << "\n";
     reconImg.TAllArray::NewtonIter(alpha,minalpha, maxalpha,
				    Elements,correctionImg,  // these are single images, the ones below exist seperately for each element
				    measuredPrjs,Hxd,reconPrj, Hx2,  // hand all elements over to determine alpha
				    dLx,dLd,gamma,backgr,weights,3);  // Hx2 as self,Hxd,Hd2,d,gamma   -> alpha

     cout << "Iteration Nr. " << i <<", Overrelaxation Factor is " << alpha << "  \n";
     reconImg.UpdateXk(correctionImg, alpha); // xk+1 = fabs(xk + alpha dk)  // also constrains x to be all positive
     
     if (ConstrainFileName!="")
       {
	 cout << " constraining image ...\n";
	 reconImg.Mul(&ConstrainArr);
       }
     
     reconImg.DSave(kflag,OFileName.c_str());
     if (CFileName!="")
       correctionImg.DSave(kflag,CFileName.c_str());   // NOT SQUARED !!
	 
     if (FFileName!="")
       to.close();
     cout << " used CPU time is " << double(clock())/double(CLOCKS_PER_SEC) << "  \n" << flush;
   } 
 reconImg.Mul(reconImg,reconImg);   // to yield the reconstructed image f
 reconImg.DSave(kflag,OFileName.c_str());
 if (hffile) fclose(hffile);
 
 return 0;
}


