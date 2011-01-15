// mulwave : multiplies complex input for every element with a exp(-ikx) wave
// thus the Fourier components are shited to the correct placed and added

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

// the wave is estimated from delta peaks in the fft, which  will be counted in the following way:
// example : 0 order and +- 1 orderX and +-1 orderY in 2D
//      0   1   2
//      3   4   5
//      6   7   8

// Thus the zero order is position 4

#include <iostream>
#include "rawarray.h"
#include "fftarray.h"
#include "parseargs.h"
#include <math.h>
#include <algorithm>
#include <string>
#include <complex>

const double Pi=2.0 * asin(1.0);
const int MaxOrders=20;

struct AnOrder {
    int X;
    int Y;
    int Elem;
    double DistSqr;
    double OrderFactor;
};

class CmpOrders {
public:
  int operator()(const AnOrder & o1, const AnOrder &o2)
  {return o1.DistSqr < o2.DistSqr;}
};

typedef complex<float> ArrayBType;
typedef TFFTArray<ArrayBType> TImgArray;

static TImgArray Img,Sum,Weight,Variance;  // a Variance of zero in the image indicates infinity
static TImgArray Tmp1,Tmp2,Tmp3;

TArray3d<float> OrderFac;
int kflag=0;

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-k] [-i inputimg] [-periodx Grid spacing x] [-periody grid spacing Y] [-nx orders along x] [-ny orders along y] [-o outputimg] \n" << flush;
  exit(-1);
}

float FFTPos(float pos, float Size)
{
  return (pos < Size/2) ? pos : pos - Size;
}

void GeneratePattern(double pi, double pj, double kxx,double kxy,double kyx,double kyy, AnOrder * Orders, int NumOrders, TImgArray & Weight)
{
  double kx,ky,Factor;
  int i,j,e,d;
  int DimX = Weight.GetSize(0);
  int DimY = Weight.GetSize(1);

  Weight.Clear();
  d=0;
  for (int num=0;num <NumOrders;num++)
      {
	cout << "\nGenerating Patterns : " << num << ", Phases ("<<pi<<","<<pj<<") \n";
	for (int y=0;y < DimY; y++)
	  for (int x=0;x < DimX; x++)
	    for (e=num; fabs(Orders[e].DistSqr-Orders[num].DistSqr) < 0.01;e++)
	      {
		i = Orders[e].X;
		j = Orders[e].Y;
		Factor = Orders[e].OrderFactor;
		kx = i*kxx + j*kyx;
		ky = i*kxy + j*kyy;
		Weight.SetValue(x,y,d,Weight.Value(x,y,d)+
				ArrayBType(Factor)*exp(ArrayBType(0,-(kx*x+ky*y+ i*pi + j*pj))));
	      }
	d++;
	num = e-1;
      }
  
}

// computes the strenght of the OTF from one order only (after having been placed to [kcx,kcy])
// r is the relative radius to the border of this order
double ComputeOrderStrength(double r, double exponent=1.3)
{
  const double a=2.29,b=2.62,c=1.35;

  if (r > 1.0) return 0.0;
  else 
    {
      double val = 1.0 - a * abs(r) + b * abs(r*r) - c * abs(r*r*r);// estimate optained from a third order fit to a simulation
      // double val = pow(1.0 -  abs(r),exponent); // estimate optained from a third order fit to a simulation
      if (val < 0.0) val = 0.0;
      return val;
    }
  // else return 1.0 -r;       // crude function to estimate the inherent weight included into the Weight image
}

// computes the strenght of the OTF as it should be in the final image
double ComputeFinalStrength(double r)
{
    if (r < 1.0)
      return cos(0.5*Pi * r);
    else
      return 0.0;  
}

void SumOrder( TImgArray & Img, double weight)  // just sum the orders up
{
  int DimX = Img.GetSize(0);
  int DimY = Img.GetSize(1);
  for (int y=0;y < DimY; y++)
    for (int x=0;x < DimX; x++)
      {
	Sum.SetValue(x,y,0,Sum.Value(x,y,0)+Img.Value(x,y,0) * ArrayBType(weight,0.0));
      }
}

// This function determines an order factor from the data, which has allready been placed.
// It is achieved by a weigthed comparison of the complex values in the actual order i,j in relation to
// the information stored in Sum and Weigths.
double DetermineOrderFac(int i, int j, double kxx,double kxy,double kyx,double kyy,
		    double kxmax, double kymax, 
		    double ignoreRadius, // radii are in relative units
		    TImgArray & Img, TImgArray & Weights, TImgArray & Variance)  // weights an order by a function
{
  int DimX = Img.GetSize(0);
  int DimY = Img.GetSize(1);
  float kx = float(DimX) * (kxx * float(i) + kyx * float(j)) / (2.0 * Pi);
  float ky = float(DimY) * (kxy * float(i) + kyy * float(j)) / (2.0 * Pi);
  float datavariance; // weight=1.0,
  double Value1,Value2,Variance1,DataStrength,OrderFac,DataSqr,OrderVar,fac,SumFac,SumWeights;

  cout << "Determining order weigths for order: i " << i << ", j " << j << "\n";

  const double eps = 1e-9;
  for (int y=0;y < DimY; y++)
    for (int x=0;x < DimX; x++)
      {
	float dx = (FFTPos(float(x),float(DimX)) - kx)/(kxmax / 2.0);  // distance to order, when shifted
	float dy = (FFTPos(float(y),float(DimY)) - ky)/(kymax / 2.0);
	float r_order = sqrt(dx*dx+dy*dy);

	float orderstr = ComputeOrderStrength(r_order);
	if (x == 0 && y == 0) cout << "Orderstrength : " << orderstr << " normed radial distance to zero : " << r_order << "\n";
	if (x == 0 && y == 0) 
	  if (r_order >= 1.0) 
	    cout << " Warning !! Radius bigger than 1.0 at central order (" << i << ", " << j << ") radius : " << r_order << "\n";

	if (ignoreRadius > 0)  // The raw variance is assumed to be one in every pixel
	  datavariance = 1.0 / (1.000000001 - exp(-r_order*r_order/ignoreRadius/ignoreRadius));
	else
	  datavariance = 1.0;


	double weight1= real(Weight.Value(x,y,0));
	if (weight1 > 0)
	  {
	    Value1 = abs(Sum.Value(x,y,0)) / weight1;
	    Variance1 = real(Variance.Value(x,y,0)) / weight1 / weight1;
	  }
	else
	  Value1 = 0;

	if (fabs(orderstr) > eps)
	  fac = 1.0 / orderstr;   // factor by which this pixel (x,y) has to be scaled to reach the final result
	else
	  fac = 1.0 / eps;

	Value2 = fac*abs(Img.Value(x,y,0));
	datavariance *= fac*fac;  // std.Dev. scales with fac, thus variance scales with fac^2
	DataStrength = abs(Value2);
	OrderFac = Value1 / Value2;
	DataSqr = DataStrength*DataStrength;
	OrderVar = Variance1 / (DataSqr*DataSqr) + datavariance* abs(Value1)*abs(Value1) / DataSqr;

	if (Value1 > 0 && Value2 > 0)
	  {
	    SumFac += OrderFac / OrderVar;
	    SumWeights += 1.0 / OrderVar;
	  }
      }
  return SumFac / SumWeights;
}

// This function weights each order according to the wanted final strenght
// The associated error in each pixel is magnified accordingly and will be used
// in the averaging procedure during computation of the sum
void WeightOrderSum(int i, int j, double kxx,double kxy,double kyx,double kyy,
		    double kxmax, double kymax, 
		    double ignoreRadius, double strength,   // radii are in relative units
		    TImgArray & Img, TImgArray & Weights, TImgArray & Variance)  // weights an order by a function
{
  int DimX = Img.GetSize(0);
  int DimY = Img.GetSize(1);
  float kx = float(DimX) * (kxx * float(i) + kyx * float(j)) / (2.0 * Pi);
  float ky = float(DimY) * (kxy * float(i) + kyy * float(j)) / (2.0 * Pi);
  float weight=1.0,
    fac=1.0, // factor by which this pixel (x,y) has to be scaled to reach a constant strength
    datavariance;


  // double strength = abs(Img.Value(0,0,0));
  cout << "i " << i << ", j " << j << ", strength " << strength << "\n";

  const double eps = 1e-9;
  for (int y=0;y < DimY; y++)
    for (int x=0;x < DimX; x++)
      {
	float dx = (FFTPos(float(x),float(DimX)) - kx)/(kxmax / 2.0);  // distance to order, when shifted
	float dy = (FFTPos(float(y),float(DimY)) - ky)/(kymax / 2.0);
	float r_order = sqrt(dx*dx+dy*dy);

	float orderstr = strength * ComputeOrderStrength(r_order);
	if (x == 0 && y == 0) cout << "Orderstrength : " << orderstr << " normed radial distance to zero : " << r_order << "\n";
	if (x == 0 && y == 0) 
	  if (r_order >= 1.0) 
	    cout << " Warning !! Radius bigger than 1.0 at central order (" << i << ", " << j << ") radius : " << r_order << "\n";

	if (ignoreRadius > 0)  // The raw variance is assumed to be one in every pixel
	  datavariance = 1.0 / (1.000000001 - exp(-r_order*r_order/ignoreRadius/ignoreRadius));
	else
	  datavariance = 1.0;

	if (fabs(orderstr) > eps)
	  fac = 1.0 / orderstr;   // factor by which this pixel (x,y) has to be scaled to reach constant strength
	else
	  fac = 1.0 / eps;

	datavariance *= fac*fac;  // std.Dev. scales with fac, thus variance scales with fac^2
	if (datavariance <= eps)  // This means the result is actually forced to zero here (weight -> inf)
	  weight = eps; // Does not really matter, since Img was forced to zero
	else
	  weight = 1.0 / datavariance;

	datavariance *= weight * weight; // after application of weight the variance will be reduced
	
	if ((strength == 0.0) || (r_order > 0.95))
	  {
	    weight = 0.0;
	    fac = 0.0;
	    datavariance = 0.0;
	  }

	if (fabs(orderstr) <= eps)
	  Img.SetValue(x,y,0,0.0);
	else
	  Img.SetValue(x,y,0,fac*Img.Value(x,y,0));
	

	Sum.SetValue(x,y,0,Sum.Value(x,y,0)+Img.Value(x,y,0) * ArrayBType(weight,0.0)); // compute sum of weigthed values
	Weights.SetValue(x,y,0,Weights.Value(x,y,0)+ArrayBType(weight,0.0));                // compute sum of weigths, contains sum of inverse variances
	Variance.SetValue(x,y,0,Variance.Value(x,y,0)+ArrayBType(datavariance,0.0));        // add all the variances (after reduction by weight)

	if (x == 0 && y == 0) cout << "Factor(0,0) : " << fac << "  weight : " << weight << "\n";
	// Img.SetValue(x,y,0,finalstr);
	// Img.SetValue(x,y,0,fac); // orderstr);
	// Img.SetValue(x,y,0,Img.Value(x,y,0)*weight);
      }
}


void GenSingleOrderWeight(TImgArray & Img, double kxmax, double kymax)
{
  int DimX = Img.GetSize(0);
  int DimY = Img.GetSize(1);

  for (int y=0;y < DimY; y++)
    for (int x=0;x < DimX; x++)
      {
	double dx = (FFTPos(x,DimX))/(kxmax / 2.0);  // distance to order, when shifted
	double dy = (FFTPos(y,DimY))/(kymax / 2.0);
	double r = sqrt(dx*dx+dy*dy);
	Img.SetValue(x,y,0,ComputeOrderStrength(r));
      }
  Img.FFTbwd();  // This is just needed for visualization
}

void CleanOrder( double ignoreRadius, TImgArray & Img)  // weights an order by a function
{
  int DimX = Img.GetSize(0);
  int DimY = Img.GetSize(1);
  double val=1.0;

  Img.FFTfwd();
  for (int y=0;y < DimY; y++)
    for (int x=0;x < DimX; x++)
      {
 	double dx = (FFTPos(x,DimX))/(ignoreRadius / 2.0);  // distance to order, when shifted
	double dy = (FFTPos(y,DimY))/(ignoreRadius / 2.0);
	double r = sqrt(dx*dx+dy*dy);
	if (ignoreRadius > 0)
	  val = 1.0 - exp(-r*r/ignoreRadius/ignoreRadius); 
	else
	  val = 1.0;
  
	Img.SetValue(x,y,0,ArrayBType(val,0.0)*Img.Value(x,y,0));
     }
  Img.FFTbwd(); 
}



void ApplyFinalWeight(TImgArray & Sum, TImgArray & Weight,double ApoRadiusX,double ApoRadiusY) 
{
  int DimX = Sum.GetSize(0);
  int DimY = Sum.GetSize(1);
  double dx,dy,sumweight,finalstr,r;

  for (int y=0;y < DimY; y++)
    for (int x=0;x < DimX; x++)
      {
	// float finalstr =  ComputeOrderStrength(r);  // Apodization of the final result using a WF apodization function
	dx = FFTPos(float(x),float(DimX))/(ApoRadiusX*DimX / 2.0);  
	dy = FFTPos(float(y),float(DimY))/(ApoRadiusY*DimY / 2.0);
	r = sqrt(dx*dx+dy*dy);
	finalstr =  ComputeFinalStrength(r);  // Apodization of the final result with the normalized radius
	sumweight = real(Weight.Value(x,y,0));

	if (sumweight != 0 && finalstr != 0)
	  Sum.SetValue(x,y,0,Sum.Value(x,y,0)/ArrayBType(sumweight)*ArrayBType(finalstr));
	else
	  Sum.SetValue(x,y,0,0.0);
      }
  Sum.FFTbwd();
}

/// This function places a separated order at the correct position in k-space
void PlaceOrder (int i, int j, double kxx,double kxy,double kyx,double kyy,
		       bool automatic, double movefirstX, double movefirstY, 
		       double Factor, TImgArray & Img)   // i : orderx, j : ordery   element -N ... N
{
  ArrayBType val,kx,ky;
  ArrayBType PhaseFac = exp(ArrayBType(0,-1.0* (i*movefirstX + j*movefirstY)));

  int DimX = Img.GetSize(0);
  int DimY = Img.GetSize(1);

  kx = ArrayBType(0,- 1.0 * (kxx * i + kyx * j));
  ky = ArrayBType(0,- 1.0 * (kxy * i + kyy * j));

  for (int y=0;y < DimY; y++)
    for (int x=0;x < DimX; x++)
      {
	val = Img.Value(x,y,0) *
	  exp(ArrayBType(y,0.0) * ky + ArrayBType(x,0.0) * kx);  // Cpx(0,1.0) is included in kx and ky
	val *= PhaseFac;
	Img.SetValue(x,y,0, val);
	// Sum.SetValue(x,y,0, Sum.Value(x,y,0)+ weight*val*float(Factor));
      }
}

int main(int argc, char *argv[])
{ 
static int INPUTSizeX=32;  
static int INPUTSizeY=32; 
static int INPUTSizeZ=32; 
 int Elements=1,ordersx=1,ordersy=1, totalordersX=3, totalordersY=3, 
   aorderx=1, aordery=1,  // number of the order to use for automatic estimation
   InputSizeX=64, InputSizeY=64;
int i,j,elem=0,doautomatic=0,doautomaticmove=0,doApoAniso=0,PatStepsX=3,PatStepsY=3;
 double Factor=1.0,periodx=4.0, periody=4.0, OrderFactor=1.0, zeroweight=1.0, // default : do not weight zero order
  angle=0.0, movefirstX=0.0,movefirstY=0.0,
  kxmax=0.5, kymax=0.5,ignoreRadius=0.05, apoRadius=0.9, apoRadiusY=0.9;

float  kxx=0,kxy=0,kyx=0,kyy=0;

string OFile = "/tmp/shifted.raw";

string OUTPUTFileName,IFileName,OFileName,OrderFileName,SumFileName,PatFileName;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=1;continue;}
   if (readArg("-auto",parg)) {doautomatic=1;continue;}
   if (readArg("-autom",parg)) {doautomaticmove=1;continue;}
   if (readArg("-aorderx", & aorderx, parg)) continue; // lowest kx-vector will be computed from this given in pixels
   if (readArg("-aordery", & aordery, parg)) continue; 
   if (readArg("-kxmax", & kxmax, parg)) continue; // maximal distance along kx to account in each order, given as fraction
   if (readArg("-kymax", & kymax, parg)) continue; 
   if (readArg("-ir", & ignoreRadius, parg)) continue; // will now be applied to all orders (including zero order)
   if (readArg("-ar", & apoRadius, parg)) continue; 
   if (readArg("-ary", & apoRadiusY, parg)) {doApoAniso=1;continue;} // if given, the apodization can be anisotropic
   if (readArg("-nx", & ordersx, parg)) continue; // orders to be used (+-)
   if (readArg("-ny", & ordersy, parg)) continue; // orders to be used (+-)
   if (readArg("-periodx", & periodx, parg)) continue; // lowest kx-vector will be computed from this given in pixels
   if (readArg("-periody", & periody, parg)) continue; 
   if (readArg("-angle", & angle, parg)) continue;   // the vectors will be rotated by this angle
   if (readArg("-mx",& movefirstX, parg)) continue;   // move the (one,0) order by this phase angle and others appropiately
   if (readArg("-my",& movefirstY,parg)) continue;    // move the (0,one) order by this phase angle (in rad)
   if (readArg("-ix",& InputSizeX, parg)) continue;
   if (readArg("-iy",& InputSizeY,parg)) continue;
   if (readArg("-i", IFileName, parg)) continue;
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-pat", PatFileName, parg)) continue;
   if (readArg("-PStepsX", & PatStepsX, parg)) continue;
   if (readArg("-PStepsY", & PatStepsY, parg)) continue;
   if (readArg("-sum", SumFileName, parg)) continue;
   if (readArg("-of", & OrderFactor, parg)) continue;
   if (readArg("-zw", & zeroweight, parg)) continue;
   if (readArg("-OrderFactors", OrderFileName, parg)) continue;
    usage(argv[0]);
  }

 if (! doApoAniso) apoRadiusY = apoRadius;
 
  if (IFileName=="") 
    { cerr << "Error : Input [-i name] is required!\n"; exit (-1);}

  if (OFileName=="") OUTPUTFileName=OFile;
    else OUTPUTFileName=OFileName;

  totalordersX=2*ordersx+1;
  totalordersY=2*ordersy+1;
  cout << "total orders in each direction is " << totalordersX*totalordersY << "\n";

  elem =0;

  ofstream to(OUTPUTFileName.c_str());
      
  if (! to )
    {
      cerr << "Couldn't open file " << OUTPUTFileName << " for writing !!\n" << flush;
      exit(-1);
    }
	  
    if (periodx != 0)
    {
  kxx = 2.0 * Pi / periodx * cos(angle * Pi / 180.0);
  kxy = 2.0 * Pi / periodx * sin(angle * Pi / 180.0);
	}
if (periody != 0)
   {
   kyx = 2.0 * Pi / periody * (- sin(angle * Pi / 180.0));
   kyy = 2.0 * Pi / periody * cos(angle * Pi / 180.0);
   }

  if (OrderFileName != "")
    {
      OrderFac.DLoad(kflag,OrderFileName.c_str(),"Float", & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,0);
      if ((INPUTSizeX != totalordersX) || (INPUTSizeY != totalordersY))
	{
	  cerr << "Error : wrong number of Orders in OrderFile, was expecting " << totalordersX << " x " << totalordersY << "\n";
	    exit (-1);
	}
    }

  int zero = (totalordersX*totalordersY)/2;
  int oneX = zero+aorderx;
  int oneY = zero+aordery*totalordersX;

  if (doautomatic)
    {  // Determines the grid spacing as well as the relative position of the first grating from the data
      // The zero component is found in slice int((totalordersx*totalordersy)/2)
      double max,Xmx,Xmy,Ymx,Ymy,dummy;

      Elements=Tmp1.DLoad(kflag,IFileName.c_str(),"Complex",
			  & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,zero);
      Weight.Resize(INPUTSizeX,INPUTSizeY,1); // will store the sum of the weights
	Tmp3.Resize(INPUTSizeX,INPUTSizeY,1); 
      // Tmp2.Resize(INPUTSizeX,INPUTSizeY,1);
	  
      if (ordersx > 0)
	{
	  Elements=Tmp2.DLoad(kflag,IFileName.c_str(),"Complex",
			      & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,oneX);

	  Tmp3.Copy(& Tmp2);
	  // cout << "ir : " << ignoreRadius << " ir2 : " <<ignoreRadius * (Tmp2.GetSize(0)+Tmp2.GetSize(1))/2.0;
	  CleanOrder( ignoreRadius* (Tmp2.GetSize(0)+Tmp2.GetSize(1))/2.0,Tmp2); 
	  
	  Tmp2.ConvMul(& Tmp1,true);
	  Tmp2.FFTfwd();  // now the cross-correlation was calculated
	  Tmp2.Magnitude();  // should not be necessary, but it is necessary

	  	  /*
	  // Normalize for sum of powerspectra
	  Tmp1.ConvMul(& Tmp1,true);
	  Tmp1.FFTfwd();  // now the cross-correlation was calculated
	  Tmp1.Magnitude();  // should not be necessary, but it is necessary
	  Tmp3.ConvMul(& Tmp3,true);
	  Tmp3.FFTfwd();  // now the cross-correlation was calculated
	  Tmp3.Magnitude();  // should not be necessary, but it is necessary
	  Tmp3.Add(& Tmp1);
	  Tmp3.Mean();
	  Tmp3.Add(abs(Tmp3.Mean())/2);
	  Tmp3.ArgDivSelf(& Tmp2);*/
	  	  
	  Tmp2.DSave(true,"/tmp/CorrelX1.kdf");
	  max = Tmp2.WrappedMaxIntensityPosition(&Xmx,&Xmy,&dummy,1,true); // range = 1, forceplanes (do not seach for shifts in Z)
	  // WrappedMax... will account for the real part only
	  
	  cout << "Determined X-shift to " << Xmx << ", " << Xmy << ", " << dummy << "\n";
	  kxx = - 2.0 * Pi * (Xmx/double(aorderx)/INPUTSizeX)  ;
	  kxy = - 2.0 * Pi * (Xmy/double(aorderx)/INPUTSizeY);
	}
      else
	{
	  kxx = 0;
	  kxy = 0;
	}

      if (ordersy > 0)
	{
	  Elements=Tmp2.DLoad(kflag,IFileName.c_str(),"Complex",
			      & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,oneY);
	  
	  CleanOrder( ignoreRadius* (Tmp2.GetSize(0)+Tmp2.GetSize(1))/2.0,Tmp2); 
	  
	  Tmp2.ConvMul(& Tmp1,true);
	  Tmp2.FFTfwd();  // now the cross-correlation was calculated
	  Tmp2.Magnitude();  // should not be necessary, but it is necessary
	  Tmp2.DSave(true,"/tmp/CorrelY.kdf");
	  max = Tmp2.WrappedMaxIntensityPosition(&Ymx,&Ymy,&dummy,1,true); // range = 1, forceplanes=true (do not seach for shifts in Z)
	  // WrappedMax... will account for the real part only
	  
	  cout << "Determined Y-shift to " << Ymx << ", " << Ymy << ", " << dummy << "\n";

	  kyx = - 2.0 * Pi * (Ymx/double(aordery)/INPUTSizeX);
	  kyy = - 2.0 * Pi * (Ymy/double(aordery)/INPUTSizeY);
	}
      else
	{
	  kyx = 0;
	  kyy = 0;
	}
    }

  if (doautomaticmove || doautomatic)
    {  // Determines the grid spacing as well as the relative position of the first grating from the data
      if (ordersx > 0)
	{
	  Elements=Tmp2.DLoad(kflag,IFileName.c_str(),"Complex",
			      & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,zero+1);
	  Sum.Resize(INPUTSizeX,INPUTSizeY,1); 

	  CleanOrder(ignoreRadius* (Tmp2.GetSize(0)+Tmp2.GetSize(1))/2.0,Tmp2); 
	  PlaceOrder (1, 0, kxx,kxy,kyx,kyy, doautomatic, 0.0, 0.0, 1.0, Tmp2); 
	  //   WeightOrderSum(1,0,kxx,kxy,kyx,kyy,  // ignore these argument, since Weights will be cleared anyway
	  // 	  kxmax * Tmp2.GetSize(0), kymax * Tmp2.GetSize(1),
	  // 	  ignoreRadius,1.0,
	  // 	  apoRadius, Tmp2,Weight); 

	  Tmp2.FFTfwd(); 
	  movefirstX = atan2(imag(Tmp2.Value(0,0,0)), real(Tmp2.Value(0,0,0)));
	  // cout << "X move components " << Tmp2.Value(0,0,0) << " DX " << movefirstX << "\n";
	}

      if (ordersy > 0)
	{
	  // cout << "Loading Element " << zero + totalordersX<< "\n";
	  Elements=Tmp2.DLoad(kflag,IFileName.c_str(),"Complex",
			 & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,zero + totalordersX);

	  CleanOrder(ignoreRadius* (Tmp2.GetSize(0)+Tmp2.GetSize(1))/2.0,Tmp2); 
	  PlaceOrder (0, 1, kxx,kxy,kyx,kyy, doautomatic, 0.0, 0.0, 1.0, Tmp2); 
	  // WeightOrderSum(0,1,kxx,kxy,kyx,kyy,  // ignore these argument, since Weights will be cleared anyway
	  // 		  kxmax * Tmp2.GetSize(0), kymax * Tmp2.GetSize(1),
	  // 	  ignoreRadius,1.0,
	  // 	  apoRadius, Tmp2,Weight); 
	  
	  Tmp2.FFTfwd();  
	  movefirstY = atan2(imag(Tmp2.Value(0,0,0)), real(Tmp2.Value(0,0,0)));
	  cout << "Y move components " << Tmp2.Value(0,0,0) << " DY " << movefirstY << "\n";
	}
    } 


  cout << "kxx " << kxx << ", kxy " << kxy << "\n";
  cout << "kyx " << kyx << ", kyy " << kyy << "\n";
  cout << "movefirstX " << movefirstX << ", movefirstY " << movefirstY << "\n";

  // Now a list with the orders and their distances is generated, sorted and traversed

  int num=0;
  AnOrder Orders[MaxOrders];
  for (j=-ordersy;j<=ordersy;j++)   // starting top left order
    for (i=-ordersx;i<=ordersx;i++)
      {
	Orders[num].X = i;
	Orders[num].Y = j;
	Orders[num].Elem = num;
	Orders[num].DistSqr = double(i)*i+double(j)*j;
	num++;
      }

  sort(Orders,&Orders[totalordersX*totalordersY],CmpOrders());

  for (num=0;num <totalordersX*totalordersY;num++)
      {
	i = Orders[num].X;
	j = Orders[num].Y;
	elem = Orders[num].Elem;
	Factor = Orders[num].OrderFactor;
	cout << "\nProcessing : " << num << ", Order ("<<i<<","<<j<<"), element "<<elem<<":\n";
  
	Elements=Img.DLoad(kflag,IFileName.c_str(),"Complex",
			   & INPUTSizeX,& INPUTSizeY,& INPUTSizeZ,elem);


	if (num==0)  // first run
	  {
	    if (kflag)
	      {
		if (kflag) WriteKhorosHeader(& to,"Generated by mulwave 2001","Complex",INPUTSizeX,INPUTSizeY,INPUTSizeZ,Elements);
		cerr << "Writing output file header \n" << flush;
	      }
	    // else nothing
	    Sum.Resize(INPUTSizeX,INPUTSizeY,1);
	    Weight.Resize(INPUTSizeX,INPUTSizeY,1); // will store the sum of the weights
	    Variance.Resize(INPUTSizeX,INPUTSizeY,1); // will store the total variance of each pixel stored in Sum

	    Sum.Clear();  // might have been set to something in autoalign procedure
	    Weight.Clear(); 
	    Variance.Clear(); 
	  }
	if (Elements != totalordersX*totalordersY)
	  {
	    cerr << "Error : wrong number of Elements to shift! Was expecting " << totalordersX*totalordersY << "\n";
	    exit (-1);
	  }

	PlaceOrder (i, j, kxx,kxy,kyx,kyy, doautomatic, movefirstX, movefirstY, Factor, Img); 
	if (apoRadius > 0.0)
	  {
	  Img.FFTfwd();
	  if (OrderFileName!="")
	    Factor = OrderFac.Value(i+ordersx,j+ordersy,0);
	  else
	      {
		if ((zeroweight == 0) && (abs(i)+abs(j) == 1))
		  Factor = 1.0;
		else if (num != 0) //  && Factor==0)
		  Factor=DetermineOrderFac(i, j, kxx,kxy,kyx,kyy,
					   kxmax * Img.GetSize(0), kymax * Img.GetSize(1),
					   ignoreRadius, // radii are in relative units
					   Img, Weight, Variance);
		else
		  Factor = 1.0;

		for (int e=num; fabs(Orders[e].DistSqr-Orders[num].DistSqr) < 0.01;e++)
		  Orders[e].OrderFactor = Factor;
	      }

	  if (num == 0) 
	    Factor*=zeroweight;  // Treat the zero order separately

	  Factor *= pow(OrderFactor,sqrt(Orders[num].DistSqr));
	  cout << "Order Factor is: "<< Factor <<"\n";

	  WeightOrderSum(i,j,kxx,kxy,kyx,kyy,
		       kxmax * Img.GetSize(0), kymax * Img.GetSize(1),
		       ignoreRadius,Factor, 
		       Img,Weight,Variance); 
	  // if (num == 1) Sum.DSave(kflag,"/tmp/Sum.kdf");
	  Img.FFTbwd();
	  }
	else
	  {
	    if (OrderFileName!="")
	      Factor = OrderFac.Value(i+ordersx,j+ordersy,0);
	    else
	      Factor = 1.0;

	    cout << "Apodization Radius was <= 0 summing all orders\n";
	    cout << "Order Factor is: "<< Factor <<"\n";
	    SumOrder(Img,Factor);
	  }

	Img.Write(& to);

   }

  if (apoRadius > 0.0)
    {
      ApplyFinalWeight(Sum,Weight,apoRadius,apoRadiusY);
      // GenSingleOrderWeight(Weight,kxmax * Img.GetSize(0), kymax * Img.GetSize(1));
      // Weight.DSave(kflag,"/tmp/Weight.kdf"); Variance.DSave(kflag,"/tmp/Variance.kdf");
    }

  if (SumFileName!="")
    Sum.DSave(kflag,SumFileName.c_str());

  if (PatFileName!="")
    {
      ofstream pto(PatFileName.c_str());
      
      if (! pto )
	{
	  cerr << "Couldn't open file " << PatFileName << " for writing !!\n" << flush;
	  exit(-1);
	}
      int d=0,e=0;
      for (int num=0;num <totalordersX*totalordersY;num++)  // count the components of equal abs(k), having equal z-behaviour
	{
	  for (e=num; fabs(Orders[e].DistSqr-Orders[num].DistSqr) < 0.01;e++)
	    ;
	  num=e-1;
	  d++;
	}
      Weight.Resize(Weight.GetSize(0),Weight.GetSize(1),d);	  
      if (kflag) WriteKhorosHeader(& pto,"Generated by mulwave 2001","Complex",Weight.GetSize(0),Weight.GetSize(1),Weight.GetSize(2),PatStepsX*PatStepsY);

      for (int pj=0;pj<PatStepsY;pj++)   // Iterate over the pattern phases
	for (int pi=0;pi<PatStepsX;pi++)
	  {
	    double PhaseX=2*M_PI*pi/PatStepsX + movefirstX;
	    double PhaseY=2*M_PI*pj/PatStepsY + movefirstY;
	    GeneratePattern(PhaseX,PhaseY, kxx, kxy, kyx,kyy, Orders,totalordersX*totalordersY , Weight);
	    Weight.Write(& pto);
	  }
      pto.close();
    }

}
