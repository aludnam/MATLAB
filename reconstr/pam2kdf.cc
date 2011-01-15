
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

#include <stdio.h>
#include <stdlib.h>
#include "khoros.h"
#include "parseargs.h"
#include "rawarray.h"

#include <iostream>
#include <string>
#include <math.h>

typedef float ArrayBType;
typedef unsigned short OrigType;


typedef TArray3d<ArrayBType>  TImgArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TArray3d<OrigType>  TOrigArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...
typedef TArray3d<double>  TResArray; // All arrays have the same FFTabel format, so that they can be compared, divided, ...

TOrigArray MyData;
TImgArray MyOutData,MyOutData2;
TResArray Results;  // For computation of plane mean and variance

int ISizeX,ISizeY,ISizeZ,MoreSlices=0;
int ROI_TLX,ROI_TLY,ROI_BRX,ROI_BRY,BgCalib ,MySensCMean=0, BgMean=0,NumFrames=1,CamNr=0,Elements=1;
int PatXSize=0,PatYSize=0,Every=0,BlankN=0;
int ROI1xs,ROI1ys,ROI1xe,ROI1ye;
int ROI2xs,ROI2ys,ROI2xe,ROI2ye;

bool ClipLimits=true,BgCor=false,SensCor=false,IsPseudoConf=false,IsHadamard=false,dosubslice=false,bleechcor=false;
bool FRAP=false;  // if true only the FRAP ROI will be evaluated
float ConstAdd=0;
float ConstAdd0=0;
float ConstAdd1=0;
int subslice=0; // offset for subslicing (when m_every > 1)

void usage(char * filename)
{
    	fprintf(stderr, "usage: %s [-k] -o outputfile -i inputfile\n",filename);
    	fflush (stderr); 
	exit(-1);
}

double sqr(double val)
{
  return val*val;
}

void GetFRAPData(char * FileName,char *  HeaderFile)
{
  Results.Resize(ISizeZ/Every,1,1);  // For mean
  int ROISizeX,ROISizeY;
  int elem, Elements;
  double Sum=0;
  char * tmp,* tmp2;
  ofstream to(FileName);


  // ROI1xs= ROI_TLX,ROI1ys=ROI_TLY,ROI1xe=ROI_BRX,ROI1ye=ROI_BRY;

  tmp=strstr(HeaderFile,"UseROI");  // Move forward to the first unused ROI, which is a FRAP ROI

  Elements = 0;
  for (;;) {  // Just count the number of FRAP ROIs
    tmp2=getVal("MyROI.TopLeft().x",tmp,& ROI2xs,true);
    cout << ROI2xs << "\n";
    // cout << tmp2;
    if (tmp2 == tmp) break;
    Elements ++;
    tmp = tmp2;
    tmp ++;  // to ensure that this one is not catched again !
  }
  cout << Elements << " ROIs found \n";

  tmp=strstr(HeaderFile,"UseROI");  // Start again now
  Results.DHeader(true,to,Elements);
  for (elem = 0; elem < Elements; elem ++)
    {
      tmp=getVal("MyROI.TopLeft().x",tmp,& ROI2xs);
      tmp=getVal("MyROI.TopLeft().y",tmp,& ROI2ys);
      tmp=getVal("MyROI.BottomRight().x",tmp,& ROI2xe);
      tmp=getVal("MyROI.BottomRight().y",tmp,& ROI2ye) ;

      ROISizeX = ROI2xe- ROI2xs;
      ROISizeY = ROI2ye- ROI2ys;
      cout << "Processing FRAP ROI # " << elem << ". ROI sizes are " << ROISizeX << " x " << ROISizeY << "\n";

    for (int z=0;z <ISizeZ; z++)   // take every section
      {
	Sum=0;
	for (int y=ROI2ys-ROI1ys;y < ROI2ye - ROI1ys; y++)
	  for (int x=ROI2xs-ROI1xs;x < ROI2xe - ROI1xs; x++)
	    {
	      Sum += MyOutData.Value(x,y,z);
	    }
	Sum /= ROISizeX*ROISizeY;
	cout << "Plane " << z << " mean : " << Sum << "\n";
	Results.SetValue(z,0,0,Sum);
      }

    Results.Write(& to);
  }
}

void ProcessPlanes(char * FileName, double Kappa)
{
  double Factor=Kappa/(1.0-Kappa);
  double Variance=0,Mean=0,Sum=0,SumSqr=0,dat1=0,dat2=0,Sum1=0,Sum2=0,SumSqr1=0,SumSqr2=0;
  double SumW=0,SumSqrW=0,SumSqrWO=0,MeanDiffW=0,MeanW=0,MeanWO=0,VarianceW=0,VarianceWO=0,SumWO=0;
  double VC1,VC2,VNC1,VNC2,VW1,VW2;
  double Mean1,Mean2,Var1,Var2,MeanDiff1,MeanDiff2=0,MeanDiff=0,MeanDiffWO=0,
    PlaneSum=0,FirstPlaneMean1=0,FirstPlaneMean2=0,FirstPlaneMean=0,LastPlaneMean=0; // FirstPlaneMean=0;

  double BleechTau=0.0;  // if bleechcor is selected, the factor will decrease exponentially
  int OtherIdent=1;  // when every 2nd is repeated, next slice is at identical position

  cout << "Kappa : " << Kappa << ", Factor : "<< Factor << "\n";

  if ((Every == 4) && (BlankN ==2))    // 2 repetitions of widefield interleaved
    {
      OtherIdent=2; 
      Results.Resize(ISizeZ/Every,12,1);  // For mean and variance
    }
  else if ((Every == 2) && (BlankN <= 1))  // Two similar measurements with or without widefields
    {
      OtherIdent=1;
      Results.Resize(ISizeZ/Every,9,1);  // For mean and variance
    }
  else   // Just evaluate the data mean and no variances
    {
      OtherIdent=0;  // just use the same plane (no stddev statistics)
      cerr << "Warning File has not exactly 2 or 4 repeats !";
      cerr << "m_Every is "<< Every << ", interpreting data slicewise and computing only means \n";
      Results.Resize(ISizeZ/Every,5,1);  // For meanC and meanNC, Sum and WF
    }


  // calculate mean of first and last planes
  double SumC=0,SumNC=0,SumLC=0,SumLNC=0;
  for (int y=0;y <ISizeY; y++)
    for (int x=0;x <ISizeX; x++)
      {
	VC1 = MyOutData.Value(x,y,0); VC2 = MyOutData.Value(x,y,0+OtherIdent);
	VNC1 = MyOutData2.Value(x,y,0); VNC2 = MyOutData2.Value(x,y,0+OtherIdent);
	SumC+=(VC1+VC2)/2.0;  // Mean of conjugate img 1 and img 2
	SumNC+=(VNC1+VNC2)/2.0;  // Mean of non-conjugate

	VC1 = MyOutData.Value(x,y,ISizeZ-Every); VC2 = MyOutData.Value(x,y,ISizeZ-Every+OtherIdent);
	VNC1 = MyOutData2.Value(x,y,ISizeZ-Every); VNC2 = MyOutData2.Value(x,y,ISizeZ-Every+OtherIdent);
	SumLC+=(VC1+VC2)/2.0;  // Mean of conjugate img 1 and img 2
	SumLNC+=(VNC1+VNC2)/2.0;  // Mean of non-conjugate
      }
  SumC /= double(ISizeX*ISizeY);   // mean of conjugate
  SumNC /= double(ISizeX*ISizeY);   // mean of non-conjugate
  SumLC /= double(ISizeX*ISizeY);   // mean of conjugate
  SumLNC /= double(ISizeX*ISizeY);   // mean of non-conjugate

  FirstPlaneMean=SumC + SumNC;
  LastPlaneMean=SumLC + SumLNC;
  FirstPlaneMean1=FirstPlaneMean * Kappa;  // Estimate from the Widefield sum
  FirstPlaneMean2=FirstPlaneMean * (1-Kappa);
  cout << "calculated FirstPlane (M1+M2)*kappa: Mean1 = " << FirstPlaneMean1 << ", Mean2 (M1+M2*(1-kappa))= " << FirstPlaneMean2 << "\n";
  cout << "FirstPlane Mean of Sum signal = " << FirstPlaneMean << "\n";

  if (bleechcor)
    {
      BleechTau = - log(LastPlaneMean/FirstPlaneMean) / (ISizeZ-1);  // Will be applied only to signal as exp(- BleechTau * z)
      cout << "Bleeching correction will be applied ! BleechTau = "<<BleechTau<<" is " << exp(- BleechTau * (ISizeZ-1))*100 << "% in total stack\n";
    }

  for (int z=0;z <ISizeZ; z+=Every)   // normally 2, but may be 4 (when Wilson interleaved)
    {Variance=0,Sum=0,SumSqr=0;Sum1=0;Sum2=0,SumSqr1=0,SumSqr2=0;MeanDiff1=0;MeanDiff2=0;MeanDiff=0;
    SumW=0;SumSqrW=0;SumSqrWO=0;MeanDiffW=0;MeanDiffWO=0;VarianceW=0;VarianceWO=0;MeanW=0;SumWO=0;MeanWO=0;
    if (bleechcor)
      PlaneSum = FirstPlaneMean*exp(-BleechTau * z);  // simulate bleeching (FWD)
    else
      PlaneSum = FirstPlaneMean;
    for (int y=0;y <ISizeY; y++)
      for (int x=0;x <ISizeX; x++)
	{
	  VC1 = MyOutData.Value(x,y,z); VC2 = MyOutData.Value(x,y,z+OtherIdent);
	  VNC1 = MyOutData2.Value(x,y,z); VNC2 = MyOutData2.Value(x,y,z+OtherIdent);

	  Sum1+=(VC1+VC2)/2.0;  // Mean of conjugate img 1 and img 2
	  SumSqr1+=sqr(VC1-VC2);  // Mean of squared difference between img1 and img2
	  MeanDiff1+=(VC1-VC2);  // Average Difference of conjugate

	  Sum2+=(VNC1+VNC2)/2.0;  // Mean of non-conjugate
	  SumSqr2+=sqr(VNC1-VNC2);  //  Mean of squared difference between img1 and img2
	  MeanDiff2+=(VNC1-VNC2);  // Average Difference of non-conjugate

	  // old : dat1=VC1 - VNC1*Factor;    // conjugate - non conjug
	  // old : dat2=VC2 - VNC2*Factor;  // now second identical data
	  // compute the weighted average
	  dat1=VC1*(1.0-Kappa) + (PlaneSum - VNC1)*Kappa;    // conjugate - non conjug
	  dat2=VC2*(1.0-Kappa) + (PlaneSum - VNC2)*Kappa;  // now second identical data
	  Sum+=(dat1+dat2)/2.0; 
	  MeanDiff+=(dat1-dat2); 
	  SumSqr += sqr(dat1-dat2);
	  // data collected below is only for the Wilson comparison
	  if (BlankN == 2)  // A Wf was taken interleaved
	    {
	      VW1 = MyOutData.Value(x,y,z+1); VW2 = MyOutData.Value(x,y,z+1+OtherIdent);
	      SumWO+=(VW1+VW2)/2.0; // just original widefield data without extra processing
	      MeanDiffWO+=(VW1-VW2); 
	      dat1=(VC1 - VW1*Kappa);  // conjugate - non conjug, scaled by two  // old : 2*
	      dat2=(VC2 - VW2*Kappa);  // now second identical data, scaled by two  // old : 2*
	      SumW+=(dat1+dat2)/2.0; 
	      MeanDiffW+=(dat1-dat2); 
	      SumSqrW += sqr(dat1-dat2);
	      SumSqrWO += sqr(VW1-VW2);   // variance of WF only
	    }
	}
    Mean1=Sum1/double(ISizeX*ISizeY);   // mean of conjugate
    Mean2=Sum2/double(ISizeX*ISizeY);   // mean of non-conjugate
    MeanDiff1/=double(ISizeX*ISizeY);   // calculate average difference
    MeanDiff2/=double(ISizeX*ISizeY);
    Var1=SumSqr1/(ISizeX*ISizeY)-sqr(MeanDiff1);
    Var2=SumSqr2/(ISizeX*ISizeY)-sqr(MeanDiff2);

    Mean = Sum/(ISizeX*ISizeY);         // mean of weighted average image
    MeanDiff/=double(ISizeX*ISizeY);
    Variance = SumSqr/(ISizeX*ISizeY)-sqr(MeanDiff);
    if (BlankN == 2)  // Only for Wilson comparison or Widefield data
      {
	MeanW = SumW/(ISizeX*ISizeY);
	MeanWO = SumWO/(ISizeX*ISizeY);
	MeanDiffW/=double(ISizeX*ISizeY);
	MeanDiffWO/=double(ISizeX*ISizeY);
	VarianceW = SumSqrW/(ISizeX*ISizeY)-sqr(MeanDiffW);
	VarianceWO = SumSqrWO/(ISizeX*ISizeY)-sqr(MeanDiffWO);
      }




    if (bleechcor)
      {
	Mean1 /= exp(- BleechTau * z);
	Mean2 /= exp(- BleechTau * z);
	Mean /= exp(- BleechTau * z);
	MeanW /= exp(- BleechTau * z);
	MeanWO /= exp(- BleechTau * z);
      }

    if (OtherIdent != 0)
      {
	Results.SetValue(z/Every,6,0,Mean1);    // The bleeching correction is not applied to this data
	Results.SetValue(z/Every,7,0,Mean2);
	Results.SetValue(z/Every,8,0,Mean1+Mean2);
      }

    if (OtherIdent != 0)   // statistic variances were collected
      {
	Results.SetValue(z/Every,0,0,Mean1);  // old : (Mean1-FirstPlaneMean1)*(1+Factor));
	Results.SetValue(z/Every,1,0,FirstPlaneMean - Mean2);  // old : (FirstPlaneMean2-Mean2)*(1+Factor));
	Results.SetValue(z/Every,2,0,Mean);
	Results.SetValue(z/Every,3,0,sqrt(Var1));  // old : (1+Factor)*sqrt(Var1)
	Results.SetValue(z/Every,4,0,sqrt(Var2));  // old : (1+Factor)*sqrt(Var1)
	Results.SetValue(z/Every,5,0,sqrt(Variance)); // old : sqrt(Variance)

	if (BlankN == 2)  // Only for Wilson comparison or Widefield data
	  {
	    // Results.SetValue(z/Every,9,0,MeanWO);    
	    Results.SetValue(z/Every,10,0,sqrt(VarianceW));
	    Results.SetValue(z/Every,11,0,sqrt(VarianceWO));  // standard deviation of widefield only
	    Results.SetValue(z/Every,9,0,MeanW);    
	  }      
      }
    else  // just the bleaching-corrected raw data means
      {
	Results.SetValue(z/Every,0,0,Mean1);  // just the mean in Plane C, but eventually bleeching corrected
	Results.SetValue(z/Every,1,0,Mean2);  // just the mean in Plane NC, but eventually bleeching corrected
	Results.SetValue(z/Every,2,0,Mean1+Mean2);  // just the mean in Plane NC, but eventually bleeching corrected
	Results.SetValue(z/Every,3,0,FirstPlaneMean - Mean2);  // Mean 2 Flipped over
	Results.SetValue(z/Every,4,0,MeanWO);   
      }
    }    // for z = ...
  // cout << "FirstPlaneMean is " << FirstPlaneMean << "\n";
  if (OtherIdent)
    cout << "Output file contains along Y : \nMean1, \n(FirstPlaneMean2-Mean2),\nMean,\nsqrt(Var1),\nsqrt(Var2),\nsqrt(Variance)\n";
  else
    cout << "Output file contains along Y : \nMean1, \nMean2,\nMean1+Mean2,\nFirstPlaneMean - Mean2\n";
  // cout << "Output file contains along Y : \n(Mean1-FirstPlaneMean1)*(1+Factor), \n(FirstPlaneMean2-Mean2)*(1+Factor),\nMean,\n(1+Factor)*sqrt(Var1),\n(1+Factor)*sqrt(Var2),\nsqrt(Variance)\n";
  Results.DSave(true,FileName);
}

void CalculateMeans(double & SensMean,double & BgMean)
{
  SensMean=0.0;
  BgMean=0.0;
  for (int y=0;y <ISizeY; y++)
  for (int x=0;x <ISizeX; x++)
    if (BgCalib)
      {
	BgMean += MyData.Value(x,y,0);  // In this case the first allways the background
	if (MySensCMean)
	  SensMean += MyData.Value(x,y,1);  // In this case the first allways the background
      }
  else
    if (MySensCMean)
      SensMean += MyData.Value(x,y,0);  // In this case the first allways the background

  BgMean /= double(ISizeX)*double(ISizeY);
  SensMean /= double(ISizeX)*double(ISizeY);
  cout << "Using newly calculated means of ROI : Bg.:"; 
  if (BgCalib && BgCor)
       cout << BgMean;
  else
    cout << "not applied";
  cout << " Sens.: ";
  if (MySensCMean && SensCor)
       cout << SensMean;
  else
    cout << "not applied";
  cout << "\n";
}

void ProcessStack(char * FileName,char * Ext, TImgArray & OutData)
{
  float dat, sens;
  double SensMean,BgMean;

  strcat(FileName,Ext);  // append this to the Filename for loading the first image

  OrigType dummy;
  int Slices=ISizeZ+MoreSlices;
  MyData.DLoad(0,FileName,TypeString(dummy),& ISizeX,& ISizeY,& Slices);

  if (BgCalib && BgCor)
	cerr << "Subtracting Background from data \n";

  if (MySensCMean && SensCor)
	cerr << "Calibrated with sensitivity map \n";

  CalculateMeans(SensMean,BgMean);

  for (int z=0;z <ISizeZ; z++)
  for (int y=0;y <ISizeY; y++)
  for (int x=0;x <ISizeX; x++)
    { dat=MyData.Value(x,y,z+MoreSlices);
    if (BgCalib && BgCor)
      {
	dat -= MyData.Value(x,y,0);  // In this case the first allways the background
      }

    if (MySensCMean && SensCor)
      {
	if (BgCalib)  // bg data is in the stack
	  if (BgCor)  // user wants to correct it
	    sens=MyData.Value(x,y,1)- BgMean; // this was wrong : -MyData.Value(x,y,0);
	  else
	    sens=MyData.Value(x,y,1);
	else sens=MyData.Value(x,y,0);

      if (sens <= 0) sens = 1;
      dat = (dat*(SensMean+ConstAdd))/(sens+ConstAdd);   // The const is just a trick to get around problems
      if (ClipLimits)
	if (dat < 0)
	  dat=0;
	else if (dat > 0x0fff) dat = 0x0fff;
      }

    if (! dosubslice)
      OutData.SetValue(x,y,z,dat);
    else
      if ((z%Every) == subslice)  // is a valid slice ?
	 OutData.SetValue(x,y,(z-subslice)/Every,dat);
    }
}

void ReadParams(char *  HeaderFile)
{
  char * tmp;
  MoreSlices=0;
  getVal("MyROI.TopLeft().x",HeaderFile,& ROI_TLX);
  getVal("MyROI.TopLeft().y",HeaderFile,& ROI_TLY);
  getVal("MyROI.BottomRight().x",HeaderFile,& ROI_BRX );
  getVal("MyROI.BottomRight().y",HeaderFile,& ROI_BRY );

  tmp=strstr(HeaderFile,"MySensCMean");
  if (tmp) 
    {
      getVal("MySensCMean",HeaderFile,& MySensCMean);
      if (MySensCMean) MoreSlices++;
    }
  else
    SensCor=false;  // make sure no sensitivity correction will be applied

  getVal("BgCalib",HeaderFile,& BgCalib);
  if (BgCalib) MoreSlices++;

  getVal("NumFrames",HeaderFile,& NumFrames);
  // getVal("CamNr",HeaderFile,& CamNr);
  getVal("BlankN",HeaderFile,& BlankN);  // Every "BlankN" a wide field is displayed
  // if (BlankN >= 1) MoreSlices++;  // unfortunately the widefield image named "Last" is also saved after bg and sensmap
  getVal("m_Every",HeaderFile,& Every);  // Every "Every" the motor is moved
  getVal("m_PatXSize",HeaderFile,& PatXSize);
  getVal("m_PatYSize",HeaderFile,& PatYSize);
  getVal("m_HadamardRandScan",HeaderFile,& IsHadamard);
  IsPseudoConf=0;  // In case it is not found in the file
  getVal("m_OtherRandScan",HeaderFile,& IsPseudoConf,true);  // ignore, if not found

  ISizeX=ROI_BRX-ROI_TLX;
  ISizeY=ROI_BRY-ROI_TLY;
  if (NumFrames == 0) NumFrames = 1;  // data was taken not in Stack mode
  ISizeZ=NumFrames;

  fprintf(stdout,"ISizeX = %d\nISizeY = %d\nISizeZ = %d\nAdditional Slices=%d\nMeanSensMap=%d\n",ISizeX ,ISizeY,ISizeZ,MoreSlices,MySensCMean);
  if (IsHadamard && (IsPseudoConf == 0))
  fprintf(stdout,"Hadamard scan,  ");
  else if (IsPseudoConf)
    fprintf(stdout,"Pseudo confocal scan,  ");
  else
    fprintf(stdout,"Dot scan scan,  ");
  fprintf(stdout,"PatXSize=%d,  PatYSize=%d\n",PatXSize,PatYSize);
}

int main(int argc, char * argv[])
{
  bool kflag=false;
  double Kappa=0;
  char * selfname = argv[0];

  FILE * textfile;
  char HFile[10000];
  char * HeaderFile=&HFile[0];
 
  bool TwoCameras=false;

  char OutName[1024],FileName[1024],PlaneMeanFile[1024];
  OutName[0]=0;FileName[0]=0,PlaneMeanFile[0]=0;

  char ** parg= & argv[1];

  while (* parg)
  {
    if (readArg("-k", parg)) {kflag=true;continue;}
    if (readArg("-add0", & ConstAdd0, parg)) {continue;}
    if (readArg("-add1", & ConstAdd1, parg)) {continue;}
    if (readArg("-subslice", & subslice, parg)) {dosubslice=true;continue;}
    if (readArg("-nl", parg)) {ClipLimits=false;continue;}
    if (readArg("-bg", parg)) {BgCor=true;continue;}
    if (readArg("-blc", parg)) {bleechcor=true;continue;}
    if (readArg("-sc", parg)) {SensCor=true;continue;}
    if (readArg("-fr", parg)) {FRAP=true;continue;}  // if selected the plane data is only extracted for the FRAP ROI
    if (readArg("-o",& OutName,parg)) {continue;}
    if (readArg("-pmf",& PlaneMeanFile,parg)) {continue;}  // calculates the mean and std.dev of each slice
    if (readArg("-kappa", & Kappa, parg)) {continue;}  // a value for kappa was stated
    if (readArg("-i",& FileName,parg)) {continue;}

    usage(selfname);
  }

  if (FileName[0]==0 || OutName[0]==0) usage(selfname);


  if ((textfile = fopen(FileName,"r")) == NULL) {
    cerr << "couldnt open datafile " << FileName << "\n";
    exit(-1);
    }
  int i;
  for (i=0;i< 10000;i++)
 	if (fread(& HeaderFile[i],1,1,textfile) != (unsigned) 1 )
	  break;

  HeaderFile[i]=0;
  fclose(textfile);

  /* printf("Lng : %d Header : \n %s",i,HeaderFile); */

  ReadParams(HeaderFile);
  ROI1xs= ROI_TLX,ROI1ys=ROI_TLY,ROI1xe=ROI_BRX,ROI1ye=ROI_BRY;

  if (! dosubslice)
    MyOutData.Resize(ISizeX,ISizeY,ISizeZ);
  else
    {
      if (subslice >= Every)
	{
	  cerr << "Error in subslicing ! Only have " << Every << " slices \n";
	  subslice = Every -1;
	}

      MyOutData.Resize(ISizeX,ISizeY,(ISizeZ+Every-1)/Every);  // round up
    }

char * tmpc=strstr(HeaderFile,"CAMERA:");
  if (tmpc!=0)
    {
      tmpc=strstr(tmpc+1,"CAMERA:");
      if (tmpc!=0)
	{
	  Elements=2;
	  HeaderFile=tmpc+1; // now examine only second Cam
	  TwoCameras=true;
	}
    }


  ofstream to(OutName);

  MyOutData.DHeader(kflag,to,Elements);
  ConstAdd=ConstAdd0;
  ProcessStack(FileName,".0",MyOutData);
  FileName[strlen(FileName)-2]=0; // clear ending

  if (FRAP)
    GetFRAPData(PlaneMeanFile,HeaderFile);

  MyOutData.Write(& to);

  if (TwoCameras)
    {
      ConstAdd=ConstAdd1;
      ReadParams(HeaderFile);

      if (! dosubslice)
	MyOutData2.Resize(ISizeX,ISizeY,ISizeZ);
      else
	MyOutData2.Resize(ISizeX,ISizeY,(ISizeZ+Every-1)/Every);  // round up

      ProcessStack(FileName,".1",MyOutData2);
      MyOutData2.Write(& to);

      if (Kappa == 0)  // compute automatically
	Kappa=1.0/(PatXSize*PatYSize);
      if (IsHadamard)
	Kappa=1.0/2.0;
      if (IsPseudoConf)
	Kappa=1.0/6.0;  // only for the case of 31 scan length

      if (PlaneMeanFile[0])
	ProcessPlanes(PlaneMeanFile,Kappa);

    }
  to.close();

}
