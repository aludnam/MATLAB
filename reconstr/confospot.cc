// This program tries to select spots in an confocal image and determines their position
// #define PRINTALL

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


static float PsfScale=1.0;  // if bigger than 1.0  -> Psf has higher resolution

#include "confospot.h"       // has all declarations needed

void usage(char * filename)
{
  cerr << "usage: " << filename << " [-k] [-mX X ...] [-i inputfile] [-p psffile] [-o outputfile] [-pX X ...] [-b background] [-n NewLimit] [-t ThroughLimit] [-c SmallestCorr] [-m MaxSpots]\n" << flush;
  cerr << "-mX X -mY Y -mZ Z :  X,Y,Z == sizes of measured projection\n";
  cerr << "-pX X Y Z : sizes of psf-file\n";
  cerr << "-k : khoros data format flag\n";
  cerr << "-vo : filename for binary vector output\n";
  exit(-1);
}

int main(int argc, char *argv[])
{  
  int i,Elements=1;
  bool kflag=false;             // KDF ?
  bool restart=false;  // if true, every frame will be restarted

  string INPUTFileName, OFileName,PSFFileName,
            SimFileName,RemainFileName,BinFileName, SpotFileName;

  double Overrelax=0.3; // 0.5; // 4.5;  // 1.0;

  int MEASUREDSizeX=256; // 128
  int MEASUREDSizeY=256; // 128
  int MEASUREDSizeZ=32;  // 25

char ** parg= & argv[1];
printCall(argc,argv);

argc = 0; // just to use it

 while (* parg)
  {
   if (readArg("-k",parg)) {kflag=true;continue;}
   if (readArg("-v",parg)) {verbose=true;continue;}
   if (readArg("-restart",parg)) {restart=true;continue;}
   if (readArg("-i", INPUTFileName, parg)) continue;
   if (readArg("-spots", SpotFileName, parg)) continue; // starting values for spot positions
   if (readArg("-o", OFileName, parg)) continue;
   if (readArg("-p", PSFFileName, parg)) continue;
   if (readArg("-vo",BinFileName, parg)) continue;
   if (readArg("-sim", SimFileName, parg)) continue;
   if (readArg("-r",  RemainFileName, parg)) continue;
   if (readArg("-region", & region, parg)) continue;

   if (readArg("-mX",& MEASUREDSizeX, parg)) continue;
   if (readArg("-mY",& MEASUREDSizeY,parg)) continue;
   if (readArg("-mZ",& MEASUREDSizeZ,parg)) continue;

   if (readArg("-pX",& PSFSizeX, parg)) continue;
   if (readArg("-pY",& PSFSizeY, parg)) continue; 
   if (readArg("-pZ",& PSFSizeZ, parg)) continue; 
   if (readArg("-pf",& PsfScale, parg)) continue; 

   if (readArg("-b",& BACKGR, parg)) continue;            // constant background value
   if (readArg("-n",& NewSpotLimit, parg)) continue;      // if divided image is below that, no spots are added anymore
   if (readArg("-t",& ThrowAwayLimit, parg)) continue;  // Spot, whos intensity is below this will be thrown away
   if (readArg("-freeze",& FreezeLimit, parg)) continue;  // Spot, whos intensity is below this will be fozen (no more movement before reactivation
   if (readArg("-identdist",& IdentDistance, parg)) continue;  // New spots are identified as old spots if these are frozen and at a distance below this value
   if (readArg("-c",& SmallestIntenCorr, parg)) continue; // A new spot is added, if corrections in all others are below this
   if (readArg("-m",& SpotNumber, parg)) continue;        // Maximal Number of Spots 
   if (readArg("-a",& IterationsAfter, parg)) continue;   // Iterations that are allways performed after all spots are localized
   if (readArg("-over",& Overrelax, parg)) continue;   // Maximal Number of Spots 

   usage(argv[0]);
  }

  if (SimFileName != "")
    PrjFile=SimFileName;

  if (RemainFileName != "")
    RemainFile= RemainFileName;
  
  ReconPrj.SetBackground(BACKGR);

  Elements=MeasuredPrj.Array.DLoad(kflag,INPUTFileName.c_str(),"Float",& MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ);
  
  if (MeasuredPrj.Array.Minimum() < 0)
    {
      cerr << "\n\nWARNING : Measured Image has values lower than zero ! Clipping them to zero !\n\n" << flush;
      MeasuredPrj.Array.ClipAt(0);
     }
  // cout << "Checkpoint 7\n" << flush;

  ValuePSFArray.DLoad(kflag,PSFFileName.c_str(),"Float",& PSFSizeX,& PSFSizeY,& PSFSizeZ);    // Load PSF from disk
  if (ValuePSFArray.Minimum() < 0)
    {
      cerr << "\n\nWARNING : PSF has values lower than zero ! Clipping them to zero !\n\n";
      ValuePSFArray.ClipAt(0);
    }

  ReconPrj.Array.Resize(MEASUREDSizeX,MEASUREDSizeY,MEASUREDSizeZ);
  ReconSpotsImg.Array.Resize(SpotDimension,SpotNumber);
  CorrectionImg.Array.Resize(SpotDimension,SpotNumber);
  OldCorr.Resize(SpotDimension,SpotNumber);
  XPSFArray.Resize(PSFSizeX,PSFSizeY,PSFSizeZ);
  YPSFArray.Resize(PSFSizeX,PSFSizeY,PSFSizeZ);
  ZPSFArray.Resize(PSFSizeX,PSFSizeY,PSFSizeZ);

  ReconPrj.SetPsfFac(PsfScale);     // enlargement of psf
  MeasuredPrj.SetPsfFac(PsfScale);  // enlargement of psf

  ReconPrj.SetClip((PSFSizeY/4.0-0.3)/PsfScale,(PSFSizeX-0.3)/2.0/PsfScale,(PSFSizeZ-0.2)/2.0/PsfScale); // each +1.0 ?
  cout << "PSF Clipping : " << (PSFSizeY/4.0-0.3)/PsfScale << "  " << (PSFSizeX-0.3)/2.0/PsfScale << "  " << (PSFSizeZ-0.2)/2.0/PsfScale << "\n";
  cout << "PSF Sizes : " << PSFSizeX << "  " << PSFSizeY << "  " << PSFSizeZ << "\n";
  cout << "PSF Scale : " << PsfScale << "\n";

  for (i=0;i<= SpotNumber;i++) 
    ActiveSpots[i] = 0;
   for (i=0;i< SpotNumber;i++) 
     ClearSpot(i);

  ActualSpotNumber=0;

  if (SpotFileName != "")
    {
      TImgArray StartSpots;
      int SpotSizeX,SpotSizeY,SpotSizeZ;

      StartSpots.DLoad(kflag,SpotFileName.c_str(),"Float",& SpotSizeX,& SpotSizeY,& SpotSizeZ);

      for (i=0;i< SpotSizeY;i++) 
	{
	  if (i > SpotNumber)
	    cerr << "Not enough spots! The number of spots in the input spotfile is bigger!\n",exit(-1);

	  if (StartSpots.Value(0,i)==0) ActiveSpots[i] = 0;
	  else 
	    {
	      ActiveSpots[i] = 1;
	      ReconSpotsImg.Array.SetValue(0,i,StartSpots.Value(0,i));  // copy spots to SpotImg.
	      ReconSpotsImg.Array.SetValue(1,i,StartSpots.Value(1,i));
	      ReconSpotsImg.Array.SetValue(2,i,StartSpots.Value(2,i));
	      ReconSpotsImg.Array.SetValue(3,i,StartSpots.Value(3,i));
	      ActualSpotNumber=i+1;
	    }
	}
      cout << "Added the " << SpotSizeY << " spots from file !\n";
    }

  // Now all arguments are evaluated and Arrays are resized

  OldCorr.Set(0.0);

  cout << "Normalize PSF  \n";

  ValuePSFArray.NormalizeToMax(1.0);             // normalizes integral, but maximum has to be normalized to one before

  XPSFArray.DeriveFrom(&ValuePSFArray,0);
  YPSFArray.DeriveFrom(&ValuePSFArray,1);
  ZPSFArray.DeriveFrom(&ValuePSFArray,2);

  ValuePSFArray.DSave(kflag,"/tmp/mypsf.raw");
  XPSFArray.DSave(kflag,"/tmp/myxpsf.raw");
  YPSFArray.DSave(kflag,"/tmp/myypsf.raw");
  ZPSFArray.DSave(kflag,"/tmp/myzpsf.raw");

  ReconSpotsImg.TakeNewComperator(& SubComperator);
  PSFIntegral=ValuePSFArray.Integral();
  ComputePSFCenter(& ValuePSFArray);

  cout << "Starting Reconstruction .. \n";

 // Now start to reconstruction process

  ReconPrj.Array.Clear();

  cout << "# Simulation :  Backgr : " << BACKGR << "\n"; // << " Clip : " << ClipVal << "\n";

  cout << "# ( IError[i], XError[i], YError[i] ), XCOMError, YCOMError \n";

  ofstream * oi = 0;
  if (BinFileName != "")
    {
            oi = new ofstream(BinFileName.c_str());
          if (! oi)
            cerr << "Error ! Couldnt open " << BinFileName << " for writing !! \n",exit(-1);
     }
	      
  ofstream * to = 0;
  if (OFileName != "")
  {
          to = new ofstream(OFileName.c_str());
         if (! to)
                cerr << "Error ! Couldnt open " << OFileName << " for writing !! \n",exit(-1);
    }
    
  ofstream * sim=0;
  if (Elements > 1 && SimFileName !="")
    {
      sim = new ofstream(SimFileName.c_str());
      if (! sim)
	cerr << "Error ! Couldnt open " << SimFileName << " for writing !! \n",exit(-1);
    }

  if (Elements > 1)
    PrjFile="";

  // MeasuredPrj.Array.DSave(false,"/tmp/start.raw");
  // ReconSpotsArray.Load("/tmp/2D/tstpos.raw");  // otherwise use old positions as startingpoints
  
  for (int el=0;el < Elements; el++)   // There are other elements present. They will be processed now
    {
      if (el > 0)
	{
           if (restart)
       {  
                   for (i=0;i< ActualSpotNumber;i++)
                      ClearSpot(i);
                    ActualSpotNumber=0;
       }

       
	  cout << "\n\nProcessing element " << el <<"\n"<<flush;
	  Elements=MeasuredPrj.Array.DLoad(kflag,INPUTFileName.c_str(),"Float",& MEASUREDSizeX,& MEASUREDSizeY,& MEASUREDSizeZ,el);
  
	  if (MeasuredPrj.Array.Minimum() < 0)
	    {
	      cerr << "\n\nWARNING : Measured Image has values lower than zero ! Clipping them to zero !\n\n";
	      MeasuredPrj.Array.ClipAt(0);
	    } 
	}
      
      if (ActualSpotNumber == 0)
	AddNewSpots();
        else
        {
	ActivateAll();  // now do Intensity-Only iterations with all spots
                for (i=0;i< IterationsBtw;i++)
            {
        	    if (IterationStep(Overrelax,0) < SmallestIntenCorr) break;  // Iterate without movement
            }
  	cout << "Intensity-Iterations finished. Iterations needed : " << i << "\n";

              FreezeDimSpots();  // freezes the spots at still fixed positions
       }

        
      do {
	for (i=0;i< IterationsBtw;i++)
	  {
#ifdef PRINTALL     
	    cout << "Iteration : " << i << " Overrelax = " << Overrelax << "\n";
#endif
	    if ((IterationStep(Overrelax) < SmallestIntenCorr)  && (i > 4)) break;
	  }
	
	cout << "Surrounding ready. Iterations needed : " << i << "\n";
	
	for (i=0;i< IterationsBtw2;i++)
	  IterationStep(Overrelax);
	
	ActivateAll();  // now do some iterations with all spots
	for (i=0;i< IterationsBtw;i++)
	  {
	    if (IterationStep(Overrelax) < SmallestIntenCorr) break;
	  }
	
	cout << "All spots ready. Iterations needed : " << i << "\n";
      } while (AddNewSpots());
      
      
      // FreezeDimSpots();
      RemoveSmallSpots();

      for (i=0;i<IterationsAfter;i++)
	{
#ifdef PRINTALL     
	  cout << "Iteration : " << i << "\n";
#endif
	  IterationStep(Overrelax);
	}


      // cout << "Spot Values : \n";
  /*    cout << "Element "<< el << ": Spot Nr., Intensity, PosX, PosY, PosZ\n";
      for (i=0;i<ActualSpotNumber;i++)
	{
	 // if (Elements > 1)
	 //  cout << "Element: " << el << " \n";
	  cout << i << "  ";
	  cout << ReconSpotsImg.Array.Value(0,i) << "  ";
	  cout << ReconSpotsImg.Array.Value(1,i) << "  ";
	  cout << ReconSpotsImg.Array.Value(2,i) << "  ";
	  cout << ReconSpotsImg.Array.Value(3,i) << "\n";
	}
      cout << "\n";
      */
      
      if (OFileName != "")
	{
	  if (el == 0)
	    (* to) << " Element Nr, Spot Nr., Intensity, PosX, PosY, PosZ\n";
	  //if (Elements > 1)
	 //   (* to)<< "Element: " << el << " \n";
	  for (i=0;i<ActualSpotNumber;i++)
	    {
	      (* to)<< el << "  ";
	      (* to)<< i << "  ";
	      (* to)<< ReconSpotsImg.Array.Value(0,i) << "  ";
	      (* to)<< ReconSpotsImg.Array.Value(1,i) << "  ";
	      (* to)<< ReconSpotsImg.Array.Value(2,i) << "  ";
	      (* to)<< ReconSpotsImg.Array.Value(3,i) << "\n";
	    }
	    (* to)<< "\n";
	}
      if (BinFileName != "")
	{
	  if (el == 0)
	    {
	      if (Elements == 1)
		SaveArray.Resize(SpotDimension,ActualSpotNumber);  // Copy valid values into array
	      else
		SaveArray.Resize(SpotDimension,SpotNumber);  // Copy valid values into array

	      SaveArray.DHeader(kflag,* oi,Elements);
	    }
	  SaveArray.Copy(& ReconSpotsImg.Array);
	  
	  SaveArray.Write(oi);
	}

      if (sim && Elements > 1)
	{
	  if (el == 0)
	    ReconPrj.Array.DHeader(kflag,* sim,Elements);

	  ReconPrj.ProjectImage(&ReconSpotsImg); // only adds into image
	  ReconPrj.Array.Write(sim);
	}
    }

  if (sim)
    sim->close();
  if (oi)
    oi->close();
  if (to)
    to->close();
    
return 0;
}
