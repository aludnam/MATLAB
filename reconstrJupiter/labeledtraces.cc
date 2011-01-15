// labeledtraces.cc : Normalizes float image by different methods

/*   This file is part of a software package written by 
     Rainer Heintzmann
     MPI bpc
     Am Fassberg 11
     37077 Goettingen, Germany
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
#include "rawarray.h"
#include "parseargs.h"

typedef float ArrayBType;
typedef TArray3d<ArrayBType>    TImgArray; 
typedef TArray3d<int>    TLabelArray;

static TImgArray TimeSeries,Out,Pixels,Stat,Histo;
static TLabelArray LabelImg;


void usage(char * filename)
{
  cerr <<  "usage: " << filename << " -i inputfile -s timeseries -o outputfile [-sX SizeX] [-sY SizeY] [-sZ SizeZ] \n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 

static int SizeX=64,SizeY=64,SizeZ=64;
static bool kflag=false,bincolor=false,allElem=true;
 int Elem=0,Elements=1;
int minvol=1,maxvol=-1;
 double endmax=1e37;

string  IFileName,SeriesFileName,StatFileName,HistoFileName,FilmFileName,OverlayFileName,  // Input for film generation
            OFileName;

 double gapthresh=0.5,  // in relative units
   mindecay=0,decay;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
    if (readArg("-sX", & SizeX, parg)) continue; 
    if (readArg("-sY", & SizeY, parg)) continue;
    if (readArg("-sZ", & SizeZ, parg)) continue;
    if (readArg("-i",  IFileName, parg)) continue;
    if (readArg("-minv", & minvol, parg)) continue;
    if (readArg("-maxv", & maxvol, parg)) continue;
    if (readArg("-endmax", & endmax, parg)) continue;
    if (readArg("-s",  SeriesFileName, parg)) continue;
    if (readArg("-stat",  StatFileName, parg)) continue;
    if (readArg("-histo",  HistoFileName, parg)) continue;
    if (readArg("-film",  FilmFileName, parg)) continue;
    if (readArg("-overlay",  OverlayFileName, parg)) continue;  // if give, this is taken as the overlay
    if (readArg("-gapthresh", & gapthresh, parg)) continue;  // defines which pixel-to-pixel decay counts as gap
    if (readArg("-mindecay", & mindecay, parg)) continue;  // only traces count that have at least this overall decay
    if (readArg("-o",  OFileName, parg)) continue;
    if (readArg("-k", parg)) {kflag=true;continue;}
    if (readArg("-bincolor", parg)) {bincolor=true;continue;}
    if (readArg("-allelem", parg)) {allElem=false;continue;}
    usage(argv[0]);
  }

  if (OFileName=="") 
    usage(argv[0]);

   ofstream to(OFileName.c_str());
   ofstream tos(StatFileName.c_str());
   ofstream toh(HistoFileName.c_str());

   double val;
   int label;

 for (Elem=0;Elem<Elements;Elem++)
   { 
       cout << "Loading Element " << Elem << "\n";
       int SizeZ2=0;
       Elements= LabelImg.DLoad(kflag,IFileName.c_str(),"Integer",& SizeX, &SizeY, & SizeZ2,Elem);
       Elements= TimeSeries.DLoad(kflag,SeriesFileName.c_str(),"Float",& SizeX, &SizeY, & SizeZ,Elem);

       int Max = LabelImg.Maximum();
       int VesicleCount = 0;

       if (Max > 1000)
         {
  	 cerr<<"Error : Too many labels ! " << Max << "\n";
                LabelImg.DSave(true,"/tmp/LabelImg.kdf");
  	 exit(-1);
         }

     if (Elem==0)
{
     Out.Resize(TimeSeries.GetSize(2),Max+1,1);
     Histo.Resize(TimeSeries.GetSize(2),6,1);  // histogram of gaps positions, histogram of max to gap, cumulative histogram, histogram of gap positions for multi-decays, dist between decays, cumulative first position for multigaps
     Stat.Resize(8,Max+1,1);  // stores : invalid-flag, maximum position, maximum value, max gap position, max gap value, num gaps, sum delta-decay, classification
     Pixels.Resize(TimeSeries.GetSize(2),Max+1,1);
}
     Out.Clear();
     Pixels.Clear();
     Histo.Clear();
     Stat.Clear();

       for (int y=0;y < LabelImg.GetSize(1);y++)
	 for (int x=0;x < LabelImg.GetSize(0);x++)
	   {
	     label = LabelImg.Value(x,y,0);
	     if (label > 0)
	       {
		 // cout << "found label # " << label << " size " << SizeZ << "\n";
	       for (int t=0;t < TimeSeries.GetSize(2);t++)
		 {
		   val = TimeSeries.Value(x,y,t);
		   * Out.Pointer(t,label,0) += val;
		   (* Pixels.Pointer(t,label,0)) ++;
		 }
	       }
	   }

    if (Elem == 0 || allElem)
      {
               for (int i=1;i <= Max;i++)
	 if (* Out.Pointer(TimeSeries.GetSize(2)-1,i,0) / Pixels.Value(TimeSeries.GetSize(2)-1,i,0) <= endmax)
	   for (int t=0;t < TimeSeries.GetSize(2);t++)
	     {
	       * Out.Pointer(t,i,0) /= Pixels.Value(t,i,0);
	       if (Pixels.Value(t,i,0) >= minvol &&
		   ((maxvol < 0) || (Pixels.Value(t,i,0) <= maxvol)))
		 {
		   Stat.SetValue(0,i,0,1.0);  // Mark this trace as valid
		   if (t == 0)
		     VesicleCount ++;
		 }
	       else
		   Stat.SetValue(0,i,0,-1.0);  // Mark this trace as invalid due to volume criterium
	     }


       int valid=0;
       int singlegaps=0;
       int havegaps=0;
       for (int i=1;i <= Max;i++)
	 {
	   double gap;
	   int gapp=0,lastgap=-1,
	       firstgap=-1,firstgapdist=-1;    // distance between first two gaps
	   double gapv = Out.Value(0,i,0)-Out.Value(1,i,0);
	   int maxp=0;
	   double maxv = Out.Value(0,i,0);
	   int maxpbg=0;    // will look for a maximum only before the first gap
	   double maxvbg = Out.Value(0,i,0);

	   for (int t=0;t < TimeSeries.GetSize(2);t++) // looks for the maximal value and the biggest gap
		{
	         val = Out.Value(t,i,0);
	         if (val > maxv)
		   maxv=val, maxp=t;
    		}

	   int numgaps=0;
	   if (gapv / maxv > gapthresh)   // check for an initial gap in the first 2 frames
	       {
		 if (firstgap < 0)
		    firstgap = 0;
		 numgaps ++;
		 lastgap=0;
	       }

	   for (int t=0;t < TimeSeries.GetSize(2);t++) // looks for the maximal value and the biggest gap
	     {
	       val = Out.Value(t,i,0);
	       if ((val > maxvbg) && (numgaps == 0))
		  maxvbg=val, maxpbg=t;
	       if (t < TimeSeries.GetSize(2)-1)
		 {
		   gap = val - Out.Value(t+1,i,0);
		   if (gap / maxv > gapthresh) // Change this  line to gap / val !?
		     if (t > lastgap+1)  // two neighboring gaps do not count as two but as one !
		       {
			 if (firstgapdist < 0 && numgaps > 0)  // this is allready the second gap!
				firstgapdist = t-firstgap;
			 if (firstgap < 0)
				firstgap = t;
			 numgaps ++;
			 lastgap=t;
		         // Out.SetValue(t,i,0,val-10.0);
		       }
		 }

	       if (gap > gapv)
		 gapv=gap, gapp=t;
	     }

	   // Position is validity flag
	   if (Out.Value(0,i,0) != 0)
	     decay = Out.Value(TimeSeries.GetSize(2)-1,i,0) / maxv;  // decay is counted as relation between maximum in trace and last value
	   else
	     decay = 1e37;
	   if ((Stat.Value(0,i,0) != -1) && (decay > mindecay))
	     Stat.SetValue(0,i,0,0);  // This will be marked as no-exocytosis vesicle

	   Stat.SetValue(1,i,0,(double) maxp);
	   Stat.SetValue(2,i,0, maxv);
	   Stat.SetValue(3,i,0,(double) gapp);
	   Stat.SetValue(4,i,0, gapv);
	   Stat.SetValue(5,i,0,(double) numgaps);
	   Stat.SetValue(6,i,0,decay);   // THis is equivalent to the sum of derivatives

	   if (Stat.Value(0,i,0) == -1) // Is invalid 
	   	Stat.SetValue(7,i,0,1);   // INVALID 
	   else if (Stat.Value(0,i,0) == 1) // Is exocytosis vesicle
	     {
	       valid ++;
	       if (numgaps >= 1) 
		 {
		   havegaps++;
		   if (numgaps == 1) 
		     {
		       singlegaps ++;
		       (* Histo.Pointer(gapp,0,0) ) += 1.0;  // Insert the maximal gap into the histogram
		       // if (gapp - maxp < 0)
			// cout << "Warning at trace " << i << " maximum "<<maxp<<" is after gap " <<gapp<<"\n";
		       // else
			 (* Histo.Pointer(gapp-maxpbg,1,0) ) += 1.0;  // Insert the gap to max distance into the histogram
		     }
		   if (numgaps > 1) 
			{
		       (* Histo.Pointer(firstgap,3,0) ) += 1.0;  // Insert the maximal gap into the histogram
			if (firstgapdist > 0)
		            (* Histo.Pointer(firstgapdist,4,0) ) += 1.0;  // Insert the maximal gap distance into the histogram
			}
		 }
		if (Stat.Value(5,i,0) == 1)  // This is a one-gap vesicle
	   		Stat.SetValue(7,i,0,2);   // One Gap
	        else  if (Stat.Value(5,i,0) > 1)  // This is a multi-gap vesicle
	   		Stat.SetValue(7,i,0,3);   // Multi Gap
		else   // no gap
	   		Stat.SetValue(7,i,0,4);   // Slow Decay
	     }
	     else  // Stat.Value == 0
	       Stat.SetValue(7,i,0,5);   // Not enough Decay but counting vesicle
	   
	 }
	     
       cout << "Statistics : \n Number of Vesicles: " << Max << ", number of counting vesicles : " << VesicleCount<< "\n";
       cout << "Relative number of exocytosis (" << double(valid) << ") vesicles : " << double(valid) / double(VesicleCount) << "\n";
       cout << "One-gap-vesicles (" << singlegaps << ") over valid vesicles (without ot with gaps): " << double(singlegaps) / double(valid) << "\n";
       cout << "Multiple-gap-vesicles over vesicles with gaps (" <<havegaps << ") : " << double(havegaps-singlegaps) / double(havegaps) << "\n";
	}  // if (Elem == 0 || allElem)   

     if (Elem == 0)
       {
	   Out.DHeader(kflag,to,Elements);
	   cout << "Statfile contains: invalid-flag (-1 volume wrong, 0 too little overall decay, 1 valid), maximum position, maximum value, max gap position, max gap value, num gaps, sum delta-decay, classification\n";
	   cout << "Histo-file contains: histogram of gaps positions, histogram of max to gap, cumulative histogram, histogram of first gap positions for multi-decays, dist between decays, cumulative first position for multigaps \n";
	   Stat.DHeader(kflag,tos,Elements);
	   Histo.DHeader(kflag,toh,Elements);
       }


     float sum=0.0,sum2=0.0;
     for (int i=0; i < Histo.GetSize(0);i++)  // generate cumulative histograms
       {
	 sum += Histo.Value(i,0,0);  // Sum up the histogram
	 (* Histo.Pointer(i,2,0) ) = sum;  // Sum up the histogram
	 sum2 += Histo.Value(i,3,0);  // Sum up the histogram
	 (* Histo.Pointer(i,5,0) ) = sum2;  // Sum up the histogram
       }

	 
     Histo.Write(& toh);
     Stat.Write(& tos);
     Out.Write(& to);

     if (FilmFileName!="")
       {
	 cout << "Genrating Film ...\n";
	 cout << "Color coding :\n purple (4) means slowly decaying counting vesicle\n ";
	 cout << " yellow (2) means counting vesicle with exactly one recognized drop (gap) in intensity\n ";
	 cout << " red (3) means counting vesicle with more than one recognized drops (gap) in intensity\n ";
	 cout << " cyan (5) means counting vesicle, but did not decayed enough and no drops were detected\n ";
	 cout << " green (1) means vesicle did not fulfill the criteria for being a counting vesicle (too small, too big)\n ";

	 TImgArray Film;
	 Film.Resize(TimeSeries.GetSize(0),TimeSeries.GetSize(1),TimeSeries.GetSize(2));
	 ofstream tof(FilmFileName.c_str());
	 Film.DHeader(kflag,tof,3);

	 if (OverlayFileName!="")
	   TimeSeries.DLoad(kflag,OverlayFileName.c_str(),"Float",& SizeX, &SizeY, & SizeZ,Elem);

	 if (! TimeSeries.SizesEqual(& Film))
	   {cerr << "Error, Overlay sizes to not equal time series!\n"; exit(-1);}

	int label;
	double value;
	 for (int e=0;e<3;e++)
	   {
	     for (int t=0;t < TimeSeries.GetSize(2);t++)
	       for (int y=0;y < TimeSeries.GetSize(1);y++)
		 for (int x=0;x < TimeSeries.GetSize(0);x++)
		   {
		     label = LabelImg.Value(x,y,0);
		     value = TimeSeries.Value(x,y,t);
		     if (bincolor)
		       value = 1.0;

		     if (label == 0) 
			{
		        if (bincolor)
				Film.SetValue(x,y,t, 0.0);  // dim black and white for outside label
			else
				Film.SetValue(x,y,t, value / 2.0);  // dim black and white for outside label
			}
		      else // inside a labeled region
	               if (Stat.Value(0,label,0) == -1)  // This is excluded due to size being too big or too small
			   {
			   if (e == 0)
				Film.SetValue(x,y,t,0);
			   else if (e == 1)
				Film.SetValue(x,y,t,value);
			   else if (e == 2)
				Film.SetValue(x,y,t,0);
			   }
		      else if (Stat.Value(0,label,0) == 1)  // This is considered exocytosis
			   {

			   if (Stat.Value(5,label,0) == 1)  // This is a one-gap vesicle
			       {
			   	if (e == 0) // yellow for one gap
					Film.SetValue(x,y,t,2.0*value);
			  	else if (e == 1)
					Film.SetValue(x,y,t,2.0*value);
				else if (e == 2)
					Film.SetValue(x,y,t,0);
			       }
			   else  if (Stat.Value(5,label,0) > 1)  // This is a multi-gap vesicle
				{
			   	if (e == 0) // yellow for one gap
					Film.SetValue(x,y,t,2.0*value);
			  	else if (e == 1)
					Film.SetValue(x,y,t,0);
				else if (e == 2)
					Film.SetValue(x,y,t,0);
			        }
			   else 
				{
			   	if (e == 0) // yellow for one gap
					Film.SetValue(x,y,t,value);
			  	else if (e == 1)
					Film.SetValue(x,y,t,0);
				else if (e == 2)
					Film.SetValue(x,y,t,value);
			        }
		       	   }
			else  // is not an exocytosis
				{
			   	if (e == 0) // yellow for one gap
					Film.SetValue(x,y,t,0);
			  	else if (e == 1)
					Film.SetValue(x,y,t,value);
				else if (e == 2)
					Film.SetValue(x,y,t,value);
			        }
		   }
	     Film.Write(& tof);  // write an element
	   }
	 tof.close();
       }
   }

 toh.close();
 tos.close();
 to.close();


}
