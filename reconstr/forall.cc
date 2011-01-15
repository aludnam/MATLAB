// forall.cc : Does some command n times inserting a number

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

#include <iostream.h>
#include "parseargs.h"
#include <stdio.h>
#include <stdlib.h>

void usage(char * filename)
{
  cerr <<  "usage: " << filename << " [-s startnumber] [-e endnumber] [-step step] [-i stepp between %d s] commands_with_%d_in_it\n" << flush;
  exit(-1);
}

int main(int argc, char *argv[])
{ 
  int start=0,end=0,step=1,i,j,increment=0;

char ** parg= & argv[1];
argc=0;  // to prevent warning

 while (* parg)
  {
   if (readArg("-s", & start, parg)) continue; 
   if (readArg("-e", & end, parg)) continue; 
   if (readArg("-i", & increment, parg)) continue; 
   if (readArg("-step", & step, parg)) continue; 
   // must be the commands now
    break;
  }

  if (! *parg) 
    usage(argv[0]);

  char command[10000];
  command[0]=0;

  for (i=start;i<=end;i+=step)
    {
      for (j=0;parg[j];j++)
	{      
	  sprintf(command,parg[j],i+1*increment,i+2*increment,i+3*increment,i+4*increment,i+5*increment,i+6*increment,i+7*increment,i+8*increment,i+9*increment,i+10*increment);
	  cout << "executing : " << command << "\n";
	  system(command);
	}
    }
}
