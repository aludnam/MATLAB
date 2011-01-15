
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

#include <stdio>
#include "vec.h"
#include "matrix.h"
#include "veccalc.h"

static VecValType XPix=(9.9206/256),
  zscalea=1; // (6.144 / 37) / XPix;

// this programm calculates z-scaling, turn-angles from matching of three given Vectors in two images
const double DEG= 360/ (2*M_PI);
Vector MXB(3),MYB(3),MZB(3),POSB(3),MXF(3),MYF(3),MZF(3),POSF(3);   // Unitvectors

Vector ma1(4,0,3), ma2(0,0,0), ma3(0,1,0);  // points which should match
Vector mb1(3,0,2.1*4), mb2(0,0,0), mb3(0,1,0);
// Vector ma1(66.586, 50.489, 27.163), ma2(238.758, 190.552, 7.293), ma3(181.416, 83.921, 13.551);  // points which should match
// Vector mb1(51.882, 42.726, 9.417),  mb2(196.751, 182.281, 31.77), mb3(146.879, 75.37, 24.435);

int main (void)
{

  CalcFromPos(ma1,ma2,ma3,mb1,mb2,mb3,POSB,MXB,MYB,MZB,POSF,MXF,MYF,MZF,zscalea);

  Vector tmp(3),tmp2(3);
  cout << "Vectors apply : \n";
  cout << "Point 1 : \n";
  tmp.copy(& POSB);
  tmp2.copy(& MXB);
  tmp2.Mul(mb1.comp(0));
  tmp.Add(& tmp2);
  tmp2.copy(& MYB);
  tmp2.Mul(mb1.comp(1));
  tmp.Add(& tmp2);
  tmp2.copy(& MZB);
  tmp2.Mul(mb1.comp(2));
  tmp.Add(& tmp2);
  tmp.show();

  tmp.copy(& POSB);
  tmp2.copy(& MXB);
  tmp2.Mul(mb2.comp(0));
  tmp.Add(& tmp2);
  tmp2.copy(& MYB);
  tmp2.Mul(mb2.comp(1));
  tmp.Add(& tmp2);
  tmp2.copy(& MZB);
  tmp2.Mul(mb2.comp(2));
  tmp.Add(& tmp2);
  tmp.show();

  tmp.copy(& POSB);
  tmp2.copy(& MXB);
  tmp2.Mul(mb3.comp(0));
  tmp.Add(& tmp2);
  tmp2.copy(& MYB);
  tmp2.Mul(mb3.comp(1));
  tmp.Add(& tmp2);
  tmp2.copy(& MZB);
  tmp2.Mul(mb3.comp(2));
  tmp.Add(& tmp2);
  tmp.show();


  return EXIT_SUCCESS;
}
