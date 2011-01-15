/* The poidev routine gives a poisson-distributed random value.
   aruments : float xm, long * idum
   xm : mean expectation
   idum : pointer to negative integer (always the same pointer)
*/

#ifndef POIDEF_H	
#define POIDEF_H

/*   This file is part of a software package written by 
     Rainer Heintzmann
     Institute of Applied Optics and Information Processing
     Albert Ueberle Strasse 3-5
     69120 Heidelberg
     Tel.: ++49 (0) 6221  549264
     No garantee, whatsoever is given about functionallaty and  savety.
     No warranty is taken for any kind of damage it may cause.
     No support for it is provided !

     THIS IS NOT FREE SOFTWARE. All rights are reserved, and the software is only
     given to a limited number of persons for evaluation purpose !
     Please do not modify and/or redistribute this file or any other files, libraries
     object files or executables of this software package !
 
    The code in this module was directly obtained from the Numerical Recipies in C book.
    The copyright stated there does therefore also apply to this code.          
*/


float ran1(long *);

float gammln(float);

float poidev(float, long *);

#endif
