#ifndef STLFIXES
#define STLFIXES

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
*/

#include <math.h>
#include <algorithm>
#include <functional>
// #include <ext/algorithm>    // for copy_n in the gnu distribution (not part of the standard)
// #include <ext/functional>

// #define M_PI            3.14159265358979323846  /* pi */
using namespace std;

// #ifdef STLFIX
template <class InputIterator, class Size, class OutputIterator>
OutputIterator copy_n(InputIterator first, Size count,
                    OutputIterator result) {
  for (;count;count--)
    *result++ = *first++;
  return result;
}

template <class InputIterator, class OutputIterator>
struct project1st : public binary_function<InputIterator,OutputIterator,InputIterator> {
	project1st() {}
  InputIterator operator () (InputIterator & first, OutputIterator & second) {
    return first;}
};


// #endif

 
template <class T>
struct conjugate : public unary_function<T,T> {
  T operator()(const T & x) const { return conj(x); }
};
     
template <class T>
struct mulconj : public binary_function<T,T,T> {
  T operator()(const T & x,const T & y) const { return x*conj(y); }
};

template <class T>
struct raisetopower : public binary_function<T,T,T> {
  T operator()(const T & x,const T & y) const { if (CastToReal(x) > 0) return pow(x,y); else return 0;}
};
          

#endif
