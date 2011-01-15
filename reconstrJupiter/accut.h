#ifndef accut_h  		// -*- C++ -*- 
#define accut_h

#include <complex>


template <class basetype>
struct Accu {
  typedef double type;  // for all simple types
};

template <>
struct Accu <complex<float> > {  // applies only if the type is a complex
  typedef complex<double> type;
};
template <>
struct Accu <complex<double> > {  // applies only if the type is a complex
  typedef complex<double> type;
};
template <>
struct Accu <complex<int> > {  // applies only if the type is a complex
  typedef complex<double> type;
};

#endif
