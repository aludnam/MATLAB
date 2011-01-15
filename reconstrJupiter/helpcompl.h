#ifndef helpcompl_h  		// -*- C++ -*- 
#define helpcompl_h

#include <complex>

/*template <class basetype>
basetype CastToReal(basetype x) {return x;}*/

float CastToReal(float x) {return x;}

double CastToReal(double x) {return x;}

int CastToReal(int x) {return x;}

template <class basetype>
basetype CastToReal(complex<basetype> x) {return real(x);}

template <class basetype>
class CastToComplex {
public:
  typedef complex<basetype> Type;
};

template <>
class CastToComplex<complex<float> > {
public:
  typedef complex<float> Type;
};

template <>
class CastToComplex<complex<double> > {
  typedef complex<double> Type;
};

template < class For, class Pred, class T>
void replace_if_else(For first, For last, Pred p, const T& val1,const T& val2)
{
	while (first != last) {
	if (p(*first)) *first = val1;
	else *first = val2;
	++first;
	}
}

/*template <class basetype>
basetype MulI(basetype x) {return x;}*/

template <class basetype>
complex<basetype> MulI(complex<basetype> val) {return val*complex<basetype> (0,1.0);}

/*template <class basetype>
basetype Conj(basetype x) {return x;}*/ 

template <class basetype>
complex<basetype> Conj(complex<basetype> x) {return conj(x);}

template <class T>
bool operator <(T val1,T val2)
{
  return CastToReal(val1) < CastToReal(val2);
}

template <class T>
bool operator >(T val1,T val2)
{
  return CastToReal(val1) > CastToReal(val2);
}
/*
template <>
bool operator <(class complex<float> val1,class complex<float> val2)
{
  return CastToReal(val1) < CastToReal(val2);
}

template <class basetype>
bool operator >(class complex<basetype> val1,class complex<basetype> val2)
{
  return CastToReal(val1) > CastToReal(val2);
}*/

template <class sourceType>
bool LowerThan(sourceType val1,sourceType val2)
{
  return CastToReal(val1) <  CastToReal(val2);
}

template <class sourceType>
bool GreaterThan(sourceType val1,sourceType val2)
{
  return CastToReal(val1) >  CastToReal(val2);
}

template <class basetype>
basetype Phasor(basetype val, CalcFloatType wiener)
{
  return val/basetype((abs(val)+wiener));
}

// Useful Predicate for STL:

template <class T> struct between : public unary_function<T,bool> {
  T Min,Max;
  between(T & VMin, T & VMax) : Min(VMin), Max(VMax) {}
  bool operator() (const T& x) const {return (x > Min) && (x <= Max);}
};

#endif
