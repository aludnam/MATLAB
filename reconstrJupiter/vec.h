#ifndef VEC_H		// -*- C++ -*- 
#define VEC_H

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
#include <stdlib.h>
#include <iostream>
#include <assert.h>

using namespace std;

typedef double VecValType;
typedef double CalcFloatType;

#define CHECK(vec,fktname) assert (vec->Dimension == Dimension) 
#define CHECKGT(vec,fktname) assert (vec->Dimension <= Dimension)
#define CHECKLT(vec,fktname) assert (vec->Dimension >= Dimension)
#define CHECK2(fktname)  assert (Dimension == 2)
#define CHECK3(fktname)  assert (Dimension == 3)
#define CHECKGED(dim,fktname) assert (Dimension >= (dim))

/*
class VecValue
{
protected:
  VecValType Value;
public:
  VecValue() {Value = 0;}
  VecValue(VecValType ItsValue) {Value = ItsValue;}
  void show(void) {cout << Value << " ";}
}; */


/** Here the structure of a vector is given. Its concept is multidimensional.  Easy constructors exist for 2 and 3D Vectors.*/
class Vector
{
protected:
  const int Dimension;
  VecValType* TheVector;
public:
  Vector(const int ItsDim);                           // yields zero Vector of dimension ItsDim
  Vector(VecValType v1,VecValType v2 );               // constructor for 2 parameters yields 2D vector 
  Vector(VecValType v1,VecValType v2, VecValType v3); // constructor for 3 parameters yields 3D vector 
  Vector(const Vector &);                             // copy constructor
  Vector & operator=(const Vector &);                          // Assignment copies the values in Vector
  ~Vector();

  void copy(const Vector *); // copy without new
  void copy(VecValType v1,VecValType v2 );
  void copy(VecValType v1,VecValType v2, VecValType v3);

  void show(void);
  void SetZero(void);
  void Set(void);
  void SetUnit(int j);
  CalcFloatType SqrDistance(Vector* ToCalc);
  CalcFloatType Distance(Vector* ToCalc);
  void Add(const Vector * pVec);
  void Sub(const Vector * pVec);
  void Mul(CalcFloatType Val);
  void Mul(const Vector * pVec);
  void MulAdd(Vector * pVec,CalcFloatType Val);        // Multiplies *PVec with Val and adds it to self 
  void Div(CalcFloatType Val);
  void Div(const Vector * pVec);
  CalcFloatType SProd(Vector * pVec);
  CalcFloatType AbsSqr(void);
  CalcFloatType Norm(void);
  CalcFloatType GaussVal(Vector * Center,CalcFloatType width, CalcFloatType hight);
  void rot90Z(void);
  void rot90X(void);
  void flipXY(void);
  void flipYZ(void);
  void flipX(void);
  void flipY(void);
  void flipZ(void);
  void flipAt(Vector * FlipVec);  // Flip at plane perpendicular to FlipVec
  void neg(void);
  void MkUnit(void);
  CalcFloatType CalcArea(void);          // Calculates the area covered by a Rectangular shaped region with edges Vec,-Vec and other two
  CalcFloatType MulLinDist(Vector * Vec2);				      
  CalcFloatType Angle(Vector * v2);
  CalcFloatType SqrDistToLine(Vector * AtPosition,Vector * ItsDirection);   // calculates the square of the distance to the line described by AtPosition and ItsDirection
                                                                                    // ItsDistance MUST be of unit length !!

  VecValType comp(int adim) const;
  void ChangeComp(int adim,VecValType NewVal);
  int ItsDimension(void) const;
  int asciiread(istream *);
};


inline void Vector::show(void)
{
  int i;
  cout << "Vector [";
  // cout.precision(8);
  for (i=0; i < Dimension; i++)
    cout << TheVector[i] << ", ";
  cout << "]\n";
}


inline void Vector::SetZero(void) // sets a Vector, which allready exists to zero
{
  int i;
  for (i=0; i < Dimension; i++)
    TheVector[i]=0.0;  
}

inline void Vector::SetUnit(int j) // sets a Vector, which allready exists to zero
{
  int i;
  for (i=0; i < Dimension; i++)
    if (i==j)
      TheVector[i]=1.0;
    else
      TheVector[i]=0.0;
}


inline Vector::Vector(const int ItsDim) : Dimension(ItsDim)
{  
//  cout << "Vector::Vector(" << ItsDim << ")\n";
  TheVector = new VecValType[Dimension];  
  SetZero();
}

/* inline Vector::Vector()
{  
//  cout << "Vector::Vector()\n";
  Dimension=0;
  TheVector=NULL;
}*/

inline void Vector::copy(const Vector * old)  // copy without new !
{
  int i;
  CHECKLT(old,"VECCOPY");
  for (i=0;i<Dimension;i++)    // old can be bigger, but only Dimensions are copied
    TheVector[i] = old->TheVector[i];
}


inline void Vector::copy(VecValType v1, VecValType v2)  // copy without new !
{
    CHECK2("2DVECCOPY");
    TheVector[0] = v1;
    TheVector[1] = v2;
}

inline void Vector::copy(VecValType v1, VecValType v2,VecValType v3)  // copy without new !
{
    CHECK3("3DVECCOPY");
    TheVector[0] = v1;
    TheVector[1] = v2;
    TheVector[2] = v3;
}


inline Vector::Vector(const Vector & old) : Dimension(old.Dimension) // copy constructor
{
  int i;
//  cout << "Vector::Vector(const Vector&)\n";

  TheVector = new VecValType[Dimension];
  for (i=0;i<Dimension;i++)
    TheVector[i] = old.TheVector[i];
//  show();
}

inline Vector &  Vector::operator=(const Vector & other)                           // Assignment copies the values in Vector
{
  int i;
   for (i=0;i<Dimension;i++)
    TheVector[i] = other.TheVector[i];
   return *this;
} 

inline Vector::~Vector()
{
//  cout << "~Vector()\n";
//  show();
  if (TheVector) delete TheVector;
}


inline Vector::Vector(VecValType v1,VecValType v2) : Dimension(2)
{
  TheVector = new VecValType[Dimension];
  TheVector[0]=v1;
  TheVector[1]=v2;
}

inline Vector::Vector(VecValType v1,VecValType v2, VecValType v3) : Dimension(3)
{
  TheVector = new VecValType[Dimension];
  TheVector[0]=v1;
  TheVector[1]=v2;
  TheVector[2]=v3;
}

inline CalcFloatType Vector::SqrDistance(Vector* ToCalc) // distance using 2-norm
{
  int i;
  CalcFloatType res=0,tmp=0;
  CHECK(ToCalc,"SqrDistance");

  for (i=0; i < Dimension; i++)
    {
      tmp = ((CalcFloatType) TheVector[i]) - ((CalcFloatType) ToCalc->TheVector[i]); 
      res += (tmp*tmp);
    }

  return res;
}

inline CalcFloatType Vector::Distance(Vector* ToCalc) // distance using 2-norm
{
  return sqrt(SqrDistance(ToCalc));
}

inline void  Vector::Add(const Vector * pVec)
{
  int i;
  CHECK(pVec,"Add");

  for (i=0; i < Dimension; i++)
    TheVector[i] += (pVec->TheVector[i]);

}

inline void  Vector::Sub(const Vector * pVec)
{
  int i;
  CHECK(pVec,"Sub");

  for (i=0; i < Dimension; i++)
    TheVector[i] -= pVec->TheVector[i];

}

inline void  Vector::Mul(CalcFloatType Val)
{
  int i;
  
  for (i=0; i < Dimension; i++)
    TheVector[i] *= Val;
}

inline void  Vector::Mul(const Vector * pVec)     // componentwise multiplication
{
  int i;
  CHECK(pVec,"Mul");
 
  for (i=0; i < Dimension; i++)
    TheVector[i] *= (pVec->TheVector[i]);
}

inline void  Vector::MulAdd(Vector * pVec,CalcFloatType Val)
{
  int i;
  CHECK(pVec,"MulAdd");
  for (i=0; i < Dimension; i++)
    TheVector[i] += (pVec->TheVector[i])*Val;
}

inline void  Vector::Div(CalcFloatType Val)
{
  int i;
  
  for (i=0; i < Dimension; i++)
    TheVector[i] /= Val;
}

inline void  Vector::Div(const Vector * pVec)     // componentwise division
{
  int i;
  CHECK(pVec,"Mul");
 
  for (i=0; i < Dimension; i++)
    TheVector[i] /= (pVec->TheVector[i]);
}


inline CalcFloatType Vector::SProd(Vector * pVec)
{
  int i;
  CalcFloatType rtv=0;
  CHECK(pVec,"SProd");

  for (i=0; i < Dimension; i++)
    rtv += TheVector[i] * (pVec->TheVector[i]);

  return rtv;
}

inline  void Vector::neg(void)
{
  int i;
  for (i=0; i < Dimension; i++)
    TheVector[i] = - TheVector[i];
}

inline  void Vector::rot90Z(void)  // rot 90 deg right around 3rd dimension
{
  VecValType tmp= TheVector[1];

  (TheVector[1]) =   - TheVector[0];
  (TheVector[0]) =   tmp;
}

inline  void Vector::rot90X(void)  // rot 90 deg right around 3rd dimension
{
  VecValType tmp= TheVector[1];

  (TheVector[1]) =   - TheVector[2];
  (TheVector[2]) =   tmp;
}

inline  void Vector::flipX(void)  // rot 90 deg around 3rd dimension
{
  (TheVector[0]) =   - TheVector[0];
}

inline  void Vector::flipY(void)  // rot 90 deg around 3rd dimension
{
  (TheVector[1]) =   - TheVector[1];
}

inline  void Vector::flipZ(void)  // rot 90 deg around 3rd dimension
{
  (TheVector[2]) =   - TheVector[2];
}

inline  void Vector::flipXY(void)  // flip X and Y coordinate
{
  VecValType tmp= TheVector[1];

  (TheVector[1]) =   TheVector[0];
  (TheVector[0]) =   tmp;
}

inline  void Vector::flipYZ(void)  // flip X and Y coordinate
{
  VecValType tmp= TheVector[1];

  (TheVector[1]) =   TheVector[2];
  (TheVector[2]) =   tmp;
}

inline  void Vector::flipAt(Vector * FlipVec)  // rot 90 deg around 3rd dimension
{
  Vector Tmp(* FlipVec);
  CHECK(FlipVec,"flipAt");
  Tmp.Mul(SProd(FlipVec));
  Tmp.Mul(-2.0);
  Add(& Tmp);
}


inline CalcFloatType Vector::AbsSqr(void)
{
  CalcFloatType abssqr=0.0;
  int i;
  
  for (i=0; i < Dimension; i++)
    abssqr += TheVector[i]*TheVector[i];
  return abssqr;
}

inline CalcFloatType Vector::Norm(void)
{
  return sqrt(AbsSqr());
}

inline CalcFloatType Vector::Angle(Vector * v2)   // calculates angle between vectors
{
  return acos(SProd(v2) / (Norm() * v2->Norm()));
}


inline CalcFloatType Vector::GaussVal(Vector * Center,CalcFloatType FWHM, CalcFloatType hight) // FWHM !
{
   CalcFloatType width = FWHM/sqrt(log(2.0))/2.0;
   Vector Tmp(* Center);
   CHECK(Center,"GaussVal");
   Tmp.Sub(this);
   return hight*exp(- Tmp.AbsSqr()/(width*width));
} 


inline CalcFloatType Vector::MulLinDist(Vector * Vec2)   // calculates the multilinear distance (1.0-XDist)*(1.0-YDist)*(1.0-ZDist) ...
{
  int i;
  CalcFloatType weight=1.0,dirweight;
  CHECK(Vec2,"MulLinDist");
  for (i=0;i<Dimension;i++)
    {
    if((dirweight = (1.0-fabs(Vec2->TheVector[i]-TheVector[i]))) > 0)
      weight *= dirweight;
    else return 0.0;
    }
  return weight;
}


inline CalcFloatType Vector::SqrDistToLine(Vector * AtPosition,Vector * ItsDirection)   // calculates the square of the distance to the line described by AtPosition and ItsDirection
                                                                                        // ItsDistance MUST be of unit length !!
{
  int i;
  CalcFloatType RelPositionComp,AbsSqr=0.0,ScalarProd=0.0;

  // AbsSqr(RelPosition) - SProd(RelPosition,ItsDirection)
  CHECKGT(AtPosition,"SqrDistToLine");
  CHECKGT(ItsDirection,"SqrDistToLine");
  for (i=0;i<AtPosition->Dimension;i++)
    {
      RelPositionComp = TheVector[i]-AtPosition->TheVector[i];
      AbsSqr += RelPositionComp*RelPositionComp;
      ScalarProd += RelPositionComp*(ItsDirection->TheVector[i]);
    }

  return AbsSqr-ScalarProd*ScalarProd;

}


inline CalcFloatType Vector::CalcArea(void)
{
  int i;
  CalcFloatType mul=2;
  for (i=0;i<Dimension;i++)
    mul *= 2;
  return AbsSqr()*mul;
}

inline void Vector::MkUnit(void)   // scale vector to unit-length
{
  Div(Norm());
}

inline  VecValType Vector::comp(int adim)
const
{
  CHECKGED(adim,"Vector::comp");
  return TheVector[adim];
}

inline  void Vector::ChangeComp(int adim,VecValType NewVal)
{
  CHECKGED(adim,"Vector::ChangeComp");
  TheVector[adim]=NewVal;
}

inline int Vector::ItsDimension(void) const
{
  return Dimension;
}

static Vector NullVec3d(0.0,0.0,0.0);
static Vector NullVec2d(0.0,0.0);

inline int even(int val)
{
  return ((val/2)*2 == val);
}

inline int uneven(int val)
{
  return ! even(val);
}

int Vector::asciiread(istream * astream)  // Read a vector from ascii-file
{
  int i;
  for (i=0;i<Dimension;i++)
    {
      (*astream) >> TheVector[i];
    }
  return 1;
}

#endif
