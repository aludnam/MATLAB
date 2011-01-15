#ifndef MATRIX_H		// -*- C++ -*- 
#define MATRIX_H

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
#include <stdio.h>
#include "vec.h"

using namespace std;

/** Here the structure of a matrix is given. Its concept is multidimensional. 
    Everything which can be done one-dimensional works this way. TheVector is the Matri. */
class Matrix : public Vector 
{
protected:
  const int DimensionX;
  const int DimensionY;
  VecValType * TmpArea;
public:
  /// yields zero Vector of dimension ItsDim
  Matrix(int ItsDimX,int ItsDimY);
  /// copy constructor
  Matrix(const Matrix &); 
  ~Matrix();
  /// copy without new
  void copy(const Matrix *);
  VecValType Value(int X, int Y) const;
  void ChangeVal(int X, int Y, VecValType val); 
  void Transp(void);
  void FlipDims(int dim1, int dim2);
  void Generate2DTurn(CalcFloatType alpha);
  void Generate3DTurn(CalcFloatType alpha,int axis);
  void Generate3DTurn(Vector * unitvec);
  void Generate3DTurn(Vector * a,Vector * b,Vector * c,Vector * d);
  void MulVec(Vector * AVectorToChange) const;
  void CheckOrthoMatrix(void) const;
  void Mul(Matrix * m1, Matrix * m2);
  void Mul(Matrix * m);
  void show(void);
};

inline Matrix::Matrix(int ItsDimX,int ItsDimY) : Vector(ItsDimX*ItsDimY), DimensionX(ItsDimX), DimensionY(ItsDimY)
{  
  // TheVector = new VecValType[Dimension]; // is allready done in Vector(old.Dimension)
  TmpArea = new VecValType[DimensionY];   // allocate matrix and a TmpArea for Multiplications
}

inline void Matrix::copy(const Matrix * old)  // copy without new !
{
  int i;
  if (old->DimensionX != DimensionX) cout << "FATAL ERROR MATRIXCOPY. NONMATCHING X DIMENSIONS \n",exit(-1);
  if (old->DimensionY != DimensionY) cout << "FATAL ERROR MATRIXCOPY. NONMATCHING Y DIMENSIONS \n",exit(-1);

  for (i=0;i<Dimension;i++)
    TheVector[i] = old->TheVector[i];
}


inline VecValType  Matrix::Value(int X, int Y)
  const  // does not change Matrix
{
  return TheVector[X+Y*DimensionX];
}

inline void  Matrix::ChangeVal(int X, int Y, VecValType val)
{
  TheVector[X+Y*DimensionX] = val;
}


inline Matrix::Matrix(const Matrix & old) : Vector(old.Dimension), DimensionX(old.DimensionX), DimensionY(old.DimensionY) // copy constructor
{
  int i;

  //  TheVector = new VecValType[Dimension]; // is allready done in Vector(old.Dimension)
  TmpArea = new VecValType[DimensionY];   // allocate matrix and a TmpArea for Multiplications
  for (i=0;i<Dimension;i++)
    TheVector[i] = old.TheVector[i];
}

inline Matrix::~Matrix() 
{
  // if (TheVector) delete TheVector;
  if (TmpArea) delete TmpArea;
};

inline void Matrix::MulVec(Vector * AVectorToChange)   // uses the TmpArea after the LastRow
  const  // does not change Matrix !
{
  CalcFloatType sum;
  int x,y;

  if (DimensionX != DimensionY) cout << "Matrix has to be square for Mukltiplication";
  for (y=0;y<DimensionY;y++)  
    {
      sum =0.0;
      for (x=0;x<DimensionX;x++)
      {
	sum += Value(x,y)*AVectorToChange->comp(x);
      }
    TmpArea[y] = sum;
    }

  for (y=0;y<DimensionY;y++)
    AVectorToChange->ChangeComp(y,TmpArea[y]);
}

inline void Matrix::Generate2DTurn(CalcFloatType alpha)
{
  if (DimensionX != 2 || DimensionY != 2) cout << "Matrix has to be a 2-dim square for Generate2DTurn\n",exit(-1);
  ChangeVal(0,0, cos(alpha));
  ChangeVal(1,0, +sin(alpha));
  ChangeVal(0,1, -sin(alpha));
  ChangeVal(1,1, cos(alpha));
}

/// dim1-axis is now dim2 axis, dim2-axis is now dim1-axis
inline void Matrix::Transp(void) 
{
  int x,y;
  VecValType tmp;
  if (DimensionX != DimensionY) cerr << "Error : Matriox hat to be square for transpose ! Has dim : "<< DimensionX << "x" << DimensionY << "\n",exit(-1);

  for (y=0;y< DimensionY;y++)
    for (x=y+1;x<DimensionX;x++)
      {
	tmp = Value(x,y);
        ChangeVal(x,y,Value(y,x));
        ChangeVal(y,x,tmp);
      }
}

/// dim1-axis is now dim2 axis, dim2-axis is now dim1-axis
inline void Matrix::FlipDims(int dim1, int dim2) 
{
  int i;
  VecValType tmp;
  for (i=0;i < DimensionX;i++)
    {
      tmp = Value(dim1,i);
      ChangeVal(dim1,i,Value(dim2,i));
      ChangeVal(dim2,i,tmp);
    }
  for (i=0;i < DimensionX;i++)
    {
      tmp = Value(i,dim1);
      ChangeVal(i,dim1,Value(i,dim2));
      ChangeVal(i,dim2,tmp);
    }
}

/// axis : 0 == X, 1 == Y, 2 == Z
inline void Matrix::Generate3DTurn(CalcFloatType alpha, int axis) 
{
  if (DimensionX != 3 || DimensionY != 3) cout << "Matrix has to be a 3-dim square for Generate3DTurn\n",exit(-1);
  SetZero();
  if (axis == 2) alpha = -alpha;  // to remain in a right handed system
  ChangeVal(0,0, cos(alpha));
  ChangeVal(1,0, +sin(alpha));
  ChangeVal(0,1, -sin(alpha));
  ChangeVal(1,1, cos(alpha));
  ChangeVal(2,2, 1);              // this is a z-axis turnmatrix now flip
  FlipDims(axis,2);
}

/// args must not be this
inline void Matrix::Mul(Matrix * m1, Matrix * m2)  // left times right
  {
    if ( this == m1 || this == m2) cerr << "Fatal error : matrix must not be its argument in matrix multiplication\n",exit(-1);

    int xt,yt,x; // target matrix
    VecValType sum;
    for (xt=0;xt < DimensionX;xt++)
      for (yt=0;yt < DimensionY;yt++)
	{
	  sum=0;
	  for (x=0;x < DimensionX;x++)
	      {
		sum += m1->Value(x,yt) * m2->Value(xt,x);
	      }
	  ChangeVal(xt,yt,sum);
	}
} 

inline void Matrix::Mul(Matrix * m)    // this times right matrix
{
  Matrix tmp(* this);
  Mul(& tmp,m);
}

/// Turn-matrix from z-axis-vector to vector vec
inline void Matrix::Generate3DTurn(Vector * unitvec) 
{
  VecValType angle_y= acos(unitvec->comp(2));                    // angle about the y-axis
  VecValType angle_z= atan2(unitvec->comp(1),unitvec->comp(0));  // angle about the z-axis
  Matrix tmp(3,3);
  Generate3DTurn(angle_z,2);
  tmp.Generate3DTurn(angle_y,1);

  // cout << " Rotate about z-axis (" << angle_z << ") : ";
  // this->show(); 
  // cout << " Rotate about y-axis (" << angle_y << ") : ";
  // tmp.show();
  Mul(& tmp);         // turns unita1 into z-direction
}

inline void Matrix::CheckOrthoMatrix(void) const // check, if matrix is unitary
{
  cout << "CheckOrthoMatrix  : Scalar Products are : ";
  cout << Value(0,0)*Value(1,0)+Value(0,1)*Value(1,1)+Value(0,2)*Value(1,2) << "\n";
  cout << Value(0,0)*Value(2,0)+Value(0,1)*Value(2,1)+Value(0,2)*Value(2,2) << "\n";
  cout << Value(0,0)*Value(0,1)+Value(1,0)*Value(1,1)+Value(2,0)*Value(2,1) << "\n";
  cout << Value(0,0)*Value(0,2)+Value(1,0)*Value(1,2)+Value(2,0)*Value(2,2) << "\n";
}


/// this matrix truns unita into unitb
inline void Matrix::Generate3DTurn(Vector * unita1, Vector * unitb1, Vector * unita2, Vector * unitb2)
  {
  if (DimensionX != 3 || DimensionY != 3) cout << "Matrix has to be a 3-dim square for Generate3DTurn\n",exit(-1);

  Vector unita2s(* unita2), unitb2s(* unitb2);
 
  Matrix rota_z(3,3),rotb_z(3,3),rotabout_z(3,3); // rotate to a to z, b to z , about z

  rota_z.Generate3DTurn(unita1);
  rotb_z.Generate3DTurn(unitb1);


  rota_z.Transp();
  rota_z.MulVec(& unita2s);
  rota_z.Transp();     // to make it again turn from z-axis to a

  rotb_z.Transp();
  rotb_z.MulVec(& unitb2s);

VecValType angle2= atan2(unita2s.comp(1),unita2s.comp(0)) - atan2(unitb2s.comp(1),unitb2s.comp(0)) ;  // angle about the z-axis
  rotabout_z.Generate3DTurn(angle2,2);

  // [rotate z to a] * [rotate about z] * [ rotate b to z]
  copy(& rota_z);        
  Mul(& rotabout_z);
  Mul(& rotb_z);  
  }

/// just show that you are there, and your contents
inline void Matrix::show(void)
{
  int x,y;
  printf("Matrix : %d x %d\n",DimensionX,DimensionY);
  for (y=0;y<DimensionY;y++)
   {
     for (x=0;x<DimensionX;x++)
      printf("%5f  ",Value(x,y));
    printf("\n");
   }
}

#endif













