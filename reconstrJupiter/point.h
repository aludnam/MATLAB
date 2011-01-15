#ifndef Point_h		// -*- C++ -*- 
#define Point_h

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


#include "vec.h"
#include "matrix.h"

/** definitions of the Point class, which is the basic class from which
   more complicated object having Space coordinates and directions are derived.
   Also the Pixel class which is derived from point is given here. It contains a direction additionally. */
class Point
{
protected:
  /// Position in Images / projections is allways meant to be of middle voxel
  Vector	Position; 

public:
// Point(){Position=0;};

/// Constructor does nothing
  Point(int Dim) : Position(Dim) {}

/// Constructor does nothing
Point(const Vector * pItsPos) : Position(* pItsPos)
{
  //  cout << "Point::Point(Vector *):\n";
}

/// copy without new
void copy(const Vector * pNewPos)
{
  Position.copy(pNewPos);
}

void show(void)
{
  cout << "Point at ";
  Position.show();
}

void SetCoordinate(const Vector * pItsPos)
{ 
  copy(pItsPos);  // Does NOT allocate new space
}

void SetPosZero()
{
  Position.SetZero();
}

/// point's squared distance to another vector.
CalcFloatType  SqrDistance(Vector * pVec)
{
  //  cout << "Point::Distance\n";
  return Position.SqrDistance(pVec);   // Calculate distance in Vector Class
}

/// point's distance to another vector.
CalcFloatType  Distance(Vector * pVec)
{
  //  cout << "Point::Distance\n";
  return Position.Distance(pVec);   // Calculate distance in Vector Class
}

/// vectorial position
void  Add(Vector * pVec){Position.Add(pVec);}
/// vectorial position
void  Sub(Vector * pVec){Position.Sub(pVec);}
/// of vectorial position
void  SProd(Vector * pVec){Position.SProd(pVec);}                           
/// vectorial position
void  SetPos(const Vector * pItsPos) {Position.copy(pItsPos);}
};


/// has a direction in addition to point
class Pixel : public Point
{
protected:
  /// the addidional direction.
  Vector  Direction;

public:
/// Constructor does nothing
  Pixel(int Dim) : Point(Dim), Direction(Dim) {}

// Pixel(){Direction=0;};
Pixel(const Vector* pItsPos, const Vector* pItsDir) : Point(pItsPos), Direction(* pItsDir)
{
   if (Direction.Norm() != 0.0) Direction.MkUnit();                   //  Vector has to be made unit-vector !
}				

Pixel(Pixel * pPix) : Point(& pPix -> Position), Direction(pPix -> Direction)  // copy cunstructor in Vector class
{
}

/// sets Dir and Pos to zero.
void SetPosDirZero(void)
{
  Position.SetZero();
  Direction.SetZero();
}

void SetDir(const Vector * pItsDir) {Direction.copy(pItsDir);}

Vector * GiveDir(void)
{
  return & Direction;
}

};

#endif
