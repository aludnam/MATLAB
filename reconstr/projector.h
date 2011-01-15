#ifndef Projector_h		// -*- C++ -*- 
#define Projector_h

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


#include <vector>
#include "point.h"
      
/** This class knows about how to intersect itself with another pixel.
 The Projector and its derived classes are the only ones, who should be informed of any
 geometrical problems. So they care about moving, turning, and so on. Derived classes also know the exact calculations
 for beam intersections. They know how to compute any weights.*/
class Projector : public Pixel	
{
private:
  void UpdatePos2(void);

protected:
 /// indicates the edge of the clipping rectangle relative to center
 
  // Vector * dirs, * dimbegins;
  vector<Vector> dirs, dimbegins;

public:
  Projector(const int Dim) : Pixel(Dim)
    {
      int i;
      if (Dim >= 2) Direction.ChangeComp(1,1.0);
      
      Vector x(Dim);
      for (i=0;i<Dim;i++)
	{
	  dirs.push_back(x);
	  dimbegins.push_back(x);
	}

      //dirs = new Vector[Dim] (Dim);
      //dimbegins = new Vector[Dim] (Dim);
      
      for (i=0;i<Dim;i++)
	{
	  dirs[i].SetUnit(i);
	  dimbegins[i].SetZero();
	}
    };

  /// Copy Constructor has to allocate space for Vectors, but not for Array
  Projector(const Projector & CopyFromPrj) : Pixel(&CopyFromPrj.Position ,& CopyFromPrj.Direction)
    {
      int Dim=CopyFromPrj.dirs[0].ItsDimension();
      
      Vector x(Dim);
      for (int i=0;i<Dim;i++)
	{
	  dirs.push_back(x);
	  dimbegins.push_back(x);
	}
      // dirs = new Vector[Dim]  (Dim);
      // dimbegins = new Vector[Dim] (Dim);
      
      // int i;
      // for (i=0;i<Dim;i++) {dirs[i].copy (& CopyFromPrj.dirs[i]); dimbegins[i].copy (& CopyFromPrj.dimbegins[i]);}
      copy(& CopyFromPrj);
    }
  
 ~Projector()
    {
      // delete dirs;
      // delete dimbegins;
    };
  

  /// just copy, without constructing.
  void copy(const Projector * APrj)
    {
      int i;
      Position.copy(& APrj->Position);
      Direction.copy(& APrj->Direction);
      for (i=0;i<dimbegins[0].ItsDimension();i++)
	{
	  dirs[i].copy (& APrj -> dirs[i]);
	  dimbegins[i].copy (& APrj -> dimbegins[i]);
	}
    } ;
  /// copy without allocating space.
  void copy(const Vector* pItsPos,
	    const Vector* pItsDir)
    {
      Position.copy(pItsPos);
      Direction.copy(pItsDir);
    };      
  
  void TakeNewChangeVec(int dim, const Vector * pVec)
    {
      dirs[dim].copy(pVec);
    };

  void InitMove()
    {
      int i;
      for (i=0;i<dimbegins[0].ItsDimension();i++)
	dimbegins[i].copy(& Position);
    };

  /// returns 1, if a new initialization of the array is needed (change of direction or width ...), otherwise 0
  int Move(int indimension)
    {
      if (indimension <= 0)          // move in normal (X) direction
	Position.Add(& dirs[0]);
      else
	{
	  int i;           // dimbegins[0] is not used at all !
	  dimbegins[indimension].Add(& dirs[indimension]);
	  Position.copy(& dimbegins[indimension]);
	  for (i=indimension;i>1;i--)
	dimbegins[i-1].copy(& dimbegins[indimension]);
	}
      return 0;         // no array initialization needed in this case
    };

  /** computes intersection of itself with one pixel.
  pPosition is Position in Image (image-coordinates)
  computes intersection of itself with one pixel. 
  Returns value between 0.0 and 1.0*Value
  This Function is to be written by the user
  Return value of 1.0 means full contribution, 0 means none. */
  // #pragma argsused
CalcFloatType ProjectYourselfFwd(void * pArr)	   
{
  return 1.0;
  pArr=0;  // to prevent warning
}

 /// computes intersection of itself with one pixel for the backward Projection.See  ProjectYourselfFwd for description
CalcFloatType ProjectYourselfBwd(void * pArr)	
{
  return 1.0;
  pArr=0;  // to prevent warning
}

  Vector * GiveDirection(void)  {return (& Direction);};
  Vector * GivePosition(void)   {return (& Position);};
  
  /// Turns direction using the matrix
  void Turn(Matrix * TurnMatrix)
    {
      TurnMatrix->MulVec(& Direction);
    }; 
  void show(void)
    {
      cout << "Projector \n";
      cout << "  Position : ";
      Position.show();
      cout << "  Direction : ";
      Direction.show();
    };

};

#endif

