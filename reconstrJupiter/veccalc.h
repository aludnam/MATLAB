#ifndef VecCalc_h		// -*- C++ -*- 
#define VecCalc_h

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

template <class TFloat>
TFloat sqr(TFloat val)
{
  return val*val;
}


Vector scalea(3), scaleb(3);

void CalcScales(const Vector & da1, const Vector & db1, const Vector & da2, const Vector & db2)
{
  double ca1,ca2,cb1,cb2,sca,scb;

  ca1 = (sqr(da1.comp(0))+sqr(da1.comp(1))-sqr(db1.comp(0))-sqr(db1.comp(1)))/sqr(db1.comp(2));
  ca2 = sqr(da1.comp(2))/sqr(db1.comp(2));

  cb1 = (sqr(da2.comp(0))+sqr(da2.comp(1))-sqr(db2.comp(0))-sqr(db2.comp(1)))/sqr(db2.comp(2));
  cb2 = sqr(da2.comp(2))/sqr(db2.comp(2));

  sca = sqrt(-(ca1-cb1) / (ca2-cb2) );
  scb = sqrt(ca1+ca2*sqr(sca));
  cout << "Scale A : " << sca << "\n";
  cout << "Scale B : " << scb << "\n";
  cout << " c : " << ca1 << "  " << cb1 << "  " << ca2 << "  " << cb2 << "\n";
}

double CalcOtherScale (const Vector & da1, const Vector & db1, double zscale)
{
 return  sqrt((sqr(da1.comp(0))+sqr(da1.comp(1))+sqr(zscale*da1.comp(2))-sqr(db1.comp(0))-sqr(db1.comp(1)))/sqr(db1.comp(2)));
}

void CalcShiftAndVecs(const Matrix & Turn,const Vector & ma1, const Vector & mb1,const Vector & ScaleA,const Vector & ScaleB, 
		      Vector & POS, Vector & MX, Vector & MY, Vector & MZ)
{
  Vector UX(1,0,0),UY(0,1,0),UZ(0,0,1);   // Unitvectors
  Vector tmp(mb1),shift(ma1);  // ma1 and mb1 are the positions that are to allign exactly !

  tmp.Mul(& ScaleB);
  Turn.MulVec(& tmp);
  tmp.Div(& ScaleA);
  shift.Sub(& tmp);   // calculate shift between ma1 and mb1
  POS.copy(& shift);  // write to POS

  UX.Mul(& ScaleB);
  Turn.MulVec(& UX);
  UX.Div(& ScaleA);
  MX.copy(& UX);

  UY.Mul(& ScaleB);
  Turn.MulVec(& UY);
  UY.Div(& ScaleA);
  MY.copy(& UY);

  UZ.Mul(& ScaleB);
  Turn.MulVec(& UZ);
  UZ.Div(& ScaleA);
  MZ.copy(& UZ);
 
}

void CalcFromPos(const Vector & ma1, const Vector & ma2, const Vector & ma3, const Vector & mb1, const Vector & mb2, const Vector & mb3,
                 Vector & POSB, Vector & MXB, Vector & MYB, Vector & MZB,Vector & POSF, Vector & MXF, Vector & MYF, Vector & MZF,
		 double ZXAspect, double ZXAspect2=0,Matrix * TurnF=0)
{
Vector da1(3),da2(3), db1(3), db2(3);
double scb1=CalcOtherScale(da1,db1,ZXAspect),scb2=CalcOtherScale(da2,db2,ZXAspect);
Vector SCA(1.0,1.0,ZXAspect),SCB(1.0,1.0,ZXAspect2),SCC(1.0,1.0,(scb1+scb2)/2.0);

  da1.copy(& ma2); da1.Sub(& ma1);         // calculate difference vectors
  da2.copy(& ma3); da2.Sub(& ma1);

  db1.copy(& mb2); db1.Sub(& mb1);
  db2.copy(& mb3); db2.Sub(& mb1);

  scalea.copy(&SCA);
  
  // CalcScales(da1,db1,da2,db2);  // still makes problems

  cout << "Z has to be scaled by : " << ZXAspect << "\n";
  if (ZXAspect2 != 0)
    scaleb.copy(&SCB);
  else
    {
      scaleb.copy(&SCC);
      cout << "Otherscale : " << scb1 << " and " << scb2 << "\n";
    }
  da1.Mul(& scalea);
  da2.Mul(& scalea);

  db1.Mul(& scaleb);
  db2.Mul(& scaleb);

  cout << "Scale : " << scaleb.comp(2) << " error is : "<< db2.Norm()-da2.Norm() << " and " << db1.Norm()-da1.Norm() << "\n";

  // angle = da1.Angle(& da2);
  // cout << "Angle between vectors is : " << angle* DEG << "\n";

Vector unita1(da1),unitb1(db1), unita2(da2),unitb2(db2);
  Matrix Turn1(3,3);

  unita1.MkUnit(), unitb1.MkUnit(),unita2.MkUnit(), unitb2.MkUnit();

  Turn1.Generate3DTurn(&unita1,&unitb1,&unita2,&unitb2);  // turns unita1 into unitb1 and unita2 into unitb2
  cout << "Matrix which transforms : \n";
  Turn1.show();

Matrix Turn2(Turn1);    // for the backprojection
  Turn2.Transp();
  cout << "backward matrix : " ;
  Turn2.show();
  cout << "scalea : " ;
  scalea.show();
  cout << " scaleb " ;
  scaleb.show();
  cout  << "\n";

  CalcShiftAndVecs(Turn2,mb1,ma1,scaleb,scalea,POSB,MXB,MYB,MZB);
  // #ifdef DEBUG_PRINT
  cout << "Backward Projecting Vectors :\n";
  cout << "shift vector for point (0,0,0) is : ";
  POSB.show();
  cout << "move vector X is : ";
  MXB.show();
  cout << "move vector Y is : ";
  MYB.show();
  cout << "move vector Z is : ";
  MZB.show();
  // #endif

  CalcShiftAndVecs(Turn1,ma1,mb1,scalea,scaleb,POSF,MXF,MYF,MZF);
  // #ifdef DEBUG_PRINT
  cout << "Forward Projecting Vectors :\n";
  cout << "shift vector for point (0,0,0) is : ";
  POSF.show();
  cout << "move vector X is : ";
  MXF.show();
  cout << "move vector Y is : ";
  MYF.show();
  cout << "move vector Z is : ";
  MZF.show();
  // #endif

  if (TurnF)
    TurnF->copy(& Turn1);
}

void ShiftCorrect(Vector * vec, double ZXAspect,double HDAspect)
{
  vec->ChangeComp(0,vec->comp(0) + (vec->comp(2)*ZXAspect)*HDAspect);    // "+" to correct for - shift
#ifdef DEBUG_PRINT
vec->show();
#endif
}

#endif
