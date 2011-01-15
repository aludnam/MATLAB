// This class offer an easy way to use simple statistic functions as mean and stddev

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

template <class NType>
class Stat {
  private :

  NType sum;
  NType SumSqr;
  int Num;
  NType max;
  NType min;
  int maxNum;
  int minNum;

  public:

  Stat(void) : sum(0.0), SumSqr(0.0), Num(0), max(0.0), min(0.0), maxNum(0), minNum(0) {};

  void Init(NType Val) // Register can also be used to init with first value
    {
      sum=Val;
      SumSqr=Val*Val;
      Num=1;
      max=Val;
      min=Val;
      maxNum=1;
      minNum=1;
    }

  void Reset(void)
    {
      sum=0.0;
      SumSqr=0.0;
      Num=0;
      min=0.0;
      max=0.0;
      maxNum=0;
      minNum=0;
    }

  void Register(NType Val)
    {
      if (Num == 0) 
	{
	  Init(Val);
	  return;
	}

      sum += Val;
      SumSqr += Val*Val;
      Num ++;
 
      if (Val > max) 
	{
	  max=Val;
	  maxNum=Num;
	}
      if (Val < min) 
	{
	  min=Val;
	  minNum=Num;
	}
    }
  
  int Events(void)
    {
      return Num;
    }

  NType Sum(void)
    {
      return sum;
    }

  NType Mean(void)
    {
      return sum/Num;
    }

  NType StdDev(void)  // for first data point NaN can result due to round of errors !!! sqrt(-eps)
    {
      return sqrt((SumSqr-sum*sum/Num)/Num);
    }

  NType Max(void)
    {
      return max;
    }

  NType Min(void)
    {
      return min;
    }

  int MaxNum(void)
    {
      return maxNum;
    }

  int MinNum(void)
    {
      return minNum;
    }

};
