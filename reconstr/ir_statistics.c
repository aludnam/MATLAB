/* 
 * Ir library of image restoration programs.
 * 
 * Copyright (C) 1998  P. J. Verveer and the Max Planc Society, Germany
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include <math.h>
#include "ir.h"

ir_Error ir_Sum(ir_Image *im, double *sum)
{
  ir_Float *p;
  int ii, s, t;
  ir_Error error = IR_OK;
  
  IRXJ(ir_Check(im));
  IRXJ(ir_ImageType(im, &t));
  
  *sum = 0.0;
  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  IRXJ(ir_ImageSize(im, &s));
  IRXJ(ir_ImageData(im, &p));
  for(ii = 0; ii < s; ii++)
    *sum += p[ii];
    
ir_error:

  return error;
}


ir_Error ir_Mean(ir_Image *im, double *mean)
{
  ir_Float *p;
  int ii, s, t;
  ir_Error error = IR_OK;
  
  IRXJ(ir_Check(im));
  IRXJ(ir_ImageType(im, &t));

  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  *mean = 0.0;
  IRXJ(ir_ImageSize(im, &s));
  IRXJ(ir_ImageData(im, &p));
  for(ii = 0; ii <  s; ii++)
    *mean += p[ii];
  *mean /= s;
  
ir_error:

  return error;
}


ir_Error ir_InProduct(ir_Image *in1, ir_Image *in2, double *ip)
{
  ir_Float *p1, *p2;
  ir_Error error = IR_OK;
  int ii, jj, even, rs, rx, x, s, d, t;
  int *dims = 0;

  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_ImageType(in1, &t));

  switch (t) {
  case IR_TYPE_IMAGE:
    IRXJ(ir_ImageSize(in1, &s));
    IRXJ(ir_ImageData(in1, &p1));
    IRXJ(ir_ImageData(in2, &p2));
    *ip = 0.0;
    for(ii = 0; ii < s; ii++)
      *ip += p1[ii] * p2[ii];
    break;
  case IR_TYPE_SPECTRUM:
    IRXJ(ir_ImageSize(in1, &s));
    IRXJ(ir_ImageData(in1, &p1));
    IRXJ(ir_ImageData(in2, &p2));
    IRXJ(ir_ImageDimensions(in1, &dims, &d));
    IRXJ(ir_ImageRealX(in1, &rx));

    if (d == 0) {
      *ip = *p1 * *p2 + *p1 * *p2;
      return error;
    }

    x = dims[0];
    even = rx == 2 * (x - 1);
    s /= x;

    *ip = 0.0;

    for(ii = 0; ii < s; ii++) {
      *ip += *p1++ * *p2++;
      *ip += *p1++ * *p2++;
      for(jj=2;jj<x;jj++) {
	*ip += 2.0 * *p1++ * *p2++;
	*ip += 2.0 * *p1++ * *p2++;
      }
      if (even) {
	*ip += *p1++ * *p2++;
	*ip += *p1++ * *p2++;
      } else {
	*ip += 2.0 * *p1++ * *p2++;
	*ip += 2.0 * *p1++ * *p2++;
      }
    }
    
    IRXJ(ir_ImageRealSize(in1, &rs));
    *ip /= rs;
    break;
  }
  
ir_error:
  if (dims) free(dims);

  return error;
}

ir_Error ir_TriProduct(ir_Image *in1, ir_Image *in2, ir_Image *in3, 
			    double *ip)
{
  ir_Float *p1, *p2, *p3;
  ir_Error error = IR_OK;
  int ii, s, t;
   
  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_CheckTwo(in1, in3));
  IRXJ(ir_ImageType(in1, &t));

  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  *ip = 0.0; 
  IRXJ(ir_ImageSize(in1, &s));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(in2, &p2));
  IRXJ(ir_ImageData(in3, &p3));
  for(ii = 0; ii < s; ii++)
    *ip += p1[ii] * p2[ii] *p3[ii];
   
ir_error:

  return error;
}


ir_Error ir_Ulnv(ir_Image *in1, ir_Image *in2, double *ulnv)
{
  ir_Float *p1, *p2;
  ir_Error error = IR_OK;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_ImageType(in1, &t));

  *ulnv = 0.0;
  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  IRXJ(ir_ImageSize(in1, &s));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(in2, &p2));
  for(ii = 0; ii < s; ii++)
    if (p2[ii] > 0.0) *ulnv += p1[ii] * log(p2[ii]);
  
ir_error:

  return error;
}
 
