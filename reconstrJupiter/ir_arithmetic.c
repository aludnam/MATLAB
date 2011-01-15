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

ir_Error ir_Add(ir_Image *in1, ir_Image *in2, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *p2, *po;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(in2, &p2));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  IRXJ(ir_ImageType(in1, &t));
  if (t == IR_TYPE_SPECTRUM) s *= 2;
  for(ii = 0; ii < s; ii++)
    po[ii] = p1[ii] + p2[ii];
  
ir_error:
  
  return error;
}

ir_Error ir_Sub(ir_Image *in1, ir_Image *in2, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *p2, *po;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(in2, &p2));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  IRXJ(ir_ImageType(in1, &t));
  if (t == IR_TYPE_SPECTRUM) s *= 2;
  for(ii = 0; ii < s; ii++)
    po[ii] = p1[ii] - p2[ii];

ir_error:
  
  return error;
}

ir_Error ir_Mul(ir_Image *in1, ir_Image *in2, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *p2, *po, r, i;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(in2, &p2));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  
  IRXJ(ir_ImageType(in1, &t));
  switch(t) {
  case IR_TYPE_IMAGE:
    for(ii = 0; ii < s; ii++)
      po[ii] = p1[ii] * p2[ii];
    break;
  case IR_TYPE_SPECTRUM:
    for(ii = 0; ii < s; ii++) {
      r = *p1 * *p2 - *(p1+1) * *(p2+1);
      i = *p1 * *(p2+1) + *(p1+1) * *p2;
      *po++ = r;
      *po++ = i;
      p1++;
      p1++;
      p2++;
      p2++;
    }
    break;
  }

ir_error:
  
  return error;
}

ir_Error ir_Div(ir_Image *in1, ir_Image *in2, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *p2, *po, r, i, rf, df;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(in2, &p2));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
    
  IRXJ(ir_ImageType(in1, &t));
  switch(t) {
  case IR_TYPE_IMAGE:
    for(ii = 0; ii < s; ii++)
      po[ii] = p1[ii] / p2[ii];
      break;
  case IR_TYPE_SPECTRUM:
    for(ii = 0; ii < s; ii++) {
      if (fabs(*p2) >= fabs(*(p2+1))) {
	rf = *(p2+1) / *p2;
	df =  *p2 + rf * *(p2+1);
	r = (*p1 + rf * *(p1+1)) / df;
	i = (*(p1+1) - rf * *p1) / df;
      } else {
	rf = *p2 / *(p2+1);
	df = *p2 * rf + *(p2+1);
	r = (*p1 * rf + *(p1+1)) / df;
	i = (*(p1+1) * rf - *p1) /df;
      }
      *po++ = r;
      *po++ = i;
      p1++;
      p1++;
      p2++;
      p2++;
    }
    break;
  }
   
ir_error:
  
  return error;
}

ir_Error ir_AddConst(ir_Image *in1, ir_Image *out, double c, double i)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  IRXJ(ir_ImageType(in1, &t));
  switch(t) {
  case IR_TYPE_IMAGE:
    for(ii = 0; ii < s; ii++)
      po[ii] = p1[ii] + c;
    break;
  case IR_TYPE_SPECTRUM:
    for(ii = 0; ii < s; ii++) {
      *po++ = *p1++ + c;
      *po++ = *p1++ + i;
    }
    break;
  }

ir_error:
  
  return error;
}

ir_Error ir_SubConst(ir_Image *in1, ir_Image *out, double c, double i)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  IRXJ(ir_ImageType(in1, &t));
  switch(t) {
  case IR_TYPE_IMAGE:
    for(ii = 0; ii < s; ii++)
      po[ii] = p1[ii] - c;
    break;
  case IR_TYPE_SPECTRUM:
    for(ii = 0; ii < s; ii++) {
      *po++ = *p1++ - c;
      *po++ = *p1++ - i;
    }
    break;
  }
    
ir_error:
  
  return error;
}

ir_Error ir_MulConst(ir_Image *in1, ir_Image *out, double re, double im)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po, r, i;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  IRXJ(ir_ImageType(in1, &t));
  switch(t) {
  case IR_TYPE_IMAGE:
    for(ii = 0; ii < s; ii++)
      po[ii] = p1[ii] * re;
    break;
  case IR_TYPE_SPECTRUM:
    for(ii = 0; ii < s; ii++) {
      r = *p1 * re - *(p1+1) * im;
      i = *p1 * im + *(p1+1) * re;
      *po++ = r;
      *po++ = i;
      p1++;
      p1++;
    }
    break;
  }
      
ir_error:
  
  return error;
}

ir_Error ir_DivConst(ir_Image *in1, ir_Image *out, double re, double im)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po, r, i, rf, df;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  IRXJ(ir_ImageType(in1, &t));
  switch(t) {
  case IR_TYPE_IMAGE:
    for(ii = 0; ii < s; ii++)
      po[ii] = p1[ii] / re;
    break;
  case IR_TYPE_SPECTRUM:
    for(ii = 0; ii < s; ii++) {
      if (fabs(re) >= fabs(im)) {
	rf = im / re;
	df =  re + rf * im;
	r = (*p1 + rf * *(p1+1)) / df;
	i = (*(p1+1) - rf * *p1) / df;
      } else {
	rf = re / im;
	df = re * rf + im;
	r = (*p1 * rf + *(p1+1)) / df;
	i = (*(p1+1) * rf - *p1) /df;
      }
      *po++ = r;
      *po++ = i;
      p1++;
      p1++;
    }
    break;
  }

ir_error:
  
  return error;
}

ir_Error ir_WeightedAdd(ir_Image *in1, ir_Image *in2, ir_Image *out, double w)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *p2, *po;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_ImageType(in1, &t));
  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(in2, &p2));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  for(ii = 0; ii < s; ii++)
    po[ii] = p1[ii] + w * p2[ii];
  
ir_error:
  
  return error;
}


ir_Error ir_WeightedMul(ir_Image *in1, ir_Image *in2, ir_Image *out, double w)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *p2, *po;
  int ii, s, t;
  
  IRXJ(ir_CheckTwo(in1, in2));
  IRXJ(ir_ImageType(in1, &t));
  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(in2, &p2));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  for(ii = 0; ii < s; ii++)
    po[ii] = p1[ii] * w * p2[ii];
     
ir_error:
  
  return error;
}


ir_Error ir_Abs(ir_Image *in1, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po;
  int ii, s, t;

  IRXJ(ir_Check(in1));
  IRXJ(ir_ImageType(in1, &t));
  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  for(ii = 0; ii < s; ii++)
    po[ii] = fabs(p1[ii]);
  
ir_error:
  
  return error;
}

ir_Error ir_Clip(ir_Image *in1, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po;
  int ii, s, t;
  
  IRXJ(ir_Check(in1));
  IRXJ(ir_ImageType(in1, &t));
  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  for(ii = 0; ii < s; ii++)
    po[ii] = p1[ii] > 0.0 ? p1[ii] : 0.0;
    
ir_error:
  
  return error;
}

ir_Error ir_Square(ir_Image *in1, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po;
  int ii, s, t, d1, d2;
  int *sdim = 0, *odim = 0;

  IRXJ(ir_Check(in1));
  IRXJ(ir_ImageType(in1, &t));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageSize(in1, &s));  
    
  IRXJ(ir_ImageType(in1, &t));
  switch (t) {
  case IR_TYPE_IMAGE:
    IRXJ(ir_CheckTwo(out, in1));
    IRXJ(ir_ImageData(out, &po));
    
    for(ii = 0; ii < s; ii++)
      po[ii] = p1[ii] * p1[ii];
    break;
  case IR_TYPE_SPECTRUM:
    IRXJ(ir_ImageType(out, &t));
    IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_DIFFERS);
    IRXJ(ir_ImageDimensions(in1, &sdim, &d1));
    IRXJ(ir_ImageDimensions(out, &odim, &d2));
    IRTJ(d1 != d2, IR_SIZES_DIFFER);
    for(ii=0; ii<d1; ii++) {
      IRTJ(sdim[ii] != odim[ii], IR_SIZES_DIFFER);
    }
    IRXJ(ir_ImageData(out, &po));
    
    for(ii = 0; ii < s; ii++) {
      *po = *p1 * *p1;
      p1++;
      *po++ += *p1 * *p1;
      p1++;
    }
    break;
  }
  
 ir_error:

  if (sdim) free(sdim);
  if (odim) free(odim);
  
  return error;
}

ir_Error ir_Sqrt(ir_Image *in1, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po;
  int ii, s, t;
  
  IRXJ(ir_ImageType(in1, &t));
  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  for(ii = 0; ii < s; ii++)
    po[ii] = p1[ii] > 0.0 ? sqrt(p1[ii]) : 0.0;
  
ir_error:
  
  return error;
}

ir_Error ir_Ln(ir_Image *in1, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *po;
  int ii, s, t;
  
  IRXJ(ir_ImageType(in1, &t));
  IRTJ(t != IR_TYPE_IMAGE, IR_TYPE_UNSUPPORTED);
  
  IRXJ(ir_CheckTwo(out, in1));
  IRXJ(ir_ImageData(in1, &p1));
  IRXJ(ir_ImageData(out, &po));
  IRXJ(ir_ImageSize(in1, &s));  
  for(ii = 0; ii < s; ii++)
    po[ii] = p1[ii] > 0.0 ? log((double) p1[ii]) : 0.0;
   
ir_error:
  
  return error;
}


ir_Error ir_MulConjugate(ir_Image *in1, ir_Image *in2, ir_Image *out)
{
  ir_Error error = IR_OK;
  ir_Float *p1, *p2, *po, r, i;
  int ii, s, t;

  IRXJ(ir_CheckTwo(in1, in2));

  IRXJ(ir_ImageType(in1, &t));
  IRXJ(ir_CheckTwo(out, in1));
  if (t != IR_TYPE_SPECTRUM) {
    IRXJ(ir_Mul(in1, in2, out));
  } else {
    IRXJ(ir_ImageData(in1, &p1));
    IRXJ(ir_ImageData(in2, &p2));
    IRXJ(ir_ImageData(out, &po));
    IRXJ(ir_ImageSize(in1, &s));
    for(ii = 0; ii < s; ii++) {
      r = *p1 * *p2 + *(p1+1) * *(p2+1);
      i = *p1 * *(p2+1) - *(p1+1) * *p2;
      *po++ = r;
      *po++ = i;
      p1++;
      p1++;
      p2++;
      p2++;
    }
  }
   
ir_error:
  
  return error;
}
