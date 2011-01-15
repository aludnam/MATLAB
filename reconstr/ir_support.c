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
#include <stdlib.h>
#include <stdio.h>
#include "ir.h"

#define TYPE  type
#define DIM   d
#define REALX realX
#define SIZES sizes
#define DATA  data

ir_Error ir_Check(ir_Image *in)
{
  ir_Error error = IR_OK;
  
  IRTJ(!in, IR_NONEXISTANT);
  IRTJ(!(in->TYPE == IR_TYPE_IMAGE || in->TYPE == IR_TYPE_SPECTRUM),
       IR_TYPE_UNSUPPORTED);

  IRTJ(in->DATA == 0, IR_EMPTY);
  
ir_error:
  
  return error;
}


ir_Error ir_CheckTwo(ir_Image *in1, ir_Image *in2)
{
  int ii;
  ir_Error error = IR_OK;

  IRTJ(!in1 || !in2, IR_NONEXISTANT);
  IRTJ(in1->TYPE != in2->TYPE, IR_TYPE_DIFFERS);
  IRTJ(!(in1->TYPE == IR_TYPE_IMAGE || in1->TYPE == IR_TYPE_SPECTRUM) ||
       !(in2->TYPE == IR_TYPE_IMAGE || in2->TYPE == IR_TYPE_SPECTRUM),
       IR_TYPE_UNSUPPORTED);

  IRTJ(in1->DATA == 0 || in2->DATA == 0, IR_EMPTY);
  IRTJ(in1->DIM != in2->DIM, IR_DIMENSIONALITY_DIFFERS);
  for (ii=0; ii < in1->DIM; ii++) {
    IRTJ(in1->SIZES[ii] != in2->SIZES[ii], IR_SIZES_DIFFER);
  }

ir_error:

  return error;
}

ir_Error ir_ImageType(ir_Image *im, int *r)
{
  *r = im->TYPE;

  return IR_OK;
}


ir_Error ir_ImageDimensionality(ir_Image *im, int *r)
{
  ir_Error error = IR_OK;
  
  IRXJ(ir_Check(im));

  *r = im->DIM;

ir_error:

  return error;
}

ir_Error ir_ImageDimensions(ir_Image *im, int **r, int *d)
{
  int t, ii;
  ir_Error error = IR_OK;

  IRXJ(ir_Check(im));

  IRXJ(ir_ImageDimensionality(im, &t));

  if (d) *d = t;
  
  *r = (int *) malloc(t * sizeof(int));
  
  for(ii=0; ii<t; ii++) (*r)[ii] = im->SIZES[ii];
  
ir_error:

  return error;
}

ir_Error ir_ImageSize(ir_Image *im, int *r)
{
  int *sizes = 0;
  int d, ii;
  ir_Error error = IR_OK;

  IRXJ(ir_Check(im));
  
  IRXJ(ir_ImageDimensions(im, &sizes, &d));
  
  *r = 1;
  for(ii=0; ii<d; ii++) *r *= sizes[ii];
  
ir_error:

  if (sizes) free(sizes);
  return error;
}

ir_Error ir_ImageData(ir_Image *im, ir_Float **r)
{
  ir_Error error = IR_OK;

  IRXJ(ir_Check(im));
  
  *r = im->DATA;

ir_error:

  return error;
}

ir_Error ir_ImageRealX(ir_Image *im, int *r)
{
  ir_Error error = IR_OK;

  IRXJ(ir_Check(im));

  *r = im->REALX;
  
ir_error:

  return error;
}

ir_Error ir_ImageRealSize(ir_Image *im, int *r)
{
  int ii, d, x;
  ir_Error error = IR_OK;
  int *sizes = 0;

  IRXJ(ir_Check(im));
  
  IRXJ(ir_ImageDimensions(im, &sizes, &d));
  
  *r = 1;
  for(ii=1; ii<d; ii++) *r *= sizes[ii];
  
  IRXJ(ir_ImageRealX(im, &x));

  *r *= x;


ir_error:
  if (sizes) free(sizes);

  return error;
}


ir_Error ir_Init(ir_Image **im)
{
  *im = (ir_Image *) malloc(sizeof(ir_Image));
  (*im)->TYPE = 0;
  (*im)->DIM = 0;
  (*im)->SIZES = 0;
  (*im)->DATA = 0;
  (*im)->REALX = 0;

  return IR_OK;
}


ir_Error ir_Destroy(ir_Image **im)
{
  ir_Error error = IR_OK;
  
  if (*im) {
    if ((*im)->DATA) free((*im)->DATA);
    if ((*im)->SIZES) free((*im)->SIZES);
    free(*im);
    *im = 0;
  }
  
  return error;
}


ir_Error ir_ImageLink(ir_Image *im, ir_Float *data, int d, int *dims, 
		      int type)
{
  int ii, *sizes = 0;
  ir_Error error = IR_OK;
  
  IRTJ(dims == 0 || d <= 0, IR_EMPTY_INPUT);

  IRTJ(!(type == IR_TYPE_IMAGE || type == IR_TYPE_SPECTRUM), 
       IR_TYPE_UNSUPPORTED);

  IRTJ(!im, IR_NONEXISTANT);

  im->TYPE = type;
  im->DIM = d;
  sizes = (int *) malloc(d * sizeof(int));
  for(ii=1; ii<d; ii++) sizes[ii] = dims[ii];
  switch(type) {
  case IR_TYPE_IMAGE:
    sizes[0] = dims[0];
    break;
  case IR_TYPE_SPECTRUM:
    sizes[0] = dims[0] / 2 + 1;
    break;
  }
  im->SIZES = sizes;
  im->DATA = data;
  im->REALX = dims[0];

ir_error:

  return error;
}


ir_Error ir_ImageUnlink(ir_Image **im)
{
  ir_Error error = IR_OK;

  if (*im) {
    if ((*im)->SIZES) free((*im)->SIZES);
    free(*im);
    *im = 0;
  }

  return error;
}


ir_Error ir_ImageChange(ir_Image  *im, int d, int *dims, int type)
{
  int ii, s, e, *osizes = 0;
  int *sizes = 0;
  ir_Error error = IR_OK;

  IRTJ(!(type == IR_TYPE_IMAGE || type == IR_TYPE_SPECTRUM), 
       IR_TYPE_UNSUPPORTED);
  IRTJ(dims  == 0 || d <= 0, IR_EMPTY_INPUT); 
  
  if (im->DATA) {
    IRXJ(ir_ImageDimensions(im, &sizes, &e));
    if (im->TYPE == type && e == d) {
      s = 1;
      switch(type) {
      default:
      case IR_TYPE_IMAGE:
	if (sizes[0] != dims[0]) s = 0;
	break;
      case IR_TYPE_SPECTRUM:
	if (sizes[0] != dims[0] / 2 + 1) s = 0;
	break;
      }
      for (ii=1; ii<d; ii++) 
	if (sizes[ii] != dims[ii]) {
	  s = 0;
	  break;
	}
      if (s) {
	goto ir_error;
      }
    }
  }
  
  osizes = (int *) malloc(d * sizeof(int));
  
  switch(type) {
  default:
  case IR_TYPE_IMAGE:
    osizes[0] = s = dims[0];
    break;
  case IR_TYPE_SPECTRUM:
    osizes[0] = s = dims[0] / 2 + 1;
    break;
  }
  for(ii=1; ii<d; ii++) {
    osizes[ii] = dims[ii] > 0 ? dims[ii] : 1;
    s *= osizes[ii];
  } 
  
  if (im->SIZES) free(im->SIZES);
  im->SIZES = osizes;
  im->DIM = d;
  im->REALX = dims[0];
  im->TYPE = type;
  
  if (type == IR_TYPE_SPECTRUM) s *= 2;
  
  if (im->DATA) free(im->DATA);
  im->DATA = (ir_Float *) malloc(s * sizeof (ir_Float));
  
  IRTJ(im->DATA == 0, IR_NO_MEMORY);

ir_error:  
  if  (sizes) free(sizes);

  return error;
}

ir_Error ir_Clone(ir_Image *im, ir_Image *other)
{
  int t, d;
  int *osizes = 0;
  ir_Error error = IR_OK;

  if (im == other) return IR_OK;

  IRXJ(ir_Check(other));
  IRXJ(ir_ImageType(other, &t));

  IRXJ(ir_ImageDimensions(other, &osizes, &d));
  if (t == IR_TYPE_SPECTRUM) IRXJ(ir_ImageRealX(other, osizes));
  IRXJ(ir_ImageChange(im, d, osizes,  t));
    
 ir_error:
  if (osizes) free(osizes);

  return error;
}



