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
#include <stdlib.h>
#include "ir.h"
#include "stdio.h"

#define IR_PI 3.14159265358979323846264
#define IR_SMALL 1e-20

extern double j0(double);
extern double j1(double);

/* the function below tabellates a number of values to integrate over. Especially the Bessel functions are tabellated for speed reasons */
static void ir_FillAL2D(double sina, double sfu, double sfv, 
			int u, int v, double h, double *Ar, double *Ai,
			double *L0, double *L1, double *L2, int ScalarTheory)
{  
  double sin2a, sinh, cosh, coshsin2a, sinha, vsinha, t, f0, f1, f2;
  int ii;

  sin2a = sina * sina;
  sinh = sin(h);
  cosh = cos(h);

  coshsin2a = cosh / sin2a;
  for (ii = 0; ii < u; ii++) {
    t = -sfu * (double) (ii - u + 1);
    Ar[ii] = cos(t * coshsin2a);
    Ai[ii] = sin(t * coshsin2a);
  }

  sinha = sinh / sina;
  if (ScalarTheory)
    f0 = sqrt(cosh) * sinh * (2.0);   /* this helps to keep the energy constant for the scalar theory */
  else
          f0 = sqrt(cosh) * sinh * (1.0 + cosh);
  f1 = sqrt(cosh) * sinh * sinh;
  f2 = sqrt(cosh) * sinh * (1.0 - cosh);
  for(ii = 0; ii < v; ii++) {
    vsinha =  sfv * sinha * (double) ii;
    *L0 = j0(vsinha);
    *L1 = j1(vsinha);
    *L2 = ii == 0 ? 0.0 : 2.0 * *L1 / vsinha - *L0;
    *L0++ *= f0;
    *L1++ *= f1;
    *L2++ *= f2;
  }
}

static ir_Error ir_Integrate2D(double sina,   /* sin(alpha) */
			       double sfu, double sfv,    /* scaling factor for u and v */
			       int u, int v,   /* points in u and v */
                               double *I0r, double *I0i, 
			       double *I1r, double *I1i, double *I2r, 
			       double *I2i, int max, double pr,
                               int ScalarTheory)
{
  ir_Error error = IR_OK;
  int sz, ii, jj, hh, kk, ll, gg;
  double *Ar = 0, *Ai = 0, *L0 = 0, *L1 = 0, *L2 = 0;
  double *tI0r = 0, *tI0i = 0, *tI1r = 0, *tI1i = 0, *tI2r = 0, *tI2i = 0;
  double h, x1, x, tmp, mp;

  sz = u * v;

  tI0r = (double *) malloc(sz * sizeof(double));
  tI0i = (double *) malloc(sz * sizeof(double));
  tI1r = (double *) malloc(sz * sizeof(double));
  tI1i = (double *) malloc(sz * sizeof(double));
  tI2r = (double *) malloc(sz * sizeof(double));
  tI2i = (double *) malloc(sz * sizeof(double));

  Ar = (double *) malloc(u * sizeof(double));
  Ai = (double *) malloc(u * sizeof(double));
  L0 = (double *) malloc(v * sizeof(double));
  L1 = (double *) malloc(v * sizeof(double));
  L2 = (double *) malloc(v * sizeof(double));

  IRTJ(!Ar || !Ai || !L0 || !L1 || !L2 || !tI0r || !tI0i || !tI1r || !tI1i || 
       !tI2r || !tI2i, IR_NO_MEMORY);

  h = asin(sina);  
  x1 = 0.0;

  for(ii = 0; ii < sz; ii++) 
    I0r[ii] = I0i[ii] = I1r[ii] = I1i[ii] = I2r[ii] = I2i[ii] = 0.0;    /* clear all I0, I1, and I2 values */

  ir_FillAL2D(sina, sfu, sfv, u, v, h, Ar, Ai, L0, L1, L2,ScalarTheory);
  ll = 0;
  for (ii = 0; ii < u; ii++) {
    for(jj = 0; jj < v; jj++) {
      I0r[ll] += 0.5 * h * Ar[ii] * L0[jj]; 
      I0i[ll] += 0.5 * h * Ai[ii] * L0[jj];
      I1r[ll] += 0.5 * h * Ar[ii] * L1[jj]; 
      I1i[ll] += 0.5 * h * Ai[ii] * L1[jj];
      I2r[ll] += 0.5 * h * Ar[ii] * L2[jj]; 
      I2i[ll] += 0.5 * h * Ai[ii] * L2[jj];
      ++ll;
    }
  }
  
  hh = 1;
  x1 -= h;
  for(kk = 0; kk < max; kk++) {
    x1 += 0.5 * h;
    x = x1 + h;
    ll = 0;
    ir_FillAL2D(sina, sfu, sfv, u, v, x, Ar, Ai, L0, L1, L2,ScalarTheory);
    for (ii = 0; ii < u; ii++) {
      for(jj = 0; jj < v; jj++) {
	tI0r[ll] = Ar[ii] * L0[jj]; 
	tI0i[ll] = Ai[ii] * L0[jj];
	tI1r[ll] = Ar[ii] * L1[jj]; 
	tI1i[ll] = Ai[ii] * L1[jj];
	tI2r[ll] = Ar[ii] * L2[jj]; 
	tI2i[ll] = Ai[ii] * L2[jj];
	++ll;
      }
    }
    for(gg = 2; gg <= hh; gg++) {
      x += h;
      ll = 0;
      ir_FillAL2D(sina, sfu, sfv, u, v, x, Ar, Ai, L0, L1, L2,ScalarTheory);
      for (ii = 0; ii < u; ii++) {
	for(jj = 0; jj < v; jj++) {
	  tI0r[ll] += Ar[ii] * L0[jj]; 
	  tI0i[ll] += Ai[ii] * L0[jj];
	  tI1r[ll] += Ar[ii] * L1[jj]; 
	  tI1i[ll] += Ai[ii] * L1[jj];
	  tI2r[ll] += Ar[ii] * L2[jj]; 
	  tI2i[ll] += Ai[ii] * L2[jj];
	  ++ll;
	}
      }
    }
 
    hh *= 2;
    h *= 0.5;   
 
    mp = 0.0;
    for(gg = 0; gg < sz; gg++) {
      tmp = I0r[gg];
      I0r[gg] = 0.5 * I0r[gg] + h * tI0r[gg]; 
      mp += tmp > IR_SMALL ? fabs(I0r[gg] - tmp) / tmp : 0.0;
      tmp = I0i[gg];
      I0i[gg] = 0.5 * I0i[gg] + h * tI0i[gg];
      mp += tmp > IR_SMALL ? fabs(I0i[gg] - tmp) / tmp : 0.0;
      tmp = I1r[gg];
      I1r[gg] = 0.5 * I1r[gg] + h * tI1r[gg];
      mp += tmp > IR_SMALL ? fabs(I1r[gg] - tmp) / tmp : 0.0;
      tmp = I1i[gg];
      I1i[gg] = 0.5 * I1i[gg] + h * tI1i[gg];
      mp += tmp > IR_SMALL ? fabs(I1i[gg] - tmp) / tmp : 0.0;
      tmp = I2r[gg];
      I2r[gg] = 0.5 * I2r[gg] + h * tI2r[gg];
      mp += tmp > IR_SMALL ? fabs(I2r[gg] - tmp) / tmp : 0.0;
      tmp = I2i[gg];
      I2i[gg] = 0.5 * I2i[gg] + h * tI2i[gg];
      mp += tmp > IR_SMALL ? fabs(I2i[gg] - tmp) / tmp : 0.0;
    }
    if (mp / sz / 6.0 < pr) 
      break;
  }

 ir_error:

  if (Ar) free(Ar);
  if (Ai) free(Ai);
  if (L0) free(L0);
  if (L1) free(L1);
  if (L2) free(L2);
  if (tI0r) free(tI0r);
  if (tI0i) free(tI0i);
  if (tI1r) free(tI1r);
  if (tI1i) free(tI1i);
  if (tI2r) free(tI2r);
  if (tI2i) free(tI2i);
  
  return error;
}

#define AbsSqr(real,imag) (((real)*(real))+((imag)*(imag)))

static ir_Error ir_WidefieldPSF2D(double *out, double sfu, double sfv, int u,
				  int v, double sina, int ms, double pr, int Pi4, double relIntensity, double relPhase, int ScalarTheory)
{
  ir_Error error = IR_OK;
  int hh, s;
  double *I0r = 0, *I0i = 0, *I1r = 0, *I1i = 0, *I2r = 0, *I2i = 0;
  double tmp_r,tmp_i;
  double a_r=sqrt(relIntensity)*cos(2*IR_PI*relPhase),a_i=sqrt(relIntensity)*sin(2*IR_PI*relPhase);

  s = u * v;

  I0r = (double *) malloc(s * sizeof(double));
  I0i = (double *) malloc(s * sizeof(double));
  I1r = (double *) malloc(s * sizeof(double));
  I1i = (double *) malloc(s * sizeof(double));
  I2r = (double *) malloc(s * sizeof(double));
  I2i = (double *) malloc(s * sizeof(double));

  IRTJ(!I0r || !I0i || !I1r || !I1i || !I2r || !I2i, IR_NO_MEMORY);

  IRXJ(ir_Integrate2D(sina, sfu, sfv, u, v, I0r, I0i, I1r, I1i, I2r, I2i,
		      ms, pr, ScalarTheory));  /* For now: allways compute vector theory for Widefield PSF */

if (ScalarTheory)  /* only account for I0 */
        {                  /* The equations are like in Nagorni2001 but with an additional factor of i before the field */
                if (! Pi4)
		  for (hh = 0; hh < s; hh++) {
    			*out++ = (I0r[hh] * I0r[hh] + I0i[hh] * I0i[hh]);                    /* |I2|^2 */
                        }
		else
		  for (hh = 0; hh < s; hh++) {
    			tmp_r = (a_r+1) * I0r[hh] - a_i*I0i[hh];
    			tmp_i = a_i*I0r[hh]+(a_r-1)*I0i[hh];
    			*out = AbsSqr(tmp_r,tmp_i);
  			}
        }
else
/* Electric field in a widefield Microscope is E = -i (I0 ,  I2 , I1 sqrt(2)/i )
 with circularly polarised light, with intensities corresponding to random polarisation as well
 See Nagorni et al. JOSA A 18 (2001) 1-13    */
  if (! Pi4)
  for (hh = 0; hh < s; hh++) {
    *out++ = (I0r[hh] * I0r[hh] + I0i[hh] * I0i[hh] +                   /* |I0|^2 */
	     2.0 * (I1r[hh] * I1r[hh] + I1i[hh] * I1i[hh]) +            /* 2*|I1|^2 */
	     I2r[hh] * I2r[hh] + I2i[hh] * I2i[hh]);                    /* |I2|^2 */
  }
  else  /* Do the 4Pi case */
  for (hh = 0; hh < s; hh++) {
/*out++ = (I0r[hh] * I0r[hh] +                  |I0|^2 
	     2.0 * I1i[hh] * I1i[hh] +          2*|I1|^2 
	     I2r[hh] * I2r[hh]);                  |I2|^2 

 The general solution extended from Nagornis paper for any complex amplitude ar and ai of the second wave is
  |a I0 + conj(aI0)|^2 +
  |a I2 - conj(I2)|^2  +
  |i a sqrt(2) I1 - i sqrt(2) conj(I1) |^2     */
    tmp_r = (a_r+1) * I0r[hh] - a_i*I0i[hh];
    tmp_i = a_i*I0r[hh]+(a_r-1)*I0i[hh];
    *out = AbsSqr(tmp_r,tmp_i);
    tmp_r = (a_r+1) * I2r[hh] - a_i*I2i[hh];
    tmp_i = a_i*I2r[hh]+(a_r-1)*I2i[hh];
    *out += AbsSqr(tmp_r,tmp_i);
    tmp_r = (a_r-1) * I1r[hh] - a_i*I1i[hh];
    tmp_i = a_i*I1r[hh]+(a_r+1)*I1i[hh];
    *out++ += 2*AbsSqr(tmp_r,tmp_i);
  }


 ir_error:

  if (I0r) free(I0r);
  if (I0i) free(I0i);
  if (I1r) free(I1r);
  if (I1i) free(I1i);
  if (I2r) free(I2r);
  if (I2i) free(I2i);

  return error;
}


#define maxs 15
#define oversampling 2

ir_Error ir_WidefieldASF(ir_Image *EX, ir_Image *EY, ir_Image *EZ, double sx, double sy, double sz,
			 double l, double na, double ri, double pr, int circPol, int ScalarTheory, int Pi4)
{
  ir_Error error = IR_OK;
  int x, y, z, x2, y2, z2, um, vm, ii, jj, kk, zi, idx, d, *dims = 0;
  double sina, sv, sfu, sfv, *sl = 0, scy, scx, scx2,scy2, v;
  double *I0r = 0, *I0i = 0, *I1r = 0, *I1i = 0, *I2r = 0, *I2i = 0;
  double cosPhi = cos (IR_PI / 4.0), cos2Phi = cos (IR_PI / 2.0), sin2Phi = sin (IR_PI / 2.0), Phi;
  double Ex1r,Ex1i,Ex2r,Ex2i,Ey1r,Ey1i,Ey2r,Ey2i,Ez1r,Ez1i,Ez2r,Ez2i, conjugate=1;

  ir_Float *px, * py, * pz;

  IRXJ(ir_Check(EX));
  IRXJ(ir_Check(EY));
  IRXJ(ir_Check(EZ));

  IRTJ(sx <= 0, IR_BAD_PARAMETER);
  IRTJ(sy <= 0, IR_BAD_PARAMETER);
  IRTJ(sz <= 0, IR_BAD_PARAMETER);
  IRTJ(l <= 0,  IR_BAD_PARAMETER);
  IRTJ(na <= 0, IR_BAD_PARAMETER);
  IRTJ(ri <= 0, IR_BAD_PARAMETER);

  IRXJ(ir_ImageDimensions(EX, &dims, &d));
  IRTJ(d > 3, IR_DIMENSIONALITY_UNSUPPORTED);
    
  x = dims[0] / 2;   /* this accounts for complex numbers */ 
  y = d > 1 ? dims[1] : 1;
  z = d > 2 ? dims[2] : 1;
    
  if (pr <= 0.0)
    pr = 0.01;

  x2 = x / 2;
  y2 = y / 2;
  z2 = z / 2;   /* due to symmetry */ 

  sina = na / ri;

  um = z2 + 1;
  sv = sqrt(sx * sx + sy * sy) / oversampling;
  vm = (int)(sqrt(sx * sx * x2 * x2 + sy * sy * y2 * y2) / sv) + oversampling;

  sfu = 2.0 * IR_PI * ri * sina * sina * sz / l;
  sfv = 2.0 * IR_PI * ri * sina * sv / l;

  I0r = (double *) malloc((um*vm) * sizeof(double));
  I0i = (double *) malloc((um*vm) * sizeof(double));
  I1r = (double *) malloc((um*vm) * sizeof(double));
  I1i = (double *) malloc((um*vm) * sizeof(double));
  I2r = (double *) malloc((um*vm) * sizeof(double));
  I2i = (double *) malloc((um*vm) * sizeof(double));

  IRTJ(!I0r || !I0i || !I1r || !I1i || !I2r || !I2i, IR_NO_MEMORY);

  /* makes the 2D I0, I1 and I2 information in u,v space */
  IRXJ(ir_Integrate2D(sina, sfu, sfv, um, vm, I0r, I0i, I1r, I1i, I2r, I2i,
		      maxs, pr,ScalarTheory));

  IRXJ(ir_ImageData(EX, &px));   /* pointer p is set to the beginning of out */
  IRXJ(ir_ImageData(EY, &py));   /* pointer p is set to the beginning of out */
  IRXJ(ir_ImageData(EZ, &pz));   /* pointer p is set to the beginning of out */
  /* interpolates the 2D PSF to a 3D volume */
  for (ii = -z2; ii < z - z2; ii++) {      /* This is the z direction than means also u direction */
  if (ii < 0)    /* accounts for symmetry in u */
    {
    zi = (z2 + ii) * vm;  
    conjugate=-1;
    }
  else
    {
    zi = (z2 - ii) * vm;
    conjugate=1;
    }
    for(jj = -y2; jj < y - y2; jj++) {     /* This is the y direction */
      scy = sy * jj;
      scy2 = scy * scy;
      for(kk = -x2; kk < x - x2; kk++) {   /* This is the x direction */
	scx = sx * kk;
	scx2 = scx * scx;
	v = sqrt(scx2 + scy2);
	idx = (int) (v / sv);
        if (! circPol)
        {
            if (v != 0)
            {
                  Phi = atan2(scy,scx);  /* first the square was accidently used here */
                  cosPhi = cos(Phi);
                  cos2Phi = cos(2.0*Phi);
                  sin2Phi = sin(2.0*Phi);
            }
            else
            {
                  Phi = 0;  /* first the square was accidently used here */
                  cosPhi = cos(Phi);
                  cos2Phi = cos(2.0*Phi);
                  sin2Phi = sin(2.0*Phi);
            }
        }
        if (ScalarTheory)  /* only account for I0 */
        {                  /* The equations are like in Nagorni2001 but with an additional factor of i before the field */
        	Ex1r = I0r[zi + idx];
	        Ex1i = conjugate*I0i[zi + idx];
        	Ex2r = I0r[zi + idx+1];
        	Ex2i = conjugate*I0i[zi + idx+1];
                if (Pi4)
                        {
                            Ex1r += Ex1r;
                            Ex1i -= Ex1i; /* this is because of the negative u with an inverted conjugate: should result to zero */
                            Ex2r += Ex2r;
                            Ex2i -= Ex2i; /* should result to zero */
                        }
        	*px++ = Ex1r + (Ex2r - Ex1r) * (v - idx * sv) / sv;   /* Linear interpolation */
        	*px++ = Ex1i + (Ex2i - Ex1i) * (v - idx * sv) / sv;   /* Linear interpolation */
        }
        else
        {
	Ex1r = I0r[zi + idx] + cos2Phi*I2r[zi + idx];
	Ex1i = conjugate*(I0i[zi + idx] + cos2Phi*I2i[zi + idx]);
	Ex2r = I0r[zi + idx+1] + cos2Phi*I2r[zi + idx+1];
	Ex2i = conjugate*(I0i[zi + idx+1] + cos2Phi*I2i[zi + idx+1]);

        Ey1r = sin2Phi*I2r[zi + idx];
	Ey1i = conjugate*(sin2Phi*I2i[zi + idx]);
	Ey2r = sin2Phi*I2r[zi + idx+1];
	Ey2i = conjugate*(sin2Phi*I2i[zi + idx+1]);

        Ez1r = conjugate*2.0*cosPhi*I1i[zi + idx];      /* the conjugation of I1 affects only the real part of Ez */
	Ez1i = -2.0*cosPhi*I1r[zi + idx];
	Ez2r = conjugate*2.0*cosPhi*I1i[zi + idx+1];
	Ez2i = -2.0*cosPhi*I1r[zi + idx+1];
        if (Pi4)
           {             /* incorporates the Y and Z directional inversion matrix (page 4 Nagorni) */
           Ex1r += Ex1r;
           Ex1i -= Ex1i; /* u -> -u : should result to zero: Phi is inverted, but cos(2 Phi)=cos(- 2 Phi) */
           Ex2r += Ex2r;
           Ex2i -= Ex2i; /* should result to zero: Phi is inverted, but cos(2 Phi)=cos(- 2 Phi) */

           Ey1r += Ey1r; /* This is because -sin2Phi = sin(-2Phi) and an additional multiplication by -1 */
           Ey1i += Ey1i; /* should result to zero */
           Ey2r += Ey2r;
           Ey2i += Ey2i; /* should result to zero */

           Ez1r -= Ez1r;  /* cos(Phi) = cos(-Phi) but a multiplication by -1 */
           Ez1i += Ez1i; /* should result to zero */
           Ez2r -= Ez2r;
           Ez2i += Ez2i; /* should result to zero */
           }
	*px++ = Ex1r + (Ex2r - Ex1r) * (v - idx * sv) / sv;   /* Linear interpolation */
	*px++ = Ex1i + (Ex2i - Ex1i) * (v - idx * sv) / sv;   /* Linear interpolation */
	*py++ = Ey1r + (Ey2r - Ey1r) * (v - idx * sv) / sv;   /* Linear interpolation */
	*py++ = Ey1i + (Ey2i - Ey1i) * (v - idx * sv) / sv;   /* Linear interpolation */
	*pz++ = Ez1r + (Ez2r - Ez1r) * (v - idx * sv) / sv;   /* Linear interpolation */
	*pz++ = Ez1i + (Ez2i - Ez1i) * (v - idx * sv) / sv;   /* Linear interpolation */
        }    
      }
    }
  }

  /* no normalization */
  
 ir_error:
  
  if (dims) free(dims);
  if (sl) free(sl);
  if (I0r) free(I0r);
  if (I0i) free(I0i);
  if (I1r) free(I1r);
  if (I1i) free(I1i);
  if (I2r) free(I2r);
  if (I2i) free(I2i);

  return error;
}

ir_Error ir_WidefieldPSF(ir_Image *out, double sx, double sy, double sz,
			 double l, double na, double ri, double pr, int Pi4, double relIntensity, double relPhase,int ScalarTheory)
{
  ir_Error error = IR_OK;
  int x, y, z, x2, y2, z2, um, vm, ii, jj, kk, zi, idx,  d, *dims = 0, size;
  double sina, sv, sfu, sfv, *sl = 0, * pos_sl=0, * neg_sl=0, scy, scx, v, sl1, sl2, sum;
  ir_Float *p;

  IRXJ(ir_Check(out));

  IRTJ(sx <= 0, IR_BAD_PARAMETER);
  IRTJ(sy <= 0, IR_BAD_PARAMETER);
  IRTJ(sz <= 0, IR_BAD_PARAMETER);
  IRTJ(l <= 0,  IR_BAD_PARAMETER);
  IRTJ(na <= 0, IR_BAD_PARAMETER);
  IRTJ(ri <= 0, IR_BAD_PARAMETER);

  IRXJ(ir_ImageDimensions(out, &dims, &d));
  IRTJ(d > 3, IR_DIMENSIONALITY_UNSUPPORTED);

  x = dims[0];
  y = d > 1 ? dims[1] : 1;
  z = d > 2 ? dims[2] : 1;

  if (pr <= 0.0)
    pr = 0.01;

  x2 = x / 2;
  y2 = y / 2;
  z2 = z / 2;

  sina = na / ri;

  um = z2 + 1;
  sv = sqrt(sx * sx + sy * sy) / oversampling;
  vm = (int)(sqrt(sx * sx * x2 * x2 + sy * sy * y2 * y2) / sv) + oversampling;

  sfu = 2.0 * IR_PI * ri * sina * sina * sz / l;
  sfv = 2.0 * IR_PI * ri * sina * sv / l;

  pos_sl = (double*) malloc(um * vm * sizeof(double));

  IRXJ(ir_WidefieldPSF2D(pos_sl, sfu, sfv, um, vm, sina, maxs, pr, Pi4,relIntensity,relPhase,ScalarTheory));

  if (Pi4)
  {
    neg_sl = (double*)malloc(um * vm * sizeof(double));
    IRXJ(ir_WidefieldPSF2D(neg_sl, sfu, sfv, um, vm, sina, maxs, pr, Pi4,relIntensity,-relPhase,ScalarTheory));
  }
  else
   neg_sl =pos_sl;


  IRXJ(ir_ImageData(out, &p));
  for (ii = -z2; ii < z - z2; ii++) {
    if (ii >= 0)
	sl =pos_sl;
    else
	sl =neg_sl;

    zi = (z2 - (ii < 0 ? -ii : ii)) * vm;
    for(jj = -y2; jj < y - y2; jj++) {
      scy = sy * jj;
      scy = scy * scy;
      for(kk = -x2; kk < x - x2; kk++) {
	scx = sx * kk;
	scx = scx * scx;
	v = sqrt(scx + scy);
	idx = (int) (v / sv);
	sl1 = sl[zi + idx];
	sl2 = sl[zi + idx + 1];
	*p++ = sl1 + (sl2 - sl1) * (v - idx * sv) / sv;
      }
    }
  }

  /* now perform a normalization */

  sum = 0.0;
  IRXJ(ir_ImageData(out, &p));
  IRXJ(ir_ImageSize(out, &size));
  for (ii = 0; ii < size; ii++)
    sum += *p++;
  IRXJ(ir_ImageData(out, &p));
  if (fabs(sum) > IR_SMALL)
    for (ii = 0; ii < size; ii++)
      *p++ /= sum;

 ir_error:

  if (dims) free(dims);
  if (sl) free(sl);

  return error;
}


static ir_Error ir_Ellipsoid(ir_Image *im, double rx, double ry)  /* Generates the Fourier-transform of a circular disc */
{
  ir_Error error = IR_OK;
  double fw, fh, r, xys, j1(double);
  int x, y, w, h, *dims = 0;
  ir_Float *p;
  double *iX=0, *iY=0;

  IRXJ(ir_ImageDimensions(im, &dims, 0));
  w = dims[0];
  h = dims[1];
  IRXJ(ir_ImageRealX(im, &x));
  fw = (double) x / rx;
  fh = (double) h / ry;
  xys = rx * ry;

  iX = (double *) malloc(w * sizeof(double));
  iY = (double *) malloc(h * sizeof(double));
  
  for(x=0; x<w; x++) iX[x] = x * x / fw / fw;
  x = 0;
  for(y=0; y<(h+1)/2; y++) iY[x++] = y * y / fh / fh;
  for(y=-h/2; y<0; y++) iY[x++] = y * y / fh / fh;

  IRXJ(ir_ImageData(im, &p));
  for(y = 0; y < h; y++) {
    for(x = 0; x < w; x++) {
      r = sqrt(iX[x] + iY[y]);
      if (r == 0.0) 
	*p++ = IR_PI * xys;
      else
	*p++ = xys * j1(2.0 * IR_PI * r) / r;
      *p++ = 0.0;
    }
  }

 ir_error:

  if (dims) free(dims);
  if (iX) free(iX);
  if (iY) free(iY);

  return error;
}


ir_Error ir_ConfocalPSF(ir_Image *out, double sx, double sy, double sz, 
			double ex, double em, double na, double ri, double d, 
			double pr, int Pi4Ex,int Pi4Em, int Normalize, int twophoton, double relExInt, double relExPhase,double relEmInt, double relEmPhase, int ScalarTheory)
{
  ir_Image *ti = 0, *ts = 0, *pinhole = 0;
  ir_Float *pi, *po;
  int x, y, z, i, j, s, e, SizeZ;
  int *dims = 0, dims2[5]={0,0,0,0,0};
  double t;
  ir_Error error = IR_OK;
  
  IRXJ(ir_Check(out));

  IRTJ(sx <= 0, IR_BAD_PARAMETER);
  IRTJ(sy <= 0, IR_BAD_PARAMETER);
  IRTJ(sz <= 0, IR_BAD_PARAMETER);
  IRTJ(ex <= 0, IR_BAD_PARAMETER);
  IRTJ(em <= 0, IR_BAD_PARAMETER);
  IRTJ(na <= 0, IR_BAD_PARAMETER);
  IRTJ(ri <= 0, IR_BAD_PARAMETER);
  IRTJ(pr <= 0.0, IR_BAD_PARAMETER);

  IRXJ(ir_ImageDimensions(out, &dims, &e));
  IRTJ(e > 3, IR_DIMENSIONALITY_UNSUPPORTED);
  SizeZ=dims[2];
  if (SizeZ==0) SizeZ=1;

  for (i=0;i<5;i++) dims2[i]=dims[i];
  if ((dims2[0] / 2)*2 != dims2[0])
    dims2[0]+=1;
  
  x = dims[0];
  y = e > 1 ? dims[1] : 1;
  z = e > 2 ? dims[2] : 1;

  IRXJ(ir_WidefieldPSF(out, sx, sy, sz, em, na, ri, pr, Pi4Em,relEmInt,relEmPhase,ScalarTheory));  /* compute emission PSF */

  if (d > 0.0) {
    if (e > 2) {
      e--;
      dims[2] = 1;
    }
    IRXJ(ir_Init(&ti));
    IRXJ(ir_Init(&ts));
    IRXJ(ir_Init(&pinhole));
    IRXJ(ir_ImageChange(ti, e, dims2, IR_TYPE_IMAGE));
    IRXJ(ir_ImageChange(ts, e, dims2, IR_TYPE_SPECTRUM));
    IRXJ(ir_ImageChange(pinhole, e, dims2, IR_TYPE_SPECTRUM));
    IRXJ(ir_Ellipsoid(pinhole, 0.5 * d / sx, 0.5 * d / sy));  /* generate a Fourier-transformed pinhole (in Fourierspace) */

    s = x * y;
    IRXJ(ir_ImageData(out, &po));  /* gets a pointer to access the image */
    for(i = 0; i < z; i++) {     
      IRXJ(ir_ImageData(ti, &pi));
      for (j = 0; j < s; )   /* Copy content in out into ti */
	{
	  *pi++ = *po++;
	  j++;
	  if (j%dims[0]==0 && dims[0] != dims2[0])
	    {
	    (*pi)=0.0;
	    pi++;  /* Skip the extra bit */
	    }
	}
      IRXJ(ir_ForwardRft(ti, ts));
      IRXJ(ir_Mul(ts, pinhole, ts));
      IRXJ(ir_InverseRft(ts, ti));
      /* IRXJ(ir_InverseRft(pinhole, ti));   Just a test */
      IRXJ(ir_ImageData(ti, &pi));
      po -= s;
      for (j = 0; j < s; )   /* Copy content in ti into out */
	{
	  *po++ = *pi++;
	  j++;
	  if (j%dims[0]==0 && dims[0] != dims2[0])
	    {
	    pi++;  /* Skip the extra bit */
	    }
	}
    }
    
    IRXJ(ir_Destroy(&ti));
    IRXJ(ir_Destroy(&ts));
    IRXJ(ir_Destroy(&pinhole));
  }

  IRXJ(ir_Init(&ti));
  IRXJ(ir_Clone(ti, out));
  IRXJ(ir_WidefieldPSF(ti, sx, sy, sz, ex, na, ri, pr, Pi4Ex,relExInt,relExPhase,ScalarTheory));  /* compute excitation PSF */
  if (twophoton)
  	IRXJ(ir_Mul(ti, ti, ti));  /* Square the Excitation PSF */

  IRXJ(ir_Mul(ti, out, out));
  IRXJ(ir_Destroy(&ti));

  if (Normalize != 0)
    {
      IRXJ(ir_Sum(out, &t));
      IRXJ(ir_DivConst(out, out, t, 0.0));  /* real, imaginary */
    }
  else
      IRXJ(ir_MulConst(out, out, (double) SizeZ, 0.0));  /* corrects for the size-effect along Z to correspond to WF calculations. real, imaginary */

ir_error:

  if (dims) free(dims);
  ir_Destroy(&ti);
  ir_Destroy(&ts);
  ir_Destroy(&pinhole);

  return error;
}

#undef IR_PI
#undef IR_SMALL

