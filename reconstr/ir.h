#ifndef IR_SUPPORT_H
#define IR_SUPPORT_H

#include <stdlib.h>

/* data type */
#ifdef DOUBLE_PRECISION
typedef double ir_Float;
#else
typedef float ir_Float;
#endif

/* errors */
typedef enum {
  IR_OK = 0,
  IR_NO_MEMORY,
  IR_NONEXISTANT,
  IR_EMPTY,
  IR_EMPTY_INPUT,
  IR_TYPE_UNSUPPORTED,
  IR_TYPE_DIFFERS,
  IR_DIMENSIONALITY_UNSUPPORTED,
  IR_DIMENSIONALITY_DIFFERS,
  IR_SIZES_DIFFER,
  IR_BAD_DIMENSION,
  IR_IDENTICAL_INPUTS,
  IR_BAD_PARAMETER,
  IR_FILTER_NOT_ODD,
  IR_FILTER_TOO_SMALL,
  IR_READ_ERROR,
  IR_WRITE_ERROR,
  IR_UNKNOWN
} ir_Error;

/* error testing macros */
#define IRTJ(a, e) if (a) {error = e; goto ir_error;}
#define IRXJ(a) if ((error = a) != IR_OK) {goto ir_error;}

/* supported types */
typedef enum {
  IR_TYPE_IMAGE = 1,    /* float data   */
  IR_TYPE_SPECTRUM      /* complex data */
} ir_ImageTypes;

/* image type */
typedef struct ir_image {
  int type;                  /* image type                  */
  int d;                     /* dimension                   */
  int realX;                 /* size of the real image in X */
  int *sizes;                /* dimensions                  */
  ir_Float *data;            /* image data                  */
} ir_Image;

/* checking */
ir_Error ir_CheckTwo(ir_Image *, ir_Image *);
ir_Error ir_Check(ir_Image *);

/* creation, clone, changesize, delete */
ir_Error ir_Init(ir_Image **);
ir_Error ir_Destroy(ir_Image **);
ir_Error ir_Clone(ir_Image *, ir_Image *);
ir_Error ir_ImageChange(ir_Image *, int, int *, int);
ir_Error ir_ImageLink(ir_Image *, ir_Float *, int, int *, int );
ir_Error ir_ImageUnlink(ir_Image **);

/* get image sizes */
ir_Error ir_ImageType(ir_Image *, int *);
ir_Error ir_ImageDimensionality(ir_Image *, int *);
ir_Error ir_ImageDimensions(ir_Image *, int **, int *);
ir_Error ir_ImageSize(ir_Image *, int *);
ir_Error ir_ImageData(ir_Image *, ir_Float **);
ir_Error ir_ImageRealX(ir_Image *, int *);
ir_Error ir_ImageRealSize(ir_Image *, int *);

/* manipulation */
ir_Error ir_Copy(ir_Image *, ir_Image *);
ir_Error ir_WrapOrigin(ir_Image *, ir_Image *);
ir_Error ir_RealPart(ir_Image *, ir_Image *);

/* statistics */
ir_Error ir_Sum(ir_Image *, double *);
ir_Error ir_Mean(ir_Image *, double *);
ir_Error ir_InProduct(ir_Image *, ir_Image *, double *);
ir_Error ir_TriProduct(ir_Image *, ir_Image *, ir_Image *, double *);
ir_Error ir_Ulnv(ir_Image *, ir_Image *, double *);

/* arithmetic */
ir_Error ir_Add(ir_Image *, ir_Image *, ir_Image *);
ir_Error ir_Sub(ir_Image *, ir_Image *, ir_Image *);
ir_Error ir_Mul(ir_Image *, ir_Image *, ir_Image *);
ir_Error ir_Div(ir_Image *, ir_Image *, ir_Image *);
ir_Error ir_AddConst(ir_Image *, ir_Image *, double, double);
ir_Error ir_SubConst(ir_Image *, ir_Image *, double, double);
ir_Error ir_MulConst(ir_Image *, ir_Image *, double, double);
ir_Error ir_DivConst(ir_Image *, ir_Image *, double, double);
ir_Error ir_MulConjugate(ir_Image *, ir_Image *, ir_Image *);
ir_Error ir_WeightedMul(ir_Image *, ir_Image *, ir_Image *, double);
ir_Error ir_WeightedAdd(ir_Image *, ir_Image *, ir_Image *, double);  
ir_Error ir_Abs(ir_Image *, ir_Image *);
ir_Error ir_Clip(ir_Image *, ir_Image *);
ir_Error ir_Square(ir_Image *, ir_Image *);
ir_Error ir_Sqrt(ir_Image *, ir_Image *);
ir_Error ir_Ln(ir_Image *, ir_Image *);

/* real fourier transform */
ir_Error ir_ForwardRft(ir_Image *, ir_Image *);
ir_Error ir_InverseRft(ir_Image *, ir_Image *);

/* filtering */
ir_Error ir_Laplace(ir_Image *, ir_Image *, double *);

/***************************************************************************/
/* restoration routines */
/***************************************************************************/

ir_Error ir_PrintStatistics(int, double, double, double, double, double);

/***************************************************************************/
ir_Error ir_TikhonovCLSParameter(ir_Image *, ir_Image *,
				 double, double, double , double *);

ir_Error ir_TikhonovGCVParameter(ir_Image *, ir_Image *,
				 double, double, double *);

ir_Error ir_TikhonovSNRParameter(ir_Image *, double, double *);

ir_Error ir_Tikhonov(ir_Image *, ir_Image *, ir_Image *, double);

ir_Error ir_EM(ir_Image *, ir_Image *, ir_Image *, double, double, int, 
	       double);
ir_Error ir_Ictm(ir_Image *, ir_Image *, ir_Image *, double, int, double);
ir_Error ir_MapGG(ir_Image *, ir_Image *, ir_Image *, double, int, double);
ir_Error ir_MapGE(ir_Image *, ir_Image *, ir_Image *, double, int, double);
ir_Error ir_MapGR(ir_Image *, ir_Image *, ir_Image *, double, double *,
		  int, double);
ir_Error ir_MapPG(ir_Image *, ir_Image *, ir_Image *, double, double, int, 
		  double);
ir_Error ir_MapPE(ir_Image *, ir_Image *, ir_Image *, double, double, int, 
		  double);
ir_Error ir_MapPR(ir_Image *, ir_Image *, ir_Image *, double, double, 
		  double *, int, double);

/***************************************************************************/
/* Functions for PSF generation for fluorescence imaging                   */
/***************************************************************************/

ir_Error ir_WidefieldPSF(ir_Image *, double, double, double, double,
			 double,double, double, int,double,double,int);
ir_Error ir_WidefieldASF(ir_Image *,ir_Image *,ir_Image *, double, double, double, double,
			 double,double, double, int, int, int);
ir_Error ir_ConfocalPSF(ir_Image *, double, double, double, double, 
			double, double, double, double, double, int,int, int,int,double,double,double,double,int);

#endif


