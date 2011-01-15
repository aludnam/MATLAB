/* get the size for Khoros files 
   Bernd Rieger, Keith Lidke & Rainer Heintzmann
   June 2004 - June 2006
*/

#include <stdio.h>
#include <string.h>
#include "mex.h"

void myexit(int arg) 
{ 
    mexErrMsgTxt("Fatal error in reading khoros file. Bailing out");
}

#include "khoros.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* format:
   [dims, datatype] = readkhoros_info('blabla')*/

/* nlhs # number left handed parameters
   plhs # left handed parameter array
   nrhs # number right handed parameters
   prhs # right handed parameter array*/
    ifstream from;
    char TypeString[50]; // will be filled in
    int ValueDim[6];
    //double dp[5];
    int begindata;
    
    char Var[2];
    string matlabtype;
    char *input_buf;
    int   buflen,status;
  
    if (nrhs != 1) 
         mexErrMsgTxt("One input required.");
    else if (nlhs !=2) 
         mexErrMsgTxt("Need two output arguments.");
    if (mxIsChar(prhs[0]) != 1)
         mexErrMsgTxt("Input must be a string.");
    if (mxGetM(prhs[0]) != 1)
      mexErrMsgTxt("Input must be a row vector.");
    /* Get the length of the input string. */
    buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
    /* Allocate memory for input and output strings. */
    input_buf = (char*) mxCalloc(buflen, sizeof(char));
  
    /* Copy the string data from prhs[0] into a C string 
       input_buf. */
    status = mxGetString(prhs[0], input_buf, buflen);
    if (status != 0) 
       mexWarnMsgTxt("Not enough space. String is truncated.");
     //printf("strlen %d,%d\n",strlen(input_buf),buflen);
      
    from.open(input_buf, ios::binary);
    if (! from.is_open())
       mexErrMsgTxt("Error in readKhoros_info: Could not open file for reading.");
        
    ValueDim[3]=ReadKhorosHeader(&from,TypeString,ValueDim[0],ValueDim[1],ValueDim[2],ValueDim[4]);
    if (ValueDim[3] < 0)
            mexErrMsgTxt("Problem reading Khoros Header (maybe this file is not in khoros format?).");

    ValueDim[5] = (int) from.tellg();
    //from.seekg(0,ios:end);
    //ValueDim[5] = ((int) from.tellg()) - ValueDim[0]*ValueDim[1]*ValueDim[2]*ValueDim[3]*ValueDim[4];
    from.close();
    
    if (strcmp(TypeString,"Unsigned Byte")==0)matlabtype="uint8";
    if (strcmp(TypeString,"Integer")==0)matlabtype="sint32";
    if (strcmp(TypeString,"Long")==0)matlabtype="Doesn't support 64 bit int!";
    if (strcmp(TypeString,"Float")==0)matlabtype="float"; 
    if (strcmp(TypeString,"Double")==0)matlabtype="double"; 
    if (strcmp(TypeString,"Complex")==0)matlabtype="scomplex"; 
    if (strcmp(TypeString,"Unsigned Short")==0)matlabtype="uint16";
    if (strcmp(TypeString,"Double Complex")==0)matlabtype="dcomplex";
   
      //printf("%s\n",TypeString);
    plhs[1] = mxCreateString(matlabtype.c_str());
    plhs[0] = mxCreateDoubleMatrix(6,1, mxREAL);
    double * dp = mxGetPr(plhs[0]);
    for (int ii=0;ii<6;ii++){
      dp[ii] = (double) ValueDim[ii];
      //printf("ValueDim %d\n",ValueDim[ii]);
      //printf("ValueDim %f\n",dp[ii]);
    }
    mxFree(input_buf);
    return;
}
