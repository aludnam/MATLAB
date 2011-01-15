/* Compares a number of simulated Gaussian Spots with experimental data

 * Examle : a=exp(-((xx(20,20)-2).^2+(yy(20,20)-1.2).^2)/20)
           [fitted,params]=FitDataNDFast([0 1 0; 1 0 0],a)
 */

#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

double dosim(double * params, int paramSize, int numparams, double * gparams, int * pos) // int numgprams, 
{
    double result=gparams[0],sq,ssq;
    //printf("PosX : %d, y %d\n",pos[0],pos[1]);
    for (int n=0;n<numparams;n++)
    {
        ssq=0;
        for (int d=0;d<paramSize-1;d++)
        {sq=(pos[d]-params[paramSize*n+d+1]);
        //printf("Param %d %d: %g\n",d,n,params[paramSize*n+d+1]);
         ssq +=sq*sq;}
        result += params[paramSize*n]*exp(-ssq/gparams[1]);
        //printf("Intensity : %g\n",params[paramSize*n]);
    }
    //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
    return result;
}

// this function computes the meas square error comparing data with simulation
// params is the array of spot positions
// paramsg are the global parameters
double do_mse(double * mydata, int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double *resdiff)  // int numgparams, 
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    //for (int y=0;y<numparams;y++) 
    //   for (int x=0;x<paramsize;x++) 
    //      printf("MSE: Param %d, %d : %g\n",x,y,params[x+y*paramsize]);
    double result=0,tmp;
    int i=0,pos[3];
    for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
        for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
            for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
                {
                tmp = dosim(params, paramsize, numparams, gparams, pos); // numgparams, 
                if (res != 0)
                    res[i] = tmp; // save simulation
                tmp=mydata[i]-tmp;
                if (resdiff != 0)
                    resdiff[i] = tmp; // save difference                
                result += tmp*tmp;
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
                i++;
            }
    return result;
}

double do_idiv(double * mydata, int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double * resdiff)  // int numgparams, 
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    double result=0,tmp;
    int i=0,pos[3];
    for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
        for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
            for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
                {
                tmp = dosim(params, paramsize, numparams, gparams, pos); // numgparams, 
                if (res != 0)
                    res[i] = tmp; // save idiv image
                if (mydata[i] !=0)
                    tmp=mydata[i]*log(mydata[i]/tmp)-(mydata[i]-tmp);  // i-divergence with sterling's approximation
                else
                    tmp=-(mydata[i]-tmp);
                if (resdiff != 0)
                    resdiff[i] = tmp; // save difference                
                result += tmp;
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
                i++;
            }
    return result;
}

double do_fidiv(double * mydata, int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double * resdiff)  // int numgparams, 
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    double result=0,tmp;
    int i=0,pos[3];
    for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
        for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
            for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
                {
                tmp = dosim(params, paramsize, numparams, gparams, pos); // numgparams, 
                if (res != 0)
                    res[i] = tmp; // save idiv image
                tmp=tmp-mydata[i]*log(tmp);  // fast version omitting constants
                if (resdiff != 0)
                    resdiff[i] = tmp; // save difference                
                result += tmp;
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
                i++;
            }
    return result;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* format:
   writekhoros_info('blabla',dims,type)
   MultiGaussMSE(Data)  to prepare with experimental data
   MultiGaussMSE(Matrix with parameters, Global parameters) each row containing intensity and positions */

/* nlhs # number left handed parameters
   plhs # left handed parameter array
   nrhs # number right handed parameters
   prhs # right handed parameter array*/
    double *data=0;  // Let's hope this one survives the individual calls
    double * params, * gparams,* fixedparams, result, * res=0, * resdiff=0;
    char * input_buf;
    int PsizeX,PsizeY, buflen; // GPsize, sizes[5],
    static double * mydata=0;
    static double * cparams=0;
    static int allocatedCparams=0;
    static int sizes[100],nd=1,totalsize=1,numdims;   // To hell with ppl who use more than 100 dimensions!
    enum method {mse, idiv,fidiv};
    static method mymethod=mse;
    
    if (nrhs != 1 && nrhs != 3) 
    {
        printf ("Nr of parameters: %d\n",nrhs);
         mexErrMsgTxt("1 or 3 inputs required");
    }
    //if (mxIsChar(prhs[0]) != 1)
    //     mexErrMsgTxt("Input must be a string.");


    if (nrhs > 1)  // archieve the data by saving the pointer and remember which method to use
    {
    if (nrhs < 3)  
            mexErrMsgTxt("When submitting data, three arguments are required: data, method-string, dimensions!");
        
    /* Get the length of the input string. */
        const int *sz;
        sz = mxGetDimensions(prhs[0]);
        nd = mxGetNumberOfDimensions(prhs[0]);
        totalsize=1;
        for (int ii=0;ii<100;ii++)
        {
            sizes[ii]=1;
        }

        for (int ii=0;ii<nd;ii++)
        {
           sizes[ii]=sz[ii];  // save these values
           // printf("dim %d, size %d\n",ii,sizes[ii]);
           totalsize *= sz[ii];
        }
        data = mxGetPr(prhs[0]);
        if (mydata != 0) free(mydata);
        mydata=(double *) calloc(totalsize,sizeof(double));
        for (int i=0;i<totalsize;i++) 
            mydata[i]=data[i];
        if (nlhs != 0)
            mexErrMsgTxt("When submitting data, no output is returned!");
        // printf("Bufflen is %d x %d, pointer is %x, copied to %x\n",dataSizeX,dataSizeY,data,mydata);
    /* Get the length of the input string. */
    buflen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
    
    /* Allocate memory for input and output strings. */
    input_buf = (char*) mxCalloc(buflen, sizeof(char));
    
    /* Copy the string data from prhs[0] into a C string */
    mxGetString(prhs[1], input_buf, buflen);

    if (strcmp(input_buf,"mse") == 0)
        mymethod = mse;
    else if (strcmp(input_buf,"idiv") == 0)
        mymethod = idiv;
    else if (strcmp(input_buf,"fidiv") == 0)
        mymethod = fidiv;
    else
    {
        printf("Requested method was %s\n",input_buf);
        mexErrMsgTxt("Invalid method. Valid methods are : 'mse', 'idiv' and 'fidiv'");
    }
    /* Copy the string data from prhs[0] into a C string */
    numdims=(int) (* mxGetPr(prhs[2]));

    mxFree(input_buf);
    }
    else
    {
        PsizeX = mxGetM(prhs[0]);
        if (mxGetN(prhs[0]) > 1)
            mexErrMsgTxt("All parameters should be in a single vector.");
        if ((PsizeX-2) % (numdims +1)!=0)
            mexErrMsgTxt("Number of parameters (2 globals and rest) does not match with number of dimensions to fit.");
        else
        {
            PsizeX = numdims+1;
            PsizeY = (mxGetM(prhs[0])-2)/numdims;
        }
        gparams = mxGetPr(prhs[0]);
        params = gparams+2;  // omit the two globals
    	if (allocatedCparams< PsizeY)
        {
            if (cparams != 0)
                free(cparams);
            cparams=(double *) calloc(PsizeY+100,sizeof(double));
            allocatedCparams=PsizeY+100;
        }

        //if (mxGetM(prhs[1]) != 1)
        //    mexErrMsgTxt("Global parameters must be a row vector.");
    

        //GPsize = mxGetN(prhs[1]);
        //fixedparams = mxGetPr(prhs[1]);
        //printf("ParamsizeX: %d, Y %d\n",PsizeX,PsizeY);
        //printf("dataSizeX: %d, Y %d\n",dataSizeX,dataSizeY);
        if (mydata != 0) {
            //for (int i=0;i<dataSizeX*dataSizeY;i++) 
            //    printf("%d: %g\n",i,mydata[i]);
            //for (int y=0;y<PsizeY;y++) 
            //  for (int x=0;x<PsizeX;x++) 
            //    printf("Param %d, %d : %g\n",x,y,params[x+y*PsizeX]);

            for (int d=0;d<PsizeY;d++)   // to account for center is zero
            for (int n=0;n<PsizeX;n++)   // to account for center is zero
                if (n>0)
                    cparams[d*PsizeX+n] = params[d*PsizeX+n] + floor((double) sizes[n-1]/2);  
                else
                    cparams[d*PsizeX+n] = params[d*PsizeX+n];  

            if (nlhs >= 2)
            {
                plhs[1] = mxCreateNumericArray(nd, sizes, mxDOUBLE_CLASS, mxREAL);
                //plhs[1] = mxCreateDoubleMatrix(sizes,nd, mxREAL);
                res = mxGetPr(plhs[1]);
            }
            if (nlhs >= 3)
            {
                plhs[2] = mxCreateNumericArray(nd, sizes, mxDOUBLE_CLASS, mxREAL);
                //plhs[1] = mxCreateDoubleMatrix(sizes,nd, mxREAL);
                resdiff = mxGetPr(plhs[2]);
            }
            //for (int y=0;y<PsizeY;y++) 
            //  for (int x=0;x<PsizeX;x++) 
            //    printf("Param %d, %d : %g\n",x,y,params[x+y*PsizeX]);

            switch (mymethod)
            {
                case mse:
                    result=do_mse(mydata, sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff); // GPsize,
                break;
                case idiv:
                    result=do_idiv(mydata, sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff); // GPsize,
                break;
                case fidiv:
                    result=do_fidiv(mydata, sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff); // GPsize,
                break;
                default:
                    mexErrMsgTxt("Undefined method. Valid methods are : 'mse' and 'idiv'");
            }
                if (nlhs >= 1)
            {
                plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
                double * dp = mxGetPr(plhs[0]);
                (* dp) = result;
            }
            // printf("Result %g\n",result);
        }
        else
            mexErrMsgTxt("Please provide just the data matrix first");
    }
    
    return;
}
