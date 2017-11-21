/*

  * This is a mex interface to compute the expanded feature
  * vectors.  This function is called in the learning and
  * inference stages by the objective functions and getMu.m

  */
#include <math.h>
#include "mex.h"

/* the input arguments */
#define FEA_IN            prhs[0]

/* the output array, this is the expanded feature vector*/
#define FEA_OUT        plhs[0]


void mexFunction(int nlhs, mxArray *plhs[ ], 
			int nrhs, const mxArray *prhs[ ])
{

	int nInFea;   /* the number of input features */
	int i,j,k;       /* counting index */
	int nExpFea;   /* number of outputs */
	double normSum;  /* used for normalizing */

	double *expFea;
	double *inFea;

	/* Check for sufficient input arguments */
	if (nrhs != 1) 
	{ 
		mexErrMsgTxt("Three input arguments required."); 
    } 

	nInFea = mxGetM(FEA_IN); 
	nExpFea = 1 + nInFea + nInFea*(nInFea+1)/2;

	/* Create a array for the return argument */ 
    FEA_OUT = mxCreateDoubleMatrix(nExpFea, 1, mxREAL); 
    
	/* get the pointers to work with */
	expFea = mxGetPr(FEA_OUT);
    
	inFea = mxGetPr(FEA_IN); 

	k=0;
	expFea[k] = 1.0;
	normSum = 1.0;
	for( k=1; k<nInFea+1; k++)
	{
		expFea[k] = inFea[k-1];
		normSum = normSum + expFea[k]*expFea[k];
	}
	for( i=0; i<nInFea; i++)
	{
		for(j=i; j<nInFea; j++)
		{
			expFea[k] = inFea[i]*inFea[j];
			normSum = normSum +expFea[k]*expFea[k];
			k++;
		}
	}
/*
	normSum = sqrt( normSum );
	
	for( k=0; k<nExpFea; k++)
		expFea[k] = expFea[k]/normSum;

	expFea[0]=1.0;
	*/

	
}
