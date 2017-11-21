#include <math.h>
#include "mex.h"

/* Mean Field Inference with Complex Numbers on General Graphs
 *
 * Usage:
 * [nodeBel niter edgeBel F] = MF_GeneralComplex_C(logEdgePot,nodePot,maxIter,
        optTol,starEdge_V,starEdge_E);
 *
 * Input:
 * nodePot(n,k) - Potential at node n for state k
 * EdgePot(k1,k2,e) - Potential on edge e for states k1 and k2
 * maxIter - Maximum number of iterations
 * opTol - Optimality Tolerance
 * starEdge_V - Vertex vector in star edge representation
 * starEdge_E - Edge vector in star edge representation
 *
 * Output:
 * nodeBel(n,k) - Belief at node i for state k
 * niter - Number of iterations
 * edgeBel(k1,k1,e) - Pairwise Beliefs
 * F - Gibbs Mean Field Free Energy
 *
 */

int DEBUG = 0;

/* Function Declarations */
void checkInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs!=6)
        mexErrMsgTxt("MF_General_C requires SIX Inputs");
    if (mxIsChar(prhs[0])||mxIsClass(prhs[0], "sparse")||!mxIsComplex(prhs[0])||mxIsChar(prhs[1])||mxIsClass(prhs[1], "sparse")||!mxIsComplex(prhs[1]))
        mexErrMsgTxt("Inputs must be COMPLEX, full, and nonstring");
    if (mxGetNumberOfDimensions(prhs[0])!=3)
        mexErrMsgTxt("Edge evidence must be a three dimensional array");
    if (mxGetNumberOfDimensions(prhs[1])!=2)
        mexErrMsgTxt("local evidence must be a two dimensional array");
    if (!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32"))
        mexErrMsgTxt("edge numberings must be int32");
}


void printVectorInt(int* X,int nElem)
{
    int i;
    printf("< ");
    for(i = 0; i < nElem; i++) {
        printf("%d ",X[i]);}
    printf(">\n");
}
void printVectorDouble(double* X,int nElem)
{
    int i;
    printf("< ");
    for(i = 0; i < nElem; i++) {
        printf("%f ",X[i]);}
    printf(">\n");
}

void printVectorDoubleComplex(double* X,double* XI,int nElem)
{
    int i;
    printf("< ");
    for(i = 0; i < nElem; i++) {
        printf("%f+i%f ",X[i],XI[i]);}
    printf(">\n");
}

void printMatrixDouble(double* X,int nRows, int nCols)
{
    int i,j;
    
    for(i = 0; i < nRows; i++) {
        printf("< ");
        for(j = 0; j < nCols; j++) {
            printf("%f ",X[i+nRows*j]);}
        printf(">\n");}
}

void printMatrixDoubleComplex(double* X,double* XI,int nRows, int nCols)
{
    int i,j;
    
    for(i = 0; i < nRows; i++) {
        printf("< ");
        for(j = 0; j < nCols; j++) {
            printf("%f+%fi ",X[i+nRows*j],XI[i+nRows*j]);}
        printf(">\n");}
}

void print2DArrayDouble(double** X,int nRows,int nCols)
{
    int i,j;
    
    for(i=0;i < nRows;i++)
    {
        printf("< ");
        for(j=0;j<nCols;j++)
        {
            printf("%f ",X[i][j]);
        }
        printf(">\n");
    }
}

void print2DArrayDoubleComplex(double** X,double** XI,int nRows,int nCols)
{
    int i,j;
    
    for(i=0;i < nRows;i++)
    {
        printf("< ");
        for(j=0;j<nCols;j++)
        {
            printf("%f+%fi ",X[i][j],XI[i][j]);
        }
        printf(">\n");
    }
}

void print3DMatrixDouble(double*X,int nRows,int nCols,int height)
{
    
    int i,j,k;
    
    for(k = 0; k < height; k++)
    {
        printf("(%d,:,:) = \n",k);
        for(i = 0; i < nRows; i++)
        {
            printf("< ");
            for(j = 0; j < nCols; j++)
            {
                printf("%f ",X[i+nRows*(j+nCols*k)]);
            }
            printf(">\n");
        }
    }
}

int edgeNum(int i,int j,int* starV,int* starE)
{
    int E;
    
    if(DEBUG)printf("Looking for edge between %d and %d\n",i,j);
    
    for(E=starV[i];E<starV[i+1];E++)
    {
        if(DEBUG)printf("Checking %d: Value = %d vs. %d\n",E,starE[E],j);
        if(j==starE[E])
            return E;
    }
    mexErrMsgTxt("Could not find Edge Number!");
    return -1;
}

double absoluteDif(double n1,double n2)
{
    if (n1 > n2)
        return n1-n2;
    else return n2-n1;
}

/* return the real part of the multiplication (x+yi)*(u+vi) */
double multiplyR(double x, double y, double u, double v)
{
    return x*u - y*v;
}

/* return the complex part of the multiplication (x+yi)*(u+vi) */
double multiplyI(double x, double y, double u, double v)
{
    return x*v + y*u;
}

/* return the real part of the division (x+yi)/(u+vi) */
double divideR(double x, double y, double u, double v)
{
    return (x*u + y*v)/(u*u + v*v);
}

/* return the complex part of the division (x+yi)/(u+vi) */
double divideI(double x, double y, double u, double v)
{
    /*printf("[x,y,u,v] = [%.3f, %.3f, %.3f, %.3f], result = %.3f\n",x,y,u,v,(-x*v+y*u)/(u*u+v*v));*/
    return (-x*v+y*u)/(u*u + v*v);
}

/* return the absolute value (squared) of (x+yi) */
double absSquared(double x, double y)
{
    return x*x+y*y;
}

double logR(double x, double y)
{
    return log(sqrt(x*x+y*y));
}

double logI(double x, double y)
{
    return atan2(y,x);
}

double expR(double x, double y)
{
    return exp(x)*cos(y);
}

double expI(double x, double y)
{
    return exp(x)*sin(y);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Declare Variables */
    
    int
    i,j,k_i,k_j,n,nbrsInd,iter,E,
    nNodes,nEdges,nStates,maxIter,converged,nNbrs,
    lhs1_dims[2],lhs2_dims[2],lhs3_dims[3],
    *starV,*starE,
    *nbrs;
    
    /* Real Coefficent Variables */
    double
    optTol,sum,b,temp,temp2,temp3,U1,U2,S1,
    *nodePot,*edgePot,
    *nodeBel,*edgeBel,*nIter,*F,
    **old_bel,**new_bel;
    
    /* Complex Coefficient Variables */
    double
    sumI,bI,tempI,temp2I,temp3I,U1I,U2I,S1I,
    *nodePotI,*edgePotI,
    *nodeBelI,*edgeBelI,
    **old_belI,**new_belI,
    *FI;
    
    /* Check Input */
    
    checkInput(nlhs,plhs,nrhs,prhs);
    
    
    /* Get Input Pointers */
    
    nodePot = mxGetPr(prhs[1]);
    edgePot = mxGetPr(prhs[0]);
    maxIter = mxGetPr(prhs[2])[0];
    optTol = mxGetPr(prhs[3])[0];
    starV = mxGetPr(prhs[4]);
    starE = mxGetPr(prhs[5]);
    nNodes = mxGetDimensions(prhs[1])[0];
    nEdges = mxGetDimensions(prhs[0])[2];
    nStates = mxGetDimensions(prhs[1])[1];
    
    nodePotI = mxGetPi(prhs[1]);
    edgePotI = mxGetPi(prhs[0]);
    
    /* Avoid indexing confusion by decrementing the Star arrays */
    
    for(i = 0; i <= nNodes;i++)
        starV[i]=starV[i]-1;
    for(j = 0; j < nEdges;j++)
        starE[j]=starE[j]-1;
    
    
    
    /* Set-up Output Arrays */
    
    lhs1_dims[0] = nNodes;
    lhs1_dims[1] = nStates;
    lhs2_dims[0] = 1;
    lhs2_dims[1] = 1;
    lhs3_dims[0] = nStates;
    lhs3_dims[1] = nStates;
    lhs3_dims[2] = nEdges;
    plhs[0] = mxCreateNumericArray(2,lhs1_dims,mxDOUBLE_CLASS,mxCOMPLEX);
    plhs[1] = mxCreateNumericArray(2,lhs2_dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericArray(3,lhs3_dims,mxDOUBLE_CLASS,mxCOMPLEX);
    plhs[3] = mxCreateNumericArray(2,lhs2_dims,mxDOUBLE_CLASS,mxCOMPLEX);
    nodeBel = mxGetPr(plhs[0]);
    nodeBelI = mxGetPi(plhs[0]);
    nIter =   mxGetPr(plhs[1]);
    edgeBel = mxGetPr(plhs[2]);
    edgeBelI= mxGetPi(plhs[2]);
    F       = mxGetPr(plhs[3]);
    FI      = mxGetPi(plhs[3]);
    
    
    
    /* Allocate Memory for Auxiliary Arrays */
    
    old_bel = mxCalloc(nNodes,sizeof(double*));
    new_bel = mxCalloc(nNodes,sizeof(double*));
    old_belI = mxCalloc(nNodes,sizeof(double*));
    new_belI = mxCalloc(nNodes,sizeof(double*));
    for(i=0;i < nNodes;i++) {
        old_bel[i] = mxCalloc(nStates,sizeof(double));
        new_bel[i] = mxCalloc(nStates,sizeof(double));
        old_belI[i] = mxCalloc(nStates,sizeof(double));
        new_belI[i] = mxCalloc(nStates,sizeof(double));
    }
    
    
    /* Initialize Variables */
    
    for(i = 0; i < nNodes; i++)
    {
        sum = 0;
        sumI= 0;
        
    /* Set Initial Beliefs to Potentials */
        for(j = 0; j < nStates; j++)
        {
            old_bel[i][j]       = nodePot[i+nNodes*j];
            old_belI[i][j]       = nodePotI[i+nNodes*j];
            sum = sum+old_bel[i][j];
            sumI= sumI+old_belI[i][j];
        }
        
    /* Normalized Initial Beliefs */
        for(j = 0; j < nStates; j++)
        {
            temp = divideR(old_bel[i][j],old_belI[i][j],sum,sumI);
            tempI= divideI(old_bel[i][j],old_belI[i][j],sum,sumI);
            old_bel[i][j]       = temp;
            new_bel[i][j]       = temp;
            old_belI[i][j]       = tempI;
            new_belI[i][j]       = tempI;
        }
        
    }
    
    
    iter = 0;
    converged = 0;
    
    /* Now do the work... */
    
    while(!converged && iter < maxIter)
    {
        if(DEBUG)printf("Iteration %.0f\n",*nIter);
        
        for(i = 0; i < nNodes; i++)
        {
            if(DEBUG)printf("Working on Node %d of %d\n",i,nNodes);
            
            nbrsInd = starV[i];
            nbrs = &starE[starV[i]];
            nNbrs = starV[i+1]-starV[i];
            
            sum = 0;
            sumI= 0;
            for(k_i = 0; k_i < nStates; k_i++)
            {
                if(DEBUG)printf("Processing State %d of %d\n",k_i,nStates);
                b = 0;
                bI= 0;
                for(n = 0; n < nNbrs; n++)
                {
                    j = nbrs[n];
                    E = nbrsInd+n;
                    if(DEBUG)
                    {
                        printf("Looking at neighbor %d of %d: %d\n",j,n,nNbrs);
                        printf("Edge Index = %d\n",nbrsInd+n);
                    }
                    
                    for(k_j = 0; k_j < nStates; k_j++)
                    {
                        if(DEBUG)printf("Processing Neighbor State %d of %d\n",k_j,nStates);
                        temp = logR(edgePot[k_i+nStates*(k_j+nStates*E)],edgePotI[k_i+nStates*(k_j+nStates*E)]);
                        tempI= logI(edgePot[k_i+nStates*(k_j+nStates*E)],edgePotI[k_i+nStates*(k_j+nStates*E)]);
                        temp2 = multiplyR(new_bel[j][k_j],new_belI[j][k_j],temp,tempI);
                        temp2I = multiplyI(new_bel[j][k_j],new_belI[j][k_j],temp,tempI);
                        
                        b = b + temp2;
                        bI= bI+ temp2I;
                        
                        if(DEBUG)printf("NB = %.3f, log(EP) = %.3f, B = %.3f\n",new_bel[j][k_j],temp,b);
                        
                    }
                }
                
                if(DEBUG)printf("b = %.3f, bI = %.3f\n",b,bI);
                temp = expR(b,bI);
                tempI= expI(b,bI);
                if(DEBUG)printf("e(B) = %.3f\n",temp);
                new_bel[i][k_i] =  multiplyR(nodePot[i+k_i*nNodes],nodePotI[i+k_i*nNodes],temp,tempI);
                new_belI[i][k_i] = multiplyI(nodePot[i+k_i*nNodes],nodePotI[i+k_i*nNodes],temp,tempI);
                if(DEBUG)printf("nb = %.3f\n",new_bel[i][k_i]);
                sum = sum+new_bel[i][k_i];
                sumI= sumI+new_belI[i][k_i];
                
            }
            
            if(DEBUG)printf("NOrmalizing, sum = %.3f + %.3fi\n",sum,sumI);
            /* Normalize */
            for(k_i=0;k_i<nStates;k_i++) {
                temp =  divideR(new_bel[i][k_i],new_belI[i][k_i],sum,sumI);
                tempI = divideI(new_bel[i][k_i],new_belI[i][k_i],sum,sumI);
                new_bel[i][k_i] = temp;
                new_belI[i][k_i]= tempI;
            }
            if(DEBUG)printf("NB = (%.3f, %.3f)\n",new_bel[i][0],new_bel[i][1]);
            
        }
        
        converged = 1;
        sum = 0;
        sumI= 0;
        for(i = 0; i < nNodes; i++) {
            for(j = 0; j < nStates; j++) {
                if(DEBUG)printf("NB = %f, OB = %f, Err = %f\n",new_bel[i][j],old_bel[i][j],absoluteDif( new_bel[i][j],old_bel[i][j]));
                
                sum = sum + absoluteDif(new_bel[i][j],old_bel[i][j]);
                sumI= sumI+ absoluteDif(new_belI[i][j],old_belI[i][j]);
                
                if (absSquared(sum,sumI) >= optTol*optTol)
                    converged = 0;
            }
        }
        
        if(DEBUG)printf("Error = %.3f\n",sum);
        
        
        if(DEBUG)printf("Converged = %d\n",converged);
        
        iter++;
        
        for(i = 0; i < nNodes;i++) {
            for(j = 0; j < nStates; j++) {
                old_bel[i][j] = new_bel[i][j];
                old_belI[i][j] = new_belI[i][j];}}
        
        
        
    }
    
    
    if(converged==1)
	{
		if(DEBUG)printf("MF: Converged in %d iterations\n",iter);
	}
    else
		printf("MF: Stopped after maxIter = %d iterations\n",iter);
    
    /* Assign Output Args */
    
    for(i = 0; i < nNodes; i++) {for(j = 0; j < nStates; j++) {
        nodeBel[i + nNodes*j] = new_bel[i][j];
        nodeBelI[i + nNodes*j] = new_belI[i][j];}}
    nIter[0] = (double)iter;
    
    /* Compute Edge Bels */
    
    for(i = 0; i < nNodes;i++) {
        nbrsInd = starV[i];
        nbrs = &starE[starV[i]];
        nNbrs = starV[i+1]-starV[i];
        for(n = 0; n < nNbrs; n++)
        {
            j = nbrs[n];
            E = nbrsInd+n;
            
            for(k_i = 0; k_i < nStates; k_i++)
            {
                for(k_j = 0; k_j < nStates; k_j++)
                {
                    temp = multiplyR(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i],nodeBel[j+nNodes*k_j],nodeBelI[j+nNodes*k_j]);
                    tempI = multiplyI(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i],nodeBel[j+nNodes*k_j],nodeBelI[j+nNodes*k_j]);
                    edgeBel[k_i + nStates*(k_j + nStates*E)] = temp;
                    edgeBelI[k_i+ nStates*(k_j + nStates*E)] = tempI;
                }
            }
        }
        
    }
    
    /* Compute Gibbs Mean Field Free Energy */
    
    U1=0;U1I=0;S1=0;S1I=0;U2=0;U2I=0;
    for(i = 0; i < nNodes;i++) {
        
        /* Local Mean-Field Average Energy Term */
        
        for(k_i = 0; k_i < nStates; k_i++) {
            temp = logR(nodePot[i+nNodes*k_i],nodePotI[i+nNodes*k_i]);
            tempI= logI(nodePot[i+nNodes*k_i],nodePotI[i+nNodes*k_i]);
            temp2= multiplyR(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i],temp,tempI);
            temp2I= multiplyI(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i],temp,tempI);
            U1 = U1 +temp2;
            U1I= U1I+temp2I;
        }
        
        /* Mean-Field Entropy Term */
        for(k_i = 0; k_i < nStates; k_i++) {
            if(nodeBel[i+nNodes*k_i] + nodeBelI[i+nNodes*k_i] >= 0.000001) {
                temp = logR(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i]);
                tempI= logI(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i]);
                temp2= multiplyR(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i],temp,tempI);
                temp2I= multiplyI(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i],temp,tempI);
                S1 = S1 + temp2;
                S1I= S1I+ temp2I;
            }
        }
        
        /* For all neighbors */
        
        nbrsInd = starV[i];
        nbrs = &starE[starV[i]];
        nNbrs = starV[i+1]-starV[i];
        for(n = 0; n < nNbrs; n++)
        {
            j = nbrs[n];
            
            if (i > j)
            {
                E = nbrsInd+n;
                for(k_i = 0; k_i < nStates; k_i++)
                {
                    for(k_j = 0; k_j < nStates; k_j++)
                    {
                        temp = multiplyR(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i],nodeBel[j+nNodes*k_j],nodeBelI[j+nNodes*k_j]);
                        tempI = multiplyI(nodeBel[i+nNodes*k_i],nodeBelI[i+nNodes*k_i],nodeBel[j+nNodes*k_j],nodeBelI[j+nNodes*k_j]);
                        temp2 = logR(edgePot[k_i+nStates*(k_j+nStates*E)],edgePotI[k_i+nStates*(k_j+nStates*E)]);
                        temp2I= logI(edgePot[k_i+nStates*(k_j+nStates*E)],edgePotI[k_i+nStates*(k_j+nStates*E)]);
                        temp3 = multiplyR(temp,tempI,temp2,temp2I);
                        temp3I = multiplyI(temp,tempI,temp2,temp2I);
                        
                        U2 = U2+temp3;
                        U2I= U2I+temp3I;
                        
                    }
                }
            }
        }
        
        
    }
    
    *F = -U2-U1+S1;
    *FI= -U2I-U1I+S1I;
    if(DEBUG)printf("U1 = %.3f, U2 = %.3f, S1 = %.3f, F = %.3f\n",U1,U2,S1,*F);
    
    
    
    
    /* Free Memory */
    
    for(i=0;i < nNodes;i++) {
        mxFree(old_bel[i]);
        mxFree(new_bel[i]);
        mxFree(old_belI[i]);
        mxFree(new_belI[i]);
    }
    
    mxFree(old_bel);
    mxFree(new_bel);
    mxFree(old_belI);
    mxFree(new_belI);
    
}

