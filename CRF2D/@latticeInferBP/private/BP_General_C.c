#include <math.h>
#include "mex.h"

/* Sum-Product Belief Propagation on General Graphs
 *
 * Usage:
 * [nodeBel niter nodeMsgs] = BP_General_C(edgePot,nodePot,maximize,maxIter,
        optTol,starEdge_V,starEdge_E);
 *
 * Input:
 * nodePot(n,k) - Potential at node n for state k
 * edgePot(k1,k2,e) - Potential on edge e for states k1 and k2
 * maxmize - Should be 0, since this code only supports sum-product
 * maxIter - Maximum number of iterations
 * opTol - Optimality Tolerance
 * starEdge_V - Vertex vector in star edge representation
 * starEdge_E - Edge vector in star edge representation
 *
 * Output:
 * nodeBel(n,k) - Belief at node i for state k
 * niter - Number of iterations
 * nodeMsgs(e,k) - Message for state k at edge e
 *
 */

int DEBUG = 0;

/* Function Declarations */
void checkInput(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs!=7)
        mexErrMsgTxt("BP_General_C requires SEVEN Inputs");
    if (mxIsChar(prhs[0])||mxIsClass(prhs[0], "sparse")||mxIsComplex(prhs[0])||mxIsChar(prhs[1])||mxIsClass(prhs[1], "sparse")||mxIsComplex(prhs[1]))
        mexErrMsgTxt("Inputs must be real, full, and nonstring");
    if (mxGetNumberOfDimensions(prhs[0])!=3)
        mexErrMsgTxt("Edge evidence must be a three dimensional array");
    if (mxGetNumberOfDimensions(prhs[1])!=2)
        mexErrMsgTxt("local evidence must be a two dimensional array");
    if (!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32"))
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

void printMatrixDouble(double* X,int nRows, int nCols)
{
    int i,j;
    
    for(i = 0; i < nRows; i++) {
        printf("< ");
        for(j = 0; j < nCols; j++) {
            printf("%f ",X[i+nRows*j]);}
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
    return 0;
}

double absoluteDif(double n1,double n2)
{
    if (n1 > n2)
        return n1-n2;
    else return n2-n1;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Declare Variables */
    
    int
    i,j,k,ii,jj,iii,nbrsInd,nbrsInd2,E,E2,iter,n,n2,k_i,k_j,
    nNodes,nEdges,nStates,maximize,maxIter,converged,nNbrs,nNbrs2,
    lhs1_dims[2],lhs2_dims[2],lhs3_dims[2],lhs4_dims[3],
    *starV,*starE,
    *nbrs,*nbrs2;
    
    double
    optTol,sum,temp2,temp3,
    *nodePot,*edgePot,
    *nodeBel,*nIter,*msgs,*edgeBel,
    *pot_ij,*temp,*newm,
    **prod_of_msgs,**old_bel,**new_bel,**old_msg,**new_msg;
    
    
    /* Check Input */
    
    checkInput(nlhs,plhs,nrhs,prhs);
    
    
    /* Get Input Pointers */
    
    nodePot = mxGetPr(prhs[1]);
    edgePot = mxGetPr(prhs[0]);
    maximize = mxGetPr(prhs[2])[0];
    maxIter = mxGetPr(prhs[3])[0];
    optTol = mxGetPr(prhs[4])[0];
    starV = mxGetPr(prhs[5]);
    starE = mxGetPr(prhs[6]);
    nNodes = mxGetDimensions(prhs[1])[0];
    nEdges = mxGetDimensions(prhs[0])[2];
    nStates = mxGetDimensions(prhs[1])[1];
    
    /* Avoid indexing confusion by decrementing the Star arrays */
    
    for(i = 0; i <= nNodes;i++)
        starV[i]=starV[i]-1;
    for(j = 0; j < nEdges;j++)
        starE[j]=starE[j]-1;
    
    /******************* START DEBUG ************************/
    
    /* Output Input to Test if it worked */
    if (DEBUG == 1)
    {
        printf("nNodes = %d\n",nNodes);
        printf("nEdges = %d\n",nEdges);
        printf("nStates = %d\n",nStates);
        printf("Maximize = %d\n",maximize);
        printf("MaxIter = %d\n",maxIter);
        printf("OptTol  = %f\n",optTol);
        printf("NodePot:\n");printMatrixDouble(nodePot,nNodes,nStates);
        printf("EdgePot:\n");print3DMatrixDouble(edgePot,nStates,nStates,nEdges);
        printf("starV:\n");printVectorInt(starV,nNodes);
        printf("starE:\n");printVectorInt(starE,nEdges);
    }
    
    /******************** END DEBUG ***********************/
    
    
    /* Set-up Output Arrays */
    
    lhs1_dims[0] = nNodes;
    lhs1_dims[1] = nStates;
    lhs2_dims[0] = 1;
    lhs2_dims[1] = 1;
    lhs3_dims[0] = nEdges;
    lhs3_dims[1] = nStates;
    lhs4_dims[0] = nStates;
    lhs4_dims[1] = nStates;
    lhs4_dims[2] = nEdges;
    plhs[0] = mxCreateNumericArray(2,lhs1_dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(2,lhs2_dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericArray(2,lhs3_dims,mxDOUBLE_CLASS,mxREAL);
    plhs[3] = mxCreateNumericArray(3,lhs4_dims,mxDOUBLE_CLASS,mxREAL);
    nodeBel = mxGetPr(plhs[0]);
    nIter =   mxGetPr(plhs[1]);
    msgs =    mxGetPr(plhs[2]);
    edgeBel = mxGetPr(plhs[3]);
    
    
    /* Allocate Memory for Auxiliary Arrays */
    
    prod_of_msgs = mxCalloc(nNodes,sizeof(double*));
    old_bel = mxCalloc(nNodes,sizeof(double*));
    new_bel = mxCalloc(nNodes,sizeof(double*));
    old_msg = mxCalloc(nEdges,sizeof(double*));
    new_msg = mxCalloc(nEdges,sizeof(double*));
    newm = mxCalloc(nStates,sizeof(double));
    pot_ij = mxCalloc(nStates*nStates,sizeof(double));
    temp = mxCalloc(nStates,sizeof(double));
    for(i=0;i < nNodes;i++) {
        prod_of_msgs[i] = mxCalloc(nStates,sizeof(double));
        old_bel[i] = mxCalloc(nStates,sizeof(double));
        new_bel[i] = mxCalloc(nStates,sizeof(double));
    }
    for(i=0;i < nEdges;i++) {
        old_msg[i] = mxCalloc(nStates,sizeof(double));
        new_msg[i] = mxCalloc(nStates,sizeof(double));
    }
    
    
    /* Initialize Variables */
    
    for(i = 0; i < nNodes; i++)
    {
        for(j = 0; j < nStates; j++)
        {
            prod_of_msgs[i][j]  = nodePot[i+nNodes*j];
            old_bel[i][j]       = nodePot[i+nNodes*j];
        }
    }
    for(i = 0; i < nEdges; i++) {
        for(j = 0; j < nStates; j++) {
            old_msg[i][j] = 1.0/(double)nStates;}}
    
    iter = 0;
    converged = 0;
    
    
    
    
    /******************* START DEBUG ************************/
    
    /* Test Initalization */
    if (DEBUG == 1)
    {
        printf("prod_of_msgs:\n");
        print2DArrayDouble(prod_of_msgs,nNodes,nStates);
        printf("old_bel:\n");
        print2DArrayDouble(old_bel,nNodes,nStates);
        printf("old_msg:\n");
        print2DArrayDouble(old_msg,nEdges,nStates);
    }
    
    /******************** END DEBUG ***********************/
    
    
    
    
    /* Now do the work... */
    
    while(!converged && iter < maxIter)
    {
        if(DEBUG)printf("Iteration %.0f\n",*nIter);
        
        for(i = 0; i < nNodes; i++)
        {
            nbrsInd = starV[i];
            nbrs = &starE[starV[i]];
            nNbrs = starV[i+1]-starV[i];
            
            if(DEBUG)
                printf("Node %d has %d neighbors, indexed from %d\n",i,nNbrs,starV[i]);
            
            for(j = 0; j < nNbrs; j++)
            {
                if(DEBUG)
                {
                    printf("Processing neighbor %d: node = %d\n",j,nbrs[j]);
                    printf("Edge Index = %d\n",nbrsInd+j);
                }
                if (i < nbrs[j]) {
                    for(ii = 0; ii < nStates; ii++) {
                        for(jj = 0; jj < nStates; jj++) {
                            pot_ij[ii+nStates*jj]=edgePot[ii+nStates*(jj+nStates*(nbrsInd+j))];}}}
                else
                {
                    for(ii = 0; ii < nStates; ii++) {
                        for(jj = 0; jj < nStates; jj++) {
                            pot_ij[jj+nStates*ii]=edgePot[ii+nStates*(jj+nStates*(nbrsInd+j))];}}}
                if(DEBUG)
                {
                    printf("pot_ij for j = %d:\n",j);
                    printMatrixDouble(pot_ij,nStates,nStates);
                }
                for(ii = 0; ii < nStates; ii++) {
                    temp[ii] = nodePot[i+ii*nNodes];}
                
                for(k = 0; k < nNbrs; k++) {
                    if(k==j)
                        continue;
                    if(DEBUG)printf("Multiplying Evidence by Message from neighbor %d\n",nbrs[k]);
                    E = edgeNum(nbrs[k],i,starV,starE);
                    for(ii = 0; ii < nStates; ii++) {
                        temp[ii] = temp[ii]*old_msg[E][ii];}
                    
                }
                if(DEBUG)
                {
                    printf("temp:\n");
                    printVectorDouble(temp,nStates);
                }
                
                if (maximize)
                    mexErrMsgTxt("BP_General_C only supports sum-product\n");
                else
                {
                    
                    for(ii=0;ii<nStates;ii++) {
                        newm[ii]=0;
                        for(iii=0;iii<nStates;iii++) {
                            newm[ii]+=pot_ij[iii+nStates*ii]*temp[iii];
                        }
                    }
                    
                    if(DEBUG)
                    {
                        printf("newm:\n");
                        printVectorDouble(newm,nStates);
                    }
                    
                    sum = 0;
                    for(ii=0;ii<nStates;ii++) {
                        sum+=newm[ii];}
                    
                    for(ii=0;ii<nStates;ii++) {
                        new_msg[nbrsInd+j][ii] = newm[ii]/sum;}
                    
                }
                
                
            }
            
        }
        
        if(DEBUG)
        {
            printf("new_msg:\n");
            print2DArrayDouble(new_msg,nEdges,nStates);
        }
        
        for(i = 0;i < nNodes; i++)
        {
            nbrsInd = starV[i];
            nbrs = &starE[starV[i]];
            nNbrs = starV[i+1]-starV[i];
            
            for(j = 0; j < nStates; j++) {
                prod_of_msgs[i][j]=nodePot[i+j*nNodes];}
            
            
            for(j = 0; j < nNbrs; j++) {
                E = edgeNum(nbrs[j],i,starV,starE);
                for(k = 0; k < nStates;k++) {
                    prod_of_msgs[i][k] = prod_of_msgs[i][k]*new_msg[E][k];
                }
            }
            
            if(DEBUG)
            {
                printf("prod_of_msgs:\n");
                printVectorDouble(prod_of_msgs[i],nStates);
            }
            
            sum = 0;
            for(j = 0; j < nStates; j++) {
                sum+=prod_of_msgs[i][j];}
            
            for(j = 0; j < nStates; j++) {
                new_bel[i][j] = prod_of_msgs[i][j]/sum;}
        }
        
        if(DEBUG)
        {
            printf("new_bel:\n");
            print2DArrayDouble(new_bel,nNodes,nStates);
        }
        
        converged = 1;
        for(i = 0; i < nNodes; i++) {
            for(j = 0; j < nStates; j++) {
                if(DEBUG)printf("NB = %f, OB = %f, Err = %f\n",new_bel[i][j],old_bel[i][j],absoluteDif( new_bel[i][j],old_bel[i][j]));
                if(absoluteDif( new_bel[i][j],old_bel[i][j]) > optTol)
                    converged = 0;
            }
        }
        
        if(DEBUG)printf("Converged = %d\n",converged);
        
        iter++;
        
        for(i = 0; i < nNodes;i++) {
            for(j = 0; j < nStates; j++) {
                old_bel[i][j] = new_bel[i][j];}}
        for(i = 0; i < nEdges;i++) {
            for(j = 0; j < nStates; j++) {
                old_msg[i][j] = new_msg[i][j];}}
        
    }
    
    if(converged==1)
	{
		if(DEBUG)		printf("BP: Converged in %d iterations\n",iter);
	}
    else
		printf("BP: Stopped after maxIter = %d iterations\n",iter);//Ã»ÓÐÊÕÁ²
    
    
    /* Assign Output Args */
    
    for(i = 0; i < nNodes; i++) {for(j = 0; j < nStates; j++) {
        nodeBel[i + nNodes*j] = new_bel[i][j];}}
    nIter[0] = (double)iter;
    for(i = 0; i < nEdges; i++) {for(j = 0; j < nStates; j++) {
        msgs[i + nEdges*j] = new_msg[i][j];}}
    
    
    /* Compute Edge Bel (divide version */
    
    /*for(i = 0; i < nNodes;i++) {
        nbrsInd = starV[i];
        nbrs = &starE[starV[i]];
        nNbrs = starV[i+1]-starV[i];
        for(n = 0; n < nNbrs; n++)
        {
            j = nbrs[n];
            E = nbrsInd+n;
            E2 = edgeNum(j,i,starV,starE);
            
            sum = 0;
            for(k_i = 0; k_i < nStates; k_i++)
            {
                for(k_j = 0; k_j < nStates; k_j++)
                {
                    if (msgs[E2+nEdges*k_i] != 0)
                        temp2 = nodeBel[i+nNodes*k_i]/msgs[E2+nEdges*k_i];
                    else temp2 = 0;
                    if (msgs[E+nEdges*k_j] != 0)
                        temp3 = nodeBel[j+nNodes*k_j]/msgs[E+nEdges*k_j];
                    else temp3 = 0;
                    edgeBel[k_i+nStates*(k_j+nStates*E)]=temp2*temp3*edgePot[k_i+nStates*(k_j + nStates*E)];
                    sum = sum+edgeBel[k_i+nStates*(k_j+nStates*E)];
                }
            }
            
            for(k_i = 0; k_i < nStates; k_i++)
            {
                for(k_j = 0; k_j < nStates; k_j++)
                {
                    edgeBel[k_i+nStates*(k_j+nStates*E)]=edgeBel[k_i+nStates*(k_j+nStates*E)]/sum;
                }
            }
            
            
        }
        
    }*/
    
    
    
    /* Compute Edge Bel (multiply version) */
    
    for(i = 0; i < nNodes;i++) {
        nbrsInd = starV[i];
        nbrs = &starE[starV[i]];
        nNbrs = starV[i+1]-starV[i];
        for(n = 0; n < nNbrs; n++)
        {
            j = nbrs[n];
            E = nbrsInd+n;
            
            sum = 0;
            for(k_i = 0; k_i < nStates; k_i++)
            {
                for(k_j = 0; k_j < nStates; k_j++)
                {
                    
                    /* Product of Local Potential with Edge Potential */
                    
                    temp2 = edgePot[k_i+nStates*(k_j+nStates*E)]*nodePot[i+nNodes*k_i]*nodePot[j+nNodes*k_j];
                    
                    /* Times messages from neighbors of i except j */
                    
                     nbrsInd2 = starV[i];
                     nbrs2 = &starE[starV[i]];
                     nNbrs2 = starV[i+1]-starV[i];
                     
                     for(n2 = 0; n2 < nNbrs2; n2++)
                        {
                        k = nbrs2[n2];
                        
                        if (k != j) {
                            E2 = edgeNum(k,i,starV,starE);
                            temp2 = temp2*msgs[E2 + nEdges*k_i];
                        }
                            
                     }
                     
                     /* Times messages from neighbors of j except i */
                     
                     nbrsInd2 = starV[j];
                     nbrs2 = &starE[starV[j]];
                     nNbrs2 = starV[j+1]-starV[j];
                     
                     for(n2 = 0; n2 < nNbrs2; n2++)
                        {
                        k = nbrs2[n2];
                        
                        if (k != i) {
                            E2 = edgeNum(k,j,starV,starE);
                            temp2 = temp2*msgs[E2 + nEdges*k_j];
                        }
                            
                     }
                     
                     edgeBel[k_i + nStates*(k_j + nStates*E)] = temp2;
                     sum = sum+temp2;
                }
                
            }
            
            for(k_i = 0; k_i < nStates; k_i++)
            {
                for(k_j = 0; k_j < nStates; k_j++)
                {
                    edgeBel[k_i+nStates*(k_j+nStates*E)]=edgeBel[k_i+nStates*(k_j+nStates*E)]/sum;
                }
            }
        }
    }
    
    
    /* Free Memory */
    
    for(i=0;i < nNodes;i++) {
        mxFree(prod_of_msgs[i]);
        mxFree(old_bel[i]);
        mxFree(new_bel[i]);
    }
    for(i=0;i < nEdges;i++) {
        mxFree(old_msg[i]);
        mxFree(new_msg[i]);
    }
    
    mxFree(prod_of_msgs);
    mxFree(old_bel);
    mxFree(new_bel);
    mxFree(old_msg);
    mxFree(new_msg);
    mxFree(newm);
    mxFree(pot_ij);
    mxFree(temp);
}

