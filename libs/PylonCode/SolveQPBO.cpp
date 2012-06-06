#include "mex.h"
#include "QPBO.h"
#include <time.h>

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

// Instantiations

template class QPBO<int>;
template class QPBO<float>;
template class QPBO<double>;

template <> 
	inline void QPBO<int>::get_type_information(char*& type_name, char*& type_format)
{
	type_name = "int";
	type_format = "d";
}

template <> 
	inline void QPBO<float>::get_type_information(char*& type_name, char*& type_format)
{
	type_name = "float";
	type_format = "f";
}

template <> 
	inline void QPBO<double>::get_type_information(char*& type_name, char*& type_format)
{
	type_name = "double";
	type_format = "Lf";
}


//mexify: mex SolveQPBO.cpp qpbo_postprocessing.cpp qpbo_extra.cpp qpbo_maxflow.cpp qpbo.cpp
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    if(nrhs != 3)
        mexErrMsgTxt("3 input parameters required!");
    if(mxGetM(prhs[0]) != 2 || mxGetM(prhs[1]) != 2 || mxGetM(prhs[2]) != 4 || mxGetN(prhs[1]) != mxGetN(prhs[2]))
        mexErrMsgTxt("Inconsistent input size!");
            
    int nNodes = mxGetN(prhs[0]);
    int nEdges = mxGetN(prhs[1]);
    
    double *nodeData = mxGetPr(prhs[0]);
    double *nodeEdge = mxGetPr(prhs[1]);    
    double *edgeData = mxGetPr(prhs[2]);    
    
    QPBO<double> *qpbo = new QPBO<double>(nNodes,nEdges);
    qpbo->AddNode(nNodes);
    
    int i;
    for(i = 0; i < nNodes; i++)
        qpbo->AddUnaryTerm(i, nodeData[2*i], nodeData[2*i+1]);
    for(i = 0; i < nEdges; i++)
        qpbo->AddPairwiseTerm(int(nodeEdge[2*i])-1, int(nodeEdge[2*i+1])-1,
            edgeData[4*i], edgeData[4*i+1], edgeData[4*i+2], edgeData[4*i+3]);
    
    clock_t start = clock();
	qpbo->Solve();
    qpbo->ComputeWeakPersistencies();
    clock_t end = clock();
    //mexPrintf("QPBO took %f seconds\n", double(end-start)/CLOCKS_PER_SEC);

    double *ptr;
    if(nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        ptr = mxGetPr(plhs[1]);
        *ptr = qpbo->ComputeTwiceEnergy()*0.5;
    } 
    
  
    plhs[0] = mxCreateDoubleMatrix(1, nNodes, mxREAL);
    ptr = mxGetPr(plhs[0]);    
    
    for(int i = 0; i < qpbo->GetNodeNum(); i++)
        ptr[i] =  qpbo->GetLabel(i);
    
    
	delete qpbo;
}
