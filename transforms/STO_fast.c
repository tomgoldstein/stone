/*				STO_fast.c  by Tom Goldstein
 *   This code implements the fast Sum-To-One Transform
 * 
 */

#include <math.h>
#include <stdlib.h>
#include "mex.h"


/*********MEX Interface********/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double* u;
    double* end;
    double *first, *second, *third, *fourth;
    double *start, *stop;
    int gap;
    int m,n;
    double a,b,c,d;
    
    u = mxGetPr(prhs[0]);
    
    if(!mxIsClass(prhs[0],"double")){
        mexErrMsgTxt("Data must be of type double!!!\n");
    }
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if(n>1 && m>1){
        mexErrMsgTxt("Only vectors allowed.  No Matrices.\n");
    }
     
    if(n==1)
        n=m;
    
    end = u+n;
    
    /* This funky block checks that the length is a power of 4 */
    if ( n & (n - 1)  
     || (n & 0xAAAAAAAB ) ){
        mexErrMsgTxt("Vector length must be power of 4.\n");
    }
    
    
   
  
    /*Do the First level transform.  This is done in its own loop for efficiency*/
    first = u;
    second = u+1;
    third = u+2;
    fourth = u+3;
    for(;first<end; first+=4, second+=4, third+=4, fourth+=4){
        a = (*first)+(*second)+(*third)-(*fourth);
        b = (*first)+(*second)-(*third)+(*fourth);
        c = (*first)-(*second)+(*third)+(*fourth);
        d = -(*first)+(*second)+(*third)+(*fourth);
        *first = a*.5;
        *second = b*.5;
        *third = c*.5;
        *fourth = d*.5;
    }
    
    /* Do every other level of the transform */
    for(gap = 4; gap<n;gap*=4){
        first = u;
        second = u+gap;
        third = u+2*gap;
        fourth = u+3*gap;
        for(start = u, stop=u+gap; start<end; start+=gap*4, stop+=gap*4){
            for(;first<stop; first++, second++, third++, fourth++){
                a = (*first)+(*second)+(*third)-(*fourth);
                b = (*first)+(*second)-(*third)+(*fourth);
                c = (*first)-(*second)+(*third)+(*fourth);
                d = -(*first)+(*second)+(*third)+(*fourth);
                *first = a*.5;
                *second = b*.5;
                *third = c*.5;
                *fourth = d*.5;
            }
            first+=3*gap;
            second+=3*gap;
            third+=3*gap;
            fourth+=3*gap;
        }
    }
    
    
    return;
}


    

