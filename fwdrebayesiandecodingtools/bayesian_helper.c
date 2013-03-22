#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{


  /* cryptic variables */
  long p,c,t;
  long i,j,k;
  mxArray *out;
  double *posterior;
  double *logRateMapsPtr, *spkCountsPtr;
  long idx1, idx2, idx3;
  int *pjc, *pir;
  int nelements;
  /* no error checking */

  /* get array sizes */
  p = mxGetM( prhs[0] );
  c = mxGetN( prhs[0] );
  t = mxGetN( prhs[1] );

  logRateMapsPtr = mxGetPr( prhs[0] );  // log(rateMaps)  2D array, rateMap(:)-by-nTimeBins
  spkCountsPtr = mxGetPr( prhs[1] );  // spikecounts  nCells - by- nTimeBins 

  /* create output array */
  out = mxCreateDoubleMatrix( p, t , mxREAL ); // returns pointer to the created mxArray
  posterior = mxGetPr( out );  // ?? out and posterior point to the same mxArray, 
  
  for (i=0;i<p*t;i++) {
    posterior[i] = 1;
  }

  if (mxIsSparse( prhs[1] )) {

    idx2 = -1;
   
    pjc = mxGetJc( prhs[1] );
    pir = mxGetIr( prhs[1] );

    for (j=0;j<t;j++) {

      nelements = pjc[j+1]-pjc[j];
      idx1 = j*p;

      for (i=0;i<nelements;i++) {

		idx2++;
		idx3 = pir[idx2]*p;

		for (k=0;k<p;k++) {
	  
		posterior[idx1+k] *= pow( logRateMapsPtr[idx3+k], spkCountsPtr[idx2] );

		}
      }
    }
  } else {
    
    for (i=0; i<c; i++) {
      
      idx3 = i*p;
      
      for(j=0;j<t;j++) {
	
		idx2 = j*c+i;
		if (spkCountsPtr[idx2]!=0) {

			idx1 = j*p;
	  
			for(k=0;k<p;k++) {
				posterior[idx1+k] += spkCountsPtr[idx2]*logRateMapsPtr[idx3+k];
			}
	  
		}
	
      }
      
    }

  }

  plhs[0] = out;
	
}
