//MAC OS X
//#include <vecLib/vBLAS.h>
//#include <veclib/clapack.h>

// May need other includes for other OS 
#include "chromabase.h"
extern "C" {
int zheev_(char *,char *,
	    int *, DComplex *, int *, Double *, DComplex *, 
	    int *, Double *, int *);

int  zgeqrf_(int *, int *,
	     DComplex *, int *,
	     DComplex *, DComplex *,
	     int *, int *) ;

int zunmqr_(char *, char *,
	    int *, int *,
	    int *, DComplex *,
	    int *, DComplex *,  DComplex *,
	    int *, DComplex *,
	    int *, int *);
}

