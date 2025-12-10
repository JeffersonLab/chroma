// $Id: qdp-lapack_biortho.h,v 1.1 2009-10-21 20:50:56 kostas Exp $
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include "qdp-lapack_Complex.h"
#include "qdp-lapack_numerical.h"
#include "qdp-lapack_numerical_private.h"

void biortho_local(Complex_Z *VL, int ldvl, Complex_Z *VR, int ldvr, int m, int ni, int nf, int ktimes);

void biortho_global_C(Complex_C *VL, int ldvl, Complex_C *VR, int ldvr, int m, int ni, int nf, int ktimes, void *params);


void biortho_global_Z(Complex_Z *VL, int ldvl, Complex_Z *VR, int ldvr, int m, int ni, int nf, int ktimes, void *params);
