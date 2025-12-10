#include "qdp-lapack_config_internal.h"
#ifdef USE_QMP
#include <qmp.h>

#ifdef __cplusplus
#include <cstdio>
#define fprintf if(QMP_get_node_number()==0) std::fprintf
#define printf  if(QMP_get_node_number()==0) std::printf
#else
#define fprintf if(QMP_get_node_number()==0) fprintf
#define printf  if(QMP_get_node_number()==0) printf
#endif
#endif

void wrap_cgemv(char *transa, int *m, int *n, Complex_C *alpha, Complex_C *a,
   int *lda, Complex_C *x, int *incx, Complex_C *beta, Complex_C *y, int *incy,
   Complex_C *work, void *);
Complex_C wrap_cdot(int *n, Complex_C *x, int *incx, Complex_C *y, int *incy, void *);
void wrap_cgemm(char *transa, char *transb, int *m, int *n, int *k,
   Complex_C *alpha, Complex_C *a, int *lda, Complex_C *b, int *ldb,
   Complex_C *beta, Complex_C *c, int *ldc,
   Complex_C *work, void *params);

void wrap_zdot_small(Complex_Z *zdot_r, int *n, Complex_Z *x, int *incx, Complex_Z *y, int *incy);
void wrap_cdot_small(Complex_C *zdot_r, int *n, Complex_C *x, int *incx, Complex_C *y, int *incy);

Complex_C wrap_zsum_cdot(int *n, Complex_C *x, int *incx, Complex_C *y, int *incy, void *);
Complex_C zsum_cdot(int *n, Complex_C *x, int *incx, Complex_C *y, int *incy);

void wrap_zgemv(char *transa, int *m, int *n, Complex_Z *alpha, Complex_Z *a,
   int *lda, Complex_Z *x, int *incx, Complex_Z *beta, Complex_Z *y, int *incy,
   Complex_Z *work,  void *);
Complex_Z wrap_zdot(int *n, Complex_Z *x, int *incx, Complex_Z *y, int *incy, void *);

void globalSumDouble(void *sendBuf, void *recvBuf, int *count, void *params);
void globalSumFloat(void *sendBuf, void *recvBuf, int *count, void *params);
