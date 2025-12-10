/* $Id: qdp-lapack_Complex.h,v 1.3 2009-10-21 20:50:56 kostas Exp $ */
/* Slight modification of the SuperLU header file
 * to avoid conflicts with f2c and g2c libraries 
 * A. Stathopoulos, Feb 7, 2008                     */
/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */
#ifndef __COMPLEX_C_H
#define __COMPLEX_C_H
#ifdef __cplusplus
extern "C" {
#endif 
/*------------------------------------------------------------------------
 * Single Precision complex
 *-----------------------------------------------------------------------*/
typedef struct { 
	float r, i; 
} Complex_C;

/* Macro definitions */

/* Complex Addition c = a + b */
#define c_add_primme(c, a, b) { (c).r = (a).r + (b).r; \
			 (c).i = (a).i + (b).i; }

/* Complex Subtraction c = a - b */
#define c_sub_primme(c, a, b) { (c).r = (a).r - (b).r; \
			 (c).i = (a).i - (b).i; }

/* Complex-Double Multiplication */
#define cd_mult_primme(c, a, b) { (c).r = (a).r * (b); \
                           (c).i = (a).i * (b); }

/* Complex-Complex Multiplication */
#define cc_mult_primme(c, a, b) { \
	double cr, ci; \
    	cr = (a).r * (b).r - (a).i * (b).i; \
    	ci = (a).i * (b).r + (a).r * (b).i; \
    	(c).r = cr; \
    	(c).i = ci; \
    }

/* Complex equality testing */
#define c_eq_primme(a, b)  ( (a).r == (b).r && (a).i == (b).i )


/* Prototypes for functions in dcomplex.c */
void   c_mul_primme(Complex_C *c, Complex_C *a, Complex_C *b);
void   c_div_primme(Complex_C *, Complex_C *, Complex_C *);
float  c_abs_primme(Complex_C);     /* exact sqrt(r^2+i^2) */
float  c_abs1_primme(Complex_C);    /* approximate  |r|+|i| */
void   c_exp_primme(Complex_C *, Complex_C *);
void   s_cnjg_primme(Complex_C *, Complex_C *);
float  s_imag_primme(Complex_C *);

/*------------------------------------------------------------------------
 * Double Precision complex
 *-----------------------------------------------------------------------*/
typedef struct { 
	double r, i; 
} Complex_Z;

/* Macro definitions */

/* Complex Addition c = a + b */
#define z_add_primme(c, a, b) { (c).r = (a).r + (b).r; \
			 (c).i = (a).i + (b).i; }

/* Complex Subtraction c = a - b */
#define z_sub_primme(c, a, b) { (c).r = (a).r - (b).r; \
			 (c).i = (a).i - (b).i; }

/* Complex-Double Multiplication */
#define zd_mult_primme(c, a, b) { (c).r = (a).r * (b); \
                           (c).i = (a).i * (b); }

/* Complex-Complex Multiplication */
#define zz_mult_primme(c, a, b) { \
	double cr, ci; \
    	cr = (a).r * (b).r - (a).i * (b).i; \
    	ci = (a).i * (b).r + (a).r * (b).i; \
    	(c).r = cr; \
    	(c).i = ci; \
    }

/* Complex equality testing */
#define z_eq_primme(a, b)  ( (a).r == (b).r && (a).i == (b).i )


/* Prototypes for functions in dcomplex.c */
void   z_mul_primme(Complex_Z *c, Complex_Z *a, Complex_Z *b);
void   z_div_primme(Complex_Z *, Complex_Z *, Complex_Z *);
double z_abs_primme(Complex_Z);     /* exact sqrt(r^2+i^2) */
double z_abs1_primme(Complex_Z);    /* approximate  |r|+|i| */
void   z_exp_primme(Complex_Z *, Complex_Z *);
void   d_cnjg_primme(Complex_Z *r, Complex_Z *z);
double d_imag_primme(Complex_Z *);
#ifdef __cplusplus
};
#endif 
#endif
