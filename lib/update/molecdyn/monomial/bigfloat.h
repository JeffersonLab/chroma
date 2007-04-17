// -*- C++ -*-
// $Id: bigfloat.h,v 3.1 2007-04-17 03:13:04 edwards Exp $
/*! \file
 *  \brief Remez algorithm for finding nth roots
 */

#ifndef __bigfloat_h__
#define __bigfloat_h__

#include "chromabase.h"
#include <gmp.h>

namespace Chroma
{

  //! Bigfloat
  /*! @ingroup monomial
   *
   * Simple C++ wrapper for multiprecision datatype used for Remez
   * algorithm
   *
   */
  class bigfloat 
  {
  public:
    bigfloat() { mpf_init(x); }
    bigfloat(const bigfloat& y) { mpf_init_set(x, y.x); }
    bigfloat(const unsigned long u) { mpf_init_set_ui(x, u); }
    bigfloat(const long i) { mpf_init_set_si(x, i); }
    bigfloat(const int i) {mpf_init_set_si(x,(long)i);}
    bigfloat(const double d) { mpf_init_set_d(x, toDouble(d)); }
    bigfloat(const Real32& d) { mpf_init_set_d(x, toDouble(d)); }
    bigfloat(const Real64& d) { mpf_init_set_d(x, toDouble(d)); }
    bigfloat(const char *str) { mpf_init_set_str(x, str, 10); }
    ~bigfloat(void) { mpf_clear(x); }
    operator const double (void) const { return (double)(mpf_get_d(x)); }
    static void setDefaultPrecision(unsigned long dprec) {
      unsigned long bprec =  (unsigned long)(3.321928094 * (double)dprec);
      mpf_set_default_prec(bprec);
    }
    
    void setPrecision(unsigned long dprec) {
      unsigned long bprec =  (unsigned long)(3.321928094 * (double)dprec);
      mpf_set_prec(x,bprec);
    }
  
    unsigned long getPrecision(void) const { return mpf_get_prec(x); }

    unsigned long getDefaultPrecision(void) const { return mpf_get_default_prec(); }

    bigfloat& operator=(const bigfloat& y) {
      mpf_set(x, y.x); 
      return *this;
    }

    bigfloat& operator=(const unsigned long y) { 
      mpf_set_ui(x, y);
      return *this; 
    }
  
    bigfloat& operator=(const signed long y) {
      mpf_set_si(x, y); 
      return *this;
    }
  
    bigfloat& operator=(const double y) {
      mpf_set_d(x, toDouble(y)); 
      return *this;
    }

    bigfloat& operator=(const Real32& y) {
      mpf_set_d(x, toDouble(y)); 
      return *this;
    }

    bigfloat& operator=(const Real64& y) {
      mpf_set_d(x, toDouble(y)); 
      return *this;
    }

    size_t write(void);
    size_t read(void);

    /* Arithmetic Functions */

    bigfloat& operator+=(const bigfloat& y) { return *this = *this + y; }
    bigfloat& operator-=(const bigfloat& y) { return *this = *this - y; }
    bigfloat& operator*=(const bigfloat& y) { return *this = *this * y; }
    bigfloat& operator/=(const bigfloat& y) { return *this = *this / y; }

    friend bigfloat operator+(const bigfloat& x, const bigfloat& y) {
      bigfloat a;
      mpf_add(a.x,x.x,y.x);
      return a;
    }

    friend bigfloat operator+(const bigfloat& x, const unsigned long y) {
      bigfloat a;
      mpf_add_ui(a.x,x.x,y);
      return a;
    }

    friend bigfloat operator-(const bigfloat& x, const bigfloat& y) {
      bigfloat a;
      mpf_sub(a.x,x.x,y.x);
      return a;
    }
  
    friend bigfloat operator-(const unsigned long x, const bigfloat& y) {
      bigfloat a;
      mpf_ui_sub(a.x,x,y.x);
      return a;
    }
  
    friend bigfloat operator-(const bigfloat& x, const unsigned long y) {
      bigfloat a;
      mpf_sub_ui(a.x,x.x,y);
      return a;
    }

    friend bigfloat operator-(const bigfloat& x) {
      bigfloat a;
      mpf_neg(a.x,x.x);
      return a;
    }

    friend bigfloat operator*(const bigfloat& x, const bigfloat& y) {
      bigfloat a;
      mpf_mul(a.x,x.x,y.x);
      return a;
    }

    friend bigfloat operator*(const bigfloat& x, const unsigned long y) {
      bigfloat a;
      mpf_mul_ui(a.x,x.x,y);
      return a;
    }

    friend bigfloat operator/(const bigfloat& x, const bigfloat& y){
      bigfloat a;
      mpf_div(a.x,x.x,y.x);
      return a;
    }

    friend bigfloat operator/(const unsigned long x, const bigfloat& y){
      bigfloat a;
      mpf_ui_div(a.x,x,y.x);
      return a;
    }

    friend bigfloat operator/(const bigfloat& x, const unsigned long y){
      bigfloat a;
      mpf_div_ui(a.x,x.x,y);
      return a;
    }

    friend bigfloat sqrt_bf(const bigfloat& x){
      bigfloat a;
      mpf_sqrt(a.x,x.x);
      return a;
    }

    friend bigfloat sqrt_bf(const unsigned long x){
      bigfloat a;
      mpf_sqrt_ui(a.x,x);
      return a;
    }

    friend bigfloat abs_bf(const bigfloat& x){
      bigfloat a;
      mpf_abs(a.x,x.x);
      return a;
    }

    friend bigfloat pow_bf(const bigfloat& a, long power) {
      bigfloat b;
      mpf_pow_ui(b.x,a.x,power);
      return b;
    }

    /* Comparison Functions */

    friend int operator>(const bigfloat& x, const bigfloat& y) {
      int test;
      test = mpf_cmp(x.x,y.x);
      if (test > 0) return 1;
      else return 0;
    }

    friend int operator<(const bigfloat& x, const bigfloat& y) {
      int test;
      test = mpf_cmp(x.x,y.x);
      if (test < 0) return 1;
      else return 0;
    }

    friend int sgn(const bigfloat&);

    /* Miscellaneous Functions */

    friend bigfloat& random(void);

  private:

    mpf_t x;

  };

}  // end namespace Chroma

#endif


