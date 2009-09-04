#include <cstdio>
using namespace std;


inline
void ord_ib_stupdates_kernel_real32(int lo, int hi, int my_id, ib_stupdate_arg<REAL32>* a)
{
  int atom=a->atom;
  int low = lo*atom;
  int len = atom*(hi - lo);

  REAL32 a_r = a->a_r;
  REAL32 a_i = a->a_i;

  REAL32* r = &(a->r[low]);
  REAL32* u = &(a->u[low]);
  REAL32* v = &(a->v[low]);
  REAL32* q = &(a->q[low]);
  REAL32* r0 = &(a->r0[low]);
  REAL32* f0 = &(a->f0[low]);
  REAL32* s =  &(a->s[low]);
  REAL32* t =  &(a->t[low]);
  REAL64* norm_array = &(a->norm_space[12*my_id]);

  // Caller zeroed norm_space
  if( len % 4 == 0) { 
    for(int count = 0; count < len; count+=4) { 
      // First need s = r - alpha*v
      
      // Real Part -(+,-)
      s[count] = r[count] - a_r*v[count];
      s[count] += a_i*v[count+1];
      // Imag Part -(+,+)
      s[count+1] = r[count+1] - a_r*v[count+1];
      s[count+1] -= a_i*v[count];
      
      // Real Part -(+,-)
      s[count+2] = r[count+2] - a_r*v[count+2];
      s[count+2] += a_i*v[count+3];
      // Imag Part -(+,+)
      s[count+3] = r[count+3] - a_r*v[count+3];
      s[count+3] -= a_i*v[count+2];
      
      // Now inner products with s. 
      // ** phi=(r0,s) 
      norm_array[0] += r0[count]*s[count];
      norm_array[0] += r0[count+1]*s[count+1];
      
      norm_array[0] += r0[count+2]*s[count+2];
      norm_array[0] += r0[count+3]*s[count+3];
      
      norm_array[1] += r0[count]*s[count+1];
      norm_array[1] -= r0[count+1]*s[count];
      
      norm_array[1] += r0[count+2]*s[count+3];
      norm_array[1] -= r0[count+3]*s[count+2];
      
      // ** gamma =(f0,s) 
      norm_array[2] += f0[count]*s[count];
      norm_array[2] += f0[count+1]*s[count+1];
      
      norm_array[2] += f0[count+2]*s[count+2];
      norm_array[2] += f0[count+3]*s[count+3];
      
      norm_array[3] += f0[count]*s[count+1];
      norm_array[3] -= f0[count+1]*s[count];
      
      
      norm_array[3] += f0[count+2]*s[count+3];
      norm_array[3] -= f0[count+3]*s[count+2];
      
      
      // T update: t = u - alpha q
      // Real Part -(+,-)
      t[count]  = u[count] - a_r*q[count];
      t[count] += a_i*q[count+1];
      // Imag Part -(+,+)
      t[count+1] = u[count+1] - a_r*q[count+1];
      t[count+1] -= a_i*q[count];
      
      // Real Part -(+,-)
      t[count+2] = u[count+2] - a_r*q[count+2];
      t[count+2] += a_i*q[count+3];
      // Imag Part -(+,+)
      t[count+3] = u[count+3] - a_r*q[count+3];
      t[count+3] -= a_i*q[count+2];
      
      
      // Now inner products with q. 
      // ** pi=(r0,q) 
      norm_array[4] += r0[count]*q[count];
      norm_array[4] += r0[count+1]*q[count+1];
      
      norm_array[4] += r0[count+2]*q[count+2];
      norm_array[4] += r0[count+3]*q[count+3];
      
      
      norm_array[5] += r0[count]*q[count+1];
      norm_array[5] -= r0[count+1]*q[count];
      
      norm_array[5] += r0[count+2]*q[count+3];
      norm_array[5] -= r0[count+3]*q[count+2];
      
      // Now inner products with t. 
      // ** eta=(f0,t) 
      norm_array[6] += f0[count]*t[count];
      norm_array[6] += f0[count+1]*t[count+1];
      
      norm_array[6] += f0[count+2]*t[count+2];
      norm_array[6] += f0[count+3]*t[count+3];
      
      norm_array[7] += f0[count]*t[count+1];
      norm_array[7] -= f0[count+1]*t[count];
      
      norm_array[7] += f0[count+2]*t[count+3];
      norm_array[7] -= f0[count+3]*t[count+2];
      
      // ** theta=(t,s) 
      norm_array[8] += t[count]*s[count];
      norm_array[8] += t[count+1]*s[count+1];
      
      norm_array[8] += t[count+2]*s[count+2];
      norm_array[8] += t[count+3]*s[count+3];
      
      
      norm_array[9] += t[count]*s[count+1];
      norm_array[9] -= t[count+1]*s[count];
      
      norm_array[9] += t[count+2]*s[count+3];
      norm_array[9] -= t[count+3]*s[count+2];
      
      // ** kappa = || t ||^2
      norm_array[10] += t[count]*t[count];
      norm_array[10] += t[count+1]*t[count+1];
      
      norm_array[10] += t[count+2]*t[count+2];
      norm_array[10] += t[count+3]*t[count+3];
      
      // ** rnorm = || r ||^2
      norm_array[11] += r[count]*r[count];
      norm_array[11] += r[count+1]*r[count+1];
      
      norm_array[11] += r[count+2]*r[count+2];
      norm_array[11] += r[count+3]*r[count+3];
      
      
    }
  }
  else { 
    QDPIO::cout << "ord_ib_stubdates_kernel_generic.h: len not divisibe by 4" << endl;
    QDP_abort(1);
  }

}

inline
void ord_ib_stupdates_kernel_real64(int lo, int hi, int my_id, ib_stupdate_arg<REAL64>* a)
{
  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi-lo);

  REAL64 a_r = a->a_r;
  REAL64 a_i = a->a_i;

  REAL64* r = &(a->r[low]);
  REAL64* u = &(a->u[low]);
  REAL64* v = &(a->v[low]);
  REAL64* q = &(a->q[low]);
  REAL64* r0 = &(a->r0[low]);
  REAL64* f0 = &(a->f0[low]);
  REAL64* s = &(a->s[low]);
  REAL64* t = &(a->t[low]);
  REAL64* norm_array = &(a->norm_space[12*my_id]);
  
  // Caller zeroed norm_space
  
  if (len % 2 == 0 ) {
    for(int count = 0; count < len; count+=2) { 
      // First need s = r - alpha*v
      
      // Real Part -(+,-)
      s[count] = r[count] - a_r*v[count];
      s[count] += a_i*v[count+1];
      // Imag Part -(+,+)
      s[count+1] = r[count+1] - a_r*v[count+1];
      s[count+1] -= a_i*v[count];
      
      
      // Now inner products with s. 
      // ** phi=(r0,s) 
      norm_array[0] += r0[count]*s[count];
      norm_array[0] += r0[count+1]*s[count+1];
      
      norm_array[1] += r0[count]*s[count+1];
      norm_array[1] -= r0[count+1]*s[count];
      
      // ** gamma=(f0,s) 
      norm_array[2] += f0[count]*s[count];
      norm_array[2] += f0[count+1]*s[count+1];
      
      norm_array[3] += f0[count]*s[count+1];
      norm_array[3] -= f0[count+1]*s[count];
      
      
      // T update: t = u - alpha q
      // Real Part -(+,-)
      t[count]  = u[count] - a_r*q[count];
      t[count] += a_i*q[count+1];
      // Imag Part -(+,+)
      t[count+1] = u[count+1] - a_r*q[count+1];
      t[count+1] -= a_i*q[count];
      
      
      // Now inner products with q. 
      // ** pi=(r0,q) 
      norm_array[4] += r0[count]*q[count];
      norm_array[4] += r0[count+1]*q[count+1];
      
      norm_array[5] += r0[count]*q[count+1];
      norm_array[5] -= r0[count+1]*q[count];
      
      
      // Now inner products with t. 
      // ** eta=(f0,t) 
      norm_array[6] += f0[count]*t[count];
      norm_array[6] += f0[count+1]*t[count+1];
      
      norm_array[7] -= f0[count+1]*t[count];
      norm_array[7] += f0[count]*t[count+1];
      
      
      // ** theta=(t,s) 
      norm_array[8] += t[count]*s[count];
      norm_array[8] += t[count+1]*s[count+1];
      
      norm_array[9] += t[count]*s[count+1];
      norm_array[9] -= t[count+1]*s[count];
      
      
      // ** kappa = || t ||^2
      norm_array[10] += t[count]*t[count];
      norm_array[10] += t[count+1]*t[count+1];
      
      // ** rnorm = || r ||^2
      norm_array[11] += r[count]*r[count];
      norm_array[11] += r[count+1]*r[count+1];
      
    }
  }
  else { 
    QDPIO::cout << "ord_bi_stubdates_kernel_generic.h: len not divisible by 2" <<endl;
    QDP_abort(1);
  }
  
}
