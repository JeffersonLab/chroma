#include <xmmintrin.h>
#include <emmintrin.h>

inline
void ord_yxpaymabz_kernel(int lo, int hi, int my_id, ord_yxpaymabz_arg* a)
{
  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi - lo);

  REAL32* x_ptr = &(a->x_ptr[low]);
  REAL32* y_ptr = &(a->y_ptr[low]);
  REAL32* z_ptr = &(a->z_ptr[low]);
  
  REAL32 a_re = a->a_re;
  REAL32 a_im = a->a_im;
  REAL32 b_re = a->b_re;
  REAL32 b_im = a->b_im;
  
  typedef union {
    float v[4];
    __m128 vec;
  } Vec4;

  __m128 av_re = _mm_set_ps(a_re, a_re, a_re, a_re);
  __m128 av_im = _mm_set_ps(a_im, -a_im, a_im,-a_im);
  __m128 bv_re = _mm_set_ps(b_re, b_re, b_re, b_re);
  __m128 bv_im = _mm_set_ps(-b_im, b_im, -b_im, b_im);
  
  if( len % 4 == 0) {
    for(int count = 0; count < len; count+=4) { 
      // Load x, y, z
      __m128 xv = _mm_load_ps(&x_ptr[count]);
      __m128 yv = _mm_load_ps(&y_ptr[count]);
      __m128 zv = _mm_load_ps(&z_ptr[count]);
      
      // Step1: t = y - b*z
      //            = yv - bv_re zv + bv_im zv2 
      
      // zv  = [ re | im | re | im ]
      // zv2 = [ im | re | im | re ]
      __m128 zv2 = zv;
      zv2 = _mm_shuffle_ps(zv2,zv2,0xb1);
      
      
      // tmp = y - bv_re zv + bv_im zv2
      __m128 t1 = _mm_mul_ps(bv_re, zv);
      __m128 t2 = _mm_sub_ps(yv,t1);
      
      __m128 t3 = _mm_mul_ps(bv_im, zv2);
      t2 = _mm_add_ps(t2,t3);
      
      
      // zv2 holds t2
      zv2 = t2;
      zv2 = _mm_shuffle_ps(zv2,zv2,0xb1);
      // yv= x + av_re t2 + av_im zv2
      
      t1 = _mm_mul_ps(av_re, t2);
      yv = _mm_add_ps(xv, t1);
      t3 = _mm_mul_ps(av_im, zv2);
      yv = _mm_add_ps(yv, t3);
      
      _mm_store_ps(&y_ptr[count], yv);
      
      
    }
  }
  else { 
    QDPIO::cout << "ord_yxpaymabz_kernel_sse.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
  
}
