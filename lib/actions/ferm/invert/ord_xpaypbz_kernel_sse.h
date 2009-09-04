#include <xmmintrin.h>
#include <emmintrin.h>

inline
void ord_xpaypbz_kernel(int lo, int hi, int my_id, ord_xpaypbz_arg* a)
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

  __m128 bv_re = _mm_set_ps(b_re, b_re, b_re, b_re);
  __m128 bv_im = _mm_set_ps(b_im, -b_im, b_im, -b_im);

  __m128 av_re = _mm_set_ps(a_re, a_re, a_re, a_re);
  __m128 av_im = _mm_set_ps(a_im, -a_im, a_im, -a_im);

  if( len % 4 == 0) { 
    for(int count = 0; count < len; count+=4) { 
      __m128 xv = _mm_load_ps(&x_ptr[count]);
      __m128 yv = _mm_load_ps(&y_ptr[count]);
      __m128 zv = _mm_load_ps(&z_ptr[count]);
      
      // yv = yv_4  yv_3  yv_2  yv_1 = im | re | im | re
      // t1 = yv_3  yv_4  yv_1  yv_2 = re | im | re | im 
      
      __m128 t1 = _mm_shuffle_ps(yv,yv, 0xb1);
      __m128 t2 = _mm_mul_ps(av_re, yv);
      __m128 t3 = _mm_add_ps(xv,t2);
      __m128 t4 = _mm_mul_ps(av_im, t1);
      xv = _mm_add_ps(t3, t4);
      
      t1 = _mm_shuffle_ps(zv,zv,0xb1);
      t2 = _mm_mul_ps(bv_re, zv);
      t3 = _mm_add_ps(xv, t2);
      t4 = _mm_mul_ps(bv_im, t1);
      xv = _mm_add_ps(t3,t4);
      
      _mm_store_ps(&x_ptr[count], xv);
      // _mm_stream_ps(&x_ptr[count], xv);
      
#if 0
      
      REAL32 tmp_re, tmp_im;
      REAL32 tmp_re2, tmp_im2;
      
      tmp_re  = x_ptr[count] + a_re*y_ptr[count];
      tmp_re -=  a_im*y_ptr[count+1];
      tmp_im  = x_ptr[count+1] + a_re*y_ptr[count+1];
      tmp_im +=  a_im*y_ptr[count];
      
      
      
      tmp_re2  = x_ptr[count+2] + a_re*y_ptr[count+2];
      tmp_re2 -=  a_im*y_ptr[count+3];
      tmp_im2  = x_ptr[count+3] + a_re*y_ptr[count+3];
      tmp_im2 +=  a_im*y_ptr[count+2];
      
      
      x_ptr[count]   = tmp_re + b_re*z_ptr[count] ;
      x_ptr[count]  -=  b_im*z_ptr[count+1];
      x_ptr[count+1] = tmp_im + b_re*z_ptr[count+1];
      x_ptr[count+1]+= b_im *z_ptr[count]; 
      
      
      x_ptr[count+2]   = tmp_re2 + b_re*z_ptr[count+2] ;
      x_ptr[count+2]  -=  b_im*z_ptr[count+3];
      x_ptr[count+3] = tmp_im2 + b_re*z_ptr[count+3];
      x_ptr[count+3]+= b_im *z_ptr[count+2]; 
#endif
      
    }
  }
  else { 
    QDPIO::cout << "ord_xpaypbz_kernel_sse.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
}
