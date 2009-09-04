#include <xmmintrin.h>
#include <emmintrin.h>

inline
void ord_norm2x_cdotxy_kernel(int lo, int hi, int my_id, ord_norm2x_cdotxy_arg* a)
{
  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi-lo);

  REAL32* x_ptr = &(a->x_ptr[low]);
  REAL32* y_ptr = &(a->y_ptr[low]);
  REAL64 norm_array[3] = {0,0,0};

  __m128d sum = _mm_set_pd((double)0,(double)0);
  __m128d dotprod = _mm_set_pd((double)0,(double)0);
  __m128d mask = _mm_set_pd((double)-1, (double)1);

  if( atom % 4 == 0) { 
    for(int count = 0; count < len; count+=4) { 
      
      // 
      __m128 xv = _mm_load_ps(&x_ptr[count]);
      __m128 yv = _mm_load_ps(&y_ptr[count]);
      
      __m128d xlow,xhi;
      __m128d ylow,yhi;
      
      // Convert low 2 in X to be doubles in xlow
      xlow = _mm_cvtps_pd(xv);
      
      // Swap low and high pairs 
      xv = _mm_shuffle_ps( xv,xv, 0x4e);
      
      // Take hi part of x into xhi
      xhi = _mm_cvtps_pd(xv);
      
      // Now samething with y
      ylow = _mm_cvtps_pd(yv);
      yv = _mm_shuffle_ps( yv,yv, 0x4e);
      yhi = _mm_cvtps_pd(yv);
      
      sum = _mm_add_pd(sum,_mm_mul_pd(xlow,xlow));
      sum = _mm_add_pd(sum,_mm_mul_pd(xhi,xhi));
      
      
      // norm
      // norm_array[0] += x_ptr[count]*x_ptr[count];
      // norm_array[1] += x_ptr[count+1]*x_ptr[count+1];
      // norm_array[0] += x_ptr[count+2]*x_ptr[count+2];
      // norm_array[1] += x_ptr[count+3]*x_ptr[count+3];
      
      // Need tmp to hold (y0, y_0)
      //    dotprod[0] += x_ptr[count]*y_ptr[count];
      //    dotprod[1] -= x_ptr[count+1]*y_ptr[count];
      
      __m128d t1 = _mm_shuffle_pd(ylow, ylow, 0x0);
      __m128d t2 = _mm_mul_pd(mask, t1);
      __m128d t3 = _mm_mul_pd(xlow, t2);
      dotprod = _mm_add_pd(dotprod, t3);
      
      //dotprod[0] += x_ptr[count+1]*y_ptr[count+1];
      //dotprod[1] += x_ptr[count]*y_ptr[count+1];
      t1 = _mm_shuffle_pd(ylow, ylow, 0x3);
      t2 = _mm_shuffle_pd(xlow, xlow, 0x1);
      t3 = _mm_mul_pd(t1,t2);
      dotprod = _mm_add_pd(dotprod,t3);
      
      
      // dotprod[0] += x_ptr[count+2]*y_ptr[count+2];
      // dotprod[1] -= x_ptr[count+3]*y_ptr[count+2]
      t1 = _mm_shuffle_pd(yhi,yhi,0x0);
      t2 = _mm_mul_pd(mask,t1);
      t3 = _mm_mul_pd(xhi,t2);
      dotprod = _mm_add_pd(dotprod,t3);
      
      // dotprod[0] += x_ptr[count+3]*y_ptr[count+3];    
      // dotprod[1] += x_ptr[count+2]*y_ptr[count+3];
      t1 = _mm_shuffle_pd(yhi,yhi, 0x3);
      t2 = _mm_shuffle_pd(xhi,xhi,0x1);
      t3 = _mm_mul_pd(t1, t2);
      dotprod = _mm_add_pd(dotprod,t3);
    }
    
    
    a->norm_space[3*my_id]=((double *)&sum)[0] + ((double *)&sum)[1];
    a->norm_space[3*my_id+1]=((double *)&dotprod)[0];
    a->norm_space[3*my_id+2]=((double *)&dotprod)[1];
  }
  else { 
    QDPIO::cout << "ord_norm2x_cdotxy_kernel_sse.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
}
