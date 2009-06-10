#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>


void ord_xymz_normx_kernel(int lo, int hi, int my_id, ord_xymz_normx_arg* a)
{
  REAL64* x_ptr;
  REAL64* y_ptr;
  REAL64* z_ptr;
  REAL64 norm=0;
  
  x_ptr = &(a->x_ptr[lo]);
  y_ptr = &(a->y_ptr[lo]);
  z_ptr = &(a->z_ptr[lo]);
  
  int len = hi-lo;

  __m128d norm_vec1 = _mm_set_pd((double)0,(double)0);
  __m128d norm_vec2 = _mm_set_pd((double)0,(double)0);
  __m128d norm_vec3 = _mm_set_pd((double)0,(double)0);
  __m128d norm_vec4 = _mm_set_pd((double)0,(double)0);

  for(int count = 0; count < len; count+=8) { 
    
    __m128d xvec1,xvec2,xvec3,xvec4;
    __m128d yvec1,yvec2,yvec3,yvec4;
    __m128d zvec1,zvec2,zvec3,zvec4;

    yvec1 = _mm_load_pd(&y_ptr[count]);
    zvec1 = _mm_load_pd(&z_ptr[count]);

    xvec1 = _mm_sub_pd(yvec1,zvec1);
    _mm_store_pd(&x_ptr[count], xvec1);


    yvec2 = _mm_load_pd(&y_ptr[count+2]);
    zvec2 = _mm_load_pd(&z_ptr[count+2]);
    xvec2 = _mm_sub_pd(yvec2,zvec2);
    _mm_store_pd(&x_ptr[count+2], xvec2); 

    yvec3 = _mm_load_pd(&y_ptr[count+4]);
    zvec3 = _mm_load_pd(&z_ptr[count+4]);
    xvec3 = _mm_sub_pd(yvec3,zvec3);
    _mm_store_pd(&x_ptr[count+4], xvec3);

    yvec4 = _mm_load_pd(&y_ptr[count+6]);
    zvec4 = _mm_load_pd(&z_ptr[count+6]);
    xvec4 = _mm_sub_pd(yvec4,zvec4);
    _mm_store_pd(&x_ptr[count+6], xvec4);    

    yvec1 = _mm_mul_pd(xvec1,xvec1);
    norm_vec1 = _mm_add_pd(norm_vec1,yvec1);

    yvec2 = _mm_mul_pd(xvec2,xvec2);
    norm_vec2 = _mm_add_pd(norm_vec2,yvec2);

    yvec3 = _mm_mul_pd(xvec3,xvec3);
    norm_vec3 = _mm_add_pd(norm_vec3,yvec3);

    yvec4 = _mm_mul_pd(xvec4,xvec4);
    norm_vec4 = _mm_add_pd(norm_vec4,yvec4);

  }
  norm_vec1 = _mm_add_pd(norm_vec1, norm_vec2);
  norm_vec3 = _mm_add_pd(norm_vec3, norm_vec4);
  norm_vec1 = _mm_add_pd(norm_vec1, norm_vec3);

  a->norm_ptr[my_id] = ((double *)&norm_vec1)[0] + ((double *)&norm_vec1)[1];


}
