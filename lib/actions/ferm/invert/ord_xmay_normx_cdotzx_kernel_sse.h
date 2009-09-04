inline
void ord_xmay_normx_cdotzx_kernel(int lo, int hi, int my_id, ord_xmay_normx_cdotzx_arg *a) 
{
  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi-lo);

  REAL32* x_ptr=&(a->x_ptr[low]);
  REAL32* y_ptr=&(a->y_ptr[low]);
  REAL32* z_ptr=&(a->z_ptr[low]);
  REAL32 a_re = a->a_re;
  REAL32 a_im = a->a_im;
  

  __m128d sum = _mm_set_pd((double)0,(double)0);
  __m128d dotprod = _mm_set_pd((double)0,(double)0);
  __m128d mask = _mm_set_pd((double)-1, (double)1);
  
  __m128 av_re = _mm_set_ps( a_re, a_re,  a_re, a_re);
  __m128 av_im = _mm_set_ps(-a_im, a_im, -a_im, a_im);
  
  if( len % 4 == 0 ) { 
    for(int count = 0; count < len; count+=4) { 
      __m128 xv = _mm_load_ps(&x_ptr[count]);
      __m128 yv = _mm_load_ps(&y_ptr[count]);
      __m128 zv = _mm_load_ps(&z_ptr[count]);
      
      
      //x_ptr[count] -= a_re*y_ptr[count];
      //x_ptr[count+1] -= a_re*y_ptr[count+1];
      //x_ptr[count+2] -= a_re*y_ptr[count+2];
      //x_ptr[count+3] -= a_re*y_ptr[count+3];
      
      __m128 t1 = _mm_mul_ps(av_re, yv);
      xv = _mm_sub_ps(xv,t1);
      
      // x_ptr[count] += a_im*y_ptr[count+1];
      // x_ptr[count+1] -= a_im*y_ptr[count];
      // x_ptr[count+2] += a_im*y_ptr[count+3];
      // x_ptr[count+3] -= a_im*y_ptr[count+2];
      
      t1 = _mm_shuffle_ps(yv,yv,0xb1);
      __m128 t2 = _mm_mul_ps(av_im, t1);
      xv  = _mm_add_ps(xv, t2);
      
      _mm_store_ps(&x_ptr[count], xv);
      
      
      __m128d xlow,xhi;
      __m128d zlow,zhi;
      
      xlow = _mm_cvtps_pd(xv);
      
      // Swap low and high pairs 
      xv = _mm_shuffle_ps( xv,xv, 0x4e);
      
      // Take hi part of x into xhi
      xhi = _mm_cvtps_pd(xv);
      
      // Now samething with y
      zlow = _mm_cvtps_pd(zv);
      zv = _mm_shuffle_ps( zv,zv, 0x4e);
      zhi = _mm_cvtps_pd(zv);
      
      //norm_array[0] += x_ptr[count]*x_ptr[count];
      //    norm_array[0] += x_ptr[count+1]*x_ptr[count+1];
      //norm_array[0] += x_ptr[count+2]*x_ptr[count+2];
      //norm_array[0] += x_ptr[count+3]*x_ptr[count+3];
      
      sum = _mm_add_pd(sum,_mm_mul_pd(xlow,xlow));
      sum = _mm_add_pd(sum,_mm_mul_pd(xhi,xhi));
      
      
      //    dotprod[0] += z_ptr[count]*x_ptr[count];
      //    dotprod[1] -= z_ptr[count+1]*x_ptr[count];
      
      __m128d t1d = _mm_shuffle_pd(xlow, xlow, 0x0);
      __m128d t2d = _mm_mul_pd(mask, t1d);
      __m128d t3d = _mm_mul_pd(zlow, t2d);
      
      dotprod = _mm_add_pd(dotprod, t3d);
      
      
      //dotprod[0] += z_ptr[count+1]*x_ptr[count+1];
      //dotprod[1] += z_ptr[count]*x_ptr[count+1];    
      t1d = _mm_shuffle_pd(xlow, xlow, 0x3);
      t2d = _mm_shuffle_pd(zlow, zlow, 0x1);
      t3d = _mm_mul_pd(t1d,t2d);
      dotprod = _mm_add_pd(dotprod,t3d);
      
      
      // dotprod[0] += z_ptr[count+2]*x_ptr[count+2];
      // dotprod[1] -= z_ptr[count+3]*x_ptr[count+2];
      t1d = _mm_shuffle_pd(xhi,xhi,0x0);
      t2d = _mm_mul_pd(mask,t1d);
      t3d = _mm_mul_pd(zhi,t2d);
      dotprod = _mm_add_pd(dotprod,t3d);
      
      // dotprod[0] += z_ptr[count+3]*x_ptr[count+3];
      // dotprod[1] += z_ptr[count+2]*x_ptr[count+3];
      t1d = _mm_shuffle_pd(xhi,xhi, 0x3);
      t2d = _mm_shuffle_pd(zhi,zhi,0x1);
      t3d = _mm_mul_pd(t1d, t2d);
      dotprod = _mm_add_pd(dotprod,t3d);
      
    } 
    a->norm_space[3*my_id]=((double *)&sum)[0] + ((double *)&sum)[1];
    a->norm_space[3*my_id+1]=((double *)&dotprod)[0];
    a->norm_space[3*my_id+2]=((double *)&dotprod)[1];
  }
  else { 
    QDPIO::cout << " ord_xmay_normx_cdotzx_kernel_sse.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
}

