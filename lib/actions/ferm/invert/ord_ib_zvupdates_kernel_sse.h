// 32 BIT Version: Use vector length of 4 for easy vectorization.
// This is guaranteed good for LatticeDiracFermions

inline
void ord_ib_zvupdates_kernel_real32(int lo, int hi, int my_id, ib_zvupdates_arg<REAL32>* a)
{

  int atom=a->atom;
  int low = atom*lo;
  int len = atom*(hi - lo);

  REAL32* r = &(a->r_ptr[low]);
  REAL32* z = &(a->z_ptr[low]);
  REAL32* v = &(a->v_ptr[low]);
  REAL32* u = &(a->u_ptr[low]);
  REAL32* q = &(a->q_ptr[low]);

  REAL32 a_re = a->alpha_re;
  REAL32 a_im = a->alpha_im;

  REAL32 arb_re = a->alpha_rat_beta_re;
  REAL32 arb_im = a->alpha_rat_beta_im;

  REAL32 ad_re = a->alpha_delta_re;
  REAL32 ad_im = a->alpha_delta_im;
  
  REAL32 b_re = a->beta_re;
  REAL32 b_im = a->beta_im;

  REAL32 d_re = a->delta_re;
  REAL32 d_im = a->delta_im;
  



  __m128 ztmp, vtmp;

  __m128 zvec, rvec, vvec, uvec, qvec, tmpshuf1, tmpshuf2, tmpshuf3,tmpshuf4;
  const __m128 arb_re_vec = _mm_set_ps(arb_re,arb_re,arb_re,arb_re);
  const __m128 arb_im_vec = _mm_set_ps(arb_im,-arb_im,arb_im,-arb_im);
  const __m128 a_re_vec = _mm_set_ps(a_re,a_re,a_re,a_re);
  const __m128 a_im_vec = _mm_set_ps(a_im,-a_im,a_im,-a_im);
  const __m128 ad_re_vec = _mm_set_ps(ad_re,ad_re,ad_re,ad_re);
  const __m128 ad_im_vec = _mm_set_ps(ad_im,-ad_im,ad_im,-ad_im);
  const __m128 b_re_vec = _mm_set_ps(b_re,b_re,b_re,b_re);
  const __m128 b_im_vec = _mm_set_ps(b_im,-b_im,b_im,-b_im);
  const __m128 d_re_vec = _mm_set_ps(d_re,d_re,d_re,d_re);
  const __m128 d_im_vec = _mm_set_ps(d_im,-d_im,d_im,-d_im);

  if( len % 4 == 0 ) { 
    for(int count = 0; count < len; count+=4) { 
      
      /* ztmp = { z[count], z[count+1], z[count+2], z[count+3] } */
      ztmp = _mm_load_ps(&z[count]);
      
      /* rvec = { r[count], r[count+1], r[count+2], r[count+3] } */
      rvec = _mm_load_ps(&r[count]);
      
      /* vtmp = { v[count], v[count+1], v[count+2], v[count+3] } */
      vtmp = _mm_load_ps(&v[count]);
      
      /* uvec = { u[count], u[count+1], u[count+2], u[count+3] } */
      uvec = _mm_load_ps(&u[count]);
      
      /* qvec = { q[count], q[count+1], q[count+2], q[count+3] } */
      qvec = _mm_load_ps(&q[count]);
      
      /* tmpshuf1 =  { z[count+1], z[count], z[count+3], z[count+2] } */
      tmpshuf1 = _mm_shuffle_ps(ztmp,ztmp, 0xb1);
      
      /* tmpshuf2 =  { r[count+1], r[count], r[count+3], r[count+2] } */
      tmpshuf2 = _mm_shuffle_ps(rvec,rvec, 0xb1);
      
      /* tmpshuf3 =  { v[count+1], v[count], v[count+3], v[count+2] } */
      tmpshuf3 = _mm_shuffle_ps(vtmp,vtmp, 0xb1);
      
      /* tmpshuf3 =  { q[count+1], q[count], q[count+3], q[count+2] } */
      tmpshuf4 = _mm_shuffle_ps(qvec,qvec, 0xb1);
      
      /*******    z = (alpha_n/alpha_n-1)*beta z   *******/
      /* 
       * arb =(alpha_n/alpha_n-1)*beta
       *
       * ztmp     = { z[count], z[count+1], z[count+2], z[count+3] }
       * tmpshuf1 = { z[count+1], z[count], z[count+3], z[count+2] }
       * arb_re_vec = { arb_re, arb_re, arb_re, arb_re }
       * arb_im_vec = {-arb_im, arb_im,-arb_im, arb_im }
       *
       * z[count]    = arb_re * ztmp[0] - arb_im * ztmp[1];
       * z[count+1]  = arb_re * ztmp[1] + arb_im * ztmp[0];  
       * z[count+2]  = arb_re * ztmp[2] - arb_im * ztmp[3];
       * z[count+3]  = arb_re * ztmp[3] + arb_im * ztmp[2];  
       */
      zvec = _mm_mul_ps(arb_re_vec,ztmp);
      zvec = _mm_add_ps(zvec, _mm_mul_ps(arb_im_vec,tmpshuf1));
      
      /********** z += alpha * r ***********/
      /*
       * a_re_vec = [ a_re,a_re,a_re,a_re ] 
       
       z[count  ] += a_re * r[count];
       z[count+1] += a_re * r[count+1];
       z[count+2] += a_re * r[count+2];
       z[count+3] += a_re * r[count+3];
      */
      zvec = _mm_add_ps(zvec, _mm_mul_ps(a_re_vec,rvec));
      
      /*
       * a_im_vec = [ -a_im, a_im, -a_im, a_im ] 
       z[count  ] -= a_im * r[count+1];
       z[count+1] += a_im * r[count];
       z[count+2] -= a_im * r[count+3];
       z[count+3] += a_im * r[count+2];
      */
      
      zvec = _mm_add_ps(zvec, _mm_mul_ps(a_im_vec,tmpshuf2));
      
      /***** z -= (alpha*delta) v **********/
      /*
	z[count  ] -= ad_re * v[count] ;
	z[count+1] -= ad_re * v[count+1];
	z[count+2] -= ad_re * v[count+2] ;
	z[count+3] -= ad_re * v[count+3];
      */
      
      zvec = _mm_sub_ps(zvec, _mm_mul_ps(ad_re_vec, vtmp));
      
      /*
	z[count  ] += ad_im * v[count+1];
	z[count+1] -= ad_im * v[count];
	z[count+2] += ad_im * v[count+3];
	z[count+3] -= ad_im * v[count+2];
      */
      zvec = _mm_sub_ps(zvec, _mm_mul_ps(ad_im_vec, tmpshuf3));
      
      _mm_store_ps(&z[count],zvec);
      
      /*************  v = u + b v ****************/
      /*
	v[count]   = u[count]   + b_re*vtmp[0] - b_im*vtmp[1];
	v[count+1] = u[count+1] + b_re*vtmp[1] + b_im*vtmp[0];
	v[count+2] = u[count+2] + b_re*vtmp[2] - b_im*vtmp[3];
	v[count+3] = u[count+3] + b_re*vtmp[3] + b_im*vtmp[2];
      */
      vvec = _mm_add_ps(uvec, _mm_mul_ps(b_re_vec,vtmp));
      vvec = _mm_add_ps(vvec, _mm_mul_ps(b_im_vec,tmpshuf3));
      
      /***************** v -= d*q ******************/
      /*
	v[count]   -= d_re*q[count];
	v[count+1] -= d_re*q[count+1];
	v[count+2] -= d_re*q[count+2];
	v[count+3] -= d_re*q[count+3];
      */
      vvec = _mm_sub_ps(vvec, _mm_mul_ps(d_re_vec,qvec));
      /*
	v[count]   += d_im*q[count+1];
	v[count+1] -= d_im*q[count];
	v[count+2] += d_im*q[count+3];
	v[count+3] -= d_im*q[count+2];
      */
      vvec = _mm_sub_ps(vvec, _mm_mul_ps(d_im_vec,tmpshuf4));
      
      _mm_store_ps(&v[count],vvec);
      
      
    }  
  }
  else { 
    QDPIO::cout << "ord_ib_zvupdates_sse.h: len not divisible by 4" << endl;
    QDP_abort(1);
    
    
  }
}

// 64 BIT Version: Use vector length of 2 for easy vectorization.
// This is guaranteed good for LatticeDiracFermions


inline
void ord_ib_zvupdates_kernel_real64(int lo, int hi, int my_id, ib_zvupdates_arg<REAL64>* a)
{

  int atom=a->atom;
  int low = atom*lo;
  int len = atom*(hi - lo);

  REAL64* r = &(a->r_ptr[low]);
  REAL64* z = &(a->z_ptr[low]);
  REAL64* v = &(a->v_ptr[low]);
  REAL64* u = &(a->u_ptr[low]);
  REAL64* q = &(a->q_ptr[low]);

  REAL64 a_re = a->alpha_re;
  REAL64 a_im = a->alpha_im;

  REAL64 arb_re = a->alpha_rat_beta_re;
  REAL64 arb_im = a->alpha_rat_beta_im;

  REAL64 ad_re = a->alpha_delta_re;
  REAL64 ad_im = a->alpha_delta_im;
  
  REAL64 b_re = a->beta_re;
  REAL64 b_im = a->beta_im;

  REAL64 d_re = a->delta_re;
  REAL64 d_im = a->delta_im;
  

  __m128d ztmp, vtmp;

  __m128d zvec, rvec, vvec, uvec, qvec, tmpshuf1, tmpshuf2, tmpshuf3,tmpshuf4;
  const __m128d arb_re_vec = _mm_set_pd(arb_re,arb_re);
  const __m128d arb_im_vec = _mm_set_pd(arb_im,-arb_im);
  const __m128d a_re_vec = _mm_set_pd(a_re,a_re);
  const __m128d a_im_vec = _mm_set_pd(a_im,-a_im);
  const __m128d ad_re_vec = _mm_set_pd(ad_re,ad_re);
  const __m128d ad_im_vec = _mm_set_pd(ad_im,-ad_im);
  const __m128d b_re_vec = _mm_set_pd(b_re,b_re);
  const __m128d b_im_vec = _mm_set_pd(b_im,-b_im);
  const __m128d d_re_vec = _mm_set_pd(d_re,d_re);
  const __m128d d_im_vec = _mm_set_pd(d_im,-d_im);

  if( len % 2 == 0 ) { 
  
    for(int count = 0; count < len; count+=2) { 
      ztmp = _mm_load_pd(&z[count]);
      vtmp = _mm_load_pd(&v[count]);
      rvec = _mm_load_pd(&r[count]);
      uvec = _mm_load_pd(&u[count]);
      qvec = _mm_load_pd(&q[count]);
      
      
      tmpshuf1= _mm_shuffle_pd(ztmp,ztmp,0x1);
      tmpshuf2 = _mm_shuffle_pd(rvec,rvec,0x1);
      tmpshuf3 = _mm_shuffle_pd(vtmp,vtmp,0x1);
      tmpshuf4 = _mm_shuffle_pd(qvec,qvec,0x1);
      
      /* z = (alpha_n/alpha_n-1)*beta z */
      /*
	ztmp[0] = z[count];
	ztmp[1] = z[count+1];
	z[count]    = arb_re * ztmp[0] - arb_im * ztmp[1];
	z[count+1]  = arb_re * ztmp[1] + arb_im * ztmp[0];  
	
      */
      
      zvec = _mm_mul_pd(arb_re_vec, ztmp);
      zvec = _mm_add_pd(zvec, _mm_mul_pd(arb_im_vec, tmpshuf1));
      
      
      
      /* z += alpha*r */
      /*
	z[count  ] += a_re * r[count];
	z[count+1] += a_re * r[count+1];
	
	z[count  ] -= a_im * r[count+1];
	z[count+1] += a_im * r[count];
	
	
      */
      zvec = _mm_add_pd(zvec,_mm_mul_pd(a_re_vec,rvec));
      zvec = _mm_add_pd(zvec,_mm_mul_pd(a_im_vec,tmpshuf2));
      
      /* z -= alpha*delta*v */
      /*
	z[count  ] -= ad_re * v[count] ;
	z[count+1] -= ad_re * v[count+1];
	
	z[count  ] += ad_im * v[count+1];
	z[count+1] -= ad_im * v[count];
      */
      zvec = _mm_sub_pd(zvec, _mm_mul_pd(ad_re_vec,vtmp));
      zvec = _mm_sub_pd(zvec, _mm_mul_pd(ad_im_vec, tmpshuf3));
      _mm_store_pd(&z[count], zvec);
      
      
      
      
      /* v = u + b*v */
      /*
	v[count]   = u[count]   + b_re*vtmp[0] - b_im*vtmp[1];
	v[count+1] = u[count+1] + b_re*vtmp[1] + b_im*vtmp[0];
	
      */
      vvec = _mm_add_pd( uvec, _mm_mul_pd(b_re_vec,vtmp));
      vvec = _mm_add_pd( vvec, _mm_mul_pd(b_im_vec,tmpshuf3));
      
      /*
	v[count]   -= d_re*q[count];
	v[count+1] -= d_re*q[count+1];
      */
      vvec = _mm_sub_pd( vvec, _mm_mul_pd(d_re_vec,qvec));
      
      /*
	v[count]   += d_im*q[count+1];
	v[count+1] -= d_im*q[count];
      */
      vvec = _mm_sub_pd( vvec, _mm_mul_pd(d_im_vec, tmpshuf4));
      _mm_store_pd(&v[count],vvec);
    }
  }
  else { 
    QDPIO::cout << "ord_ib_zvupdates_sse.h: len not divisible by 2" << endl;
    QDP_abort(1);
  }
  
}

