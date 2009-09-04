#include <xmmintrin.h>
#include <emmintrin.h>

inline
void ord_ib_stupdates_kernel_real32(int lo, int hi, int my_id, ib_stupdate_arg<REAL32>* a)
{
  REAL32 a_r = a->a_r;
  REAL32 a_i = a->a_i;
  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi-lo);

  REAL32* r = &(a->r[low]);
  REAL32* u = &(a->u[low]);
  REAL32* v = &(a->v[low]);
  REAL32* q = &(a->q[low]);
  REAL32* r0 = &(a->r0[low]);
  REAL32* f0 = &(a->f0[low]);
  REAL32* s = &(a->s[low]);
  REAL32* t = &(a->t[low]);
  REAL64* norm_array = &(a->norm_space[12*my_id]);

  // Caller zeroed norm_space
  
  __m128 svec, rvec, vvec, tmpshuf1;
  __m128 qvec, uvec, tvec, tmpshuf2; 
  __m128 ar_vec = _mm_set_ps(a_r,a_r,a_r,a_r);
  __m128 ai_vec = _mm_set_ps(a_i,-a_i,a_i,-a_i);
  __m128 lvec;

  __m128d dotprod;
  __m128d llo,lhi;
  __m128d slo,shi;
  __m128d tlo,thi;
  __m128d qlo,qhi;
  __m128d r0lo,r0hi;
  __m128d f0lo,f0hi;
  __m128d mask = _mm_set_pd((double)-1,(double)1);
  __m128d t1,t2;

  if( len % 4 == 0 ) { 
    for(int count = 0; count < len; count+=4) { 
      
      vvec = _mm_load_ps(&v[count]);
      qvec = _mm_load_ps(&q[count]);
      
      rvec = _mm_load_ps(&r[count]);
      uvec = _mm_load_ps(&u[count]);
      tmpshuf1 = _mm_shuffle_ps(vvec,vvec,0xb1);
      tmpshuf2 = _mm_shuffle_ps(qvec, qvec,0xb1);
      
      // First need s = r - alpha*v    
      /*
	s[count  ] = r[count  ] - a_r*v[count  ];
	s[count+1] = r[count+1] - a_r*v[count+1];
	s[count+2] = r[count+2] - a_r*v[count+2];
	s[count+3] = r[count+3] - a_r*v[count+3];
      */
      svec = _mm_sub_ps(rvec, _mm_mul_ps(ar_vec,vvec));
      
      /*
	s[count] += a_i*v[count+1];
	s[count+1] -= a_i*v[count];
	s[count+2] += a_i*v[count+3];
	s[count+3] -= a_i*v[count+2];
      */
      svec = _mm_sub_ps(svec, _mm_mul_ps(ai_vec, tmpshuf1));
      
      // Second want t = u - alqha q
      tvec = _mm_sub_ps(uvec, _mm_mul_ps(ar_vec,qvec));
      tvec = _mm_sub_ps(tvec, _mm_mul_ps(ai_vec,tmpshuf2));
      
      _mm_store_ps(&s[count],svec);
      _mm_store_ps(&t[count],tvec);
      
      
      slo = _mm_cvtps_pd(svec);
      svec = _mm_shuffle_ps(svec,svec,0x4e);
      shi = _mm_cvtps_pd(svec);
      
      qlo = _mm_cvtps_pd(qvec);
      qvec = _mm_shuffle_ps(qvec,qvec,0x4e);
      qhi = _mm_cvtps_pd(qvec);
      
      tlo = _mm_cvtps_pd(tvec);
      tvec = _mm_shuffle_ps(tvec,tvec,0x4e);
      thi = _mm_cvtps_pd(tvec);
      
#if 0
      // ** phi=(r0,s)
      norm_array[0] += r0[count]*s[count];
      norm_array[0] += r0[count+1]*s[count+1];
      norm_array[0] += r0[count+2]*s[count+2];
      norm_array[0] += r0[count+3]*s[count+3];
      
      norm_array[1] += r0[count]*s[count+1];
      norm_array[1] -= r0[count+1]*s[count];
      norm_array[1] += r0[count+2]*s[count+3];
      norm_array[1] -= r0[count+3]*s[count+2];
#else
      
      // ** phi=(r0,s)
      
      /* Left vector = r0  */
      lvec = _mm_load_ps(&r0[count]);
      
      /* Load dotprod accumulated so far */
      dotprod = _mm_load_pd(&norm_array[0]);
      
      /* Turn left vector into double precision */
      r0lo = _mm_cvtps_pd(lvec);
      lvec = _mm_shuffle_ps(lvec,lvec,0x4e);
      r0hi = _mm_cvtps_pd(lvec);
      
      
      
      /* So now:  llo = [ r0[count+1], r0[count]   ]
	 lhi = [ r0[count+3], r0[count+2] ]  
	 rlo = [ s[count+1],  s[count]    ]
	 rhi = [ s[count+3],  s[count+2]  ]
	 
	 na_vec holds the current partial sum */
      
      /*    dotprod[0] += r0lo[0]*slo[0];         */
      /*    dotprod[1] -= r0lo[1]*slo[0];       */
      t1 = _mm_shuffle_pd(slo,slo,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(r0lo, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += r0lo[1]*slo[1]   */
      /*    dotprod[1] += r0lo[0]*slo[1];    */
      t1 = _mm_shuffle_pd(slo,slo,0x3);
      t2 = _mm_shuffle_pd(r0lo,r0lo,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      /*    dotprod[0] += r0hi[0]*shi[0];         */
      /*    dotprod[1] -= r0hi[1]*shi[0];       */
      t1 = _mm_shuffle_pd(shi,shi,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(r0hi, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += r0hi[1]*shi[1]   */
      /*    dotprod[1] += r0hi[0]*shi[1];    */
      t1 = _mm_shuffle_pd(shi,shi,0x3);
      t2 = _mm_shuffle_pd(r0hi,r0hi,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      _mm_store_pd(&norm_array[0],dotprod);
      
#endif
      
      // ** phi=(f0,s) 
      
      /* rlo and rhi are still svec as before
	 no need to turn them into double precision */
#if 0
      norm_array[2] += f0[count]*s[count];
      norm_array[2] += f0[count+1]*s[count+1];
      norm_array[2] += f0[count+2]*s[count+2];
      norm_array[2] += f0[count+3]*s[count+3];
      
      norm_array[3] += f0[count]*s[count+1];
      norm_array[3] -= f0[count+1]*s[count];
      norm_array[3] += f0[count+2]*s[count+3];
      norm_array[3] -= f0[count+3]*s[count+2];
#else 
      
      // ** phi=(f0,s)
      
      /* Left vector = f0  */
      lvec = _mm_load_ps(&f0[count]);
      
      /* Load dotprod accumulated so far */
      dotprod = _mm_load_pd(&norm_array[2]);
      
      /* Turn left vector into double precision */
      f0lo = _mm_cvtps_pd(lvec);
      lvec = _mm_shuffle_ps(lvec,lvec,0x4e);
      f0hi = _mm_cvtps_pd(lvec);
      
      
      
      /* So now:  llo = [ f0[count+1], f0[count]   ]
	 lhi = [ f0[count+3], f0[count+2] ]  
	 rlo = [ s[count+1],  s[count]    ]
	 rhi = [ s[count+3],  s[count+2]  ]
	 
	 dotprod holds the current partial sum */
      
      /*    dotprod[0] += f0lo[0]*slo[0];         */
      /*    dotprod[1] -= f0lo[1]*slo[0];       */
      t1 = _mm_shuffle_pd(slo,slo,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(f0lo, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += f0lo[1]*slo[1]   */
      /*    dotprod[1] += f0lo[0]*slo[1];    */
      t1 = _mm_shuffle_pd(slo,slo,0x3);
      t2 = _mm_shuffle_pd(f0lo,f0lo,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      /*    dotprod[0] += f0hi[0]*shi[0];         */
      /*    dotprod[1] -= f0hi[1]*shi[0];       */
      t1 = _mm_shuffle_pd(shi,shi,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(f0hi, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += f0hi[1]*shi[1]   */
      /*    dotprod[1] += f0hi[0]*shi[1];    */
      t1 = _mm_shuffle_pd(shi,shi,0x3);
      t2 = _mm_shuffle_pd(f0hi,f0hi,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      _mm_store_pd(&norm_array[2],dotprod);
#endif
      
      
#if 0
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
#else 
      // ** pi=(r0,q)
      
      /* Load dotprod accumulated so far */
      dotprod = _mm_load_pd(&norm_array[4]);
      
      /* dotprod holds the current partial sum */
      
      /*    dotprod[0] += r0lo[0]*qlo[0];         */
      /*    dotprod[1] -= r0lo[1]*qlo[0];       */
      t1 = _mm_shuffle_pd(qlo,qlo,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(r0lo, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += r0lo[1]*qlo[1]   */
      /*    dotprod[1] += r0lo[0]*qlo[1];    */
      t1 = _mm_shuffle_pd(qlo,qlo,0x3);
      t2 = _mm_shuffle_pd(r0lo,r0lo,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      /*    dotprod[0] += r0hi[0]*qhi[0];         */
      /*    dotprod[1] -= r0hi[1]*qhi[0];       */
      t1 = _mm_shuffle_pd(qhi,qhi,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(r0hi, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += r0hi[1]*qhi[1]   */
      /*    dotprod[1] += r0hi[0]*qhi[1];    */
      t1 = _mm_shuffle_pd(qhi,qhi,0x3);
      t2 = _mm_shuffle_pd(r0hi,r0hi,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      _mm_store_pd(&norm_array[4],dotprod);
#endif
      
      
#if 0
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
      
#else
      // ** eta=(f0,t)
      
      /* Load dotprod accumulated so far */
      dotprod = _mm_load_pd(&norm_array[6]);
      
      /* dotprod holds the current partial sum */
      
      /*    dotprod[0] += f0lo[0]*tlo[0];         */
      /*    dotprod[1] -= f0lo[1]*tlo[0];       */
      t1 = _mm_shuffle_pd(tlo,tlo,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(f0lo, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += f0lo[1]*tlo[1]   */
      /*    dotprod[1] += f0lo[0]*tlo[1];    */
      t1 = _mm_shuffle_pd(tlo,tlo,0x3);
      t2 = _mm_shuffle_pd(f0lo,f0lo,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      /*    dotprod[0] += f0hi[0]*thi[0];         */
      /*    dotprod[1] -= f0hi[1]*thi[0];       */
      t1 = _mm_shuffle_pd(thi,thi,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(f0hi, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += f0hi[1]*thi[1]   */
      /*    dotprod[1] += f0hi[0]*thi[1];    */
      t1 = _mm_shuffle_pd(thi,thi,0x3);
      t2 = _mm_shuffle_pd(f0hi,f0hi,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      _mm_store_pd(&norm_array[6],dotprod);
      
#endif
      
#if 0
      // ** theta=(t,s) 
      norm_array[8] += t[count]*s[count];
      norm_array[8] += t[count+1]*s[count+1];
      norm_array[8] += t[count+2]*s[count+2];
      norm_array[8] += t[count+3]*s[count+3];
      
      norm_array[9] += t[count]*s[count+1];
      norm_array[9] -= t[count+1]*s[count];
      norm_array[9] += t[count+2]*s[count+3];
      norm_array[9] -= t[count+3]*s[count+2];
#else
      // ** theta=(t,s)
      
      /* Load dotprod accumulated so far */
      dotprod = _mm_load_pd(&norm_array[8]);
      
      /* dotprod holds the current partial sum */
      
      /*    dotprod[0] += tlo[0]*slo[0];         */
      /*    dotprod[1] -= tlo[1]*slo[0];       */
      t1 = _mm_shuffle_pd(slo,slo,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(tlo, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += tlo[1]*slo[1]   */
      /*    dotprod[1] += tlo[0]*slo[1];    */
      t1 = _mm_shuffle_pd(slo,slo,0x3);
      t2 = _mm_shuffle_pd(tlo,tlo,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      /*    dotprod[0] += thi[0]*shi[0];         */
      /*    dotprod[1] -= thi[1]*shi[0];       */
      t1 = _mm_shuffle_pd(shi,shi,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(thi, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += thi[1]*shi[1]   */
      /*    dotprod[1] += thi[0]*shi[1];    */
      t1 = _mm_shuffle_pd(shi,shi,0x3);
      t2 = _mm_shuffle_pd(thi,thi,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      _mm_store_pd(&norm_array[8],dotprod);
#endif
      
#if 0
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
#else
      dotprod = _mm_load_pd(&norm_array[10]);
      
      r0lo = _mm_cvtps_pd(rvec);
      rvec = _mm_shuffle_ps(rvec,rvec,0x4e);
      r0hi = _mm_cvtps_pd(rvec);
      
      // [ti,tr] & [ri,rr] -> [tr, rr] 
      f0lo = _mm_shuffle_pd(tlo,r0lo, 0x0);
      
      // [ti,tr] & [ri,rr] -> [ti, ri] 
      f0hi = _mm_shuffle_pd(tlo,r0lo, 0x3);
      
      // [ tnorm, rnorm ] += [ tr^2, rr^2 ]
      // [ tnorm, rnorm ] += [ ti^2, ri^2 ]
      dotprod = _mm_add_pd(dotprod,_mm_mul_pd(f0lo,f0lo));
      dotprod = _mm_add_pd(dotprod,_mm_mul_pd(f0hi,f0hi));
      
      // [ti,tr] & [ri,rr] -> [tr, rr] 
      f0lo = _mm_shuffle_pd(thi,r0hi, 0x0);
      
      // [ti,tr] & [ri,rr] -> [ti, ri] 
      f0hi = _mm_shuffle_pd(thi,r0hi, 0x3);
      
      // [ tnorm, rnorm ] += [ tr^2, rr^2 ]
      // [ tnorm, rnorm ] += [ ti^2, ri^2 ]
      dotprod = _mm_add_pd(dotprod,_mm_mul_pd(f0lo,f0lo));
      dotprod = _mm_add_pd(dotprod,_mm_mul_pd(f0hi,f0hi));
      
      _mm_store_pd(&norm_array[10],dotprod);
#endif
    }
  }
  else { 
    QDPIO::cout << "ord_ib_stupdates_kernel_sse.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
}

inline
void ord_ib_stupdates_kernel_real64(int lo, int hi, int my_id, ib_stupdate_arg<REAL64>* a)
{
  REAL64 a_r = a->a_r;
  REAL64 a_i = a->a_i;
  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi-lo);

  REAL64* r = &(a->r[low]);
  REAL64* u = &(a->u[low]);
  REAL64* v = &(a->v[low]);
  REAL64* q = &(a->q[low]);
  REAL64* r0 = &(a->r0[low]);
  REAL64* f0 = &(a->f0[low]);
  REAL64* s = &(a->s[low]);
  REAL64* t = &(a->t[low]);
  REAL64* norm_array = &(a->norm_space[12*my_id]);
  
  __m128d svec, rvec, vvec;
  __m128d qvec, uvec, tvec;
  __m128d ar_vec = _mm_set_pd(a_r,a_r);
  __m128d ai_vec = _mm_set_pd(a_i,-a_i);

  __m128d dotprod;
  __m128d sv;
  __m128d tv;
  __m128d qv;
  __m128d r0v;
  __m128d f0v;
  __m128d mask = _mm_set_pd((double)-1,(double)1);
  __m128d t1,t2;

  // Caller zeroed norm_space
  if( len % 2 == 0){
    for(int count = 0; count < len; count+=2) { 
      // First need s = r - alpha*v
      rvec = _mm_load_pd(&r[count]);
      vvec = _mm_load_pd(&v[count]);
      t1 = _mm_shuffle_pd(vvec,vvec,0x1);
      
      uvec = _mm_load_pd(&u[count]);
      qvec = _mm_load_pd(&q[count]);
      t2 = _mm_shuffle_pd(qvec,qvec,0x1);
      
      svec = _mm_sub_pd(rvec, _mm_mul_pd(ar_vec, vvec));
      r0v = _mm_load_pd(&r0[count]);   
      svec = _mm_sub_pd(svec, _mm_mul_pd(ai_vec, t1));
      
      
      tvec = _mm_sub_pd(uvec, _mm_mul_pd(ar_vec, qvec));
      f0v = _mm_load_pd(&f0[count]);
      tvec = _mm_sub_pd(tvec, _mm_mul_pd(ai_vec, t2));
      
      _mm_store_pd(&s[count],svec);
      _mm_store_pd(&t[count],tvec);
      
      
      
      // Now inner products with s. 
      // ** phi=(r0,s) 
#if 0
      norm_array[0] += r0[count]*s[count];
      norm_array[0] += r0[count+1]*s[count+1];
      
      norm_array[1] += r0[count]*s[count+1];
      norm_array[1] -= r0[count+1]*s[count];
#else
      
      dotprod = _mm_load_pd(&norm_array[0]);
      
      /* dotprod holds the current partial sum */
      
      /*    dotprod[0] += r0[0]*s[0];         */
      /*    dotprod[1] -= r0[1]*s[0];       */
      t1 = _mm_shuffle_pd(svec,svec,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(r0v, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += r0[1]*s[1]   */
      /*    dotprod[1] += r0[0]*s[1];    */
      t1 = _mm_shuffle_pd(svec,svec,0x3);
      t2 = _mm_shuffle_pd(r0v,r0v,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      _mm_store_pd(&norm_array[0],dotprod);
#endif
      
#if 0
      // ** phi=(f0,s) 
      norm_array[2] += f0[count]*s[count];
      norm_array[2] += f0[count+1]*s[count+1];
      
      norm_array[3] += f0[count]*s[count+1];
      norm_array[3] -= f0[count+1]*s[count];
#else
      
      dotprod = _mm_load_pd(&norm_array[2]);
      
      /* dotprod holds the current partial sum */
      
      /*    dotprod[0] += r0[0]*s[0];         */
      /*    dotprod[1] -= r0[1]*s[0];       */
      t1 = _mm_shuffle_pd(svec,svec,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(f0v, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += r0[1]*s[1]   */
      /*    dotprod[1] += r0[0]*s[1];    */
      t1 = _mm_shuffle_pd(svec,svec,0x3);
      t2 = _mm_shuffle_pd(f0v,f0v,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      _mm_store_pd(&norm_array[2],dotprod);
#endif
      
      
#if 0
      // Now inner products with q. 
      // ** pi=(r0,q) 
      norm_array[4] += r0[count]*q[count];
      norm_array[4] += r0[count+1]*q[count+1];
      
      norm_array[5] += r0[count]*q[count+1];
      norm_array[5] -= r0[count+1]*q[count];
#else
      dotprod = _mm_load_pd(&norm_array[4]);
      
      /* dotprod holds the current partial sum */
      
      /*    dotprod[0] += r0[0]*q[0];         */
      /*    dotprod[1] -= r0[1]*q[0];       */
      t1 = _mm_shuffle_pd(qvec,qvec,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(r0v, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += r0[1]*q[1]   */
      /*    dotprod[1] += r0[0]*q[1];    */
      t1 = _mm_shuffle_pd(qvec,qvec,0x3);
      t2 = _mm_shuffle_pd(r0v,r0v,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      _mm_store_pd(&norm_array[4],dotprod);
      
#endif
      
      
#if 0
      // Now inner products with t. 
      // ** eta=(f0,t) 
      norm_array[6] += f0[count]*t[count];
      norm_array[6] += f0[count+1]*t[count+1];
      
      norm_array[7] -= f0[count+1]*t[count];
      norm_array[7] += f0[count]*t[count+1];
#else
      dotprod = _mm_load_pd(&norm_array[6]);
      
      /* dotprod holds the current partial sum */
      
      /*    dotprod[0] += f0[0]*t[0];         */
      /*    dotprod[1] -= f0[1]*t[0];       */
      t1 = _mm_shuffle_pd(tvec,tvec,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(f0v, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += f0[1]*t[1]   */
      /*    dotprod[1] += f0[0]*t[1];    */
      t1 = _mm_shuffle_pd(tvec,tvec,0x3);
      t2 = _mm_shuffle_pd(f0v,f0v,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      _mm_store_pd(&norm_array[6],dotprod);
#endif
      
      
#if 0
      // ** theta=(t,s) 
      norm_array[8] += t[count]*s[count];
      norm_array[8] += t[count+1]*s[count+1];
      
      norm_array[9] += t[count]*s[count+1];
      norm_array[9] -= t[count+1]*s[count];
#else
      dotprod = _mm_load_pd(&norm_array[8]);
      
      /* dotprod holds the current partial sum */
      
      /*    dotprod[0] += t[0]*s[0];         */
      /*    dotprod[1] -= t[1]*s[0];       */
      t1 = _mm_shuffle_pd(svec,svec,0x0);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(tvec, _mm_mul_pd(mask,t1)));
      
      /*    dotprod[0] += t[1]*s[1]   */
      /*    dotprod[1] += t[0]*s[1];    */
      t1 = _mm_shuffle_pd(svec,svec,0x3);
      t2 = _mm_shuffle_pd(tvec,tvec,0x1);
      dotprod = _mm_add_pd(dotprod, _mm_mul_pd(t2,t1));
      
      _mm_store_pd(&norm_array[8],dotprod);
#endif
      
#if 0
      // ** kappa = || t ||^2
      norm_array[10] += t[count]*t[count];
      norm_array[10] += t[count+1]*t[count+1];
      
      // ** rnorm = || r ||^2
      norm_array[11] += r[count]*r[count];
      norm_array[11] += r[count+1]*r[count+1];
#else
      dotprod = _mm_load_pd(&norm_array[10]);
      
      // [ti,tr] & [ri,rr] -> [tr, rr] 
      f0v = _mm_shuffle_pd(tvec,rvec, 0x0);
      
      // [ti,tr] & [ri,rr] -> [ti, ri] 
      r0v = _mm_shuffle_pd(tvec,rvec, 0x3);
      
      // [ tnorm, rnorm ] += [ tr^2, rr^2 ]
      // [ tnorm, rnorm ] += [ ti^2, ri^2 ]
      dotprod = _mm_add_pd(dotprod,_mm_mul_pd(f0v,f0v));
      dotprod = _mm_add_pd(dotprod,_mm_mul_pd(r0v,r0v));
      
      
      _mm_store_pd(&norm_array[10],dotprod);
      
#endif
      
      
      
    }
  }
  else { 
    QDPIO::cout << "ord_ib_stubdates_kernel_sse.h: len not divisible by 2" << endl;
    QDP_abort(1);
  }
  
}
