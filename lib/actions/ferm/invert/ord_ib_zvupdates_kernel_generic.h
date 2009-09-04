// 32 BIT Version: Use vector length of 4 for easy vectorization.
// This is guaranteed good for LatticeDiracFermions

inline
void ord_ib_zvupdates_kernel_real32(int lo, int hi, int my_id, ib_zvupdates_arg<REAL32>* a)
{

  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi-lo);

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
  
  REAL32 ztmp[4];
  REAL32 vtmp[4];

  if( len % 4 == 0 ) { 
    for(int count = 0; count < len; count+=4) { 
      
      /* z = (alpha_n/alpha_n-1)*beta z */
      ztmp[0] = z[count];
      ztmp[1] = z[count+1];
      ztmp[2] = z[count+2];
      ztmp[3] = z[count+3];
      
      vtmp[0] = v[count];
      vtmp[1] = v[count+1];
      vtmp[2] = v[count+2];
      vtmp[3] = v[count+3];
      
      z[count]    = arb_re * ztmp[0] - arb_im * ztmp[1];
      z[count+1]  = arb_re * ztmp[1] + arb_im * ztmp[0];  
      z[count+2]  = arb_re * ztmp[2] - arb_im * ztmp[3];
      z[count+3]  = arb_re * ztmp[3] + arb_im * ztmp[2];  
      
      
      /* z += alpha*r */
      z[count  ] += a_re * r[count];
      z[count+1] += a_re * r[count+1];
      z[count+2] += a_re * r[count+2];
      z[count+3] += a_re * r[count+3];
      
      z[count  ] -= a_im * r[count+1];
      z[count+1] += a_im * r[count];
      z[count+2] -= a_im * r[count+3];
      z[count+3] += a_im * r[count+2];
      
      
      /* z -= alpha*delta*v */
      z[count  ] -= ad_re * v[count] ;
      z[count+1] -= ad_re * v[count+1];
      z[count+2] -= ad_re * v[count+2] ;
      z[count+3] -= ad_re * v[count+3];
      
      z[count  ] += ad_im * v[count+1];
      z[count+1] -= ad_im * v[count];
      z[count+2] += ad_im * v[count+3];
      z[count+3] -= ad_im * v[count+2];
      
      
      v[count]   = u[count]   + b_re*vtmp[0] - b_im*vtmp[1];
      v[count+1] = u[count+1] + b_re*vtmp[1] + b_im*vtmp[0];
      v[count+2] = u[count+2] + b_re*vtmp[2] - b_im*vtmp[3];
      v[count+3] = u[count+3] + b_re*vtmp[3] + b_im*vtmp[2];
      
      v[count]   -= d_re*q[count];
      v[count+1] -= d_re*q[count+1];
      v[count+2] -= d_re*q[count+2];
      v[count+3] -= d_re*q[count+3];
      
      v[count]   += d_im*q[count+1];
      v[count+1] -= d_im*q[count];
      v[count+2] += d_im*q[count+3];
      v[count+3] -= d_im*q[count+2];
      
    }  
  }
  else { 
    QDPIO::cout << "ord_ib_zvupdates_kernel_generic.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
}

// 64 BIT Version: Use vector length of 2 for easy vectorization.
// This is guaranteed good for LatticeDiracFermions


inline
void ord_ib_zvupdates_kernel_real64(int lo, int hi, int my_id, ib_zvupdates_arg<REAL64>* a)
{

  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi-lo);

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

  REAL64 ztmp[2];
  REAL64 vtmp[2];

  if( len % 2 == 0) { 
    for(int count = 0; count < len; count+=2) { 
      
      /* z = (alpha_n/alpha_n-1)*beta z */
      ztmp[0] = z[count];
      ztmp[1] = z[count+1];
      vtmp[0] = v[count];
      vtmp[1] = v[count+1];
      
      z[count]    = arb_re * ztmp[0] - arb_im * ztmp[1];
      z[count+1]  = arb_re * ztmp[1] + arb_im * ztmp[0];  
      
      /* z += alpha*r */
      z[count  ] += a_re * r[count];
      z[count+1] += a_re * r[count+1];
      
      z[count  ] -= a_im * r[count+1];
      z[count+1] += a_im * r[count];
      
      /* z -= alpha*delta*v */
      z[count  ] -= ad_re * v[count] ;
      z[count+1] -= ad_re * v[count+1];
      
      z[count  ] += ad_im * v[count+1];
      z[count+1] -= ad_im * v[count];
      
      v[count]   = u[count]   + b_re*vtmp[0] - b_im*vtmp[1];
      v[count+1] = u[count+1] + b_re*vtmp[1] + b_im*vtmp[0];
      
      v[count]   -= d_re*q[count];
      v[count+1] -= d_re*q[count+1];
      
      v[count]   += d_im*q[count+1];
      v[count+1] -= d_im*q[count];
      
      
    }  
  }
  else { 
    QDPIO::cout << "ord_ib_zvupdates_kernel_generic.h: len not divisible by 2"<<endl;
    QDP_abort(1);
  }
}

