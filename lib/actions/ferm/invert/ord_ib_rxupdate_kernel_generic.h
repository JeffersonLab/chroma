// 32 BIT Version. Use vector length of 4 (guaranteed OK for LatticeDirac Fermion)
// for easy vectorization (with suitably good compiler).

inline
void ord_ib_rxupdate_kernel_real32(int lo, int hi, int my_id, ib_rxupdate_arg<REAL32>* a)
{

  int atom=a->atom;
  int low=atom*lo;
  int len = atom*(hi - lo);

  REAL32* s = &(a->s_ptr[low]);
  REAL32* t = &(a->t_ptr[low]);
  REAL32* z = &(a->z_ptr[low]);
  REAL32* r = &(a->r_ptr[low]);
  REAL32* x = &(a->x_ptr[low]);

  REAL32 om_re = a->omega_re;
  REAL32 om_im = a->omega_im;

  if( len % 4 == 0 ) {
    for(int count = 0; count < len; count+=4) { 
      
      
      r[count]   = s[count]   - om_re*t[count]   + om_im*t[count+1];
      r[count+1] = s[count+1] - om_re*t[count+1] - om_im*t[count];
      r[count+2] = s[count+2] - om_re*t[count+2] + om_im*t[count+3];
      r[count+3] = s[count+3] - om_re*t[count+3] - om_im*t[count+2];
      
      
      x[count]   += om_re*s[count]   - om_im*s[count+1] + z[count];
      x[count+1] += om_re*s[count+1] + om_im*s[count]   + z[count+1];
      x[count+2] += om_re*s[count+2] - om_im*s[count+3] + z[count+2];
      x[count+3] += om_re*s[count+3] + om_im*s[count+2] + z[count+3];
    }  
  }
  else { 
    QDPIO::cout << "ord_ib_rxupdate_kernel_generic.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
}

// 64 BIT Version. Use vector length of 2 (guaranteed OK for Complex Numbers)
// for easy vectorization (with suitably good compiler). Can do function 
// overloading so no need to change kernel name

inline
void ord_ib_rxupdate_kernel_real64(int lo, int hi, int my_id, ib_rxupdate_arg<REAL64>* a)
{

  int atom=a->atom;
  int low = atom*lo;
  int len =atom*(hi - lo);

  REAL64* s = &(a->s_ptr[low]);
  REAL64* t = &(a->t_ptr[low]);
  REAL64* z = &(a->z_ptr[low]);
  REAL64* r = &(a->r_ptr[low]);
  REAL64* x = &(a->x_ptr[low]);

  REAL64 om_re = a->omega_re;
  REAL64 om_im = a->omega_im;

  if( len % 2 == 0) { 
    for(int count = 0; count < len; count+=2) { 
      r[count]   = s[count]   - om_re*t[count]   + om_im*t[count+1];
      r[count+1] = s[count+1] - om_re*t[count+1] - om_im*t[count];
      
      x[count]   += om_re*s[count]   - om_im*s[count+1] + z[count];
      x[count+1] += om_re*s[count+1] + om_im*s[count]   + z[count+1];
      
    }  
  }
  else { 
    QDPIO::cout << "ord_ib_rxupdate_kernel_generic.h: len not divisible by 2"<<endl;
    QDP_abort(1);
  }
}

