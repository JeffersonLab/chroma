inline
void ord_xpaypbz_kernel(int lo, int hi, int my_id, ord_xpaypbz_arg* a)
{
  int atom = a-> atom;
  int low = atom*lo;
  int len = atom*(hi - lo);

  REAL32* x_ptr = &(a->x_ptr[low]);
  REAL32* y_ptr = &(a->y_ptr[low]);
  REAL32* z_ptr = &(a->z_ptr[low]);
  
  REAL32 a_re = a->a_re;
  REAL32 a_im = a->a_im;
  REAL32 b_re = a->b_re;
  REAL32 b_im = a->b_im;
  

  if( len % 4 == 0 ) { 
    for(int count = 0; count < len; count+=4) { 
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
      
    }
  }
  else { 
    QDPIO::cout << "ord_xpaypbz_kernel_generic.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
}
