inline
void ord_cxmayf_kernel(int lo, int hi, int my_id, ord_cxmayf_arg* arg)
{
  int atom=arg->atom;
  int low = atom*lo;
  int len = atom*(hi - lo);

  REAL32* x_ptr = &(arg->x_ptr[low]);
  REAL32* y_ptr = &(arg->y_ptr[low]);

  REAL32 a_re = arg->a_re;
  REAL32 a_im = arg->a_im;

  if( len%4 == 0) {
    for(int count=0; count < len; count+=4) { 
      x_ptr[count]   = x_ptr[count]     - a_re*y_ptr[count];
      x_ptr[count]   +=                               a_im*y_ptr[count+1];
      
      x_ptr[count+1] = x_ptr[count+1] - a_re*y_ptr[count+1];
      x_ptr[count+1] -=                               a_im*y_ptr[count];
      
      x_ptr[count+2]   = x_ptr[count+2]     - a_re*y_ptr[count+2];
      x_ptr[count+2]   +=                               a_im*y_ptr[count+3];
      x_ptr[count+3] = x_ptr[count+3] - a_re*y_ptr[count+3]; 
      x_ptr[count+3] -=                               a_im*y_ptr[count+2];
      
    }
  }
  else { 
    QDPIO::cout << "ord_cxmayf_kernel: len not divisible by 4" << endl;
    QDP_abort(1);
  }

}
