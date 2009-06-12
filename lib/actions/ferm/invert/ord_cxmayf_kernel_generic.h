void ord_cxmayf_kernel(int lo, int hi, int my_id, ord_cxmayf_arg* arg)
{
  REAL32* x_ptr = &(arg->x_ptr[lo]);
  REAL32* y_ptr = &(arg->y_ptr[lo]);

  REAL32 a_re = arg->a_re;
  REAL32 a_im = arg->a_im;

  int len = hi - lo;
  
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
