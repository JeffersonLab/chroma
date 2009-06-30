inline
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
  for(int count = 0; count < len; count+=4) { 
    
    x_ptr[count] = y_ptr[count] - z_ptr[count];
    x_ptr[count+1] = y_ptr[count+1] - z_ptr[count+1];
    x_ptr[count+2] = y_ptr[count+2] - z_ptr[count+2];
    x_ptr[count+3] = y_ptr[count+3] - z_ptr[count+3];
    
    norm += x_ptr[count]*x_ptr[count];
    norm += x_ptr[count+1]*x_ptr[count+1];
    norm += x_ptr[count+2]*x_ptr[count+2];
    norm += x_ptr[count+3]*x_ptr[count+3];
  }
  
  a->norm_ptr[my_id] = norm;
}
