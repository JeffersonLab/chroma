inline
void ord_norm2x_cdotxy_kernel(int lo, int hi, int my_id, ord_norm2x_cdotxy_arg* a)
{
  int atom = a->atom;
  int low = atom*lo;
  int len = atom*(hi-lo);

  REAL32* x_ptr = &(a->x_ptr[low]);
  REAL32* y_ptr = &(a->y_ptr[low]);
  REAL64 norm_array[3] = {0,0,0};

  if( atom % 4 == 0) { 
    for(int count = 0; count < len; count+=4) { 
      // norm
      norm_array[0] += x_ptr[count]*x_ptr[count];
      norm_array[0] += x_ptr[count+1]*x_ptr[count+1];
      norm_array[0] += x_ptr[count+2]*x_ptr[count+2];
      norm_array[0] += x_ptr[count+3]*x_ptr[count+3];
      
      // cdot re
      norm_array[1] += x_ptr[count]*y_ptr[count];
      norm_array[1] += x_ptr[count+1]*y_ptr[count+1];
      norm_array[1] += x_ptr[count+2]*y_ptr[count+2];
      norm_array[1] += x_ptr[count+3]*y_ptr[count+3];
      
      // cdot im
      norm_array[2] += x_ptr[count]*y_ptr[count+1];
      norm_array[2] -= x_ptr[count+1]*y_ptr[count];
      norm_array[2] += x_ptr[count+2]*y_ptr[count+3];
      norm_array[2] -= x_ptr[count+3]*y_ptr[count+2];
      
    }
    a->norm_space[3*my_id]=norm_array[0];
    a->norm_space[3*my_id+1]=norm_array[1];
    a->norm_space[3*my_id+2]=norm_array[2];
  }
  else { 
    QDPIO::cout << "ord_norm2x_cdotxy_kernel_generic.h: len not divisible by 4" << endl;
    QDP_abort(1);
  }
}
