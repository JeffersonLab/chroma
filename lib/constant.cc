#include "chroma_config.h"
#include "chromabase.h"

namespace Chroma  
{
  namespace
  {
    __chroma_constant* constant_data;
    bool constant_data_allocated = false;
  }
  
  const __chroma_constant& constant()
  {
    if (!constant_data)
      {
	constant_data = new __chroma_constant;
	constant_data_allocated = true;
	
	constant_data->twopi = 6.283185307179586476925286;

#if BASE_PRECISION == 32
	constant_data->fuzz = 1.0e-5;
#elif BASE_PRECISION == 64
	constant_data->fuzz = 1.0e-10;
#endif
      }

    return *constant_data;
  }


  void constant_destroy()
  {
    if (constant_data_allocated)
      delete constant_data;
    constant_data_allocated = false;
  }

}
