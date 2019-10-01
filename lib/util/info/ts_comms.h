#include <vector>
#include <string>

namespace Chroma 
{

  int ts_comms_setup( std::vector<std::string> fifos , size_t ts_size );
  void ts_comms_done();
  
  void ts_comms_send( int harom_instance , int val );
  int ts_comms_recv( int harom_instance );
  void* ts_comms_get_shm(int ts);
    
}

