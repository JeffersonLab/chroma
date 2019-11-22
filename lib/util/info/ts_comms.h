#include <vector>
#include <string>

namespace Chroma 
{

  int ts_comms_setup( std::vector<std::string> fifos , size_t ts_size );
  void ts_comms_done();
  
  void ts_comms_send( int harom_instance , int val );
  void ts_comms_send( int harom_instance , const std::string& str );
  
  int ts_comms_recv( int harom_instance );
  std::string ts_comms_recv_str( int ts );

  void* ts_comms_get_shm(int ts);
    
}

