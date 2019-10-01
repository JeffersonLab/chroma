#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <stdio.h>
#include <string.h>

#include <string>
#include <iostream>
#include <vector>
#include <cassert>

#include "ts_comms.h"

namespace Chroma 
{

  namespace
  {
    struct ts_comms_t
    {
      std::string fifo_send_name;
      std::string fifo_recv_name;
    
      int         fifo_send_fd;
      int         fifo_recv_fd;
    
      std::string shm_name;
      int         shm_fd;
    
      void*       ts_buf;
    };

    std::vector< ts_comms_t > ts_comms;

    size_t ts_size_;
  }


  void* ts_comms_get_shm(int ts)
  {
    assert( ts >= 0 && ts < ts_comms.size() );
    return ts_comms[ts].ts_buf;
  }
  

  void ts_comms_send( int ts , int val )
  {
    write( ts_comms[ts].fifo_send_fd , &val , sizeof(int) );
  }

  int ts_comms_recv( int ts )
  {
    int val;
    int r = read( ts_comms[ts].fifo_recv_fd , &val, sizeof(int) );
    if ( r != sizeof(int) ) {
      std::cout << "fifo read less than the size of an integer\n";
      exit(1);
    }
    return val;
  }


    void ts_comms_done()
    {
      for ( int i = 0 ; i < ts_comms.size() ; ++i )
	{
	  std::cout << "closing comms for TS " << i << "\n";

	  close(  ts_comms[i].fifo_send_fd );
	  close(  ts_comms[i].fifo_recv_fd );
	  unlink( ts_comms[i].fifo_send_name.c_str() );
	  unlink( ts_comms[i].fifo_recv_name.c_str() );

	  munmap( ts_comms[i].ts_buf , ts_size_ );
	  shm_unlink( ts_comms[i].shm_name.c_str() );
	}
    }


  int ts_comms_setup( std::vector<std::string> fifos , size_t ts_size )
  {
    ts_size_ = ts_size;
    
    ts_comms.resize( fifos.size() );
    
    for ( int i = 0 ; i < fifos.size() ; ++i )
      {
	std::cout << "creating comms for TS " << i << "\n";
	    
	ts_comms[i].fifo_send_name = fifos[i] + "_th";
	ts_comms[i].fifo_recv_name = fifos[i] + "_tc";

	if( access( ts_comms[i].fifo_send_name.c_str() , F_OK ) != -1 ) {
	  std::cout << "File exists: " << ts_comms[i].fifo_send_name.c_str() << "\n";
	  if ( remove(ts_comms[i].fifo_send_name.c_str() ) ) {
	      std::cout << "..and could not remove the file. giving up\n";
	      return -1;
	    }
	}
	
	if( access( ts_comms[i].fifo_recv_name.c_str() , F_OK ) != -1 ) {
	  std::cout << "File exists: " << ts_comms[i].fifo_recv_name.c_str() << "\n";
	  if ( remove(ts_comms[i].fifo_recv_name.c_str() ) ) {
	      std::cout << "..and could not remove the file. giving up\n";
	      return -1;
	    }
	}

	
	if (mkfifo( ts_comms[i].fifo_send_name.c_str() , 0666 )) {
	  std::cout << "error creating send fifo\n";
	  return -1;
	}
	std::cout << "send fifo created\n";


	if (mkfifo( ts_comms[i].fifo_recv_name.c_str() , 0666 )) {
	  std::cout << "error creating recv fifo\n";
	  return -1;
	}
	std::cout << "recv fifo created\n";
	

	std::cout << "Waiting for harom to connect on " << ts_comms[i].fifo_send_name.c_str() << "\n";
	if ((ts_comms[i].fifo_send_fd = open( ts_comms[i].fifo_send_name.c_str() , O_WRONLY)) < 0) {
	  std::cout << "error opening send fifo\n";
	  return -1;
	}
	std::cout << "ts_comms[i].fifo_send_fd = " << ts_comms[i].fifo_send_fd << "\n";
	    
	std::cout << "Waiting for harom to connect on " << ts_comms[i].fifo_recv_name.c_str() << "\n";
	while ((ts_comms[i].fifo_recv_fd = open( ts_comms[i].fifo_recv_name.c_str() , O_RDONLY)) < 0);
	std::cout << "fifo successfully setup for TS " << i << "\n";

	//
	// Sending each harom its assigned number used later to open the correct shared memory file
	//
	ts_comms_send( i , i );

	ts_comms[i].shm_name = "shm_ts_" + std::to_string(i);

	ts_comms[i].shm_fd = shm_open( ts_comms[i].shm_name.c_str() , O_CREAT | O_RDWR , 0666 );
	if (ts_comms[i].shm_fd == -1) {
	  std::cout << "Failed to open shared memory file " << ts_comms[i].shm_name << "\n";
	  return -1;
	}

	ftruncate( ts_comms[i].shm_fd , ts_size );

	ts_comms_send( i , -1 ); // used for sync

	ts_comms[i].ts_buf = mmap(NULL, ts_size, PROT_READ | PROT_WRITE, MAP_SHARED, ts_comms[i].shm_fd, 0);
	if ( ts_comms[i].ts_buf == MAP_FAILED ) {
	  std::cout << "Map failed for TS = " << i << "\n";
	  return -1;
	}
	
	std::cout << "shared memory sucessfully mapped.\n";


      }


  }

}
