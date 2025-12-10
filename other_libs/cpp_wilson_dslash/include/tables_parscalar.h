#ifndef TABLES_PARSCALAR_H
#define TABLES_PARSCALAR_H

#include <qmp.h>
#include "allocate.h"

namespace CPlusPlusWilsonDslash {

  template<typename HalfSpinor, int Nd>
  class DslashTables {
  public:
    DslashTables(int subgrid[]);
    ~DslashTables();

    // Accessors
    HalfSpinor* getChi1() {
      return chi1;
    }

    HalfSpinor* getChi2() {
      return chi2;
    }

    HalfSpinor*** getRecvBufptr() {
      return (HalfSpinor***)recv_bufptr;
    }

    HalfSpinor*** getSendBufptr() {
      return (HalfSpinor***)send_bufptr;
    }
    
    // Communications
    //
    inline
    void startReceives() {
      /* Prepost all receives */
      if (total_comm > 0) {

	// Use QMP Harness 
	if (QMP_start(recv_all_mh[0]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	  QMP_abort(1);
	}
	
	if (QMP_start(recv_all_mh[1]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	  QMP_abort(1);
	}
      }
    }

    inline void finishReceiveFromForward() 
    {  
      /* Finish all forward receives */
      if (total_comm > 0 ) { 
	if (QMP_wait(recv_all_mh[1]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	  QMP_abort(1);
	}
      }
    }

    inline void finishReceiveFromBack() 
    { 
      if( total_comm > 0 ) { 
	/* Finish all forward receives */
	if (QMP_wait(recv_all_mh[0]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	  QMP_abort(1);
	}
      }
    }
	

    inline void startSendBack() 
    { 
      if(total_comm > 0) {
	if (QMP_start(send_all_mh[1]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	  QMP_abort(1);
	}
      }
    }
    
    inline void startSendForward() 
    {
      if(total_comm > 0) {
	if (QMP_start(send_all_mh[0]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	  QMP_abort(1);
	}
      }
    }

    inline void finishSendBack() {
      if( total_comm > 0 ) {
	/* Finish all sends */
	if (QMP_wait(send_all_mh[1]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	  QMP_abort(1);
	}
      }
    }

    inline void finishSendForward() {
      if( total_comm > 0 ) {
	/* Finish all sends */
	if (QMP_wait(send_all_mh[0]) != QMP_SUCCESS) {
	  QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	  QMP_abort(1);
	}
      }
    }
  private:

    static void* xchi;
    
    HalfSpinor *chi1;
    HalfSpinor *chi2;

    HalfSpinor* recv_bufptr[2][Nd];
    HalfSpinor* send_bufptr[2][Nd];

    HalfSpinor* send_bufs;
    HalfSpinor* recv_bufs;

    QMP_msgmem_t send_msg[2][Nd];
    QMP_msgmem_t recv_msg[2][Nd];

    QMP_msghandle_t send_mh[2][Nd];
    QMP_msghandle_t recv_mh[2][Nd];

    QMP_msghandle_t send_all_mh[Nd];
    QMP_msghandle_t recv_all_mh[Nd];

    int total_comm;

   
  };

  template<typename HalfSpinor, int Nd>
    void* DslashTables<HalfSpinor, Nd>::xchi = nullptr;

  template<typename HalfSpinor, int Nd>
    DslashTables<HalfSpinor,Nd>::~DslashTables() 
    {
      /* Memory/comms handles */
      if (total_comm > 0) {
	
	
	for(int i=0; i < 2; i++) { 
	  /* If we collapsed the handles -- free the collapsed handles */
	  QMP_free_msghandle(send_all_mh[i]);
	  QMP_free_msghandle(recv_all_mh[i]);
	  
	  /* Free the msgmem structures */
	  for(int mu=0; mu < total_comm; mu++) { 
	    QMP_free_msgmem(send_msg[i][mu]);
	    QMP_free_msgmem(recv_msg[i][mu]);
	  }
	}
	
      }
      
      /* Free all space - 4 spinors and actual comms buffers  */
      /* Simon: we have no mechanism to check if there are still other
       * instances using our static xchi, so we do not free it. Typically
       * it is used till the end of the program anyways */
      
      /* Free the shift table itself */
    }

  template<typename HalfSpinor, int Nd>
    DslashTables<HalfSpinor, Nd>::DslashTables(int subgrid[]) 
    {
      struct BufTable { 
	unsigned int dir;
	unsigned int offset;
	unsigned int size;
	unsigned int pad;
      };

      /* Get the dimensions of the machine */
      const int *machine_size = QMP_get_logical_dimensions();
    
      /* Check we are in 4D */
      if (QMP_get_logical_number_of_dimensions() != 4) {
	QMP_error("init_sse_su3dslash: number of logical dimensions does not match problem");
	QMP_abort(1);
      }
      int sx = subgrid[0];
      int sy = subgrid[1];
      int sz = subgrid[2];
      int st = subgrid[3];

      /* Compute the number of bounds */
      int nbound[4];
      nbound[0]=(sy*sz*st)/2;
      nbound[1]=(sx*sz*st)/2;
      nbound[2]=(sx*sy*st)/2;
      nbound[3]=(sx*sy*sz)/2;
      int subgrid_vol_cb = sx*sy*sz*st/2;
      int subgrid_vol = sx*sy*sz*st;
      BufTable recv[2][4];

      int offset = 0;
      int num=0;
      int pad;
      
      for(int i=0; i < 2; i++) { 

	num = 0;
	/* i=0 => recv from forward/send backward */
	/* i=1 => recv from backward/send forward */
	/* Fill out the buffer descriptor structure */
	for(int mu=0; mu < Nd; ++mu) {
	  if( machine_size[mu] > 1 ) { 
	    
	    recv[i][num].dir = mu;
	    recv[i][num].offset = offset;
	    recv[i][num].size = nbound[mu]*sizeof(HalfSpinor);
	    
	    
	    /* Cache line align the next buffer */
	    if ( (offset % Cache::CacheLineSize) != 0 ) { 
	      pad = Cache::CacheLineSize - (offset % Cache::CacheLineSize);
	    }
	    else { 
	      pad = 0;
	    }
	    
	    
	    recv[i][num].pad = pad;
	    offset += recv[i][num].size + pad;
	    num++;
	  }
	}
      }
      /*** ABOVE: by now offset should 
	   i) Be big enough to cover the receive buffers
	   ii) Be cache line padded
	   iii) Be set padded assuming that the comms start on a set boundary 
	   
	   You could do two of these to cover send and recv buffers */
      
      
      /* Now for the chi spinors 
	 This is the size of one of the Chi-s either forward or backward.
	 The factor of 4 is the 4 Mu directions
         The second factor of 4 is for the 4 types.*/
      int chisize = sizeof(HalfSpinor)*subgrid_vol_cb*4*Nd;
      
      /* Total amount: 2 x offset -- for the comms.
	 2 x chisize -- for the half spinor temps (2 cb's)
	 10*CacheCacheLine - 5 lines of padding between
	 comms bufs and chi1, and chi1 and chi2 */
      int total_allocate = 2*chisize+2*offset;
      
      if(xchi == nullptr) {
        if ((xchi = CPlusPlusWilsonDslash::alloc(total_allocate)) == nullptr) {
          QMP_error("init_wnxtsu3dslash: could not initialize xchi1");
          QMP_abort(1);
        }
      }
      
      /* Get the aligned pointer out. This is the start of our memory */
      unsigned char* chi = (unsigned char *)xchi;
      unsigned char* send_bufs = chi + offset;
      
      /* Put pointers to the send and receive buffers */
      for(int i=0; i < 2; i++) { 
	for(int mu=0; mu < num; mu++) { 
	  recv_bufptr[i][mu] = (HalfSpinor *)(chi + recv[i][mu].offset + recv[i][mu].pad);
	  send_bufptr[i][mu] = (HalfSpinor *)(send_bufs + recv[i][mu].offset + recv[i][mu].pad);
	}
      }
      
      /* Chi 1 should be after the send bufs */
      /* Should be padded. and aligned */
      chi1 = (HalfSpinor *)((unsigned char *)send_bufs + offset);
      
      /* Strictly speaking I shouldn't be writing into chi2 
	 while working on chi1 and vice versa. So I don't want
	 to worry about false sharing here. For now just Pad to a
	 line
      */
      if( (chisize % Cache::CacheLineSize) != 0 ) { 
	pad = Cache::CacheLineSize - (chisize%Cache::CacheLineSize);
      }
      else { 
	pad = 0;
      }
      chi2 = (HalfSpinor *)((unsigned char *)chi1 + chisize+pad);
      
      /* Now we can set up the QMP isms... */
      for(int i=0; i < 2; i++) { 
	for(int mu=0; mu < num; mu++) { 
	  recv_msg[i][mu] = QMP_declare_msgmem(recv_bufptr[i][mu], recv[i][mu].size);
	  send_msg[i][mu] = QMP_declare_msgmem(send_bufptr[i][mu], recv[i][mu].size);
	  if( i == 0 ) { 
	    /* Recv from forward, send backward pair */
	    recv_mh[i][mu]= QMP_declare_receive_relative(recv_msg[i][mu], recv[i][mu].dir, +1, 0);
	    send_mh[i][mu]= QMP_declare_send_relative(send_msg[i][mu], recv[i][mu].dir, -1, 0);
	  }
	  else { 
	    /* Recv from backwards, send forward pair */
	    recv_mh[i][mu]= QMP_declare_receive_relative(recv_msg[i][mu], recv[i][mu].dir, -1, 0);
	    send_mh[i][mu]= QMP_declare_send_relative(send_msg[i][mu], recv[i][mu].dir, +1, 0);
	  }
	}
      }
      
      /* Combine the messages */
      if (num > 0) {
	for(int i=0; i<2; i++) { 
	  send_all_mh[i] = QMP_declare_multiple(send_mh[i], num);
	  recv_all_mh[i] = QMP_declare_multiple(recv_mh[i], num);
	}
      }
      
      total_comm = num;
    }



}

#endif
