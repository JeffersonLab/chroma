// -*- C++ -*-
// $Id: lwldslash_w.h,v 3.3 2009/11/14 20:01:46 eneil Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_llvm_h__
#define __lwldslash_llvm_h__

#include "state.h"
#include "io/aniso_io.h"
#include "actions/ferm/linop/lwldslash_base_w.h"

extern "C" {
#include "dslash_sig_0.h"
#include "dslash_sig_1.h"
}

namespace Chroma 
{ 
  //! General Wilson-Dirac dslash
  /*!
   * \ingroup linop
   *
   * DSLASH
   *
   * This routine is specific to Wilson fermions!
   *
   * Description:
   *
   * This routine applies the operator D' to Psi, putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \
   *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
   *	       /    mu			  mu
   *	       ---
   *	       mu=0
   *
   *	             Nd-1
   *	             ---
   *	             \    +
   *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *	             /    mu			   mu
   *	             ---
   *	             mu=0
   *
   */

  template<typename T, typename P, typename Q> 
  class LLVMWilsonDslashT : public WilsonDslashBase<T, P, Q>
  {
  public:

    //! Empty constructor. Must use create later
    LLVMWilsonDslashT();

    //! Full constructor
    LLVMWilsonDslashT(Handle< FermState<T,P,Q> > state);

    //! Full constructor with anisotropy
    LLVMWilsonDslashT(Handle< FermState<T,P,Q> > state,
		    const AnisoParam_t& aniso_);

    //! Full constructor with general coefficients
    LLVMWilsonDslashT(Handle< FermState<T,P,Q> > state,
		    const multi1d<Real>& coeffs_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > state);

    //! Creation routine with anisotropy
    void create(Handle< FermState<T,P,Q> > state,
		const AnisoParam_t& aniso_);

    //! Full constructor with general coefficients
    void create(Handle< FermState<T,P,Q> > state, 
		const multi1d<Real>& coeffs_);

    //! No real need for cleanup here
    ~LLVMWilsonDslashT() { comms_free(); }

    /**
     * Apply a dslash
     *
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT vector               (Read) 
     *
     * \return The output of applying dslash on psi
     */
    void apply (T& chi, const T& psi, enum PlusMinus isign, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

  protected:
    //! Get the anisotropy parameters
    const multi1d<Real>& getCoeffs() const {return coeffs;}

  private:
    void comms_setup();
    void setup();
    void comms_free();
    void comms_wait() const;
    void comms_send_receive(int i) const;

    multi1d<Real> coeffs;  /*!< Nd array of coefficients of terms in the action */
    Handle< FermBC<T,P,Q> >  fbc;
    Q  u;

    int innerCount[2];
    int faceCount[2];
    const int *innerSites[2];
    const int *faceSites[2];
    
    struct comms_t {
      bool do_comms;
      int snd_nd;
      int snd_sz;
      int rcv_nd;
      int rcv_sz;
    };
    struct comms_t comms[8];

    double * send_buf[8];
    double * recv_buf[8];
    QMP_msgmem_t msg[8][2];
    QMP_msghandle_t mh_a[8][2], mh[8];
    QMP_mem_t *send_buf_mem[8];
    QMP_mem_t *recv_buf_mem[8];    

    int comm_thread;
  };

  //! General Wilson-Dirac dslash
  /*! \ingroup linop
   * DSLASH
   *
   * This routine is specific to Wilson fermions!
   *
   * Description:
   *
   * This routine applies the operator D' to Psi, putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \
   *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
   *	       /    mu			  mu
   *	       ---
   *	       mu=0
   *
   *	             Nd-1
   *	             ---
   *	             \    +
   *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *	             /    mu			   mu
   *	             ---
   *	             mu=0
   *
   */

  //! Empty constructor
  template<typename T, typename P, typename Q>
  LLVMWilsonDslashT<T,P,Q>::LLVMWilsonDslashT() {
    for (int i=0;i<8;i++)
      comms[i].do_comms = false;
  }
  
  //! Full constructor
  template<typename T, typename P, typename Q>
  LLVMWilsonDslashT<T,P,Q>::LLVMWilsonDslashT(Handle< FermState<T,P,Q> > state)
  {
    create(state);
  }
  
  //! Full constructor with anisotropy
  template<typename T, typename P, typename Q>
  LLVMWilsonDslashT<T,P,Q>::LLVMWilsonDslashT(Handle< FermState<T,P,Q> > state,
				   const AnisoParam_t& aniso_) 
  {
    create(state, aniso_);
  }

  //! Full constructor with general coefficients
  template<typename T, typename P, typename Q>
  LLVMWilsonDslashT<T,P,Q>::LLVMWilsonDslashT(Handle< FermState<T,P,Q> > state,
				   const multi1d<Real>& coeffs_)
  {
    create(state, coeffs_);
  }

  //! Creation routine
  template<typename T, typename P, typename Q>
  void LLVMWilsonDslashT<T,P,Q>::create(Handle< FermState<T,P,Q> > state)
  {
    multi1d<Real> cf(Nd);
    cf = 1.0;
    create(state, cf);
  }

  //! Creation routine with anisotropy
  template<typename T, typename P, typename Q>
  void LLVMWilsonDslashT<T,P,Q>::create(Handle< FermState<T,P,Q> > state,
			       const AnisoParam_t& anisoParam) 
  {
    START_CODE();

    create(state, makeFermCoeffs(anisoParam));

    END_CODE();
  }



  template<typename T, typename P, typename Q>
  void LLVMWilsonDslashT<T,P,Q>::comms_setup() 
  {
    for (int i = 0 ; i < 8 ; i++ )
      if (comms[i].do_comms)
	{
	  int dstnum = comms[i].snd_sz;
	  int srcnum = comms[i].rcv_sz;
	  int dstnode = comms[i].snd_nd;
	  int srcnode = comms[i].rcv_nd;
    
	  send_buf_mem[i] = QMP_allocate_aligned_memory(dstnum,QDP_ALIGNMENT_SIZE, (QMP_MEM_COMMS|QMP_MEM_FAST) ); // packed data to send
	  if( send_buf_mem[i] == 0x0 ) { 
	    send_buf_mem[i] = QMP_allocate_aligned_memory(dstnum, QDP_ALIGNMENT_SIZE, QMP_MEM_COMMS);
	    if( send_buf_mem[i] == 0x0 ) { 
	      QDP_error_exit("Unable to allocate send_buf_mem\n");
	    }
	  }
	  recv_buf_mem[i] = QMP_allocate_aligned_memory(srcnum,QDP_ALIGNMENT_SIZE, (QMP_MEM_COMMS|QMP_MEM_FAST)); // packed receive data
	  if( recv_buf_mem[i] == 0x0 ) { 
	    recv_buf_mem[i] = QMP_allocate_aligned_memory(srcnum, QDP_ALIGNMENT_SIZE, QMP_MEM_COMMS); 
	    if( recv_buf_mem[i] == 0x0 ) { 
	      QDP_error_exit("Unable to allocate recv_buf_mem\n");
	    }
	  }
	  send_buf[i]=(double*)QMP_get_memory_pointer(send_buf_mem[i]);
	  recv_buf[i]=(double*)QMP_get_memory_pointer(recv_buf_mem[i]);


	  msg[i][0] = QMP_declare_msgmem( recv_buf[i] , srcnum );

	  if( msg[i][0] == (QMP_msgmem_t)NULL ) { 
	    QDP_error_exit("QMP_declare_msgmem for msg[0] failed in Map::operator()\n");
	  }

	  msg[i][1] = QMP_declare_msgmem( send_buf[i] , dstnum );

	  if( msg[i][1] == (QMP_msgmem_t)NULL ) {
	    QDP_error_exit("QMP_declare_msgmem for msg[1] failed in Map::operator()\n");
	  }

	  mh_a[i][0] = QMP_declare_receive_from(msg[i][0], srcnode, 0);
	  if( mh_a[i][0] == (QMP_msghandle_t)NULL ) { 
	    QDP_error_exit("QMP_declare_receive_from for mh_a[0] failed in Map::operator()\n");
	  }

	  mh_a[i][1] = QMP_declare_send_to(msg[i][1], dstnode , 0);
	  if( mh_a[i][1] == (QMP_msghandle_t)NULL ) {
	    QDP_error_exit("QMP_declare_send_to for mh_a[1] failed in Map::operator()\n");
	  }

	  mh[i] = QMP_declare_multiple(mh_a[i], 2);
	  if( mh[i] == (QMP_msghandle_t)NULL ) { 
	    QDP_error_exit("QMP_declare_multiple for mh failed in Map::operator()\n");
	  }
	}
  }


  template<typename T, typename P, typename Q>
  void LLVMWilsonDslashT<T,P,Q>::comms_send_receive(int i) const {
    QMP_status_t err;
    if ((err = QMP_start(mh[i])) != QMP_SUCCESS)
      QDP_error_exit(QMP_error_string(err));
  }


  template<typename T, typename P, typename Q>
  void LLVMWilsonDslashT<T,P,Q>::comms_wait() const {
    QMP_status_t err;
    for (int i = 0 ; i < 8 ; i++ )
      if (comms[i].do_comms)
	if ((err = QMP_wait(mh[i])) != QMP_SUCCESS)
	  QDP_error_exit(QMP_error_string(err));
    
#if QDP_DEBUG >= 3
    QDP_info("Map: calling free msgs");
#endif
  }


  template<typename T, typename P, typename Q>
  void LLVMWilsonDslashT<T,P,Q>::comms_free() {
    for (int i = 0 ; i < 8 ; i++ )
      if (comms[i].do_comms)
	{
	  QMP_free_msghandle(mh[i]);
	  // QMP_free_msghandle(mh_a[1]);
	  // QMP_free_msghandle(mh_a[0]);
	  QMP_free_msgmem(msg[i][1]);
	  QMP_free_msgmem(msg[i][0]);
	  
	  QMP_free_memory(recv_buf_mem[i]);
	  QMP_free_memory(send_buf_mem[i]);
	}
  }


  template<typename T, typename P, typename Q>
  void LLVMWilsonDslashT<T,P,Q>::setup()
  {
    //QDPIO::cout << "Setting up LLVM Wilson Dslash\n";

    // Get communication footprint
    int offnode_maps=0;
    int comm_no=-1;
    for( int dir = 3 ; dir >= 0 ; --dir ) {
      for( int isign = -1 ; isign <= +1 ; isign += 2 ) {
	comm_no++;

	const Map& map = shift.getMap(isign,dir);

	//QDPIO::cout << "dir = " << dir << "     isign = " << isign << "   id = " << id <<  "\n";

	if (map.hasOffnode()) 
	  {
	    offnode_maps |= map.getId();

	    int dstnum = shift.getMap(isign,dir).get_destnodes_num()[rb[0].getId()][0]*sizeof(double)*12; // we are sending half-spinors
	    int srcnum = shift.getMap(isign,dir).get_srcenodes_num()[rb[0].getId()][0]*sizeof(double)*12;

	    // QDPIO::cout << "receive buffer size = " << dstnum << "\n";
	    // QDPIO::cout << "send buffer size    = " << srcnum << "\n";

	    int dstnode = shift.getMap(isign,dir).get_destnodes()[0];
	    int srcnode = shift.getMap(isign,dir).get_srcenodes()[0];
	  
	    // QDPIO::cout << "send to      = " << dstnode << "\n";
	    // QDPIO::cout << "receive from = " << srcnode << "\n";

	    comms[comm_no].do_comms=true;
	    comms[comm_no].snd_nd=dstnode;
	    comms[comm_no].snd_sz=dstnum;
	    comms[comm_no].rcv_nd=srcnode;
	    comms[comm_no].rcv_sz=srcnum;
	  }
	else
	  {
	    comms[comm_no].do_comms=false;
	  }
      }
    }

    comms_setup();

    // QDPIO::cout << "comms footprint = " << offnode_maps << "\n";

    for(int cb=0 ; cb<2 ; ++cb) {
      if (offnode_maps > 0) 
	{      
	  innerCount[cb] = MasterMap::Instance().getCountInner(rb[cb],offnode_maps);
	  faceCount[cb]  = MasterMap::Instance().getCountFace(rb[cb],offnode_maps);
	  innerSites[cb] = MasterMap::Instance().getInnerSites(rb[cb],offnode_maps).slice();
	  faceSites[cb]  = MasterMap::Instance().getFaceSites(rb[cb],offnode_maps).slice();
	} 
      else
	{
	  innerCount[cb] = rb[cb].numSiteTable();
	  faceCount[cb]  = 0;
	  innerSites[cb] = rb[cb].siteTable().slice();
	  faceSites[cb]  = NULL;
	}
      // QDPIO::cout << "inner rb[" << cb << "] = " << innerCount[cb] << "   "
      // 		  << "face  rb[" << cb << "] = " << faceCount[cb] << "\n";
    }

    comm_thread = omp_get_max_threads() >=4 ? 3 : 0;
    //QDPIO::cout << "using comm thread " << comm_thread << "\n";
  } // setup




  //! Full constructor with general coefficients
  template<typename T, typename P, typename Q>
  void LLVMWilsonDslashT<T,P,Q>::create(Handle< FermState<T,P,Q> > state,
			       const multi1d<Real>& coeffs_)
  {
    // Save a copy of the aniso params original fields and with aniso folded in
    coeffs = coeffs_;

    // Save a copy of the fermbc
    fbc = state->getFermBC();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "LLVMWilsonDslash: error: fbc is null" << endl;
      QDP_abort(1);
    }

    u.resize(Nd);

    // Fold in anisotropy
    for(int mu=0; mu < u.size(); ++mu) {
      u[mu] = (state->getLinks())[mu];
    }
  
    // Rescale the u fields by the anisotropy
    for(int mu=0; mu < u.size(); ++mu)
    {
      u[mu] *= coeffs[mu];
    }

    setup();
  }




  //! General Wilson-Dirac dslash
  /*! \ingroup linop
   * Wilson dslash
   *
   * Arguments:
   *
   *  \param chi	      Result				                (Write)
   *  \param psi	      Pseudofermion field				(Read)
   *  \param isign      D'^dag or D' ( MINUS | PLUS ) resp.		(Read)
   *  \param cb	      Checkerboard of OUTPUT vector			(Read) 
   */
  template<typename T, typename P, typename Q>
  void 
  LLVMWilsonDslashT<T,P,Q>::apply (T& chi, const T& psi, 
			  enum PlusMinus isign, int cb) const
  {
    double* psi_ptr = (double*)psi.getFjit();
    double* chi_ptr = (double*)chi.getFjit();
    double* u0_ptr = (double*)u[0].getFjit();
    double* u1_ptr = (double*)u[1].getFjit();
    double* u2_ptr = (double*)u[2].getFjit();
    double* u3_ptr = (double*)u[3].getFjit();

    switch (isign)
      {
      case PLUS:
#pragma omp parallel default(shared) 
	{
	  int threads_num = omp_get_num_threads();
	  int myId = omp_get_thread_num();

	  if (comms[0].do_comms)
	    {
	      {
		const Map& map = shift.getMap(-1,3);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_M_3_0(low,high,0,true,0,soffset_slice,send_buf[0],psi_ptr,u3_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(0);

	      {
		const Map& map = shift.getMap(+1,3);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_P_3_0(low,high,0,true,0,soffset_slice,send_buf[1],psi_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(1);
	    }

	  if (comms[2].do_comms) 
	    {
	      {
		const Map& map = shift.getMap(-1,2);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_M_2_0(low,high,0,true,0,soffset_slice,send_buf[2],psi_ptr,u2_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(2);

	      {
		const Map& map = shift.getMap(+1,2);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_P_2_0(low,high,0,true,0,soffset_slice,send_buf[3],psi_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(3);
	    }

	  if (comms[4].do_comms) 
	    {
	      {
		const Map& map = shift.getMap(-1,1);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_M_1_0(low,high,0,true,0,soffset_slice,send_buf[4],psi_ptr,u1_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(4);

	      {
		const Map& map = shift.getMap(+1,1);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_P_1_0(low,high,0,true,0,soffset_slice,send_buf[5],psi_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(5);
	    }

	  if (comms[6].do_comms) 
	    {
	      {
		const Map& map = shift.getMap(-1,0);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_M_0_0(low,high,0,true,0,soffset_slice,send_buf[6],psi_ptr,u0_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(6);

	      {
		const Map& map = shift.getMap(+1,0);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_P_0_0(low,high,0,true,0,soffset_slice,send_buf[7],psi_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(7);
	    }

	  {
	    int low = innerCount[cb]*myId/threads_num;
	    int high = innerCount[cb]*(myId+1)/threads_num;

	    func_dslash_____0( low , high , 0 , false , 0 , innerSites[cb] , chi_ptr ,
			       shift.getMap(-1,3).goffset(rb[cb]).slice(), NULL, psi_ptr, u3_ptr ,
			       shift.getMap(+1,3).goffset(rb[cb]).slice(), NULL, psi_ptr, u3_ptr ,
			       shift.getMap(-1,2).goffset(rb[cb]).slice(), NULL, psi_ptr, u2_ptr ,
			       shift.getMap(+1,2).goffset(rb[cb]).slice(), NULL, psi_ptr, u2_ptr ,
			       shift.getMap(-1,1).goffset(rb[cb]).slice(), NULL, psi_ptr, u1_ptr ,
			       shift.getMap(+1,1).goffset(rb[cb]).slice(), NULL, psi_ptr, u1_ptr ,
			       shift.getMap(-1,0).goffset(rb[cb]).slice(), NULL, psi_ptr, u0_ptr ,
			       shift.getMap(+1,0).goffset(rb[cb]).slice(), NULL, psi_ptr, u0_ptr );
	  }

	  if (myId == comm_thread)
	    comms_wait();
#pragma omp barrier

	  {	
	    int low = faceCount[cb]*myId/threads_num;
	    int high = faceCount[cb]*(myId+1)/threads_num;

	    func_dslash_____0( low , high , 0 , false , 0 , faceSites[cb] , chi_ptr ,
			       shift.getMap(-1,3).goffset(rb[cb]).slice(),recv_buf[0],psi_ptr,u3_ptr,
			       shift.getMap(+1,3).goffset(rb[cb]).slice(),recv_buf[1],psi_ptr,u3_ptr,
			       shift.getMap(-1,2).goffset(rb[cb]).slice(),recv_buf[2],psi_ptr,u2_ptr,
			       shift.getMap(+1,2).goffset(rb[cb]).slice(),recv_buf[3],psi_ptr,u2_ptr,
			       shift.getMap(-1,1).goffset(rb[cb]).slice(),recv_buf[4],psi_ptr,u1_ptr,
			       shift.getMap(+1,1).goffset(rb[cb]).slice(),recv_buf[5],psi_ptr,u1_ptr,
			       shift.getMap(-1,0).goffset(rb[cb]).slice(),recv_buf[6],psi_ptr,u0_ptr,
			       shift.getMap(+1,0).goffset(rb[cb]).slice(),recv_buf[7],psi_ptr,u0_ptr);

	  }

	}

	break;

      case MINUS:

#pragma omp parallel default(shared) 
	{
	  int threads_num = omp_get_num_threads();
	  int myId = omp_get_thread_num();

	  if (comms[0].do_comms)
	    {
	      {
		const Map& map = shift.getMap(-1,3);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_M_3_1(low,high,0,true,0,soffset_slice,send_buf[0],psi_ptr,u3_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(0);

	      {
		const Map& map = shift.getMap(+1,3);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_P_3_1(low,high,0,true,0,soffset_slice,send_buf[1],psi_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(1);
	    }

	  if (comms[2].do_comms) 
	    {
	      {
		const Map& map = shift.getMap(-1,2);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_M_2_1(low,high,0,true,0,soffset_slice,send_buf[2],psi_ptr,u2_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(2);

	      {
		const Map& map = shift.getMap(+1,2);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_P_2_1(low,high,0,true,0,soffset_slice,send_buf[3],psi_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(3);
	    }

	  if (comms[4].do_comms) 
	    {
	      {
		const Map& map = shift.getMap(-1,1);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_M_1_1(low,high,0,true,0,soffset_slice,send_buf[4],psi_ptr,u1_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(4);

	      {
		const Map& map = shift.getMap(+1,1);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_P_1_1(low,high,0,true,0,soffset_slice,send_buf[5],psi_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(5);
	    }

	  if (comms[6].do_comms) 
	    {
	      {
		const Map& map = shift.getMap(-1,0);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_M_0_1(low,high,0,true,0,soffset_slice,send_buf[6],psi_ptr,u0_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(6);

	      {
		const Map& map = shift.getMap(+1,0);
		const int* soffset_slice = map.soffset(rb[cb]).slice();
		int soffset_num = map.soffset(rb[cb]).size();

		int low = soffset_num*myId/threads_num;
		int high = soffset_num*(myId+1)/threads_num;

		func_gather_P_0_1(low,high,0,true,0,soffset_slice,send_buf[7],psi_ptr);
	      }

#pragma omp barrier
	      if (myId == comm_thread)
		comms_send_receive(7);
	    }

	  {
	    int low = innerCount[cb]*myId/threads_num;
	    int high = innerCount[cb]*(myId+1)/threads_num;

	    func_dslash_____1( low , high , 0 , false , 0 , innerSites[cb] , chi_ptr ,
			       shift.getMap(-1,3).goffset(rb[cb]).slice(), NULL, psi_ptr, u3_ptr ,
			       shift.getMap(+1,3).goffset(rb[cb]).slice(), NULL, psi_ptr, u3_ptr ,
			       shift.getMap(-1,2).goffset(rb[cb]).slice(), NULL, psi_ptr, u2_ptr ,
			       shift.getMap(+1,2).goffset(rb[cb]).slice(), NULL, psi_ptr, u2_ptr ,
			       shift.getMap(-1,1).goffset(rb[cb]).slice(), NULL, psi_ptr, u1_ptr ,
			       shift.getMap(+1,1).goffset(rb[cb]).slice(), NULL, psi_ptr, u1_ptr ,
			       shift.getMap(-1,0).goffset(rb[cb]).slice(), NULL, psi_ptr, u0_ptr ,
			       shift.getMap(+1,0).goffset(rb[cb]).slice(), NULL, psi_ptr, u0_ptr );
	  }

	  if (myId == comm_thread)
	    comms_wait();
#pragma omp barrier

	  {	
	    int low = faceCount[cb]*myId/threads_num;
	    int high = faceCount[cb]*(myId+1)/threads_num;

	    func_dslash_____1( low , high , 0 , false , 0 , faceSites[cb] , chi_ptr ,
			       shift.getMap(-1,3).goffset(rb[cb]).slice(),recv_buf[0],psi_ptr,u3_ptr,
			       shift.getMap(+1,3).goffset(rb[cb]).slice(),recv_buf[1],psi_ptr,u3_ptr,
			       shift.getMap(-1,2).goffset(rb[cb]).slice(),recv_buf[2],psi_ptr,u2_ptr,
			       shift.getMap(+1,2).goffset(rb[cb]).slice(),recv_buf[3],psi_ptr,u2_ptr,
			       shift.getMap(-1,1).goffset(rb[cb]).slice(),recv_buf[4],psi_ptr,u1_ptr,
			       shift.getMap(+1,1).goffset(rb[cb]).slice(),recv_buf[5],psi_ptr,u1_ptr,
			       shift.getMap(-1,0).goffset(rb[cb]).slice(),recv_buf[6],psi_ptr,u0_ptr,
			       shift.getMap(+1,0).goffset(rb[cb]).slice(),recv_buf[7],psi_ptr,u0_ptr);

	  }

	}

      break;
    }

    LLVMWilsonDslashT<T,P,Q>::getFermBC().modifyF(chi, QDP::rb[cb]);
  }





  // typedef LLVMWilsonDslashT<LatticeFermion,
  // 			   multi1d<LatticeColorMatrix>,
  // 			   multi1d<LatticeColorMatrix> > LLVMWilsonDslash;


  // typedef LLVMWilsonDslashT<LatticeFermionF,
  // 			   multi1d<LatticeColorMatrixF>,
  // 			   multi1d<LatticeColorMatrixF> > LLVMWilsonDslashF;

  typedef LLVMWilsonDslashT<LatticeFermionD,
			   multi1d<LatticeColorMatrixD>,
			   multi1d<LatticeColorMatrixD> > LLVMWilsonDslashD;

} // End Namespace Chroma


#endif
