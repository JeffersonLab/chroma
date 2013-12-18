// -*- C++ -*-
// $Id: syssolver_mdagm_OPTeigcg.h,v 3.2 2009-06-02 15:56:40 bjoo Exp $
/*! \file
 *  \brief Solve a M^dag*M*psi=chi linear system by EigCG
 */

#ifndef __syssolver_mdagm_OPTeigcg_h__
#define __syssolver_mdagm_OPTeigcg_h__
#include "chroma_config.h"

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "lmdagm.h"
#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_OPTeigcg_params.h"
#include "actions/ferm/invert/containers.h"

#include "util/info/unique_id.h"

namespace Chroma
{

  //! Eigenvector accelerated CG system solver namespace
  namespace MdagMSysSolverOptEigCGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2 with eigenvectors
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMSysSolverOptEigCG : public MdagMSystemSolver<T>
  {
  public:

    //! Write out an OptEigInfo Type                         
    void QIOWriteOptEvecs(){
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      // A shorthand for the object                                      
      const LinAlg::OptEigInfo& obj =
	TheNamedObjMap::Instance().getData<LinAlg::OptEigInfo>(invParam.eigen_id);

      // File XML                                            
      XMLBufferWriter file_xml;
      push(file_xml, "OptEigInfo");
      write(file_xml, "id", uniqueId());
      write(file_xml, "N", obj.N);
      write(file_xml, "ncurEvals", obj.ncurEvals);
      write(file_xml, "restartTol", obj.restartTol);
      write(file_xml, "lde", obj.lde);
      write(file_xml, "ldh", obj.evals.size());
      pop(file_xml);

      // Open file                                  
      QDPFileWriter to(file_xml,
		       invParam.file.file_name,
		       invParam.file.file_volfmt,
		       QDPIO_SERIAL,QDPIO_OPEN);

      for(int v(0);v<obj.ncurEvals;v++){
	LatticeFermion lf ;
	obj.CvToLatFerm(lf,subset(),v);
	XMLBufferWriter record_xml;
	push(record_xml, "EigenVector");
	write(record_xml,"no",v);
	pop(record_xml);
	write(to, record_xml, lf);
      }
      
      {
	XMLBufferWriter record_xml;
	push(record_xml, "EigenValues");
	pop(record_xml);
	write(to, record_xml, obj.evals);
      }
      {
	XMLBufferWriter record_xml;
	push(record_xml, "H");
	pop(record_xml);
	write(to, record_xml, obj.H);
      }
      {
	XMLBufferWriter record_xml;
	push(record_xml, "HU");
	pop(record_xml);
	write(to, record_xml, obj.HU);
      }

      // Close                                                 
      close(to);
      swatch.stop();
      QDPIO::cout<<" QIOWriteOptEvecs: Time to write evecs= "
		 << swatch.getTimeInSeconds() <<" secs "<<endl ;
    }

    //-----------------------------------------------------------------------
    //! Read a OptEigInfo Type                                           
    void QIOReadOptEvecs()
    {

      StopWatch swatch;
      swatch.reset();
      swatch.start();

      // File XML                                      
      XMLReader file_xml;

      // Open file     
      QDPFileReader to(file_xml,invParam.file.file_name,QDPIO_SERIAL);

      // A shorthand for the object         
      LinAlg::OptEigInfo& obj =
	TheNamedObjMap::Instance().getData<LinAlg::OptEigInfo>(invParam.eigen_id);

      XMLReader record_xml;
      //TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);

      int ldh,lde,N ;
      float restartTol ;
      read(file_xml, "/OptEigInfo/ldh", ldh);
      read(file_xml, "/OptEigInfo/lde", lde);
      read(file_xml, "/OptEigInfo/N",    N);
      read(file_xml, "/OptEigInfo/ncurEvals", obj.ncurEvals);
      read(file_xml, "/OptEigInfo/restartTol", restartTol);

      QDPIO::cout <<__func__ << " : Reading object with following properties"<<endl ; 
      QDPIO::cout <<__func__ << " :   lde= " << lde << endl;
      QDPIO::cout <<__func__ << " :   ldh= " << ldh << endl;
      QDPIO::cout <<__func__ << " :   N  = " << N   << endl;

      QDPIO::cout << __func__ << " :   ncurEvals  = " << obj.ncurEvals<< endl;
      QDPIO::cout << __func__ << " :   restartTol = " << restartTol<< endl;
      
      if(obj.evals.size()< obj.ncurEvals){
	QDPIO::cerr<<__func__<< " : ldh of the current object is not large enough to hold the vectors" ;
	QDP_abort(1);
      }
      
      for(int v(0);v<obj.ncurEvals;v++){
        LatticeFermion lf ;
	XMLReader record_xml;
	read(to, record_xml, lf);
	obj.CvToEigCGvec(lf, subset(), v) ;
      }
      //this is just fine
      {
	XMLReader record_xml;
	multi1d<Real>    evals(ldh) ;
	multi1d<Complex> H(ldh*ldh)     ;
	multi1d<Complex> HU(ldh*ldh)    ;
	read(to, record_xml, evals);
	read(to, record_xml, H);
	read(to, record_xml, HU);
	if(ldh<=obj.evals.size()){
	  for(int i(0);i<ldh;i++)
	    obj.evals[i] = evals[i] ;
	  for(int i(0);i<H.size();i++){
	    obj.H[i] = H[i] ;
	    obj.HU[i] = HU[i] ;
	  }
	}
	else{
	  QDPIO::cerr<<__func__<< " : ldh of the current object is not large enough to hold the Cholesky factors" ;
	  QDP_abort(1);
	}
      }

    // Close                             
    close(to);

    swatch.stop();
    QDPIO::cout<<" QIOReadOptEvecs: Time to read evecs= "
	       << swatch.getTimeInSeconds() <<" secs "<<endl ;
  }

    //! Constructor
    /*!
     * \param M_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    MdagMSysSolverOptEigCG(Handle< LinearOperator<T> > A_,
			   const SysSolverOptEigCGParams& invParam_) : 
      MdagM(new MdagMLinOp<T>(A_)), A(A_), invParam(invParam_) 
      {
#ifndef QDP_IS_QDPJIT
	numMatvecs = 0 ;
	// NEED to grab the eignvectors from the named buffer here
	if (! TheNamedObjMap::Instance().check(invParam.eigen_id))
	{
	  TheNamedObjMap::Instance().create< LinAlg::OptEigInfo >(invParam.eigen_id);
	  LinAlg::OptEigInfo& EigInfo = 
	    TheNamedObjMap::Instance().getData< LinAlg::OptEigInfo >(invParam.eigen_id);
	  int N = Layout::sitesOnNode()*Nc*Ns ;
	  int VectorSpaceSize =  Nc*Ns*(A->subset()).numSiteTable();
	  EigInfo.init(invParam.Neig_max, N, VectorSpaceSize) ;
	  EigInfo.restartTol =  invParam.restartTol.elem().elem().elem().elem();
	  if(invParam.file.read){
	    QDPIO::cout<<"MdagMSysSolverOptEigCG : reading evecs from disk"<<endl ;
	    QIOReadOptEvecs() ;
	  }
	}
#endif
      }

    //! Destructor is automatic
    ~MdagMSysSolverOptEigCG()
      {
	if(invParam.file.write){
	  QDPIO::cout<<"MdagMSysSolverOptEigCG : writing evecs to disk"<<endl ;
	  QIOWriteOptEvecs() ;
	}
	if (invParam.cleanUpEvecs)
	{
	  TheNamedObjMap::Instance().erase(invParam.eigen_id);
	}
      }

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     *
     * Definitions supplied in the correspond .cc file
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const;


    //! Solve the linear system starting with a chrono guess 
    /*! 
     * \param psi solution (Write)
     * \param chi source   (Read)
     * \param predictor   a chronological predictor (Read)
     * \return syssolver results
     */

    SystemSolverResults_t operator()(T& psi, const T& chi, 
				     AbsChronologicalPredictor4D<T>& predictor) const 
    {
      
      START_CODE();

      // This solver uses InvCG2, so A is just the matrix.
      // I need to predict with A^\dagger A
      {
	Handle< LinearOperator<T> > MdagM( new MdagMLinOp<T>(A) );
	predictor(psi, (*MdagM), chi);
      }
      // Do solve
      SystemSolverResults_t res=(*this)(psi,chi);

      // Store result
      predictor.newVector(psi);
      END_CODE();
      return res;
    }


  private:

    // Hide default constructor
    MdagMSysSolverOptEigCG() {}
    int numMatvecs ;
    Handle< LinearOperator<T> > MdagM;
    Handle< LinearOperator<T> > A;
    SysSolverOptEigCGParams invParam;
  };

} // End namespace



#endif 

