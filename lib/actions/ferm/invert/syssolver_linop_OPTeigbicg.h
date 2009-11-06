// -*- C++ -*-
/*! \file
 *  \brief Solve a M*psi=chi linear system by EigBiCG
 */

#ifndef __syssolver_OPTeigbicg_h__
#define __syssolver_OPTeigbicg_h__
#include "chroma_config.h"

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_OPTeigbicg_params.h"
#include "actions/ferm/invert/containers.h"

#include "util/info/unique_id.h"

namespace Chroma
{

  //! Eigenvector accelerated biCG system solver namespace
  namespace LinOpSysSolverOptEigBiCGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by biCG with eigenvectors
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverOptEigBiCG : public LinOpSystemSolver<T>
  {
  public:

    //! Write out an OptEigInfo Type                         
    void QIOWriteOptEvecs(){
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      // A shorthand for the object                                      
      //const LinAlg::OptEigBiInfo<WordType<T>::Type_t>& obj =
      //	TheNamedObjMap::Instance().getData<LinAlg::OptEigBiInfo<WordType<T>::Type_t> >(invParam.eigen_id);

      const LinAlg::OptEigBiInfo<REAL>& obj =
      	TheNamedObjMap::Instance().getData<LinAlg::OptEigBiInfo<REAL> >(invParam.eigen_id);

      // File XML                                            
      XMLBufferWriter file_xml;
      push(file_xml, "OptEigBiInfo");
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
	obj.CvToLatFerm(lf,subset(),v,"L");//left
	{
	  XMLBufferWriter record_xml;
	  push(record_xml, "LeftEigenVector");
	  write(record_xml,"no",v);
	  pop(record_xml);
	  write(to, record_xml, lf);
	}
	obj.CvToLatFerm(lf,subset(),v,"R");//Right
	{
	  XMLBufferWriter record_xml;
	  push(record_xml, "RightEigenVector");
	  write(record_xml,"no",v);
	  pop(record_xml);
	  write(to, record_xml, lf);
	}
      }
      
      {
	XMLBufferWriter record_xml;
	push(record_xml, "EigenValues");
	pop(record_xml);
	//write(to, record_xml, obj.evals);
	
	multi1d<Complex> foo(obj.evals.size()) ;
	for(int i(0);i<obj.evals.size();i++)
	  foo[i].elem().elem().elem() = obj.evals[i] ;
	write(to, record_xml, foo);
      }
      {
	XMLBufferWriter record_xml;
	push(record_xml, "H");
	pop(record_xml);
	//write(to, record_xml, obj.H);
	multi1d<Complex> foo(obj.H.size()) ;
	for(int i(0);i<obj.H.size();i++)
	  foo[i].elem().elem().elem() = obj.H[i] ;
	write(to, record_xml, foo);
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
      LinAlg::OptEigBiInfo<REAL>& obj =
	TheNamedObjMap::Instance().getData<LinAlg::OptEigBiInfo<REAL> >(invParam.eigen_id);

      XMLReader record_xml;
      //TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);

      int ldh,lde,N ;
      float restartTol ;
      read(file_xml, "/OptEigBiInfo/ldh", ldh);
      read(file_xml, "/OptEigBiInfo/lde", lde);
      read(file_xml, "/OptEigBiInfo/N",    N);
      read(file_xml, "/OptEigBiInfo/ncurEvals", obj.ncurEvals);
      read(file_xml, "/OptEigBiInfo/restartTol", restartTol);

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
	obj.CvToEigCGvec(lf, subset(), v,"L") ;
	read(to, record_xml, lf);
	obj.CvToEigCGvec(lf, subset(), v,"R") ;
      }
      //this is just fine
      {
	XMLReader record_xml;
	multi1d<Complex>    evals(ldh) ;
	multi1d<Complex> H(ldh*ldh)     ;
	read(to, record_xml, evals);
	read(to, record_xml, H);
	if(ldh<=obj.evals.size()){
	  for(int i(0);i<ldh;i++){
	    obj.evals[i] = evals[i].elem().elem().elem() ;
	  }
	  for(int i(0);i<H.size();i++){
	    obj.H[i] = H[i].elem().elem().elem() ;
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
    LinOpSysSolverOptEigBiCG(Handle< LinearOperator<T> > A_,
			     const SysSolverOptEigBiCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {
	numMatvecs = 0 ;
	// NEED to grab the eignvectors from the named buffer here
	if (! TheNamedObjMap::Instance().check(invParam.eigen_id))
	{
	  TheNamedObjMap::Instance().create< LinAlg::OptEigBiInfo<REAL> >(invParam.eigen_id);
	  LinAlg::OptEigBiInfo<REAL>& EigInfo = 
	    TheNamedObjMap::Instance().getData< LinAlg::OptEigBiInfo<REAL> >(invParam.eigen_id);
	  int N = Layout::sitesOnNode()*Nc*Ns ;
	  int VectorSpaceSize =  Nc*Ns*(A->subset()).numSiteTable();
	  EigInfo.init(invParam.Neig_max, N, VectorSpaceSize) ;
	  EigInfo.restartTol =  invParam.restartTol.elem().elem().elem().elem();
	  if(invParam.file.read){
	    QDPIO::cout<<"LinOpSysSolverOptEigBiCG : reading evecs from disk"<<endl ;
	    QIOReadOptEvecs() ;
	  }
	}
      }

    //! Destructor is automatic
    ~LinOpSysSolverOptEigBiCG()
      {
	if(invParam.file.write){
	  QDPIO::cout<<"LinOpSysSolverOptEigBiCG : writing evecs to disk"<<endl ;
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


    /*********** NO CHRONO FOR BICGSTAB ***********
    //! Solve the linear system starting with a chrono guess 
    /*! 
     * \param psi solution (Write)
     * \param chi source   (Read)
     * \param predictor   a chronological predictor (Read)
     * \return syssolver results
     */
    /*********** NO CHRONO FOR BICGSTAB ***********
    SystemSolverResults_t operator()(T& psi, const T& chi, 
				     AbsChronologicalPredictor4D<T>& predictor) const 
    {
      
      START_CODE();

      // I need to predict with A
      {
	predictor(psi, A, chi);
      }
      // Do solve
      SystemSolverResults_t res=(*this)(psi,chi);

      // Store result
      predictor.newVector(psi);
      END_CODE();
      return res;
    }
    **********************/


  private:

    // Hide default constructor
    LinOpSysSolverOptEigBiCG() {}
    int numMatvecs ;
    Handle< LinearOperator<T> > A;
    SysSolverOptEigBiCGParams invParam;
  };

} // End namespace



#endif 

