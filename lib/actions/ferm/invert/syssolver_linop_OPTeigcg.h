// -*- C++ -*-
// $Id: syssolver_linop_OPTeigcg.h,v 1.6 2008-12-15 05:02:06 kostas Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_OPTeigcg_h__
#define __syssolver_linop_OPTeigcg_h__

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
  namespace LinOpSysSolverOptEigCGEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2 with eigenvectors
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverOptEigCG : public LinOpSystemSolver<T>
  {
  public:

    //! Write out an OptEigInfo Type                                                                           
    void QIOWriteOptEvecs(){
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

    }

    //-----------------------------------------------------------------------                                  
    //! Read a OptEigInfo Type                                                                                 
    void QIOReadOptEvecs()
    {
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
      
      //the following loop still needs work...
      for(int v(0);v<obj.ncurEvals;v++){
        LatticeFermion lf ;
	XMLReader record_xml;
	read(to, record_xml, lf);
      }
      //this is just fine
      {
	XMLReader record_xml;
	read(to, record_xml, obj.evals);
	read(to, record_xml, obj.H);
	read(to, record_xml, obj.HU);
      }

    // Close                                                                                                 
    close(to);
  }

    //! Constructor
    /*!
     * \param M_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    LinOpSysSolverOptEigCG(Handle< LinearOperator<T> > A_,
			const SysSolverOptEigCGParams& invParam_) : 
      MdagM(new MdagMLinOp<T>(A_)), A(A_), invParam(invParam_) 
      {
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
	    QDPIO::cout<<"LinOpSysSolverOptEigCG : reading evecs from disk"<<endl ;
	  }
	}
      }

    //! Destructor is automatic
    ~LinOpSysSolverOptEigCG()
      {
	if(invParam.file.write){
	  QDPIO::cout<<"LinOpSysSolverOptEigCG : writing evecs to disk"<<endl ;
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

  private:

    // Hide default constructor
    LinOpSysSolverOptEigCG() {}
    int numMatvecs ;
    Handle< LinearOperator<T> > MdagM;
    Handle< LinearOperator<T> > A;
    SysSolverOptEigCGParams invParam;
  };

} // End namespace



#endif 

