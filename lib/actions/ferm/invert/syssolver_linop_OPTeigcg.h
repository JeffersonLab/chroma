// -*- C++ -*-
// $Id: syssolver_linop_OPTeigcg.h,v 1.2 2008-04-01 04:02:28 kostas Exp $
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
//#include "actions/ferm/invert/containers.h"

namespace Chroma
{

  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverOptEigCGEnv
  {
    //! Name to be used
    extern const std::string name;

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
	}
	esize = invParam.esize*EigInfo.N ;
	NcNs = Nc*Ns ;
      }

    //! Destructor is automatic
    ~LinOpSysSolverOptEigCG()
      {
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
    //The MatVec Interface

    T XX ;
    T YY ;
    
    void MatrixMatvec(void *x, void *y, void *params) {
      //Works only in single precision CHROMA
      RComplex<float> *px = (RComplex<float> *) x;
      RComplex<float> *py = (RComplex<float> *) y;

      //XX.getF() = (Complex *) x ; //AN DOULEPSEI AUTO NA ME FTUSEIS
      //YY.getF() = (Complex *) y ; //AN DOULEPSEI AUTO NA ME FTUSEIS

      //Alliws kanoume copy
      //copy x into XX
      if((A->subset()).hasOrderedRep()){
	int count=0 ;
	//can be done with ccopy for speed...
	for(int i=s.start(); i <= s.end(); i++)
	  for(int s(0);s<Ns;s++)
	    for(int c(0);c<Nc;c++){
	      XX.elem(i).elem(s).elem(c) = *(px+count);
	      count++;
	    }
      }
      else{
	int i ;
	const int *tab = (A->subset()).siteTable().slice();
	int count=0;
	for(int x=0; x < (A->subset()).numSiteTable(); ++x){
	  i = tab[x] ;
	  for(int s(0);s<Ns;s++)
	    for(int c(0);c<Nc;c++){
	      XX.elem(i).elem(s).elem(c) = *(px+count);
	      count++;
	    }
	}
      }

      MdagM(YY,XX,PLUS) ;

      //copy back..
      if((A->subset()).hasOrderedRep()){
	int count=0 ;
	//can be done with ccopy for speed...
	for(int i=s.start(); i <= s.end(); i++)
	  for(int s(0);s<Ns;s++)
	    for(int c(0);c<Nc;c++){
	      *(px+count) = XX.elem(i).elem(s).elem(c) ;
	      count++;
	    }
      }
      else{
	int count=0;
	for(int x=0; x < (A->subset()).numSiteTable(); ++x){
	  i = tab[x] ;
	  for(int s(0);s<Ns;s++)
	    for(int c(0);c<Nc;c++){
	      *(py+count) = YY.elem(i).elem(s).elem(c) ;
	      count++;
	    }
	}
      }

      numMatvecs++;

    }
    
    int NcNs ; // This is Nc*Ns might help speed up the copy
      
    int esize ;

    // Hide default constructor
    LinOpSysSolverOptEigCG() {}
    int numMatvecs ;
    Handle< LinearOperator<T> > MdagM;
    Handle< LinearOperator<T> > A;
    SysSolverOptEigCGParams invParam;
  };

} // End namespace

#endif 

