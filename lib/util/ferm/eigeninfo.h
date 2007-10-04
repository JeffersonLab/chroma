// -*- C++ -*-
// $Id: eigeninfo.h,v 3.1 2007-10-04 13:53:57 edwards Exp $
/*! \file
 *  \brief Hold eigenvalues and eigenvectors
 */

#ifndef EIGENSTATE_H
#define EIGENSTATE_H


namespace Chroma
{
 
  //! Hold eigenvalues and eigenvectors
  /*!
   * \ingroup ferm
   */
  class EigenInfo
  {
  public:
    EigenInfo() {}

    EigenInfo(multi1d<Real>& eval, Real& l, multi1d<LatticeFermion>& evec): evalues(eval), largest(l), evectors(evec) {
      if (eval.size() != evec.size() ) {
	QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	QDP_abort(1);
      }
    }

    ~EigenInfo() {}
      
    EigenInfo(const EigenInfo& e): evalues(e.evalues), largest(e.largest), evectors(e.evectors) {
      if (e.evalues.size() != e.evectors.size() ) {
	QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	QDP_abort(1);
      }
    }


    const multi1d<Real>& getEvalues() const { 
      return evalues;
	
    }

    multi1d<Real>& getEvalues() {
      return evalues;
    }


    const Real& getLargest() const { 
      return largest;
    }

    Real& getLargest() {
      return largest;
    }

    const multi1d<LatticeFermion>& getEvectors() const { 
      return evectors;
    }

    multi1d<LatticeFermion>& getEvectors() {
      return evectors;
    }

  private:
    multi1d<Real> evalues;
    multi1d<LatticeFermion> evectors;
    Real largest;
  };
}

















#endif
