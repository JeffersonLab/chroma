// -*- C++ -*-
// $Id: eigeninfo.h,v 3.2 2008-06-17 15:36:58 edwards Exp $
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
  template<typename T>
  class EigenInfo
  {
  public:
    //! Partial constructor
    EigenInfo() {}

    //! Full constructor
    EigenInfo(const multi1d<Real>& eval, const Real& l, const multi1d<T>& evec) 
      : evalues(eval), largest(l), evectors(evec) {
      if (eval.size() != evec.size() ) {
	QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	QDP_abort(1);
      }
    }

    //! Destructor
    ~EigenInfo() {}
      
    //! Copy constructor
    EigenInfo(const EigenInfo& e): evalues(e.evalues), largest(e.largest), evectors(e.evectors) {
      if (e.evalues.size() != e.evectors.size() ) {
	QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	QDP_abort(1);
      }
    }

    //! Getter
    const multi1d<Real>& getEvalues() const {return evalues;}
    //! Setter
    multi1d<Real>& getEvalues() {return evalues;}

    //! Getter
    const Real& getLargest() const {return largest;}
    //! Setter
    Real& getLargest() {return largest;}

    //! Getter
    const multi1d<T>& getEvectors() const {return evectors;}
    //! Setter
    multi1d<T>& getEvectors() {return evectors;}

  private:
    multi1d<Real> evalues;
    multi1d<T> evectors;
    Real largest;
  };
}

















#endif
