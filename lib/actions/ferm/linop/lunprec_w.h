// -*- C++ -*-

#ifndef __lunprec_h__
#define __lunprec_h__

#include "handle.h"
#include "linearop.h"
#include "eoprec_linop.h"

namespace Chroma 
{ 
  //! Gamma(5) hermitian linear operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * This operator scales its input operator
   */
  template<typename T,typename P,typename Q>
  class Lunprec : public LinearOperator<T>
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    Lunprec(LinearOperator<T>* p) : A(p) {}

    //! Copy pointer (one more owner)
    Lunprec(Handle< LinearOperator<T> > p): A(p) {}

    //! Destructor
    ~Lunprec() {}

    //! Subset comes from underlying operator
    inline const Subset& subset() const {return A->subset();}

    //! Apply the operator onto a source std::vector
    /*! For this operator, the sign is ignored */
    inline void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {

      try { 
	EvenOddPrecLinearOperator<T,P,Q>& A_eo=
	  dynamic_cast<EvenOddPrecLinearOperator<T,P,Q>&>(*A);
	
	A_eo.unprecLinOp(chi, psi, isign);
      }
      catch(std::bad_cast& bc) {
	QDPIO::cout << "Couldnt cast Linop to Even Odd Operator" << std::endl;
	QDP_abort(1);
      }
      catch(...) { 
	QDPIO::cout << "Unknown exception occurred in lunprec::operator()" << std:: endl;
	QDP_abort(1);
      }

    }
    
  private:
    Handle< LinearOperator<T> > A;
  };


} // End Namespace Chroma


#endif
