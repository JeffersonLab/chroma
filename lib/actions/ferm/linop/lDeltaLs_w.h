// -*- C++ -*-
// $Id: lDeltaLs_w.h,v 1.1 2004-10-15 16:37:19 bjoo Exp $
/*! \file
 *   Apply DeltaLs = (1/4)(1-eps^2) to a fermion field.

 */

#ifndef __ldeltaLs_w_h__
#define __ldeltaLs_w_h__

#include "linearop.h"
#include "fermact.h" 


using namespace QDP;
using namespace Chroma; 

namespace Chroma {

class lDeltaLs: public LinearOperator<LatticeFermion>
{
public:
  //! Creation routine
  /*!
   * \ingroup linop
   *
   */
  lDeltaLs(Handle< const LinearOperator<LatticeFermion> > g5eps_ ) :
    g5eps(g5eps_) {}

  //! Destructor is automatic
  ~lDeltaLs() {}
 
  //! Only defined on the entire lattice
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

private:
  Handle<const LinearOperator<LatticeFermion> > g5eps;
};

}
#endif
