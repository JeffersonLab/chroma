// -*- C++ -*-
// $Id: lDeltaLs_w.h,v 1.2 2004-12-09 03:58:03 edwards Exp $
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
  lDeltaLs(Handle< const LinearOperator<LatticeFermion> > D_ ) :
    D(D_) {}

  //! Destructor is automatic
  ~lDeltaLs() {}
 
  //! Only defined on the entire lattice
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

private:
  Handle<const LinearOperator<LatticeFermion> > D;
};

}

using namespace Chroma;

#endif
