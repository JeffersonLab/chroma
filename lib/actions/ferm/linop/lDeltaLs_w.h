// -*- C++ -*-
// $Id: lDeltaLs_w.h,v 3.2 2007-02-22 21:11:46 bjoo Exp $
/*! \file
 *   Apply DeltaLs = (1/4)(1-eps^2) to a fermion field.

 */

#ifndef __ldeltaLs_w_h__
#define __ldeltaLs_w_h__

#include "linearop.h"
#include "handle.h"

namespace Chroma 
{

  //! GW Defect operator
  /*! \ingroup linop */
  class lDeltaLs: public LinearOperator<LatticeFermion>
  {
  public:
    //! Creation routine
    lDeltaLs(Handle< LinearOperator<LatticeFermion> > D_ ) :
      D(D_) {}

    //! Destructor is automatic
    ~lDeltaLs() {}
 
    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

  private:
    Handle< LinearOperator<LatticeFermion> > D;
  };

}


#endif
