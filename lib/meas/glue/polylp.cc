// $Id: polylp.cc,v 1.1 2003-10-08 18:40:57 edwards Exp $
/*! \file
 *  \brief Calculate the global normalized sum of the Polyakov loop
 */

#include "chromabase.h"
#include "meas/glue/mesplq.h"


using namespace QDP;

//! Return the value of the average plaquette normalized to 1
/*!
 * \ingroup glue
 *
 * \param u          gauge field (Read)
 * \param poly_loop  Polyakov loop average in direction mu (Write) 
 * \param mu         direction of Polyakov loop (Read)
 */

void polylp(const multi1d<LatticeColorMatrix>& u, DComplex& poly_loop, int mu)
{
  LatticeColorMatrix u_dble;
  LatticeColorMatrix poly;

  START_CODE("polylp");
        
  // Initial Polyakov loop
  poly = u[mu];

  for(int n = 1; n < Layout::lattSize()[mu]; ++n)  // run over all links in mu dir
  {    
    LatticeColorMatrix tmp = shift(poly, FORWARD, mu);
    poly = u[mu] * tmp;
  }

  /* Take the trace and sum up */
  poly_loop = sum(trace(poly)) / Double(Nc*Layout::vol());

  END_CODE("polylp");
}
