// -*- C++ -*-
// $Id: integrator_shared.h,v 3.1 2006-11-20 19:15:02 bjoo Exp $

#ifndef __INTEGRATOR_SHARED_H__
#define __INTEGRATOR_SHARED_H__

#include "chromabase.h"
#include "handle.h"
#include <string>
#include "update/molecdyn/monomial/abs_monomial.h"
#include "meas/inline/io/named_objmap.h"
#include "update/molecdyn/integrator/md_integrator_factory.h"

using namespace std;
using namespace QDP;
 
namespace Chroma { 

  namespace IntegratorShared {
    struct MonomialPair {
	Handle< Monomial< multi1d<LatticeColorMatrix>, 
	                  multi1d<LatticeColorMatrix> > > mon;
	std::string id;
    };
    
    //! A routine to bind Monomial IDs to an array of Monomial Handles
    void bindMonomials(const multi1d<std::string>& monomial_ids,
		       multi1d< MonomialPair >& monomials);

    //! A routine to create a sub integrator from a generic piece of subintegrator XML
    AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
			    multi1d<LatticeColorMatrix> 
                          >* 
    createSubIntegrator(const std::string& subintegrator_xml);
    
  }

}
#endif

