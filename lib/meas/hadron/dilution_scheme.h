// -*- C++ -*-
// $Id: dilution_scheme.h,v 1.8 2008-04-21 03:19:35 edwards Exp $
/*! \file
 *  \brief Dilution Schemes
 */

#ifndef __dilution_scheme_h__
#define __dilution_scheme_h__

#include "chromabase.h"

namespace Chroma
{
  //! Abstract dilution scheme
  /*! @ingroup hadron
   *
   * Supports creation of (abstract) dilution schemes used in 
   * stochastic sources and solutions
   */
  template<typename T>
  class DilutionScheme
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~DilutionScheme() {}
 
    //! The decay direction
    virtual int getDecayDir() const = 0;

    //! The seed identifies this quark
    virtual const Seed& getSeed() const = 0;

    virtual int getT0(int t0) const = 0 ;
		
    virtual int getDilSize(int t0) const = 0 ;

    virtual int getNumTimeSlices() const = 0;
	
    virtual Real getKappa() const = 0;

    virtual std::string getCfgInfo() const = 0;

    virtual std::string getPropHeader(int t0, int dil) const = 0; 
    
    virtual std::string getSourceHeader(int t0, int dil) const = 0; 

    //! Return the diluted source vector
    /*! MAYBE THIS SHOULD BE A CONST REFERENCE?? PROBABLY NO */
    virtual T dilutedSource(int t0, int dil ) const = 0;
    
    //! Return the solution vector corresponding to the diluted source
    /*! MAYBE THIS SHOULD BE A CONST REFERENCE?? POSSIBLY YES */
    virtual T dilutedSolution(int t0, int dil ) const = 0;

  };

} // namespace Chroma


#endif
