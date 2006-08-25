// -*- C++ -*-
// $Id: stag_phases_s.h,v 3.1 2006-08-25 23:46:37 edwards Exp $

#ifndef STAG_PHASES_H
#define STAG_PHASES_H

#include <chromabase.h>

namespace Chroma 
{
  namespace StagPhases 
  {
    // These are the K-S Phases (aka alpha in Eduardo's notation)
    /*! \ingroup gauge */
    class alphaClass 
    {
    private:
      static multi1d<LatticeInteger> phases;

      // initP should be automatically initialised to 0
      // I don't know how I can force this
      static bool initP;
      static void init();
    public:
      static const LatticeInteger& alpha(const int dim);
    };


    // These are the Mesonic Phases (a.k.a beta in Eduardo's notation)
    /*! \ingroup gauge */
    class betaClass 
    {
    private:
      static multi1d<LatticeInteger> phases;

      // initP should be automatically initialised to 0
      // I don't know how I can force this
      static bool initP;
      static void init();
    public:
      static const LatticeInteger& beta(const int dim);
    };

    static const LatticeInteger& alpha(const int dim) {
      return alphaClass::alpha(dim);
    }

    static const LatticeInteger& beta(const int dim) {
      return betaClass::beta(dim);
    }

  }

}  // end namespace Chroma


#endif
