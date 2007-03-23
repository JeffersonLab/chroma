// -*- C++ -*-
// $Id: lcm_toplevel_integrator.h,v 3.4 2007-03-23 20:21:39 bjoo Exp $

#ifndef __LCM_TOPLEVEL_INTEGRATOR_H_
#define __LCM_TOPLEVEL_INTEGRATOR_H_

#include "chromabase.h"
#include "update/molecdyn/integrator/abs_integrator.h"

using namespace QDP;

namespace Chroma {

  //! A Struct for holding pairs of IDs in the copy list
  struct IDPair {
    std::string source;
    std::string dest;
  };
  
  //! A reader for an element of the copy list
  void read(XMLReader& xml, const std::string& path, IDPair& p);

  //! A Writer for an element of the copy list
  void write(XMLWriter& xml, const std::string& path, const IDPair& p);


  //! A Structure to hold the top level parameters
  struct LCMToplevelIntegratorParams { 

    //! Construct from XML
    LCMToplevelIntegratorParams(XMLReader& xml, const std::string& path);
    Real tau0;
    multi1d<IDPair> copy_list;
    std::string integrator_xml;
    bool anisoP;
    int t_dir;
    Real xi_mom;

  };

  //! Read the Integrator Params
  void read(XMLReader& xml, const std::string& path, LCMToplevelIntegratorParams& p);

  //! Write the Integrator Params
  void write(XMLWriter& xml, const std::string& path, const LCMToplevelIntegratorParams& p);

  class LCMToplevelIntegrator : public 
  AbsMDIntegrator< multi1d<LatticeColorMatrix>,
	           multi1d<LatticeColorMatrix> >
  {
  public:

    //! Construct from parameters.
    LCMToplevelIntegrator(const LCMToplevelIntegratorParams& p);

    //! Copy Constructor
    LCMToplevelIntegrator(const LCMToplevelIntegrator& i) : params(i.params), top_integrator(i.top_integrator) {}

    //! Destructor is automagic
    ~LCMToplevelIntegrator() {} 

    //! Get the length of a trajectory
    Real getTrajLength(void) const { 
      return params.tau0;
    }

    //! Copy fields between the monomials of a copy list
    void copyFields(void) const;

    AbsComponentIntegrator< multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> >& getIntegrator(void) const { 
      return *top_integrator;
    }

  private:
    LCMToplevelIntegratorParams params;
    Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, 
				    multi1d<LatticeColorMatrix> > > top_integrator;
  };


};


#endif
