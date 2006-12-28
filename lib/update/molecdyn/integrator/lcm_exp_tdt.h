// -*- C++ -*-
// $Id: lcm_exp_tdt.h,v 3.3 2006-12-28 17:34:00 bjoo Exp $

/*! @file
 * @brief Intgrator for exp(T dt)
 *
 * A component integrator to integrate exp(T dt) for the kinetic term (momenta)
 * (ie this is a leapQ like component)
 */

#ifndef LCM_EXPTDT_H
#define LCM_EXPTDT_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"


namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatExpTdtIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatExpTdtIntegratorParams
  {
    LatColMatExpTdtIntegratorParams();
    LatColMatExpTdtIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatExpTdtIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatExpTdtIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatExpTdtIntegrator 
    : public AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatExpTdtIntegrator(int  n_steps_) : n_steps(n_steps_) {};

    // Construct from params struct and Hamiltonian
    LatColMatExpTdtIntegrator(const LatColMatExpTdtIntegratorParams& p) : 
      n_steps(p.n_steps) {}

    // Copy constructor
    LatColMatExpTdtIntegrator(const LatColMatExpTdtIntegrator& l) :
      n_steps(l.n_steps) {}

    // ! Destruction is automagic
    ~LatColMatExpTdtIntegrator(void) {};


    void operator()( AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s, 
		     const Real& traj_length) const;
   			    

    void refreshFields(AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s) const
    {}

    void resetPredictors(void) const {}

  private:
    int  n_steps;
  };

}


#endif
