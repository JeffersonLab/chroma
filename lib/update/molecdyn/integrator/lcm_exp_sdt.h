// -*- C++ -*-
// $Id: lcm_exp_sdt.h,v 3.1 2006-11-07 23:11:02 bjoo Exp $

/*! @file
 * @brief Intgrator for exp(S dt)
 *
 * A component integrator to integrate exp(Sdt) for some 
 * monomial list S (ie this is a leapP like component)
 */

#ifndef LCM_EXPSDT_H
#define LCM_EXPSDT_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"


namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatExpSdtIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatExpSdtIntegratorParams
  {
    LatColMatExpSdtIntegratorParams();
    LatColMatExpSdtIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    multi1d<int> monomial_list;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatExpSdtIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatExpSdtIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatExpSdtIntegrator 
    : public AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatExpSdtIntegrator(const int  n_steps_, 
			      const multi1d<int>& monomial_list_,
			      const Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_ ) : n_steps(n_steps_), monomial_list(monomial_list_), H_MD(H_) {};

    // Construct from params struct and Hamiltonian
    LatColMatExpSdtIntegrator(const LatColMatExpSdtIntegratorParams& p,
				   Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) : n_steps(p.n_steps), monomial_list(p.monomial_list), H_MD(H_) {}

    // Copy constructor
    LatColMatExpSdtIntegrator(const LatColMatExpSdtIntegrator& l) :
      n_steps(l.n_steps), monomial_list(l.monomial_list), H_MD(l.H_MD) {}

    // ! Destruction is automagic
    ~LatColMatExpSdtIntegrator(void) {};


    void getMonomialList(multi1d<int>& ml) const {
      ml.resize(monomial_list.size());
      for(int i=0; i < monomial_list.size(); i++) { 
	ml[i]=monomial_list[i];
      }

    }


    void operator()( AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s, 
		     const Real& traj_length) const;
   			    

  private:
    int  n_steps;
    const multi1d<int> monomial_list;
    const Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > H_MD;
    
    
  };

}


#endif
