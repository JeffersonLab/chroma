// -*- C++ -*-
// $Id: lcm_4mn5fp_recursive.h,v 3.1 2006-12-21 22:46:52 bjoo Exp $

/*! @file
 * @brief Lat Col Mat 4th order 5 force calculation minimum norm integrator
 *
 * The 4th order 5 force calculation minimum norm integrator (velocity variant)
 * from the paper by Takaishi and deForcrand (eq 24-25). Made recursive through
 * replacing the Momentum update term with a subIntegrator call.
 */

#ifndef LCM_4MN5FP_RECURSIVE_H
#define LCM_4MN5FP_RECURSIVE_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h"

namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMat4MN5FPRecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMat4MN5FPRecursiveIntegratorParams
  {
    LatColMat4MN5FPRecursiveIntegratorParams();
    LatColMat4MN5FPRecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    Real theta;
    Real rho;
    Real lambda;
    Real mu;

    multi1d<std::string> monomial_ids;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMat4MN5FPRecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMat4MN5FPRecursiveIntegratorParams& p);

  //! MD integrator interface for 4th order 5 Force Min. Norm. Integrator (Velocity variant)
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMat4MN5FPRecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:
    
    // Construct from params struct and Hamiltonian
    LatColMat4MN5FPRecursiveIntegrator(
               const LatColMat4MN5FPRecursiveIntegratorParams& p
	  )  : params(p), SubIntegrator(IntegratorShared::createSubIntegrator(p.subintegrator_xml)) {
      IntegratorShared::bindMonomials(p.monomial_ids, monomials);
    }


    // Copy constructor
    LatColMat4MN5FPRecursiveIntegrator(const LatColMat4MN5FPRecursiveIntegrator& l) :
      params(l.params), monomials(l.monomials), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMat4MN5FPRecursiveIntegrator(void) {};


    void operator()( AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s, 
		     const Real& traj_length) const;
   			    
    AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
			   multi1d<LatticeColorMatrix> >& getSubIntegrator() const {
      return (*SubIntegrator);
    }

  protected:
    //! Refresh fields in just this level
    void refreshFieldsThisLevel(AbsFieldState<multi1d<LatticeColorMatrix>,
				multi1d<LatticeColorMatrix> >& s) const {
      for(int i=0; i < monomials.size(); i++) { 
	monomials[i]->refreshInternalFields(s);
      }
    }
  private:
    LatColMat4MN5FPRecursiveIntegratorParams params;


    multi1d< Handle< Monomial< multi1d<LatticeColorMatrix>,
			       multi1d<LatticeColorMatrix> > > > monomials;

    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;
	              

    
  };

}


#endif
