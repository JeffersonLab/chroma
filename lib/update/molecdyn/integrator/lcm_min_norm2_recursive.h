// -*- C++ -*-
// $Id: lcm_min_norm2_recursive.h,v 3.1 2006-11-08 23:46:11 bjoo Exp $

/*! @file
 * @brief Intgrator for exp(S dt)
 *
 * A component integrator to integrate exp(Sdt) for some 
 * monomial list S (ie this is a leapP like component)
 */

#ifndef LCM_MIN_NORM2_RECURSIVE_H
#define LCM_MIN_NORM2_RECURSIVE_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"


namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatMinNorm2RecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatMinNorm2RecursiveIntegratorParams
  {
    LatColMatMinNorm2RecursiveIntegratorParams();
    LatColMatMinNorm2RecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    Real lambda;
    multi1d<int> monomial_list;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatMinNorm2RecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatMinNorm2RecursiveIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatMinNorm2RecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatMinNorm2RecursiveIntegrator(int  n_steps_, 
					 multi1d<int> monomial_list_,
					 Real lambda_, 
					 Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_ ,
					 Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& SubIntegrator_) : n_steps(n_steps_), monomial_list(monomial_list_), lambda(lambda_),  H_MD(H_), SubIntegrator(SubIntegrator_) {};

    // Construct from params struct and Hamiltonian
    LatColMatMinNorm2RecursiveIntegrator(
					    const LatColMatMinNorm2RecursiveIntegratorParams& p,
					    Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) : n_steps(p.n_steps), monomial_list(p.monomial_list), lambda(p.lambda), H_MD(H_), SubIntegrator(createSubIntegrator(p.subintegrator_xml, H_)) {}


    // Copy constructor
    LatColMatMinNorm2RecursiveIntegrator(const LatColMatMinNorm2RecursiveIntegrator& l) :
      n_steps(l.n_steps), monomial_list(l.monomial_list), lambda(l.lambda), H_MD(l.H_MD), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMatMinNorm2RecursiveIntegrator(void) {};


    void getMonomialList(multi1d<int>& ml) const {
      AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
	multi1d<LatticeColorMatrix> >& subInt = getSubIntegrator();

      multi1d<int> sub_monomial_list;
      subInt.getMonomialList(sub_monomial_list);
      int my_size = monomial_list.size();
      int sub_size = sub_monomial_list.size();

      ml.resize(my_size+sub_size);
      for(int i=0; i < my_size; i++) {
	ml[i] = monomial_list[i];
      }

      for(int i=0; i < sub_size; i++) { 
	ml[i+my_size] = sub_monomial_list[i];
      }

    }


    void operator()( AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s, 
		     const Real& traj_length) const;
   			    
    AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
			   multi1d<LatticeColorMatrix> >& getSubIntegrator() const {
      return (*SubIntegrator);
    }

  private:
    AbsComponentIntegrator< multi1d<LatticeColorMatrix>,
			    multi1d<LatticeColorMatrix> >* createSubIntegrator(const std::string& subintegrator_xml, 
									       Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) ;

    int  n_steps;
    const multi1d<int> monomial_list;
    Real lambda;

    const Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > H_MD;

    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;
	              

  };

}


#endif
