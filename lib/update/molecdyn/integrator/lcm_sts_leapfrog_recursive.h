// -*- C++ -*-
// $Id: lcm_sts_leapfrog_recursive.h,v 3.1 2006-11-07 23:11:05 bjoo Exp $

/*! @file
 * @brief Intgrator for exp(S dt)
 *
 * A component integrator to integrate exp(Sdt) for some 
 * monomial list S (ie this is a leapP like component)
 */

#ifndef LCM_STS_LEAPFROG_RECURSIVE_H
#define LCM_STS_LEAPFROG_RECURSIVE_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"


namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatSTSLeapfrogRecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatSTSLeapfrogRecursiveIntegratorParams
  {
    LatColMatSTSLeapfrogRecursiveIntegratorParams();
    LatColMatSTSLeapfrogRecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    multi1d<int> monomial_list;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatSTSLeapfrogRecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatSTSLeapfrogRecursiveIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatSTSLeapfrogRecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatSTSLeapfrogRecursiveIntegrator(int  n_steps_, 
			      multi1d<int> monomial_list_,
			      Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_ ,
			      Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& SubIntegrator_) : n_steps(n_steps_), monomial_list(monomial_list_), H_MD(H_), SubIntegrator(SubIntegrator_) {};

    // Construct from params struct and Hamiltonian
    LatColMatSTSLeapfrogRecursiveIntegrator(
					    const LatColMatSTSLeapfrogRecursiveIntegratorParams& p,
					    Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) : n_steps(p.n_steps), monomial_list(p.monomial_list), H_MD(H_), SubIntegrator(createSubIntegrator(p.subintegrator_xml, H_)) {}


    // Copy constructor
    LatColMatSTSLeapfrogRecursiveIntegrator(const LatColMatSTSLeapfrogRecursiveIntegrator& l) :
      n_steps(l.n_steps), monomial_list(l.monomial_list), H_MD(l.H_MD), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMatSTSLeapfrogRecursiveIntegrator(void) {};


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
			    multi1d<LatticeColorMatrix> >* createSubIntegrator(const std::string& subintegrator_xml, Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) ;

    int  n_steps;
    const multi1d<int> monomial_list;
    const Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > H_MD;

    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;
	              

    // Hmm should this really be here?
    multi1d<int> total_monomial_list;
    
  };

}


#endif
