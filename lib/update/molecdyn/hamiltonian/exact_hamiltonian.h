// -*- C++ -*-
// $Id: exact_hamiltonian.h,v 3.3 2006-11-20 19:12:54 bjoo Exp $
/*! \file
 * \brief Exact Hamiltonians
 */

#ifndef EXACT_HAMILTONIAN_H
#define EXACT_HAMILTONIAN_H

#include "chromabase.h"
#include "handle.h"

#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"

#include "io/xmllog_io.h"
#include "io/monomial_io.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{

  //! Exact Hamiltonian
  /*! @ingroup hamilton */
  class ExactLatColMatHamiltonian : public ExactAbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > 
  {
  public:
    ExactLatColMatHamiltonian(XMLReader& xml, const std::string path) 
    {
      START_CODE();

      XMLReader paramtop(xml,path);
      multi1d< Handle< 
	         ExactMonomial<
 	           multi1d<LatticeColorMatrix>,
	           multi1d<LatticeColorMatrix> 
	         > 
	> > monomial_array;

      try { 
	read(paramtop, "./Monomials", monomial_array);
      }
      catch( const std::string& e ) { 
	QDPIO::cerr << "Error Reading Monomials " << e << endl;
	QDP_abort(1);
      }
      
      QDPIO::cout << "Read " << monomial_array.size() << " monomials" << endl;
      monomials = monomial_array;
      QDPIO::cout << "Leaving function" << endl;
    
      END_CODE();
    }

    // Constructor
    ExactLatColMatHamiltonian( multi1d< Handle< ExactMonomial<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > >& monomials_ ) : monomials(monomials_) {}

    // Copy Constructor
    ExactLatColMatHamiltonian( const ExactLatColMatHamiltonian& H ) : monomials(H.monomials) {}


    //! The Kinetic Energy
    Double mesKE(const AbsFieldState<multi1d<LatticeColorMatrix>,
		 multi1d<LatticeColorMatrix> >& s) const 
    {
      START_CODE();

      // Self Description Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesKE");
	
      // may need to loop over the indices of P?
      Double KE=Double(0);
      for(int mu=0; mu < s.getP().size(); mu++) { 
	KE += norm2((s.getP())[mu]);
      }

      write(xml_out, "KE", KE);
      pop(xml_out);
    
      END_CODE();
      return KE;
    }
    
  protected:
    ExactMonomial<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& getMonomial(int i) const {
      return *(monomials[i]);
    }

    int numMonomials(void) const {
      return monomials.size();
    }

  private:
    multi1d< Handle< ExactMonomial<multi1d<LatticeColorMatrix>,
                                   multi1d<LatticeColorMatrix> > > > monomials;
  };


  //! Parameter structure for new Hamiltonian
  /*! @ingroup hamilton */
  struct NewExactHamParams { 
    
    //! Constructor 
    NewExactHamParams(XMLReader& xml, const std::string& path); 
    multi1d<std::string> monomial_ids; /*!< list of monomial IDs */

  };

  //! Read the parameters for the Hamiltonian
  void read(XMLReader& xml, const std::string& path, NewExactHamParams& p);

  //! Write the parameters for the Hamiltonian
  void write(XMLWriter& xml, const std::string& path, const NewExactHamParams& p);


  //! The Exact Hamiltonian Class - supplies internal field refreshment and energy calculations
  /*! @ingroup hamilton */
  class NewExactHam : public NewAbsHamiltonian< multi1d<LatticeColorMatrix>, 
		      multi1d<LatticeColorMatrix> >
  {
  public:

    //! Construct from a list of string monomial_ids
    NewExactHam(const multi1d<std::string>& monomial_ids_)  {
      create(monomial_ids_);
    }
   
    //! Construct from a parameter structure
    NewExactHam(const NewExactHamParams& p) {
      create(p.monomial_ids);
    }

    //! Copy constructor
    NewExactHam(const NewExactHam& H) : monomials(H.monomials) {}

    //! Destructor 
    ~NewExactHam(void) {}

    //! Internal Field Refreshment 
    void refreshInternalFields(const AbsFieldState<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >& s)  
    { 
      START_CODE();
      for(int i=0; i < monomials.size(); i++) {
	monomials[i]->refreshInternalFields(s);
      }
      END_CODE();
    }
    
    //! The Potential Energy 
    Double mesPE(const AbsFieldState< multi1d<LatticeColorMatrix>,
		        multi1d<LatticeColorMatrix> >& s) const 
    {
      START_CODE();

      // Self Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesPE");
      // Cycle through all the monomials and compute their contribution
      int num_terms = monomials.size();

      write(xml_out, "num_terms", num_terms);
      multi1d<Double> PE_terms(num_terms);
      Double PE = Double(0);
      
      // Caller writes elem rule
      push(xml_out, "PEByMonomials");
      for(int i=0; i < num_terms; i++) 
      {
	push(xml_out, "elem");

	Double tmp = monomials[i]->S(s);
	PE += tmp;
	pop(xml_out);
      }
      pop(xml_out); // PEByMonomials
      pop(xml_out); // pop(mesPE);
      
      END_CODE();
      return PE;
    }

  private:
    //! Convenience 
    typedef ExactMonomial< multi1d<LatticeColorMatrix>, 
			           multi1d<LatticeColorMatrix> >  ExactMon;

    //! This creates the hamiltonian. It is similar to the 
    void create(const multi1d<std::string>& monomial_ids);
    

    multi1d< Handle<ExactMon> >  monomials;

    
  };

}

#endif
