// -*- C++ -*-
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



  //! Parameter structure for new Hamiltonian
  /*! @ingroup hamilton */
  struct ExactHamiltonianParams { 
    
    //! Constructor 
    ExactHamiltonianParams(XMLReader& xml, const std::string& path); 
    multi1d<std::string> monomial_ids; /*!< list of monomial IDs */

  };

  //! Read the parameters for the Hamiltonian
  void read(XMLReader& xml, const std::string& path, ExactHamiltonianParams& p);

  //! Write the parameters for the Hamiltonian
  void write(XMLWriter& xml, const std::string& path, const ExactHamiltonianParams& p);


  //! The Exact Hamiltonian Class - supplies internal field refreshment and energy calculations
  /*! @ingroup hamilton */
  class ExactHamiltonian : public AbsHamiltonian< multi1d<LatticeColorMatrix>, 
		      multi1d<LatticeColorMatrix> >
  {
  public:

    //! Construct from a list of std::string monomial_ids
    ExactHamiltonian(const multi1d<std::string>& monomial_ids_)  {
      create(monomial_ids_);
    }
   
    //! Construct from a parameter structure
    ExactHamiltonian(const ExactHamiltonianParams& p) {
      create(p.monomial_ids);
    }

    //! Copy constructor
    ExactHamiltonian(const ExactHamiltonian& H) : monomials(H.monomials) {}

    //! Destructor 
    ~ExactHamiltonian(void) {}

    //! Internal Field Refreshment 
    void refreshInternalFields(const AbsFieldState<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >& s)  
    { 
      START_CODE();
      for(int i=0; i < monomials.size(); i++) {
	monomials[i]->refreshInternalFields(s);
      }
      END_CODE();
    }

    //! Internal field refreshment with OU options
    void refreshInternalFields(const AbsFieldState<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >& s,
                               const Real& traj_length,
                               const Real& gamma,
                               const InternalFieldsRefreshMode mode)
    {
      START_CODE();
      for (int i = 0; i < monomials.size(); ++i) {
        monomials[i]->refreshInternalFields(s, traj_length, gamma, mode);
      }
      END_CODE();
    }

    //! Push monomial internal fields
    void pushInternalFields(void)
    {
      for (int i = 0; i < monomials.size(); ++i) {
        if (monomials[i]->hasInternalFields()) {
          monomials[i]->pushInternalFields();
        }
      }
    }

    //! Pop monomial internal fields
    void popInternalFields(void)
    {
      for (int i = 0; i < monomials.size(); ++i) {
        if (monomials[i]->hasInternalFields()) {
          monomials[i]->popInternalFields();
        }
      }
    }

    //! Drop monomial internal-field backup
    void dropInternalFields(void)
    {
      for (int i = 0; i < monomials.size(); ++i) {
        if (monomials[i]->hasInternalFields()) {
          monomials[i]->dropInternalFields();
        }
      }
    }

    //! Save all pseudofermion internal fields into one QIO/LIME file
    bool saveInternalFields(const std::string& file,
                            QDP_volfmt_t volfmt,
                            QDP_serialparallel_t serpar) const
    {
      START_CODE();

      int num_with_internal_fields = 0;
      for (int i = 0; i < monomials.size(); ++i) {
        if (monomials[i]->hasInternalFields()) {
          ++num_with_internal_fields;
        }
      }

      if (num_with_internal_fields == 0) {
        END_CODE();
        return false;
      }

      XMLBufferWriter file_xml;
      push(file_xml, "SMDPseudoFermions");
      write(file_xml, "NumInternalFieldMonomials", num_with_internal_fields);
      pop(file_xml);

      QDPFileWriter to(file_xml, file, volfmt, serpar, QDPIO_OPEN);
      if (to.bad()) {
        QDPIO::cerr << "ExactHamiltonian::saveInternalFields: error opening " << file << std::endl;
        QDP_abort(1);
      }

      XMLBufferWriter global_meta_xml;
      push(global_meta_xml, "PseudoFermionGlobalMeta");
      write(global_meta_xml, "NumInternalFieldMonomials", num_with_internal_fields);
      pop(global_meta_xml);

      multi1d<int> global_meta(1);
      global_meta[0] = num_with_internal_fields;
      BinaryBufferWriter global_meta_bin;
      write(global_meta_bin, global_meta);
      write(to, global_meta_xml, global_meta_bin);

      for (int i = 0; i < monomials.size(); ++i) {
        if (!monomials[i]->hasInternalFields()) {
          continue;
        }
        if (!monomials[i]->saveInternalFields(to, i)) {
          QDPIO::cerr << "ExactHamiltonian::saveInternalFields: failed on monomial index " << i << std::endl;
          QDP_abort(1);
        }
      }

      close(to);
      END_CODE();
      return true;
    }

    //! Load all pseudofermion internal fields from one QIO/LIME file
    bool loadInternalFields(const std::string& file,
                            QDP_serialparallel_t serpar)
    {
      START_CODE();

      int num_with_internal_fields = 0;
      for (int i = 0; i < monomials.size(); ++i) {
        if (monomials[i]->hasInternalFields()) {
          ++num_with_internal_fields;
        }
      }

      if (num_with_internal_fields == 0) {
        END_CODE();
        return false;
      }

      XMLReader file_xml;
      QDPFileReader from(file_xml, file, serpar);
      if (from.bad()) {
        END_CODE();
        return false;
      }

      XMLReader global_meta_xml;
      multi1d<int> global_meta;
      try {
        BinaryBufferReader global_meta_bin;
        read(from, global_meta_xml, global_meta_bin);
        read(global_meta_bin, global_meta);
      }
      catch (...) {
        close(from);
        END_CODE();
        return false;
      }

      if (global_meta.size() != 1 || global_meta[0] != num_with_internal_fields) {
        close(from);
        END_CODE();
        return false;
      }

      for (int i = 0; i < monomials.size(); ++i) {
        if (!monomials[i]->hasInternalFields()) {
          continue;
        }
        if (!monomials[i]->loadInternalFields(from, i)) {
          close(from);
          END_CODE();
          return false;
        }
      }

      close(from);
      END_CODE();
      return true;
    }
 

    Double mesKE(const AbsFieldState< 
	       multi1d<LatticeColorMatrix>, 
	       multi1d<LatticeColorMatrix> > &s
	       ) const
    {
      START_CODE();

      /* Accumulate KE per site */
      multi1d<LatticeDouble> ke_per_site(Nd);

      
      /* Now add on the local Norm2 of the momenta for each link */
      for(int mu=0; mu < Nd; mu++) { 
	ke_per_site[mu] = -Double(4);
	ke_per_site[mu] += localNorm2(s.getP()[mu]);
      }

      /* Sum up the differences */
      Double KE=zero;
      for(int mu=0; mu < Nd; mu++) { 
	KE += sum(ke_per_site[mu]);
      }

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesKE");
      write(xml_out, "KE", sum(KE));
      pop(xml_out);  // pop(mesKE);


      return KE;

      END_CODE();
    }

    //! The Potential Energy 
    Double  mesPE(const AbsFieldState< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >& s) const 
    {
      START_CODE();

      // Self Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesPE");
      // Cycle through all the monomials and compute their contribution
      int num_terms = monomials.size();

      write(xml_out, "num_terms", num_terms);
      Double PE=zero;

      // Caller writes elem rule
      push(xml_out, "PEByMonomials");
      for(int i=0; i < num_terms; i++) 
      {
	push(xml_out, "elem");
	Double tmp;
	tmp=monomials[i]->S(s);
	PE += tmp;
	pop(xml_out); // elem
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
