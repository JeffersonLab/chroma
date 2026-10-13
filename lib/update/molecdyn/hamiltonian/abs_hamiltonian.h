// -*- C++ -*-
/*! \file
 * \brief Abstract Hamiltonian
 *
 * Abstract Hamiltonian
 */

#ifndef abs_hamiltonian_h
#define abs_hamiltonian_h

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "io/xmllog_io.h"

namespace Chroma 
{ 


  //! New Abstract Hamiltonian
  /*! @ingroup hamilton
   *
   * Abstraction for Hamiltonians. They can refresh Internal 
   * fields and measure energies.
   */
  template<typename P, typename Q>
  class AbsHamiltonian 
  {
  public:
    
    //! virtual descructor:
    virtual ~AbsHamiltonian() {}

    //! Refresh pseudofermsions (if any)
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& s) =0;

    //! Refresh pseudofermsions (if any) with optional OU parameters
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& s,
                                       const Real& traj_length,
                                       const Real& gamma,
                                       const InternalFieldsRefreshMode mode)
    {
      refreshInternalFields(s);
    }

    //! Push monomial internal fields onto backup stack(s)
    virtual void pushInternalFields(void) {}

    //! Pop monomial internal fields from backup stack(s)
    virtual void popInternalFields(void) {}

    //! Drop monomial internal-field backup without restoring
    virtual void dropInternalFields(void) {}

    //! Save monomial internal fields into a QIO/LIME file
    virtual bool saveInternalFields(const std::string& file,
                                    QDP_volfmt_t volfmt,
                                    QDP_serialparallel_t serpar) const
    {
      return false;
    }

    //! Load monomial internal fields from a QIO/LIME file
    virtual bool loadInternalFields(const std::string& file,
                                    QDP_serialparallel_t serpar)
    {
      return false;
    }
    
    //! Compute the energies 
    //! The total energy
    virtual void  mesE(const AbsFieldState<P,Q>& s, Double& KE, Double& PE) const 
    {
      START_CODE();

      // Self Description Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesE");

      KE = mesKE(s);
      PE = mesPE(s);

      pop(xml_out);

      END_CODE();
    }
        
    //! The Kinetic Energy
    virtual Double mesKE(const AbsFieldState<P,Q>& s) const 
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesKE");

      // Return 1/2 sum pi^2
      // may need to loop over the indices of P?
      Double KE=norm2(s.getP());

      write(xml_out, "KE", KE);
      pop(xml_out);  // pop(mesKE);
    
      END_CODE();
      return KE;
    }
    
    //! The Potential Energy 
    virtual Double mesPE(const AbsFieldState<P,Q>& s) const = 0;
   
  };

} // End namespace Chroma
#endif
