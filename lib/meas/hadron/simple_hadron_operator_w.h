// -*- C++ -*-
// $Id: simple_hadron_operator_w.h,v 1.1 2008-05-09 03:59:01 kostas Exp $
/*! \file
 *  \brief Construct simple hadron operators
 */

#ifndef __simple_hadron_operator_w_h__
#define __simple_hadron_operator_w_h__

#include "handle.h"
#include "meas/smear/quark_smearing.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleHadronOperatorEnv
  {
    extern const std::string name;
    bool registerAll();

  
     //! Construct baryon operators
  /*! @ingroup hadron
   *
   */
  template<typename T>
  class HadronOperator
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~HadronOperator() {}

    //! Construct the operator (do the contractions)
    virtual multi1d<LatticeComplex> operator()(const multi1d<T>& quarks, 
                                               enum PlusMinus isign) const = 0;

    HadronOperator(const GroupXML_t& xml){
      std::istringstream  xml_l(xml.xml);
      XMLReader  xmltop(xml_l);
      read(xmltop,"File",file); 
      read(xmltop,"Name",name); 
      read(xmltop,"Type",type); 
    }
    HadronOperator(){;}
  private:
    string type ;
    string file ;
    string name ;
  };


    //! Baryon Operator
    /*! @ingroup hadron
     *
     * Create a simple baryon
     */
    class Baryon : public HadronOperator<LatticeFermion>{
    public:
      //! Full constructor
      Baryon(const GroupXML_t& p);
      // Just hand in an identifying string that tells us the  diquark operator
      // The permutations are going to be done at post processing.

      //! Compute the operator (no permutations are done)
      multi1d<LatticeComplex> operator()(const multi1d<LatticeFermion>& quarks,
					 enum PlusMinus isign) const;


    private:
      //! Hide partial constructor
      Baryon() {}

    private:
      //      string name ;
      SpinMatrix DiqGamma ;
    };

  }  // end namespace


}  // end namespace Chroma


#endif
