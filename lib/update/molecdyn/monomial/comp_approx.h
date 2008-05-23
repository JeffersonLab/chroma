// -*- C++ -*-
// $Id: comp_approx.h,v 3.1 2008-05-23 21:31:32 edwards Exp $
/*! @file
 * @brief Components of rational approximation
 */

#ifndef __comp_approx_h__
#define __comp_approx_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{

  //----------------------------------------------------------------------------
  //! Params for each major group - action/heatbath & force
  /*! @ingroup monomial */
  struct TermApprox_t
  {
    GroupXML_t    ratApprox;     /*!< Rational approximation info for f(M^dag*M) */
    GroupXML_t    invParam;      /*!< Inverter Parameters */
  };

  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, TermApprox_t& param);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const TermApprox_t& params);


  //----------------------------------------------------------------------------
  //! Params for numerator and denominator fermion actions
  /*! @ingroup monomial */
  struct CompApprox_t
  {
    GroupXML_t    fermact;       /*!< Fermion action */
    TermApprox_t  action;        /*!< Action/heatbath params */
    TermApprox_t  force;         /*!< Force params */
  };

  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, CompApprox_t& param);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const CompApprox_t& params);


  //----------------------------------------------------------------------------
  //! Params for numerator and denominator fermion actions
  /*! @ingroup monomial */
  struct CompAction_t
  {
    GroupXML_t    fermact;       /*!< Fermion action */
  };

  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, CompAction_t& param);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const CompAction_t& params);


  //----------------------------------------------------------------------------
  //! Params for numerator and denominator fermion actions
  /*! @ingroup monomial */
  struct CompActionInv_t
  {
    GroupXML_t    fermact;       /*!< Fermion action */
    GroupXML_t    invParam;      /*!< Inverter Parameters */
  };

  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, CompActionInv_t& param);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const CompActionInv_t& params);

} //end namespace chroma

#endif
