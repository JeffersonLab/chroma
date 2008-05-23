// -*- C++ -*-
// $Id: remez_rat_approx.h,v 3.1 2008-05-23 21:31:34 edwards Exp $
/*! @file
 * @brief Remez-type rational approximation
 */

#ifndef __remez_rat_approx_h__
#define __remez_rat_approx_h__

#include "update/molecdyn/monomial/rat_approx.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup monomial */
  namespace RemezRatApproxEnv
  {
    bool registerAll();

    //! Params for Remez type rational approximation
    /*! @ingroup monomial */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      int  numPower;        /*!< Approximate x^-(numPower/denPower) */
      int  denPower;        /*!< Approximate x^-(numPower/denPower) */
      Real lowerMin;        /*!< lower bound of approximation region */
      Real upperMax;        /*!< upper bound of approximation region */
      int  degree;          /*!< degree of approximation */
      int  digitPrecision;  /*!< number of digits used for bigfloat calcs */
    };


    //! Remez type of rational approximations
    /*! @ingroup monomial */
    class RatApprox : public RationalApprox
    {
    public:
      //! Full constructor
      RatApprox(const Params& p) : params(p) {}

      //! Destructor
      ~RatApprox() {}
      
      //! Produce the partial-fraction-expansion (PFE) and its inverse (IPFE)
      void operator()(RemezCoeff_t& pfe, RemezCoeff_t& ipfe) const;

    private:
      Params  params;   /*!< remez params */
    };

  }  // end namespace

  //! Reader
  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, RemezRatApproxEnv::Params& param);

  //! Reader
  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const RemezRatApproxEnv::Params& params);

} //end namespace chroma

#endif
