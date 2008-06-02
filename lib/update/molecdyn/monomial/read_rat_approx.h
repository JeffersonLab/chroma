// -*- C++ -*-
// $Id: read_rat_approx.h,v 3.1 2008-06-02 15:42:49 bjoo Exp $
/*! @file
 * @brief Remez-type rational approximation
 */

#ifndef __read_rat_approx_h__
#define __read_rat_approx_h__

#include "update/molecdyn/monomial/rat_approx.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup monomial */
  namespace ReadRatApproxEnv
  {
    bool registerAll();

    //! Params for Remez type rational approximation
    /*! @ingroup monomial */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
      RemezCoeff_t pfe;
      RemezCoeff_t ipfe;
    };

  };

  //! Remez type of rational approximations
  /*! @ingroup monomial */
  class ReadRatApprox : public RationalApprox
  {
  public:
    //! Full constructor
    ReadRatApprox(const ReadRatApproxEnv::Params& p) : params(p) {}
    
    //! Destructor
    ~ReadRatApprox() {}
    
    //! Produce the partial-fraction-expansion (PFE) and its inverse (IPFE)
    void operator()(RemezCoeff_t& pfe, RemezCoeff_t& ipfe) const;
    
  private:
    ReadRatApproxEnv::Params  params;   /*!< remez params */
  };


  //! Reader
  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, ReadRatApproxEnv::Params& param);

  //! Reader
  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const ReadRatApproxEnv::Params& params);

} //end namespace chroma

#endif
