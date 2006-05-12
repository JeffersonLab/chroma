// -*- C++ -*-
// $Id: simple_baryon_operator_w.h,v 1.1 2006-05-12 03:38:01 edwards Exp $
/*! \file
 *  \brief Construct simple baryon operators
 */

#ifndef __simple_baryon_operator_w_h__
#define __simple_baryon_operator_w_h__

#include "meas/hadron/baryon_operator.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleBaryonOperatorEnv
  {
    extern const bool registered;

  
    //! Simple baryon operator
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    };


    //! Nucleon with Cg5
    /*! @ingroup hadron
     *
     * Create a simple nucleon
     */
    class BarNuclCg5 : public BaryonOperator<LatticeFermion>
    {
    public:
      //! Full constructor
      BarNuclCg5(const Params& p, const multi1d<LatticeColorMatrix>& u);

      //! Compute the operator
      multi1d<LatticeComplex> operator()(const LatticeFermion& quark1, 
					 const LatticeFermion& quark2, 
					 const LatticeFermion& quark3,
					 enum PlusMinus isign) const;

    private:
      //! Hide partial constructor
      BarNuclCg5() {}

    private:
      Params  params;   /*!< parameters */
      multi1d<LatticeColorMatrix> u;
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, SimpleBaryonOperatorEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const SimpleBaryonOperatorEnv::Params& param);

}  // end namespace Chroma


#endif
