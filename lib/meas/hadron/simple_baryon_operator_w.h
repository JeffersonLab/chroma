// -*- C++ -*-
// $Id: simple_baryon_operator_w.h,v 1.2 2006-05-21 04:40:21 edwards Exp $
/*! \file
 *  \brief Construct simple baryon operators
 */

#ifndef __simple_baryon_operator_w_h__
#define __simple_baryon_operator_w_h__

#include "handle.h"
#include "meas/hadron/baryon_operator.h"
#include "meas/smear/quark_smearing.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleBaryonOperatorEnv
  {
    extern const bool registered;
    extern const std::string name;

  
    //! Simple baryon operator
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string      source_quark_smearing;       /*!< xml string holding source smearing params */
      std::string      source_quark_smearing_type;  /*!< source quark smearing type name */

      std::string      sink_quark_smearing;         /*!< xml string holding sink smearing params */
      std::string      sink_quark_smearing_type;    /*!< sink quark smearing type name */

      std::string      link_smearing;               /*!< link smearing xml */
      std::string      link_smearing_type;          /*!< link smearing type name */
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

    protected:
      //! Manipulate the quark fields
      void quarkManip(multi1d<LatticeFermion>& q,
		      const LatticeFermion& q1, 
		      const LatticeFermion& q2, 
		      const LatticeFermion& q3,
		      enum PlusMinus isign) const;

    private:
      //! Hide partial constructor
      BarNuclCg5() {}

    private:
      Params  params;   /*!< parameters */
      multi1d<LatticeColorMatrix> u_smr;

      Handle< QuarkSmearing<LatticeFermion> > sourceQuarkSmearing;
      Handle< QuarkSmearing<LatticeFermion> > sinkQuarkSmearing;
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
