// -*- C++ -*-
// $Id: simple_baryon_operator_w.h,v 1.6 2007-04-11 04:55:30 juge Exp $
/*! \file
 *  \brief Construct simple baryon operators
 */

#ifndef __simple_baryon_operator_w_h__
#define __simple_baryon_operator_w_h__

#include "handle.h"
#include "meas/hadron/baryon_operator.h"
#include "meas/smear/quark_smearing.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleBaryonOperatorEnv
  {
    extern const std::string name;
    bool registerAll();

  
    //! Simple baryon operator
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      GroupXML_t       source_quark_smearing;       /*!< xml string holding smearing params */
      GroupXML_t       sink_quark_smearing;         /*!< xml string holding smearing params */
      GroupXML_t       link_smearing;               /*!< link smearing xml */
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

      multi1d<LatticeComplex> operator()(const LatticeFermion& quark1, 
					 const LatticeFermion& quark2, 
					 const LatticeFermion& quark3,
					 int* qindices,
					 enum PlusMinus isign) const
      {
	QDPIO::cerr << "Not implemented yet\n";
      }
      LatticeComplex operator()(const LatticeFermion& quark1, 
				const LatticeFermion& quark2, 
				const LatticeFermion& quark3,
				int* qindices ) const
      {
	QDPIO::cerr << "Not implemented yet\n";
      }

      //! Add this hack for the moment
      LatticeComplex operator()(const LatticeFermion& quark1, 
				const LatticeFermion& quark2, 
				const LatticeFermion& quark3) const
      {
	QDPIO::cerr << "BarNuclCg5 operator not implemented at this moment\n";
	QDP_abort(1);
      }

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
