// -*- C++ -*-
// $Id: group_baryon_operator_w.h,v 1.4 2006-05-21 04:40:21 edwards Exp $
/*! \file
 *  \brief Construct group baryon operators
 */

#ifndef __group_baryon_operator_w_h__
#define __group_baryon_operator_w_h__

#include "handle.h"
#include "meas/hadron/baryon_operator.h"
#include "meas/smear/quark_smearing.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace GroupBaryonOperatorEnv
  {
    extern const bool registered;
    extern const std::string name;

  
    //! Group baryon operator
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

      std::string operator_coeff_file;      /*!< File holding group coefficients */
      int   displacement_length;
    };


    //! Nucleon with Cg5
    /*! @ingroup hadron
     *
     * Create a group theoretical construction baryon
     */
    class GroupBaryon : public BaryonOperator<LatticeFermion>
    {
    public:
      //! Full constructor
      GroupBaryon(const Params& p, const multi1d<LatticeColorMatrix>& u);

      //! Compute the operator
      multi1d<LatticeComplex> operator()(const LatticeFermion& quark1, 
					 const LatticeFermion& quark2, 
					 const LatticeFermion& quark3,
					 enum PlusMinus isign) const;

    private:
      //! Hide partial constructor
      GroupBaryon() {}

    private:
      Params  params;   /*!< parameters */
      multi1d<LatticeColorMatrix> u_smr;
      SpinMatrix rotate_mat;

      struct CoeffTerm_t
      {
	struct QuarkTerm_t
	{
	  int  displacement;    /*!< Orig plus/minus 1-based directional displacements */

	  int  spin;            /*!< 0-based spin index */
	  int  disp_dir;        /*!< 0-based direction */
	  int  disp_len;        /*!< 0-based length */
	};

	multi1d<QuarkTerm_t>  quark;    /*!< Displacement and spin for each quark */
	Complex               coeff;    /*!< Weight on color contraction */
      };

      multi1d< multi1d< CoeffTerm_t > >  coeffs;

      Handle< QuarkSmearing<LatticeFermion> > sourceQuarkSmearing;
      Handle< QuarkSmearing<LatticeFermion> > sinkQuarkSmearing;

    protected:
      //! Construct array of maps of displacements
      void displaceQuarks(multi1d< map<int,LatticeFermion> >& disp_quarks,
			  const multi1d<LatticeFermion>& q,
			  enum PlusMinus isign) const;

      //! First displace then smear the quarks
      void displaceSmearQuarks(multi1d< map<int,LatticeFermion> >& disp_quarks,
			       const LatticeFermion& q1, 
			       const LatticeFermion& q2, 
			       const LatticeFermion& q3,
			       enum PlusMinus isign) const;

      //! First smear then displace the quarks
      void smearDisplaceQuarks(multi1d< map<int,LatticeFermion> >& disp_quarks,
			       const LatticeFermion& q1, 
			       const LatticeFermion& q2, 
			       const LatticeFermion& q3,
			       enum PlusMinus isign) const;

      //! Manipulate the quark fields
      void quarkManip(multi1d< map<int,LatticeFermion> >& disp_quarks,
		      const LatticeFermion& q1, 
		      const LatticeFermion& q2, 
		      const LatticeFermion& q3,
		      enum PlusMinus isign) const;

      //! The spin basis matrix to goto Dirac
      const SpinMatrix& rotateMat() const {return rotate_mat;}

      //! Reader
      virtual void readCoeffs(multi1d< multi1d< CoeffTerm_t > >& coef);
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, GroupBaryonOperatorEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const GroupBaryonOperatorEnv::Params& param);

}  // end namespace Chroma


#endif
