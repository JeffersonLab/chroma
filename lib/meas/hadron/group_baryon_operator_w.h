// -*- C++ -*-
// $Id: group_baryon_operator_w.h,v 1.1 2006-05-12 03:38:01 edwards Exp $
/*! \file
 *  \brief Construct group baryon operators
 */

#ifndef __group_baryon_operator_w_h__
#define __group_baryon_operator_w_h__

#include "meas/hadron/baryon_operator.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace GroupBaryonOperatorEnv
  {
    extern const bool registered;

  
    //! Group baryon operator
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string operator_coeff_file;      /*!< File holding group coefficients */
      int displacement_length;
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
      multi1d<LatticeColorMatrix> u;

      struct CoeffTerm_t
      {
	multi1d<int>  spin;
	multi1d<int>  disp_dir;
	multi1d<int>  disp_len;
	Complex       coeff;
      };
      multi1d< multi1d< CoeffTerm_t > >  coeffs;
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
