// -*- C++ -*-
// $Id: extfield_aggregate_w.h,v 1.6 2007-11-02 16:34:31 kostas Exp $
/*! \file
 *  \brief External field functions
 */

#ifndef __extfield_functions_w_h__
#define __extfield_functions_w_h__

#include "handle.h"
#include "actions/ferm/fermstates/extfield.h"


namespace Chroma 
{

  //! Name and registration
  /*! @ingroup fermstates */
  namespace ExternalFieldEnv
  {
    bool registerAll();

    //! Construct zero field
    /*!
     * \ingroup fermstates
     */
    class ZeroExternalField : public ExternalField
    {
    public:
      //! Full constructor
      ZeroExternalField() {}

      //! Return the field
      LatticeComplex operator()(int dir) const;
    };

    struct ConstantMagneticParams{
      int t_dir   ; // the time direction
      int x_dir   ; // the direction which the A_mu is turn on
      int b_dir   ; // the direction of the B field
      Real Bfield ; // The value of the B field
    } ;
    
    //! Construct constant Magnetic field
    /*!
     * \ingroup fermstates
     */
    class ConstantMagneticExternalField : public ExternalField
    {
      int t_dir ;
      int x_dir ;
      int b_dir ;
      Real Bfield ;
    public:
      //! Full constructor
      ConstantMagneticExternalField(ConstantMagneticParams p):t_dir(p.t_dir),
							      x_dir(p.x_dir),
							      b_dir(p.b_dir),
							      Bfield(p.Bfield)
      {}
      //! set time default to be Nd-1 and zero external field
      ConstantMagneticExternalField():t_dir(Nd-1),x_dir(0),b_dir(2),Bfield(0.0)
      {} 

      //! Return the field
      LatticeComplex operator()(int dummy) const;
    };

     Handle< ExternalField > reader(XMLReader& xml,  const std::string& path) ;

  }  // end namespace

  
  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, ExternalFieldEnv::ConstantMagneticParams& param);
  
  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const ExternalFieldEnv::ConstantMagneticParams& param);




#if 0
  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, DerivQuarkDisplacementEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const DerivQuarkDisplacementEnv::Params& param);
#endif

}  // end namespace Chroma

#endif
