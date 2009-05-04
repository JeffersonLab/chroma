// -*- C++ -*-
// $Id: extfield_aggregate_w.h,v 1.9 2009-05-04 17:11:50 caubin Exp $
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
      int y_dir   ; // the direction which the A_mu is turn on
      int b_dir   ; // the direction of the B field
      Real Bfield ; // The value of the B field
      bool patch  ; // do the patch ?
    } ;
    
    //! Construct constant Magnetic field
    /*!
     * \ingroup fermstates
     */
    class ConstantMagneticExternalField : public ExternalField
    {
      int t_dir ;
      int y_dir ;
      int b_dir ;
      Real Bfield ;
      bool patch ;
      int x_dir ;
    public:
      //! Full constructor
      ConstantMagneticExternalField(ConstantMagneticParams p):t_dir(p.t_dir),
							      y_dir(p.y_dir),
							      b_dir(p.b_dir),
							      Bfield(p.Bfield),
							      patch(p.patch)
      {
	for( int d(0);d<Nd;d++)
	  if((d!=t_dir)&&(d!=y_dir)&&(d!=b_dir))
	    x_dir = d ;
      }
      //! set time default to be Nd-1 and zero external field
      ConstantMagneticExternalField():t_dir(Nd-1),y_dir(1),b_dir(2),Bfield(0.0),patch(false),x_dir(0)
      {
      } 

      //! Return the field
      LatticeComplex operator()(int dummy) const;
    };

    
    struct LinearElectricParams{
      int t_dir   ; // the time direction
      int x_dir   ; // the direction which Phi is turned on, also the dir of E
      int x_src   ; // source point fo e field
      Real Efield ; // The E field is E = (Efield)*x, so
	  //this is the e field divided by the distance.
      bool patch ;
    } ;
    //! Construct Linear Electric field (E = epsilon*x)
    /*!
     * \ingroup fermstates
     */
    class LinearElectricExternalField : public ExternalField
    {
      int t_dir ;
      int x_dir ;
      int x_src   ; // source point fo e field
      Real Efield ;
      bool patch ;
    public:
      //! Full constructor
      LinearElectricExternalField(LinearElectricParams p):t_dir(p.t_dir),
							  x_dir(p.x_dir),
							  x_src(p.x_src),
							  Efield(p.Efield)
      {}
      //! set time default to be Nd-1 and zero external field
      LinearElectricExternalField():t_dir(Nd-1),x_dir(0),x_src(0),Efield(0.0)
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


  void read(XMLReader& xml, const string& path, ExternalFieldEnv::LinearElectricParams& param);
  
  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const ExternalFieldEnv::LinearElectricParams& param);


}  // end namespace Chroma

#endif
